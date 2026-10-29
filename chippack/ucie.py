"""UCIe-style die-to-die link models (standard vs advanced package).

Package and PHY parameters below are representative assumptions loosely based on
public UCIe 1.x descriptions (bump pitch ranges, lanes per module, reach) and
typical organic-substrate / silicon-bridge trace properties. They are model
inputs, not specification values.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .network import Network, abcd_to_s2, cascade, rlgc_line
from .serdes import analyze_eye, pulse_response


@dataclass
class PackageClass:
    name: str
    reach_mm: float
    bump_pitch_um: float
    lanes_per_module: int  # data lanes, one direction
    module_width_um: float  # shoreline per module (assumed)
    # per-metre trace model
    r_dc: float
    r_skin: float
    l: float
    c: float
    tan_d: float
    pad_c_ff: float  # TX/RX pad + ESD capacitance per side
    terminated: bool
    driver_ohm: float = 30.0
    rx_term_ohm: float = 50.0
    energy_fixed_pj: float = 0.15  # serializer/clocking overhead per bit (assumed)


STANDARD = PackageClass("standard (organic substrate)", reach_mm=10.0, bump_pitch_um=110.0, lanes_per_module=16,
                        module_width_um=1100.0, r_dc=8.0, r_skin=1.5e-3, l=3.3e-7, c=1.3e-10, tan_d=0.015,
                        pad_c_ff=250.0, terminated=True)
ADVANCED = PackageClass("advanced (Si bridge/interposer)", reach_mm=2.0, bump_pitch_um=45.0, lanes_per_module=64,
                        module_width_um=388.0, r_dc=2.5e4, r_skin=0.0, l=4.0e-7, c=2.0e-10, tan_d=0.005,
                        pad_c_ff=120.0, terminated=False, driver_ohm=50.0, rx_term_ohm=1e6,
                        energy_fixed_pj=0.10)


def link_transfer(pkg: PackageClass, f: np.ndarray) -> np.ndarray:
    """Voltage transfer V_rx / V_source including driver impedance, pad C and RX termination."""
    line = rlgc_line(f, pkg.reach_mm * 1e-3, r_dc=pkg.r_dc, r_skin=pkg.r_skin, l=pkg.l, g_tan=pkg.tan_d, c=pkg.c)
    w = 2 * np.pi * np.maximum(f, 1.0)
    from .network import s2_to_abcd
    a_line = s2_to_abcd(line.s, line.z0)
    ycap = 1j * w * pkg.pad_c_ff * 1e-15
    shunt = np.zeros((f.size, 2, 2), complex)
    shunt[:, 0, 0] = shunt[:, 1, 1] = 1
    shunt[:, 1, 0] = ycap
    a_tot = shunt @ a_line @ shunt
    zl = pkg.rx_term_ohm
    zs = pkg.driver_ohm
    A, B, C, D = a_tot[:, 0, 0], a_tot[:, 0, 1], a_tot[:, 1, 0], a_tot[:, 1, 1]
    # V_L / V_s for source impedance zs and load zl
    return zl / (A * zl + B + zs * (C * zl + D))


def simulate(pkg: PackageClass, data_rate: float, noise_rms: float = 5e-3, swing: float = 0.4, ber: float = 1e-15):
    f = np.linspace(10e6, max(4 * data_rate, 40e9), 2000)
    h = link_transfer(pkg, f)
    spu = 32
    _, p = pulse_response(f, h, data_rate, spu=spu, rise_ui=0.3)
    eye = analyze_eye(p, spu, noise_rms=noise_rms, ber=ber, swing=swing)
    # energy/bit: dynamic CV^2 of line + pads (unterminated) or termination current (terminated)
    c_tot = pkg.c * pkg.reach_mm * 1e-3 + 2 * pkg.pad_c_ff * 1e-15
    if pkg.terminated:
        p_static = 0.5 * (swing ** 2) / (pkg.driver_ohm + pkg.rx_term_ohm)  # ~half of time driving
        e_bit = p_static / data_rate + 0.25 * c_tot * swing ** 2
    else:
        e_bit = 0.25 * c_tot * swing ** 2  # activity factor 0.5, CV^2/2
    e_bit_pj = e_bit * 1e12 + pkg.energy_fixed_pj
    bw_density = pkg.lanes_per_module * data_rate / (pkg.module_width_um * 1e-3) / 1e9  # Gbps per mm
    return {
        "package": pkg.name, "data_rate_gtps": data_rate / 1e9, "reach_mm": pkg.reach_mm,
        "loss_at_nyq_db": float(-20 * np.log10(np.abs(np.interp(data_rate / 2, f, np.abs(h))) / np.abs(h[0]))),
        "eye_height_mv": 1e3 * eye.stat_height, "eye_width_ui": eye.stat_width_ui,
        "pd_eye_height_mv": 1e3 * eye.peak_distortion_height,
        "energy_pj_per_bit": e_bit_pj, "bw_density_gbps_per_mm": bw_density,
    }
