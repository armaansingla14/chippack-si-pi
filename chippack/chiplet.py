"""2.5D vs 3D chiplet integration trade-off model.

Combines the package models in this repository into a single comparison for a
two-chiplet system (a hot "compute" die and a cooler "cache/IO" die):

- die-to-die interconnect: UCIe-style advanced-package link for 2.5D (pulse-response
  eye + energy model) and a hybrid-bond / micro-bump vertical link for 3D,
- thermal: finite-volume solve of the full stack under a cold plate,
- warpage / stress: 2-D plane-strain FEA of the half cross-section for a cool-down,
- power delivery: PDN impedance and IR drop to the compute die.

Every number is produced by these simplified models with the assumed inputs in
``SystemSpec``; none are measured.
"""
from __future__ import annotations

from dataclasses import dataclass, field, replace

import numpy as np

from . import fea, thermal
from .pdn import PDN, optimise_decaps, target_impedance
from .ucie import ADVANCED, simulate as ucie_simulate


@dataclass
class SystemSpec:
    die_mm: float = 8.0  # both dies are square, this edge length
    die_t_um: float = 700.0
    thinned_die_t_um: float = 50.0  # bottom die thickness in a 3D stack (TSV die)
    p_compute_w: float = 90.0
    p_cache_w: float = 15.0
    package_mm: float = 24.0
    substrate_t_mm: float = 1.0
    interposer_t_um: float = 100.0
    gap_mm: float = 0.2
    vdd: float = 0.75
    ripple: float = 0.05
    d2d_rate_2p5d: float = 16e9
    hybrid_bond_pitch_um: float = 9.0
    microbump_pitch_um: float = 36.0
    d2d_signal_area_mm2: float = 4.0  # area reserved for vertical D2D signals in 3D
    d2d_signal_fraction: float = 0.25  # fraction of bonds in that area carrying signals
    rate_3d: float = 4e9  # per vertical lane
    c_bond_ff: float = 8.0  # pad + via + RX input per vertical lane (assumed)
    v_d2d: float = 0.4
    h_coldplate: float = 2.5e4
    dT_cooldown: float = -100.0  # e.g. underfill cure / reflow to room temperature
    tsv_r_mohm_total: float = 0.03  # aggregate power-TSV + bond resistance into the top die
    tsv_l_ph_total: float = 0.5


K_SI, K_SUB, K_TIM, K_LID, K_UF, K_BOND = 150.0, 20.0, 5.0, 390.0, 0.6, 2.0


def interconnect(spec: SystemSpec, style: str) -> dict:
    if style == "2.5D":
        r = ucie_simulate(ADVANCED, spec.d2d_rate_2p5d)
        shoreline_mm = spec.die_mm  # one die edge facing the neighbour
        return {"d2d_reach_mm": ADVANCED.reach_mm, "d2d_bw_tbps": r["bw_density_gbps_per_mm"] * shoreline_mm / 1e3,
                "d2d_energy_pj_bit": r["energy_pj_per_bit"], "d2d_eye_mv": r["eye_height_mv"]}
    pitch = spec.hybrid_bond_pitch_um if style == "3D-hybrid" else spec.microbump_pitch_um
    n_bonds = spec.d2d_signal_area_mm2 * 1e6 / pitch ** 2
    lanes = n_bonds * spec.d2d_signal_fraction
    c = spec.c_bond_ff * (1.0 if style == "3D-hybrid" else 2.5) * 1e-15  # micro-bumps: larger pads
    e = 0.25 * c * spec.v_d2d ** 2 * 1e12 + 0.02
    return {"d2d_reach_mm": spec.thinned_die_t_um * 1e-3, "d2d_bw_tbps": lanes * spec.rate_3d / 1e12,
            "d2d_energy_pj_bit": e, "d2d_eye_mv": None}


def _thermal(spec: SystemSpec, style: str, nx: int = 36) -> dict:
    L = spec.package_mm * 1e-3
    d = spec.die_mm * 1e-3
    c0 = L / 2
    if style == "2.5D":
        x_a = c0 - spec.gap_mm * 1e-3 / 2 - d
        x_b = c0 + spec.gap_mm * 1e-3 / 2
        y0 = c0 - d / 2
        compute = thermal.block_map(nx, nx, L, L, x_a, y0, d, d, spec.p_compute_w)
        cache = thermal.block_map(nx, nx, L, L, x_b, y0, d, d, spec.p_cache_w)
        occ = (thermal.block_map(nx, nx, L, L, x_a, y0, 2 * d + spec.gap_mm * 1e-3, d, 1.0) > 0)
        k_die = np.where(occ, K_SI, K_UF)  # mold/underfill between and around dies
        layers = [thermal.Layer("substrate", spec.substrate_t_mm * 1e-3, K_SUB, 2),
                  thermal.Layer("interposer", spec.interposer_t_um * 1e-6, K_SI, 1),
                  thermal.Layer("dies", spec.die_t_um * 1e-6, k_die, 3),
                  thermal.Layer("tim", 50e-6, K_TIM, 1),
                  thermal.Layer("lid", 1.5e-3, K_LID, 2)]
        power = {2: compute + cache}
        m = thermal.ThermalModel(L, L, nx, nx, layers, h_top=spec.h_coldplate, h_bottom=20.0,
                                 coolant=thermal.Coolant(mass_flow=0.02), power=power)
        r = m.solve()
        sel_c = compute > 0
        sel_m = cache > 0
        tdie = r.t[:, :, r.cell_layer == 2].max(axis=2)
        return {"tj_compute_c": float(tdie[sel_c].max()), "tj_cache_c": float(tdie[sel_m].max()),
                "footprint_mm2": (2 * spec.die_mm + spec.gap_mm) * spec.die_mm + 0.0,
                "energy_balance_err": r.energy_balance_error}
    # 3D: one die footprint, compute on top (next to cold plate) or bottom
    x0 = c0 - d / 2
    occ = thermal.block_map(nx, nx, L, L, x0, x0, d, d, 1.0) > 0
    k_die = np.where(occ, K_SI, K_UF)
    bond_k = np.where(occ, K_BOND if style.endswith("hybrid") else 1.0, K_UF)
    top_is_compute = not style.endswith("cache-top")
    p_top = spec.p_compute_w if top_is_compute else spec.p_cache_w
    p_bot = spec.p_cache_w if top_is_compute else spec.p_compute_w
    layers = [thermal.Layer("substrate", spec.substrate_t_mm * 1e-3, K_SUB, 2),
              thermal.Layer("bottom die", spec.thinned_die_t_um * 1e-6, k_die, 1),
              thermal.Layer("bond", (3e-6 if style.endswith("hybrid") else 25e-6), bond_k, 1),
              thermal.Layer("top die", spec.die_t_um * 1e-6, k_die, 3),
              thermal.Layer("tim", 50e-6, K_TIM, 1),
              thermal.Layer("lid", 1.5e-3, K_LID, 2)]
    power = {1: thermal.block_map(nx, nx, L, L, x0, x0, d, d, p_bot),
             3: thermal.block_map(nx, nx, L, L, x0, x0, d, d, p_top)}
    m = thermal.ThermalModel(L, L, nx, nx, layers, h_top=spec.h_coldplate, h_bottom=20.0,
                             coolant=thermal.Coolant(mass_flow=0.02), power=power)
    r = m.solve()
    t_top = r.layer_max(3)
    t_bot = r.layer_max(1)
    return {"tj_compute_c": t_top if top_is_compute else t_bot, "tj_cache_c": t_bot if top_is_compute else t_top,
            "footprint_mm2": spec.die_mm ** 2, "energy_balance_err": r.energy_balance_error}


def _warpage(spec: SystemSpec, style: str) -> dict:
    half_pkg = spec.package_mm / 2 * 1e-3
    t_sub = spec.substrate_t_mm * 1e-3
    uf = 50e-6
    blocks = [fea.Block(fea.ORGANIC_SUB, 0, half_pkg, 0, t_sub)]
    d = spec.die_mm * 1e-3
    if style == "2.5D":
        w_int = d + spec.gap_mm * 1e-3 / 2 + 0.5e-3
        ti = spec.interposer_t_um * 1e-6
        y = t_sub
        blocks += [fea.Block(fea.UNDERFILL, 0, w_int, y, y + uf),
                   fea.Block(fea.SILICON, 0, w_int, y + uf, y + uf + ti)]
        y += uf + ti
        x_die = spec.gap_mm * 1e-3 / 2
        blocks += [fea.Block(fea.UNDERFILL, 0, w_int, y, y + 20e-6),
                   fea.Block(fea.SILICON, x_die, x_die + d, y + 20e-6, y + 20e-6 + spec.die_t_um * 1e-6)]
    else:
        y = t_sub
        tb = spec.thinned_die_t_um * 1e-6
        blocks += [fea.Block(fea.UNDERFILL, 0, d / 2, y, y + uf),
                   fea.Block(fea.SILICON, 0, d / 2, y + uf, y + uf + tb)]
        y += uf + tb
        bond = 3e-6 if style.endswith("hybrid") else 25e-6
        bond_mat = fea.SILICON if style.endswith("hybrid") else fea.UNDERFILL
        blocks += [fea.Block(bond_mat, 0, d / 2, y, y + bond),
                   fea.Block(fea.SILICON, 0, d / 2, y + bond, y + bond + spec.die_t_um * 1e-6)]
    r = fea.solve(blocks, spec.dT_cooldown, nx=100, y_breaks_n=3, plane="strain")
    mats = sorted({b.material for b in blocks}, key=lambda m: [b.material for b in blocks].index(m))
    si_idx = mats.index(fea.SILICON)
    # half-symmetric model: the deflection range along the substrate bottom is the
    # centre-to-edge warpage of the full package cross-section
    return {"warpage_um": r.warpage(0) * 1e6,
            "max_si_von_mises_mpa": float(r.von_mises[r.elem_material == si_idx].max() / 1e6)}


def _pdn(spec: SystemSpec, style: str) -> dict:
    f = np.logspace(3, np.log10(2e8), 200)
    i_dc = spec.p_compute_w / spec.vdd
    zt = target_impedance(spec.vdd, spec.ripple, 0.5 * i_dc)  # 50 % load step
    # scale the default (10 A class) PDN to this current: more bumps/balls/planes in parallel
    # and on-die capacitance proportional to the load
    k = i_dc / 10.0
    d = PDN()
    base = PDN(r_vrm=d.r_vrm / k, l_vrm=d.l_vrm / k, r_bga=d.r_bga / k, l_bga=d.l_bga / k, r_pkg=d.r_pkg / k,
               l_pkg=d.l_pkg / k, c_die=d.c_die * k, r_die=d.r_die / k)
    if style == "2.5D":
        base = replace(base, r_pkg=base.r_pkg + 0.02e-3, l_pkg=base.l_pkg + 0.3e-12)  # interposer TSVs/RDL
    elif style in ("3D-hybrid", "3D-microbump"):  # compute die on top: power crosses the bottom die
        base = replace(base, r_pkg=base.r_pkg + spec.tsv_r_mohm_total * 1e-3, l_pkg=base.l_pkg + spec.tsv_l_ph_total * 1e-12)
    opt = optimise_decaps(base, f, zt)
    r_path = base.r_vrm + base.r_bga + base.r_pkg
    return {"pdn_ztarget_mohm": zt * 1e3, "pdn_zmax_mohm": opt.max_z_ohm * 1e3, "pdn_meets": opt.meets_target,
            "decap_cost": opt.cost, "ir_drop_mv": i_dc * r_path * 1e3}


def evaluate(spec: SystemSpec, style: str) -> dict:
    out = {"config": style}
    base_style = style.replace("-cache-top", "")
    out.update(interconnect(spec, base_style))
    out.update(_thermal(spec, style))
    out.update(_warpage(spec, base_style))
    out.update(_pdn(spec, style))
    return out


STYLES = ("2.5D", "3D-microbump", "3D-hybrid", "3D-hybrid-cache-top")
