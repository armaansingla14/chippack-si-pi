import numpy as np
import pytest

from chippack.network import rlgc_line
from chippack.serdes import analyze_eye, ddj_from_pulse, extrapolate_to_dc, impulse_response, pulse_response
from chippack.ucie import ADVANCED, STANDARD, link_transfer, simulate

F = np.linspace(10e6, 100e9, 4000)


def test_extrapolate_to_dc():
    f = np.arange(4, 400) * 25e6
    f2, h2 = extrapolate_to_dc(f, np.exp(-2j * np.pi * f * 1e-10))
    assert f2[0] == 0 and np.isclose(f2[1] - f2[0], 25e6) and f2.size == 400
    assert np.allclose(h2[:4], np.exp(-2j * np.pi * f2[:4] * 1e-10))


def test_ideal_channel_pulse_is_rectangle():
    _, p = pulse_response(F, np.ones(F.size, complex), 10e9, spu=32)
    assert np.isclose(p.max(), 1.0, atol=0.06)  # band-limited: small Gibbs overshoot
    eye = analyze_eye(p, 32)
    assert eye.stat_height > 1.9 and eye.isi_sum < 0.05


def test_impulse_area_equals_dc_gain():
    h = 0.5 * np.exp(-2j * np.pi * F * 50e-12)
    dt, ht = impulse_response(F, h)
    assert np.isclose(ht.sum() * dt, 0.5, atol=0.02)


@pytest.mark.parametrize("rate", [8e9, 16e9, 32e9])
def test_longer_line_smaller_eye(rate):
    def eye(length):
        n = rlgc_line(F, length)
        _, p = pulse_response(F, n.sij(2, 1), rate)
        return analyze_eye(p, 32).stat_height
    assert eye(0.05) < eye(0.01)


def test_noise_closes_eye():
    n = rlgc_line(F, 0.02)
    _, p = pulse_response(F, n.sij(2, 1), 16e9)
    a = analyze_eye(p, 32, noise_rms=0.0)
    b = analyze_eye(p, 32, noise_rms=0.02)
    assert b.stat_height < a.stat_height
    assert a.peak_distortion_height <= a.stat_height + 1e-9


def test_lower_ber_target_smaller_eye():
    n = rlgc_line(F, 0.03)
    _, p = pulse_response(F, n.sij(2, 1), 16e9)
    assert analyze_eye(p, 32, 0.01, ber=1e-15).stat_height < analyze_eye(p, 32, 0.01, ber=1e-3).stat_height


def test_ddj_grows_with_loss():
    vals = []
    for length in (0.02, 0.2):
        n = rlgc_line(F, length)
        _, p = pulse_response(F, n.sij(2, 1), 16e9)
        vals.append(ddj_from_pulse(p, 32, 6))
    assert vals[1] > vals[0] >= 0


@pytest.mark.parametrize("pkg", [STANDARD, ADVANCED])
def test_link_transfer_dc(pkg):
    h = link_transfer(pkg, np.array([1e6]))
    zl, zs = pkg.rx_term_ohm, pkg.driver_ohm
    r_line = pkg.r_dc * pkg.reach_mm * 1e-3
    assert np.isclose(abs(h[0]), zl / (zl + zs + r_line), rtol=0.02)


@pytest.mark.parametrize("pkg", [STANDARD, ADVANCED])
def test_ucie_eye_shrinks_with_rate(pkg):
    lo = simulate(pkg, 4e9)
    hi = simulate(pkg, 32e9)
    assert hi["eye_height_mv"] < lo["eye_height_mv"]
    assert hi["bw_density_gbps_per_mm"] == pytest.approx(8 * lo["bw_density_gbps_per_mm"])


def test_advanced_package_denser_and_cheaper_per_bit():
    s = simulate(STANDARD, 16e9)
    a = simulate(ADVANCED, 16e9)
    assert a["bw_density_gbps_per_mm"] > s["bw_density_gbps_per_mm"]
    assert a["energy_pj_per_bit"] < s["energy_pj_per_bit"]
