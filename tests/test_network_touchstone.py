import numpy as np
import pytest

from chippack import touchstone as ts
from chippack.network import (Network, abcd_to_s2, cascade, causality_check, coupled_line_4port, enforce_passivity,
                              mixed_mode, passivity_violation, reciprocity_error, rlgc_line, s2_to_abcd, s_to_t,
                              s_to_y, s_to_z, t_to_s, z_to_s)

F = np.linspace(10e6, 40e9, 400)


@pytest.mark.parametrize("length_mm", [1, 5, 20])
def test_line_passive_reciprocal(length_mm):
    n = rlgc_line(F, length_mm * 1e-3)
    ok, worst = passivity_violation(n)
    assert ok and worst < 1
    assert reciprocity_error(n) < 1e-12


@pytest.mark.parametrize("a,b", [(1, 2), (5, 5), (3, 17)])
def test_cascade_lengths_add(a, b):
    c = cascade(rlgc_line(F, a * 1e-3), rlgc_line(F, b * 1e-3))
    assert np.allclose(c.s, rlgc_line(F, (a + b) * 1e-3).s, atol=1e-10)


def test_cascade_4port_lengths_add():
    p = coupled_line_4port(F, 0.01).reorder([1, 3, 2, 4])
    q = coupled_line_4port(F, 0.02).reorder([1, 3, 2, 4])
    assert np.allclose(cascade(p, p).s, q.s, atol=1e-10)


def test_cascade_identity_thru():
    thru = Network(F, np.tile(np.array([[0, 1], [1, 0]], complex), (F.size, 1, 1)))
    n = rlgc_line(F, 0.01)
    assert np.allclose(cascade(thru, n).s, n.s)


def test_cascade_mismatch_raises():
    with pytest.raises(ValueError):
        cascade(rlgc_line(F, 0.01), coupled_line_4port(F, 0.01))


@pytest.mark.parametrize("z0", [25.0, 50.0, 100.0])
def test_s_z_roundtrip(z0, rng):
    s = 0.3 * (rng.standard_normal((5, 3, 3)) + 1j * rng.standard_normal((5, 3, 3)))
    assert np.allclose(z_to_s(s_to_z(s, z0), z0), s)
    assert np.allclose(s_to_y(s, z0) @ s_to_z(s, z0), np.eye(3))


def test_abcd_roundtrip():
    n = rlgc_line(F, 0.01)
    assert np.allclose(abcd_to_s2(s2_to_abcd(n.s, 50), 50), n.s)


def test_abcd_series_resistor():
    m = np.array([[[1, 50.0], [0, 1]]], complex)
    s = abcd_to_s2(m, 50.0)
    assert np.isclose(s[0, 0, 0], 1 / 3) and np.isclose(s[0, 1, 0], 2 / 3)


@pytest.mark.parametrize("n", [1, 2])
def test_t_roundtrip(n, rng):
    s = 0.2 * (rng.standard_normal((4, 2 * n, 2 * n)) + 1j * rng.standard_normal((4, 2 * n, 2 * n)))
    s[:, n:, :n] += 0.7 * np.eye(n)
    assert np.allclose(t_to_s(s_to_t(s)), s)


def test_renormalize_roundtrip():
    n = rlgc_line(F, 0.01)
    assert np.allclose(n.renormalize(40).renormalize(50).s, n.s)


def test_mixed_mode_symmetric_pair():
    p = coupled_line_4port(F, 0.02)
    mm = mixed_mode(p)
    assert np.max(np.abs(mm["cd"])) < 1e-12 and np.max(np.abs(mm["dc"])) < 1e-12
    # differential mode sees the odd-mode impedance: ~ |Sdd21| close to 1 at low frequency
    assert abs(mm["dd"][0, 1, 0]) > 0.95


def test_mixed_mode_detects_skew():
    p = coupled_line_4port(F, 0.02)
    s = p.s.copy()
    s[:, 3, 2] *= np.exp(-2j * np.pi * F * 5e-12)  # 5 ps skew on the negative leg
    s[:, 2, 3] = s[:, 3, 2]
    mm = mixed_mode(Network(F, s))
    assert np.max(np.abs(mm["cd"][:, 1, 0])) > 0.05


def test_mixed_mode_requires_4port():
    with pytest.raises(ValueError):
        mixed_mode(rlgc_line(F, 0.01))


def test_passivity_violation_detected_and_fixed():
    n = rlgc_line(F, 0.001)
    bad = Network(F, n.s * 1.2)
    ok, worst = passivity_violation(bad)
    assert not ok and worst > 1
    ok2, worst2 = passivity_violation(enforce_passivity(bad))
    assert ok2 and worst2 <= 1


@pytest.mark.parametrize("length_mm", [5, 20, 50])
def test_causal_line_passes(length_mm):
    n = rlgc_line(F, length_mm * 1e-3)
    ok, ratio = causality_check(n.sij(2, 1), F)
    assert ok, ratio


@pytest.mark.parametrize("advance_ps", [200, 500])
def test_noncausal_detected(advance_ps):
    n = rlgc_line(F, 0.005)
    h = n.sij(2, 1) * np.exp(2j * np.pi * F * advance_ps * 1e-12)
    ok, ratio = causality_check(h, F)
    assert not ok


def test_bad_shape():
    with pytest.raises(ValueError):
        Network(F, np.zeros((3, 2, 2)))


@pytest.mark.parametrize("fmt", ["RI", "MA", "DB"])
@pytest.mark.parametrize("unit", ["HZ", "MHZ", "GHZ"])
@pytest.mark.parametrize("ports", [2, 4])
def test_touchstone_roundtrip(tmp_path, fmt, unit, ports):
    n = rlgc_line(F, 0.01) if ports == 2 else coupled_line_4port(F, 0.01)
    path = str(tmp_path / f"x.s{ports}p")
    ts.write(path, n, fmt, unit, comment="test")
    m = ts.read(path)
    assert np.allclose(m.f, n.f, rtol=1e-8)
    assert np.allclose(m.s, n.s, atol=1e-9)


def test_touchstone_2port_order(tmp_path):
    path = tmp_path / "a.s2p"
    path.write_text("! test\n# GHZ S RI R 50\n1.0 0.1 0 0.9 0 0.2 0 0.3 0\n")
    n = ts.read(str(path))
    assert np.isclose(n.s[0, 1, 0], 0.9) and np.isclose(n.s[0, 0, 1], 0.2)
    assert n.f[0] == 1e9


def test_touchstone_1port_and_z0(tmp_path):
    path = tmp_path / "a.s1p"
    path.write_text("# MHZ S MA R 75\n100 0.5 90\n200 0.25 -90\n")
    n = ts.read(str(path))
    assert n.z0 == 75 and np.isclose(n.s[0, 0, 0], 0.5j)


def test_touchstone_bad_option(tmp_path):
    path = tmp_path / "a.s1p"
    path.write_text("# GHZ Y RI R 50\n1 0 0\n")
    with pytest.raises(ValueError):
        ts.read(str(path))


def test_interpolate():
    n = rlgc_line(F, 0.01)
    m = n.interpolate(F[::3])
    assert np.allclose(m.s, n.s[::3], atol=1e-9)
