import numpy as np
import pytest

from chippack.pdn import (LIBRARY, PDN, Capacitor, antiresonances, baseline_equal_mix, baseline_uniform,
                          optimise_decaps, target_impedance, worst_violation)

F = np.logspace(3, np.log10(2e8), 300)


def test_target_impedance():
    assert np.isclose(target_impedance(1.0, 0.05, 10.0), 5e-3)


@pytest.mark.parametrize("cap", LIBRARY)
def test_capacitor_srf_minimum(cap):
    f = np.logspace(4, 10, 4000)
    z = np.abs(cap.z(f))
    assert np.isclose(f[np.argmin(z)], cap.srf_hz, rtol=0.02)
    assert np.isclose(z.min(), cap.esr, rtol=0.05)


def test_pdn_low_frequency_is_resistive():
    p = PDN()
    z = p.impedance(np.array([1.0]))
    assert np.isclose(abs(z[0]), p.r_vrm + p.r_bga + p.r_pkg, rtol=0.05)


def test_pdn_high_frequency_set_by_die_cap():
    p = PDN()
    f = np.array([5e9])
    assert abs(p.impedance(f)[0]) < abs(p.r_die + 1 / (2j * np.pi * f[0] * p.c_die)) * 1.01 + 0.01


def test_adding_caps_lowers_vrm_inductive_region():
    # around 1 MHz the bare PDN is dominated by the VRM inductance resonating with the die cap
    f = np.array([1e6])
    a = np.abs(PDN().impedance(f))
    b = np.abs(PDN(caps={LIBRARY[4]: 10}).impedance(f))
    assert b[0] < a[0] / 10


def test_antiresonance_found():
    p = PDN()
    peaks = antiresonances(F, p.impedance(F))
    assert peaks and all(z > 0 for _, z in peaks)


@pytest.mark.parametrize("i_step", [3.0, 5.0])
def test_optimiser_meets_target_cheaper_than_baseline(i_step):
    zt = target_impedance(0.75, 0.05, i_step)
    opt = optimise_decaps(PDN(), F, zt)
    assert opt.meets_target and worst_violation(PDN(caps={c: opt.counts.get(c.name, 0) for c in LIBRARY}), F, zt) <= 1
    mix = baseline_equal_mix(PDN(), F, zt)
    assert mix.meets_target
    assert opt.cost <= mix.cost


def test_single_type_cannot_meet_tight_target():
    zt = target_impedance(0.75, 0.05, 5.0)
    assert not baseline_uniform(PDN(), F, zt, LIBRARY[0], max_n=50).meets_target


def test_optimiser_reports_failure_when_impossible():
    zt = 1e-6
    r = optimise_decaps(PDN(), F, zt, max_iter=20)
    assert not r.meets_target
