import pytest

from chippack.chiplet import STYLES, SystemSpec, evaluate, interconnect


@pytest.fixture(scope="module")
def results():
    return {s: evaluate(SystemSpec(), s) for s in STYLES}


def test_all_styles_evaluate(results):
    for r in results.values():
        assert r["energy_balance_err"] < 1e-9
        assert r["warpage_um"] > 0 and r["tj_compute_c"] > 25


def test_3d_smaller_footprint(results):
    assert results["3D-hybrid"]["footprint_mm2"] < results["2.5D"]["footprint_mm2"]


def test_hybrid_bond_bandwidth_beats_microbump():
    s = SystemSpec()
    assert interconnect(s, "3D-hybrid")["d2d_bw_tbps"] > interconnect(s, "3D-microbump")["d2d_bw_tbps"]


def test_compute_on_top_runs_cooler(results):
    assert results["3D-hybrid"]["tj_compute_c"] < results["3D-hybrid-cache-top"]["tj_compute_c"]


def test_more_power_hotter():
    a = evaluate(SystemSpec(p_compute_w=60), "2.5D")["tj_compute_c"]
    b = evaluate(SystemSpec(p_compute_w=120), "2.5D")["tj_compute_c"]
    assert b > a
