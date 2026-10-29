import numpy as np
import pytest

from chippack import fea, thermal


@pytest.mark.parametrize("h", [1e3, 1e4, 5e4])
@pytest.mark.parametrize("power", [1.0, 25.0])
def test_1d_stack_matches_analytic(h, power):
    L = 0.01
    layers = [thermal.Layer("si", 500e-6, 150, 4), thermal.Layer("tim", 50e-6, 5, 2), thermal.Layer("cu", 1e-3, 390, 3)]
    m = thermal.ThermalModel(L, L, 8, 8, layers, h_top=h, h_bottom=1e-9, power={0: np.full((8, 8), power / 64)})
    r = m.solve()
    A = L * L
    # heat generated uniformly in the Si layer; the peak is half a Si-layer resistance below its bottom cell
    expected_top = 25 + power / (h * A) + power * (1e-3 / 390 + 50e-6 / 5) / A
    expected_max = expected_top + power * 500e-6 / (2 * 150 * A)
    assert np.isclose(r.t_max, expected_max, rtol=1e-6)
    assert r.energy_balance_error < 1e-9


@pytest.mark.parametrize("nx", [10, 20])
def test_energy_balance_with_hotspot(nx):
    L = 0.02
    pm = thermal.block_map(nx, nx, L, L, 0.005, 0.005, 0.004, 0.004, 40.0)
    assert np.isclose(pm.sum(), 40.0)
    m = thermal.ThermalModel(L, L, nx, nx, [thermal.Layer("sub", 1e-3, 20, 2), thermal.Layer("si", 7e-4, 150, 2)],
                             h_top=1e4, h_bottom=50, power={1: pm})
    r = m.solve()
    assert r.energy_balance_error < 1e-9
    assert r.t.min() >= 25 - 1e-9


def test_coolant_energy_balance():
    L = 0.01
    m = thermal.ThermalModel(L, L, 10, 10, [thermal.Layer("si", 5e-4, 150, 2)], h_top=1e4, h_bottom=1e-9,
                             coolant=thermal.Coolant(mass_flow=0.002, t_in=30.0), power={0: np.full((10, 10), 0.2)})
    r = m.solve()
    t_out = r.t_fluid[-1].mean()
    assert np.isclose(t_out, 30 + 20.0 / (0.002 * 4180), rtol=1e-6)
    assert np.all(np.diff(r.t_fluid.mean(axis=1)) > 0)  # coolant warms along the flow


def test_higher_conductivity_lower_peak():
    L = 0.01
    pm = thermal.block_map(16, 16, L, L, 0.004, 0.004, 0.002, 0.002, 10.0)
    temps = []
    for k in (20, 150):
        m = thermal.ThermalModel(L, L, 16, 16, [thermal.Layer("x", 5e-4, k, 2)], h_top=1e4, h_bottom=1e-9,
                                 power={0: pm})
        temps.append(m.solve().t_max)
    assert temps[1] < temps[0]


def test_layer_map_conductivity():
    L = 0.01
    k = np.full((8, 8), 150.0)
    k[:, :4] = 1.0
    m = thermal.ThermalModel(L, L, 8, 8, [thermal.Layer("x", 5e-4, k, 2)], power={0: np.full((8, 8), 0.1)})
    r = m.solve()
    assert r.energy_balance_error < 1e-9


@pytest.mark.parametrize("plane", ["stress", "strain"])
@pytest.mark.parametrize("ratio", [0.5, 1.0, 2.0])
def test_bimetal_matches_timoshenko(plane, ratio):
    L, t2 = 10e-3, 0.4e-3
    t1 = ratio * t2
    r = fea.solve([fea.Block(fea.ORGANIC_SUB, 0, L, 0, t1), fea.Block(fea.SILICON, 0, L, t1, t1 + t2)], -100.0,
                  nx=160, y_breaks_n=6, plane=plane)
    k_ref = fea.timoshenko_curvature(fea.ORGANIC_SUB, t1, fea.SILICON, t2, -100.0, plane)
    assert np.isclose(r.curvature(), k_ref, rtol=0.03)


def test_same_cte_no_warpage():
    m2 = fea.Material("Si2", 50e9, fea.SILICON.nu, fea.SILICON.cte)
    r = fea.solve([fea.Block(fea.SILICON, 0, 5e-3, 0, 5e-4), fea.Block(m2, 0, 5e-3, 5e-4, 1e-3)], -100.0, nx=40,
                  plane="stress")  # plane strain would (correctly) carry sigma_zz
    assert r.warpage() < 1e-12 and r.von_mises.max() < 1e-3


def test_warpage_scales_linearly_with_dT():
    blocks = [fea.Block(fea.ORGANIC_SUB, 0, 6e-3, 0, 5e-4), fea.Block(fea.SILICON, 0, 3e-3, 5e-4, 1e-3)]
    a = fea.solve(blocks, -50.0, nx=40).warpage()
    b = fea.solve(blocks, -100.0, nx=40).warpage()
    assert np.isclose(b, 2 * a, rtol=1e-9)


def test_void_cells_and_sliver_coordinates():
    eps = 1e-19
    blocks = [fea.Block(fea.ORGANIC_SUB, 0, 6e-3, 0, 5e-4), fea.Block(fea.UNDERFILL, 0, 3e-3, 5e-4 - eps, 5.2e-4),
              fea.Block(fea.SILICON, 0, 3e-3, 5.2e-4, 1e-3)]
    r = fea.solve(blocks, -100.0, nx=30)
    assert np.any(r.elem_material < 0)  # void above the substrate outside the die
    assert len(np.unique(np.round(r.y, 12))) == r.y.size
