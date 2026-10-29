"""Accuracy of the thermal FV and FEA solvers against closed-form references.

    python experiments/solver_validation.py
"""
import numpy as np

from common import md_table, write_csv, write_md

from chippack import fea, thermal


def main():
    rows = []
    L = 0.01
    for nz in (1, 2, 4):
        layers = [thermal.Layer("si", 500e-6, 150, nz), thermal.Layer("tim", 50e-6, 5, nz),
                  thermal.Layer("cu", 1e-3, 390, nz)]
        m = thermal.ThermalModel(L, L, 8, 8, layers, h_top=1e4, h_bottom=1e-9, power={0: np.full((8, 8), 20 / 64)})
        r = m.solve()
        A = L * L
        top = 25 + 20 / (1e4 * A) + 20 * (1e-3 / 390 + 50e-6 / 5) / A
        exact_max = top + 20 * 500e-6 / (2 * 150 * A)  # analytic peak of uniformly heated Si slab
        rows.append({"check": f"1-D 3-layer stack, {nz} cells/layer", "model": r.t_max, "reference": exact_max,
                     "rel_err_pct": 100 * abs(r.t_max - exact_max) / (exact_max - 25),
                     "energy_balance_err": r.energy_balance_error})
    for nx, nyb in ((40, 2), (80, 4), (160, 6)):
        t1 = t2 = 0.4e-3
        r = fea.solve([fea.Block(fea.ORGANIC_SUB, 0, 10e-3, 0, t1), fea.Block(fea.SILICON, 0, 10e-3, t1, t1 + t2)],
                      -100.0, nx=nx, y_breaks_n=nyb, plane="stress")
        ref = fea.timoshenko_curvature(fea.ORGANIC_SUB, t1, fea.SILICON, t2, -100.0, "stress")
        rows.append({"check": f"bimetal curvature (1/m), {nx}x{2 * nyb} Q4", "model": r.curvature(), "reference": ref,
                     "rel_err_pct": 100 * abs(r.curvature() - ref) / abs(ref), "energy_balance_err": None})
    for r in rows:
        print(r)
    write_csv("solver_validation.csv", rows)
    write_md("solver_validation.md",
             "# Solver validation against closed-form references\n\n"
             "Thermal: relative error of the peak temperature rise of a uniformly heated 3-layer stack "
             "(the FV peak is at a cell centre, so it converges to the slab-centre value). FEA: curvature of a "
             "10 mm organic/Si bimetal strip after a -100 K change vs Timoshenko's formula (plane stress).\n\n"
             "Command: `python experiments/solver_validation.py`\n\n"
             + md_table(rows, None, {"model": ".5g", "reference": ".5g", "rel_err_pct": ".3f",
                                     "energy_balance_err": ".1e"}))


if __name__ == "__main__":
    main()
