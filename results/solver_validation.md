# Solver validation against closed-form references

Thermal: relative error of the peak temperature rise of a uniformly heated 3-layer stack (the FV peak is at a cell centre, so it converges to the slab-centre value). FEA: curvature of a 10 mm organic/Si bimetal strip after a -100 K change vs Timoshenko's formula (plane stress).

Command: `python experiments/solver_validation.py`

| check | model | reference | rel_err_pct | energy_balance_err |
|---|---|---|---|---|
| 1-D 3-layer stack, 1 cells/layer | 47.846 | 47.846 | 0.000 | 3.1e-14 |
| 1-D 3-layer stack, 2 cells/layer | 47.846 | 47.846 | 0.000 | 1.8e-14 |
| 1-D 3-layer stack, 4 cells/layer | 47.846 | 47.846 | 0.000 | 5.3e-14 |
| bimetal curvature (1/m), 40x4 Q4 | -1.7985 | -1.9183 | 6.244 | - |
| bimetal curvature (1/m), 80x8 Q4 | -1.8865 | -1.9183 | 1.657 | - |
| bimetal curvature (1/m), 160x12 Q4 | -1.9091 | -1.9183 | 0.479 | - |
