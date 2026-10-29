# 2.5D vs 3D chiplet trade-off (all values modeled)

Two 8 x 8 mm dies: compute 90 W, cache/IO 15 W, 24 mm organic package, cold plate h = 25000 W/m2K with a 0.02 kg/s water stream (25 C inlet), cool-down dT = -100 K for warpage (plane strain, half cross-section). 2.5D: side by side on a 100 um Si interposer, UCIe-style advanced-package link at 16 GT/s across one 8 mm die edge. 3D: face-to-back stack, 50 um thinned bottom die, 4 mm2 of vertical D2D area at 4 Gb/s per lane (9 um hybrid bond or 36 um micro-bump pitch). In '3D-hybrid' the compute die is on top (next to the cold plate), in '3D-hybrid-cache-top' it is underneath. PDN: decaps optimised for a 50 % load step; IR drop is the DC drop from VRM to die.

Command: `python experiments/tradeoff_study.py`

## Main comparison

| config | footprint_mm2 | d2d_reach_mm | d2d_bw_tbps | d2d_energy_pj_bit | tj_compute_c | tj_cache_c | warpage_um | max_si_von_mises_mpa | pdn_zmax_mohm | decap_cost | ir_drop_mv |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 2.5D | 129.6 | 2.00 | 21.1 | 0.126 | 63.5 | 49.5 | 83.7 | 110 | 0.625 | 233.5 | 14.4 |
| 3D-microbump | 64.0 | 0.05 | 3.1 | 0.021 | 70.3 | 75.3 | 47.4 | 89 | 0.624 | 250.7 | 15.6 |
| 3D-hybrid | 64.0 | 0.05 | 49.4 | 0.020 | 70.5 | 71.0 | 48.4 | 86 | 0.624 | 250.7 | 15.6 |
| 3D-hybrid-cache-top | 64.0 | 0.05 | 49.4 | 0.020 | 106.2 | 74.1 | 48.4 | 86 | 0.624 | 216.9 | 12.0 |

## Peak compute-die temperature (C) vs compute power

| p_compute_w | tj_compute_2.5D | tj_compute_3D-microbump | tj_compute_3D-hybrid | tj_compute_3D-hybrid-cache-top |
|---|---|---|---|---|
| 40.0 | 42.8 | 48.8 | 49.0 | 64.6 |
| 70.0 | 55.2 | 61.7 | 61.9 | 89.5 |
| 100.0 | 67.6 | 74.6 | 74.8 | 114.5 |
| 130.0 | 80.0 | 87.5 | 87.7 | 139.5 |

Runtime 22 s.
