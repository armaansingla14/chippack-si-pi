# PDN target impedance and decap optimisation (modeled)

Bare PDN (no decaps): max |Z| = 2285.6 mOhm at 1.62 MHz; anti-resonances at 1.62 MHz (2285.6 mOhm).
Band 1 kHz - 200 MHz, Vdd 0.75 V, 5 % ripple. Costs are relative units from the assumed part library. Baselines: n of every part type (equal mix), and the cheapest single part type (up to 200 pieces).

Command: `python experiments/pdn_study.py`

| i_step_a | z_target_mohm | opt_cost | opt_parts | opt_meets | opt_zmax_mohm | opt_bom | equal_mix_cost | equal_mix_meets | best_single_type | best_single_meets | best_single_cost |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 2 | 18.75 | 10.7 | 4 | True | 16.70 | 0402_1uF x1, 0402_4u7F x1, 0805_22uF x1, bulk_330uF x1 | 15.7 | True | 0402_4u7F | True | 235.5 |
| 3 | 12.50 | 13.9 | 7 | True | 11.96 | 0201_100nF x2, 0402_1uF x2, 0402_4u7F x1, 0805_22uF x1, bulk_330uF x1 | 15.7 | True | 0201_100nF | False | 200.0 |
| 5 | 7.50 | 26.9 | 15 | True | 7.48 | 0201_100nF x9, 0402_1uF x2, 0402_4u7F x1, 0805_22uF x1, bulk_330uF x2 | 31.4 | True | 0201_100nF | False | 200.0 |
| 7 | 5.36 | 44.5 | 26 | True | 5.35 | 0201_100nF x15, 0402_1uF x5, 0402_4u7F x1, 0805_22uF x2, bulk_330uF x3 | 1570.0 | False | 0201_100nF | False | 200.0 |
