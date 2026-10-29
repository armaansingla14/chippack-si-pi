# S-parameter workflow on modeled coupled lines

All networks are analytic even/odd-mode line models written to and read back from Touchstone (.s4p). Max round-trip error: 7.0e-13. Eyes use Sdd21, 0.8 V differential swing, 5 mV rms noise, BER 1e-12, 0.3 UI rise time, no equalisation.

Command: `python experiments/sparam_study.py`

| network | sdd21_db_8g | scd21_db_8g | sdd21_db_16g | scd21_db_16g | sdd11_max_db | passive | max_singular_value | causal | noncausal_energy_ratio | reciprocity_err | eye_h_mv_16g | eye_w_ui_16g | ddj_ui_16g | eye_h_mv_32g | eye_w_ui_32g | ddj_ui_32g |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| pkg 15 mm | -1.09 | -300.00 | -1.50 | -300.00 | -22.18 | True | 0.9883 | True | 2.5e-05 | 0.0e+00 | 1300.7 | 0.938 | 0.010 | 1328.6 | 0.938 | 0.021 |
| pcb 100 mm | -5.09 | -300.00 | -7.20 | -300.00 | -28.85 | True | 0.9417 | True | 1.1e-05 | 1.0e-15 | 665.9 | 0.844 | 0.061 | 452.4 | 0.719 | 0.122 |
| pkg+pcb+pkg | -7.24 | -300.00 | -10.21 | -300.00 | -21.49 | True | 0.9190 | True | 3.1e-05 | 5.6e-16 | 409.7 | 0.750 | 0.106 | 127.6 | 0.469 | 0.218 |
| pkg+pcb(3 ps skew)+pkg | -7.59 | -33.40 | -10.02 | -43.19 | -21.55 | True | 0.9189 | True | 3.2e-05 | 6.8e-16 | 367.0 | 0.750 | 0.112 | 115.0 | 0.438 | 0.245 |
