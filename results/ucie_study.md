# UCIe-style die-to-die link (modeled)

Standard: 10 mm organic-substrate trace, 50 ohm terminated. Advanced: 2 mm RC-dominated silicon bridge/interposer wire, unterminated. Eye from the pulse response at BER 1e-15 with 5 mV rms noise and 0.4 V swing (no equalisation). Energy/bit = termination + CV^2 activity + an assumed fixed overhead. Bandwidth density = lanes per module x rate / assumed module shoreline. Package parameters are model assumptions, not spec values.

Command: `python experiments/ucie_study.py`

| package | data_rate_gtps | reach_mm | loss_at_nyq_db | eye_height_mv | eye_width_ui | pd_eye_height_mv | energy_pj_per_bit | bw_density_gbps_per_mm |
|---|---|---|---|---|---|---|---|---|
| standard (organic substrate) | 4 | 10 | -0.01 | 422.9 | 0.91 | 493.4 | 0.472 | 58 |
| standard (organic substrate) | 8 | 10 | 0.39 | 405.6 | 0.91 | 476.7 | 0.347 | 116 |
| standard (organic substrate) | 12 | 10 | 1.20 | 355.6 | 0.94 | 426.6 | 0.305 | 175 |
| standard (organic substrate) | 16 | 10 | 0.91 | 380.5 | 0.84 | 449.4 | 0.285 | 233 |
| standard (organic substrate) | 24 | 10 | 2.38 | 328.2 | 0.78 | 392.9 | 0.264 | 349 |
| standard (organic substrate) | 32 | 10 | 1.98 | 259.6 | 0.72 | 318.9 | 0.253 | 465 |
| advanced (Si bridge/interposer) | 4 | 2 | 0.79 | 669.6 | 0.94 | 748.5 | 0.126 | 660 |
| advanced (Si bridge/interposer) | 8 | 2 | 2.41 | 557.0 | 0.91 | 636.0 | 0.126 | 1320 |
| advanced (Si bridge/interposer) | 12 | 2 | 3.97 | 441.1 | 0.88 | 517.0 | 0.126 | 1979 |
| advanced (Si bridge/interposer) | 16 | 2 | 5.16 | 334.7 | 0.81 | 416.8 | 0.126 | 2639 |
| advanced (Si bridge/interposer) | 24 | 2 | 6.55 | 201.9 | 0.62 | 279.8 | 0.126 | 3959 |
| advanced (Si bridge/interposer) | 32 | 2 | 7.34 | 94.2 | 0.47 | 160.7 | 0.126 | 5278 |
