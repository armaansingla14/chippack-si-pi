"""UCIe-style standard vs advanced package die-to-die link vs data rate (modeled).

    python experiments/ucie_study.py
"""
from common import md_table, write_csv, write_md

from chippack.ucie import ADVANCED, STANDARD, simulate


def main():
    rows = []
    for pkg in (STANDARD, ADVANCED):
        for rate in (4e9, 8e9, 12e9, 16e9, 24e9, 32e9):
            r = simulate(pkg, rate)
            rows.append(r)
            print(r)
    write_csv("ucie_study.csv", rows)
    write_md("ucie_study.md",
             "# UCIe-style die-to-die link (modeled)\n\n"
             "Standard: 10 mm organic-substrate trace, 50 ohm terminated. Advanced: 2 mm RC-dominated "
             "silicon bridge/interposer wire, unterminated. Eye from the pulse response at BER 1e-15 with "
             "5 mV rms noise and 0.4 V swing (no equalisation). Energy/bit = termination + CV^2 activity "
             "+ an assumed fixed overhead. Bandwidth density = lanes per module x rate / assumed module "
             "shoreline. Package parameters are model assumptions, not spec values.\n\n"
             "Command: `python experiments/ucie_study.py`\n\n"
             + md_table(rows, None, {"data_rate_gtps": ".0f", "reach_mm": ".0f", "loss_at_nyq_db": ".2f",
                                     "eye_height_mv": ".1f", "eye_width_ui": ".2f", "pd_eye_height_mv": ".1f",
                                     "energy_pj_per_bit": ".3f", "bw_density_gbps_per_mm": ".0f"}))


if __name__ == "__main__":
    main()
