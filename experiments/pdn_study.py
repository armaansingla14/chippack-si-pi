"""PDN decap optimisation vs simple baselines over several target impedances (modeled).

    python experiments/pdn_study.py
"""
import numpy as np

from common import md_table, write_csv, write_md

from chippack.pdn import LIBRARY, PDN, antiresonances, baseline_equal_mix, baseline_uniform, optimise_decaps, \
    target_impedance


def main():
    f = np.logspace(3, np.log10(2e8), 300)
    base = PDN()
    z0 = np.abs(base.impedance(f))
    peaks = antiresonances(f, base.impedance(f))
    rows = []
    for i_step in (2.0, 3.0, 5.0, 7.0):
        zt = target_impedance(0.75, 0.05, i_step)
        opt = optimise_decaps(base, f, zt)
        mix = baseline_equal_mix(base, f, zt)
        best_single = min((baseline_uniform(base, f, zt, c, max_n=200) for c in LIBRARY),
                          key=lambda r: (not r.meets_target, r.cost))
        rows.append({"i_step_a": i_step, "z_target_mohm": zt * 1e3,
                     "opt_cost": opt.cost, "opt_parts": sum(opt.counts.values()), "opt_meets": opt.meets_target,
                     "opt_zmax_mohm": opt.max_z_ohm * 1e3, "opt_bom": ", ".join(f"{k} x{v}" for k, v in opt.counts.items()),
                     "equal_mix_cost": mix.cost, "equal_mix_meets": mix.meets_target,
                     "best_single_type": next(iter(best_single.counts)), "best_single_meets": best_single.meets_target,
                     "best_single_cost": best_single.cost})
        print(rows[-1])
    write_csv("pdn_study.csv", rows)
    write_md("pdn_study.md",
             "# PDN target impedance and decap optimisation (modeled)\n\n"
             f"Bare PDN (no decaps): max |Z| = {z0.max() * 1e3:.1f} mOhm at {f[np.argmax(z0)] / 1e6:.2f} MHz; "
             f"anti-resonances at {', '.join(f'{p[0] / 1e6:.2f} MHz ({p[1] * 1e3:.1f} mOhm)' for p in peaks)}.\n"
             "Band 1 kHz - 200 MHz, Vdd 0.75 V, 5 % ripple. Costs are relative units from the assumed part "
             "library. Baselines: n of every part type (equal mix), and the cheapest single part type "
             "(up to 200 pieces).\n\n"
             "Command: `python experiments/pdn_study.py`\n\n"
             + md_table(rows, None, {"z_target_mohm": ".2f", "opt_cost": ".1f", "opt_zmax_mohm": ".2f",
                                     "equal_mix_cost": ".1f", "best_single_cost": ".1f", "i_step_a": ".0f"}))


if __name__ == "__main__":
    main()
