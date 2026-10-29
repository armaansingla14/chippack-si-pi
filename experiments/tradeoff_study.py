"""2.5D vs 3D chiplet integration trade-off (modeled).

    python experiments/tradeoff_study.py
"""
import time

from common import md_table, write_csv, write_md

from chippack.chiplet import STYLES, SystemSpec, evaluate


def main():
    t0 = time.time()
    spec = SystemSpec()
    rows = [evaluate(spec, s) for s in STYLES]
    for r in rows:
        print(r)
    sweep = []
    for p in (40.0, 70.0, 100.0, 130.0):
        sp = SystemSpec(p_compute_w=p)
        row = {"p_compute_w": p}
        for s in STYLES:
            row[f"tj_compute_{s}"] = evaluate(sp, s)["tj_compute_c"]
        sweep.append(row)
        print(row)
    write_csv("tradeoff_main.csv", rows)
    write_csv("tradeoff_power_sweep.csv", sweep)
    cols = ["config", "footprint_mm2", "d2d_reach_mm", "d2d_bw_tbps", "d2d_energy_pj_bit", "tj_compute_c",
            "tj_cache_c", "warpage_um", "max_si_von_mises_mpa", "pdn_zmax_mohm", "decap_cost", "ir_drop_mv"]
    fmt = {"footprint_mm2": ".1f", "d2d_reach_mm": ".2f", "d2d_bw_tbps": ".1f", "d2d_energy_pj_bit": ".3f",
           "tj_compute_c": ".1f", "tj_cache_c": ".1f", "warpage_um": ".1f", "max_si_von_mises_mpa": ".0f",
           "pdn_zmax_mohm": ".3f", "decap_cost": ".1f", "ir_drop_mv": ".1f"}
    s = spec
    text = (
        "# 2.5D vs 3D chiplet trade-off (all values modeled)\n\n"
        f"Two {s.die_mm:.0f} x {s.die_mm:.0f} mm dies: compute {s.p_compute_w:.0f} W, cache/IO {s.p_cache_w:.0f} W, "
        f"{s.package_mm:.0f} mm organic package, cold plate h = {s.h_coldplate:.0f} W/m2K with a 0.02 kg/s water "
        f"stream (25 C inlet), cool-down dT = {s.dT_cooldown:.0f} K for warpage (plane strain, half cross-section). "
        "2.5D: side by side on a 100 um Si interposer, UCIe-style advanced-package link at "
        f"{s.d2d_rate_2p5d / 1e9:.0f} GT/s across one {s.die_mm:.0f} mm die edge. 3D: face-to-back stack, "
        f"{s.thinned_die_t_um:.0f} um thinned bottom die, {s.d2d_signal_area_mm2:.0f} mm2 of vertical D2D area at "
        f"{s.rate_3d / 1e9:.0f} Gb/s per lane ({s.hybrid_bond_pitch_um:.0f} um hybrid bond or "
        f"{s.microbump_pitch_um:.0f} um micro-bump pitch). In '3D-hybrid' the compute die is on top (next to "
        "the cold plate), in '3D-hybrid-cache-top' it is underneath. PDN: decaps optimised for a 50 % "
        "load step; IR drop is the DC drop from VRM to die.\n\n"
        "Command: `python experiments/tradeoff_study.py`\n\n"
        "## Main comparison\n\n" + md_table(rows, cols, fmt)
        + "\n\n## Peak compute-die temperature (C) vs compute power\n\n"
        + md_table(sweep, None, {k: ".1f" for k in sweep[0]})
        + f"\n\nRuntime {time.time() - t0:.0f} s.\n")
    write_md("tradeoff.md", text)


if __name__ == "__main__":
    main()
