"""S-parameter workflow on modeled coupled-line channels: Touchstone write/read, cascading,
mixed-mode conversion, passivity and causality checks, and the resulting eye.

    python experiments/sparam_study.py
"""
from __future__ import annotations

import os

import numpy as np

from common import ROOT, md_table, write_csv, write_md

from chippack import touchstone as ts
from chippack.network import (Network, cascade, causality_check, coupled_line_4port, mixed_mode,
                              passivity_violation, reciprocity_error)
from chippack.serdes import analyze_eye, ddj_from_pulse, pulse_response


def main():
    f = np.arange(1, 401) * 100e6  # 100 MHz .. 40 GHz, uniform (needed for causality/time domain)
    data = os.path.join(ROOT, "data")
    os.makedirs(data, exist_ok=True)
    pkg = coupled_line_4port(f, 0.015, z_even=58, z_odd=46, loss_db_per_m_per_sqrtghz=25.0)
    pcb = coupled_line_4port(f, 0.10, z_even=55, z_odd=48, loss_db_per_m_per_sqrtghz=18.0)
    ts.write(os.path.join(data, "pkg_pair_15mm.s4p"), pkg, "RI", "GHZ", comment="modeled coupled pair, 15 mm")
    ts.write(os.path.join(data, "pcb_pair_100mm.s4p"), pcb, "RI", "GHZ", comment="modeled coupled pair, 100 mm")
    a = ts.read(os.path.join(data, "pkg_pair_15mm.s4p"))
    b = ts.read(os.path.join(data, "pcb_pair_100mm.s4p"))
    rt_err = float(max(np.abs(a.s - pkg.s).max(), np.abs(b.s - pcb.s).max()))
    # cascade pkg + pcb + pkg (reorder to left/right port grouping first)
    order = [1, 3, 2, 4]
    chan = cascade(a.reorder(order), b.reorder(order), a.reorder(order)).reorder(order)
    # skewed variant: 3 ps extra delay on one leg of the pcb section
    sk = pcb.s.copy()
    d = np.exp(-2j * np.pi * f * 3e-12)
    sk[:, 3, 2] *= d
    sk[:, 2, 3] *= d
    skewed = cascade(a.reorder(order), Network(f, sk).reorder(order), a.reorder(order)).reorder(order)
    rows = []
    for name, net in (("pkg 15 mm", a), ("pcb 100 mm", b), ("pkg+pcb+pkg", chan), ("pkg+pcb(3 ps skew)+pkg", skewed)):
        mm = mixed_mode(net)
        sdd21 = mm["dd"][:, 1, 0]
        scd21 = mm["cd"][:, 1, 0]
        sdd11 = mm["dd"][:, 0, 0]
        pas, worst_sv = passivity_violation(net)
        caus, nc = causality_check(sdd21, f)
        row = {"network": name}
        for fq in (8e9, 16e9):
            k = int(np.argmin(np.abs(f - fq)))
            row[f"sdd21_db_{int(fq / 1e9)}g"] = float(20 * np.log10(abs(sdd21[k])))
            row[f"scd21_db_{int(fq / 1e9)}g"] = float(20 * np.log10(max(abs(scd21[k]), 1e-15)))
        row["sdd11_max_db"] = float(20 * np.log10(np.abs(sdd11).max()))
        row["passive"] = pas
        row["max_singular_value"] = worst_sv
        row["causal"] = caus
        row["noncausal_energy_ratio"] = nc
        row["reciprocity_err"] = reciprocity_error(net)
        for rate in (16e9, 32e9):
            _, p = pulse_response(f, sdd21, rate, spu=32, rise_ui=0.3)
            e = analyze_eye(p, 32, noise_rms=0.005, ber=1e-12, swing=0.8)
            row[f"eye_h_mv_{int(rate / 1e9)}g"] = 1e3 * e.stat_height
            row[f"eye_w_ui_{int(rate / 1e9)}g"] = e.stat_width_ui
            row[f"ddj_ui_{int(rate / 1e9)}g"] = ddj_from_pulse(p, 32, 7)
        rows.append(row)
        print(row)
    write_csv("sparam_study.csv", rows)
    fmt = {k: ".2f" for k in rows[0] if "db" in k} | {"max_singular_value": ".4f", "noncausal_energy_ratio": ".1e",
                                                      "reciprocity_err": ".1e"}
    fmt |= {k: ".1f" for k in rows[0] if k.startswith("eye_h")} | {k: ".3f" for k in rows[0] if "_ui_" in k}
    write_md("sparam_study.md",
             "# S-parameter workflow on modeled coupled lines\n\n"
             "All networks are analytic even/odd-mode line models written to and read back from Touchstone "
             f"(.s4p). Max round-trip error: {rt_err:.1e}. Eyes use Sdd21, 0.8 V differential swing, "
             "5 mV rms noise, BER 1e-12, 0.3 UI rise time, no equalisation.\n\n"
             "Command: `python experiments/sparam_study.py`\n\n" + md_table(rows, None, fmt))


if __name__ == "__main__":
    main()
