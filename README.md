# chippack-si-pi

Models for advanced packaging: signal integrity for SerDes and UCIe-style die-to-die links,
S-parameter handling, power-delivery impedance, a finite-volume thermal solver, a 2-D
thermo-mechanical FEA solver for warpage and stress, and a 2.5D vs 3D chiplet trade-off study
that ties them together.

**Everything here is modeled.** The channels, package stack-ups, capacitor library and
material properties are analytic or assumed inputs. They are not measured or extracted from
a real design. The numbers below show what these models produce; they are not silicon or
lab data.

## Features

- **S-parameters** (`network.py`, `touchstone.py`)
  - Touchstone 1.x read/write for `.s1p` through `.s4p` and up (RI/MA/DB, Hz-GHz, reference
    impedance, comments, the 2-port column-order rule).
  - S <-> Z/Y, 2-port S <-> ABCD, 2N-port S <-> T, renormalisation, port reordering and
    interpolation.
  - Cascading of 2-ports and 4-ports.
  - Single-ended to mixed-mode conversion (Sdd/Sdc/Scd/Scc).
  - Passivity check (max singular value) and simple enforcement, reciprocity error, and a
    causality check that measures non-causal impulse energy.
  - Analytic RLGC line and even/odd-mode coupled-pair generators.
- **SerDes / eye analysis** (`serdes.py`): pulse response from S21 (DC extrapolation,
  windowing), peak-distortion eye, statistical eye at a target BER (exact ISI PDF + Gaussian
  noise), eye width at that BER, and data-dependent jitter from edge-pattern crossings.
- **UCIe-style D2D** (`ucie.py`): standard (organic, terminated) and advanced (Si
  bridge/interposer, unterminated, RC-dominated) package classes. Each includes driver
  impedance, pad capacitance and termination, and reports eye, loss, energy/bit and
  shoreline bandwidth density.
- **PDN** (`pdn.py`): VRM/board/package/die ladder impedance, target impedance, an
  anti-resonance finder, and a greedy add-then-prune decap optimiser compared against two
  baselines (equal mix of every part, and a single part type).
- **Thermal** (`thermal.py`): 3-D steady finite-volume conduction over stacked layers, with
  per-cell conductivity maps, power maps, a convective bottom, and a top surface that is
  either a fixed ambient or a 1-D advecting coolant stream (the "CFD-lite" part).
- **FEA** (`fea.py`): 2-D bilinear Q4 elements, plane stress or plane strain, rectangular
  material blocks with void regions, and a uniform temperature-change load. Outputs are
  warpage, curvature, stresses and von Mises. Checked against Timoshenko's bimetal formula.
- **Chiplet trade-off** (`chiplet.py`): 2.5D (side-by-side on an interposer) vs 3D (micro-bump
  or hybrid-bond stack, compute die on top or bottom). Compares D2D bandwidth/energy,
  junction temperatures, warpage/stress, PDN and IR drop.

## Layout

```
chippack/     network, touchstone, serdes, ucie, pdn, thermal, fea, chiplet
experiments/  sparam_study, ucie_study, pdn_study, solver_validation, tradeoff_study
results/      Markdown/CSV tables written by the experiments
data/         modeled .s4p files written by sparam_study.py
tests/        pytest suite (103 tests)
```

## Running

```
pip install -r requirements.txt
python -m pytest -q
cd experiments
python tradeoff_study.py      # ~20 s
python sparam_study.py
python ucie_study.py
python pdn_study.py
python solver_validation.py
```

## Results (modeled)

### 2.5D vs 3D trade-off

`python experiments/tradeoff_study.py`: two 8 x 8 mm dies (90 W compute, 15 W cache/IO), a
24 mm organic package, a water cold plate, and a -100 K cool-down for warpage. Every input is
listed in `SystemSpec`; full text is in `results/tradeoff.md`.

| config | footprint (mm2) | D2D BW (Tb/s) | D2D energy (pJ/bit) | Tj compute (C) | Tj cache (C) | warpage (um) | max Si von Mises (MPa) | PDN Zmax (mOhm) | decap cost | IR drop (mV) |
|---|---|---|---|---|---|---|---|---|---|---|
| 2.5D | 129.6 | 21.1 | 0.126 | 63.5 | 49.5 | 83.7 | 110 | 0.625 | 233.5 | 14.4 |
| 3D-microbump | 64.0 | 3.1 | 0.021 | 70.3 | 75.3 | 47.4 | 89 | 0.624 | 250.7 | 15.6 |
| 3D-hybrid | 64.0 | 49.4 | 0.020 | 70.5 | 71.0 | 48.4 | 86 | 0.624 | 250.7 | 15.6 |
| 3D-hybrid-cache-top | 64.0 | 49.4 | 0.020 | 106.2 | 74.1 | 48.4 | 86 | 0.624 | 216.9 | 12.0 |

Peak compute-die temperature (C) vs compute power:

| compute W | 2.5D | 3D-microbump | 3D-hybrid | 3D-hybrid-cache-top |
|---|---|---|---|---|
| 40 | 42.8 | 48.8 | 49.0 | 64.6 |
| 70 | 55.2 | 61.7 | 61.9 | 89.5 |
| 100 | 67.6 | 74.6 | 74.8 | 114.5 |
| 130 | 80.0 | 87.5 | 87.7 | 139.5 |

With these assumptions:

- 3D halves the footprint and cuts D2D energy per bit by about 6x compared with the 2.5D UCIe-style link.
- Hybrid bonding gives about 16x the micro-bump bandwidth for the same D2D area.
- 3D runs the compute die 7 C hotter with compute on top, and 43 C hotter with the cache
  stacked on top of it.
- 3D gives lower package warpage, because the interposer is gone and the silicon is smaller.
- The PDN numbers differ only slightly between styles. The power-TSV and interposer path
  values are assumptions.

### UCIe-style D2D link

`python experiments/ucie_study.py` (BER 1e-15, 5 mV noise, 0.4 V swing, no EQ):

| package | GT/s | eye height (mV) | eye width (UI) | energy (pJ/bit) | BW density (Gb/s/mm) |
|---|---|---|---|---|---|
| standard, 10 mm | 16 | 380.5 | 0.84 | 0.285 | 233 |
| standard, 10 mm | 32 | 259.6 | 0.72 | 0.253 | 465 |
| advanced, 2 mm | 16 | 334.7 | 0.81 | 0.126 | 2639 |
| advanced, 2 mm | 32 | 94.2 | 0.47 | 0.126 | 5278 |

The RC-limited advanced-package wire closes faster with data rate. The terminated standard
link keeps a larger eye but costs more energy per bit and has about 11x lower shoreline density.

### S-parameter workflow

`python experiments/sparam_study.py` cascades a 15 mm package pair, a 100 mm PCB pair and a
second 15 mm package pair, all modeled and round-tripped through `.s4p` files (max error 7e-13).
Every network passed the passivity, causality and reciprocity checks. A 3 ps intra-pair skew in
the PCB section raised Scd21 at 8 GHz from numerically zero to -33.4 dB. It also reduced the
32 Gb/s Sdd21 eye from 127.6 to 115.0 mV (0.8 V differential swing, no EQ).

### PDN decap optimisation

`python experiments/pdn_study.py` (0.75 V, 5 % ripple, 1 kHz-200 MHz):

| load step (A) | Z target (mOhm) | optimiser cost (parts) | equal-mix baseline cost | single-type baseline |
|---|---|---|---|---|
| 2 | 18.75 | 10.7 (4) | 15.7 | 235.5 (0402 4.7 uF) |
| 3 | 12.50 | 13.9 (7) | 15.7 | fails within 200 parts |
| 5 | 7.50 | 26.9 (15) | 31.4 | fails |
| 7 | 5.36 | 44.5 (26) | fails (100 of each) | fails |

Costs are relative units from the assumed part library.

### Solver validation

`python experiments/solver_validation.py`:

- **Thermal**: the FV solver reproduces the analytic peak temperature of a uniformly heated
  3-layer stack to round-off. This is expected, because cell-centred FV is exact for that 1-D
  quadratic profile. Energy balance closes to about 1e-14.
- **FEA**: bimetal curvature error vs Timoshenko is 6.2 % / 1.7 % / 0.48 % on 40x4 / 80x8 /
  160x12 meshes.

## Limitations

- All geometry, materials, package and PHY parameters are representative assumptions. The
  UCIe-style numbers (bump pitch, lanes per module, module width) are loosely based on public
  descriptions, not taken from the specification.
- The thermal model is steady-state conduction with lumped convection and a 1-D coolant energy
  balance. It does not solve the flow field, and it has no temperature-dependent properties or
  spreading-resistance calibration.
- The FEA is 2-D, linear-elastic and small-deflection, with uniform temperature change and no
  viscoelastic underfill/mold or solder creep. Warpage numbers are cross-section estimates.
- The PDN is a lumped ladder with no spatial plane model. The decap optimiser is greedy and does
  not guarantee the minimum cost.
- The eye analysis is NRZ, linear and time-invariant, with no equalisation or crosstalk
  aggressors. The DDJ estimate ignores bits after the edge.

## License

MIT, see `LICENSE`.
