"""Power-delivery-network impedance model, target impedance and decap optimisation.

Lumped model seen from the die bumps (all elements in parallel to ground unless
noted):  VRM (R + L, with the PCB/package spreading path in series) || board
bulk caps || package caps || on-die capacitance, separated by package R/L.

    die --[R_pkg + L_pkg]-- package node --[R_bga + L_bga]-- board node -- VRM
           |                        |                          |
         C_die + R_die           pkg decaps                board decaps
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass(frozen=True)
class Capacitor:
    name: str
    c: float  # F
    esr: float  # ohm
    esl: float  # H
    cost: float = 1.0
    location: str = "pkg"  # pkg | board

    def z(self, f: np.ndarray) -> np.ndarray:
        w = 2 * np.pi * f
        return self.esr + 1j * w * self.esl + 1 / (1j * w * self.c)

    @property
    def srf_hz(self) -> float:
        return 1 / (2 * np.pi * np.sqrt(self.esl * self.c))


# Representative MLCC-like parts (values are typical-order assumptions, not a vendor library)
LIBRARY = [
    Capacitor("0201_100nF", 100e-9, 0.030, 120e-12, 1.0, "pkg"),
    Capacitor("0402_1uF", 1e-6, 0.012, 200e-12, 1.2, "pkg"),
    Capacitor("0402_4u7F", 4.7e-6, 0.008, 250e-12, 1.5, "pkg"),
    Capacitor("LICC_2uF", 2.2e-6, 0.010, 50e-12, 4.0, "pkg"),
    Capacitor("0805_22uF", 22e-6, 0.004, 500e-12, 2.0, "board"),
    Capacitor("bulk_330uF", 330e-6, 0.008, 1.5e-9, 6.0, "board"),
]


@dataclass
class PDN:
    r_vrm: float = 0.5e-3
    l_vrm: float = 10e-9  # effective VRM output inductance (loop bandwidth folded in)
    r_bga: float = 0.2e-3
    l_bga: float = 40e-12
    r_pkg: float = 0.5e-3
    l_pkg: float = 10e-12
    c_die: float = 1e-6
    r_die: float = 1.0e-3
    caps: dict[Capacitor, int] = field(default_factory=dict)

    def impedance(self, f: np.ndarray) -> np.ndarray:
        w = 2 * np.pi * f
        z_vrm = self.r_vrm + 1j * w * self.l_vrm
        y_board = 1 / z_vrm
        y_pkg = np.zeros_like(f, dtype=complex)
        for cap, n in self.caps.items():
            if n <= 0:
                continue
            y = n / cap.z(f)
            if cap.location == "board":
                y_board = y_board + y
            else:
                y_pkg = y_pkg + y
        z_board = 1 / y_board
        z_pkg_node = 1 / (y_pkg + 1 / (z_board + self.r_bga + 1j * w * self.l_bga))
        z_die_cap = self.r_die + 1 / (1j * w * self.c_die)
        return 1 / (1 / (z_pkg_node + self.r_pkg + 1j * w * self.l_pkg) + 1 / z_die_cap)


def target_impedance(vdd: float, ripple_frac: float, i_transient: float) -> float:
    """Z_target = Vdd * allowed ripple / transient current step."""
    return vdd * ripple_frac / i_transient


def worst_violation(pdn: PDN, f: np.ndarray, z_target: float) -> float:
    return float(np.max(np.abs(pdn.impedance(f)) / z_target))


@dataclass
class DecapResult:
    counts: dict[str, int]
    cost: float
    max_z_ohm: float
    meets_target: bool
    iterations: int


def optimise_decaps(base: PDN, f: np.ndarray, z_target: float, library=LIBRARY, max_per_type: int = 64,
                    max_iter: int = 400) -> DecapResult:
    """Greedy add (best penalty reduction per unit cost) then prune redundant parts.

    The greedy step minimises a smooth penalty, the mean squared excess of
    |Z|/Z_target above 1 over the band (a pure worst-case objective stalls because
    adding one part often just moves an anti-resonance). Parts are added until the
    worst case meets the target or the budget runs out, then removed one at a time
    (most expensive first) while the target still holds.
    """
    counts = {c: 0 for c in library}
    it = 0

    def evaluate(cn):
        p = PDN(**{**base.__dict__, "caps": cn})
        r = np.abs(p.impedance(f)) / z_target
        return float(np.mean(np.maximum(r - 1, 0) ** 2)), float(r.max())

    pen, cur = evaluate(counts)
    while cur > 1.0 and it < max_iter:
        it += 1
        best = None
        for cap in library:
            if counts[cap] >= max_per_type:
                continue
            trial = dict(counts)
            trial[cap] += 1
            tp, tw = evaluate(trial)
            gain = (pen - tp) / cap.cost
            if best is None or gain > best[0]:
                best = (gain, cap, tp, tw)
        if best is None or best[0] <= 0:
            break
        counts[best[1]] += 1
        pen, cur = best[2], best[3]

    def score(cn):
        return evaluate(cn)[1]
    # prune: try removing parts from most to least expensive
    changed = True
    while changed and cur <= 1.0:
        changed = False
        for cap in sorted(library, key=lambda c: -c.cost):
            if counts[cap] == 0:
                continue
            trial = dict(counts)
            trial[cap] -= 1
            s = score(trial)
            if s <= 1.0:
                counts, cur, changed = trial, s, True
    cost = sum(c.cost * n for c, n in counts.items())
    return DecapResult({c.name: n for c, n in counts.items() if n}, cost, cur * z_target, cur <= 1.0, it)


def baseline_uniform(base: PDN, f: np.ndarray, z_target: float, cap: Capacitor, max_n: int = 400) -> DecapResult:
    """Baseline: use only one capacitor type and increase the count until the target is met."""
    for n in range(0, max_n + 1):
        p = PDN(**{**base.__dict__, "caps": {cap: n}})
        v = worst_violation(p, f, z_target)
        if v <= 1.0:
            return DecapResult({cap.name: n}, n * cap.cost, v * z_target, True, n)
    return DecapResult({cap.name: max_n}, max_n * cap.cost, v * z_target, False, max_n)


def antiresonances(f: np.ndarray, z: np.ndarray) -> list[tuple[float, float]]:
    m = np.abs(z)
    i = np.nonzero((m[1:-1] > m[:-2]) & (m[1:-1] > m[2:]))[0] + 1
    return [(float(f[k]), float(m[k])) for k in i]


def baseline_equal_mix(base: PDN, f: np.ndarray, z_target: float, library=LIBRARY, max_n: int = 100) -> DecapResult:
    """Baseline: n of every capacitor type, increasing n until the target is met."""
    v = np.inf
    for n in range(0, max_n + 1):
        p = PDN(**{**base.__dict__, "caps": {c: n for c in library}})
        v = worst_violation(p, f, z_target)
        if v <= 1.0:
            break
    return DecapResult({c.name: n for c in library}, n * sum(c.cost for c in library), v * z_target, v <= 1.0, n)
