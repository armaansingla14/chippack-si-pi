"""Touchstone v1.x reader/writer (.s1p .. .s4p and beyond).

Supports frequency units Hz/kHz/MHz/GHz, formats RI/MA/DB, the ``R`` reference
impedance, comment lines and the 2-port column-order quirk (S11 S21 S12 S22).
"""
from __future__ import annotations

import os
import re

import numpy as np

from .network import Network

_UNITS = {"HZ": 1.0, "KHZ": 1e3, "MHZ": 1e6, "GHZ": 1e9}


def _ports_from_name(path: str) -> int:
    m = re.search(r"\.s(\d+)p$", path.lower())
    if not m:
        raise ValueError(f"cannot infer port count from {path!r}")
    return int(m.group(1))


def read(path: str, n_ports: int | None = None) -> Network:
    n = n_ports or _ports_from_name(path)
    unit, fmt, z0 = 1e9, "MA", 50.0
    values: list[float] = []
    with open(path) as fh:
        for line in fh:
            line = line.split("!", 1)[0].strip()
            if not line:
                continue
            if line.startswith("#"):
                tok = line[1:].upper().split()
                i = 0
                while i < len(tok):
                    t = tok[i]
                    if t in _UNITS:
                        unit = _UNITS[t]
                    elif t in ("RI", "MA", "DB"):
                        fmt = t
                    elif t == "R":
                        z0 = float(tok[i + 1])
                        i += 1
                    elif t != "S":
                        raise ValueError(f"unsupported option {t!r} (only S-parameters are handled)")
                    i += 1
                continue
            if line.startswith("["):
                raise ValueError("Touchstone 2.0 keywords are not supported")
            values.extend(float(x) for x in line.split())
    per = 1 + 2 * n * n
    if len(values) % per:
        raise ValueError(f"data length {len(values)} not a multiple of {per}")
    arr = np.array(values).reshape(-1, per)
    f = arr[:, 0] * unit
    a, b = arr[:, 1::2], arr[:, 2::2]
    if fmt == "RI":
        v = a + 1j * b
    elif fmt == "MA":
        v = a * np.exp(1j * np.deg2rad(b))
    else:
        v = 10 ** (a / 20) * np.exp(1j * np.deg2rad(b))
    s = v.reshape(-1, n, n)
    if n == 2:
        s = np.swapaxes(s, 1, 2)  # file order S11 S21 S12 S22 is column-major for 2-ports
    return Network(f, s, z0, os.path.basename(path))


def write(path: str, net: Network, fmt: str = "RI", unit: str = "GHZ", comment: str = "") -> None:
    fmt, unit = fmt.upper(), unit.upper()
    n = net.n_ports
    s = np.swapaxes(net.s, 1, 2) if n == 2 else net.s
    flat = s.reshape(s.shape[0], -1)
    with open(path, "w") as fh:
        for c in comment.splitlines():
            fh.write(f"! {c}\n")
        fh.write(f"# {unit} S {fmt} R {net.z0:g}\n")
        for fi, row in zip(net.f, flat):
            if fmt == "RI":
                pairs = np.stack([row.real, row.imag], 1)
            elif fmt == "MA":
                pairs = np.stack([np.abs(row), np.rad2deg(np.angle(row))], 1)
            elif fmt == "DB":
                pairs = np.stack([20 * np.log10(np.maximum(np.abs(row), 1e-300)), np.rad2deg(np.angle(row))], 1)
            else:
                raise ValueError(fmt)
            nums = [f"{fi / _UNITS[unit]:.9g}"] + [f"{x:.12g}" for x in pairs.reshape(-1)]
            # Touchstone 1.x: at most 4 complex pairs per line for n >= 3
            if n <= 2:
                fh.write(" ".join(nums) + "\n")
            else:
                head, rest = nums[0], nums[1:]
                for k in range(0, len(rest), 2 * n):
                    chunk = rest[k:k + 2 * n]
                    fh.write((head + " " if k == 0 else "  ") + " ".join(chunk) + "\n")
