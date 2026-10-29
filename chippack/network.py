"""N-port network container and conversions (S, Z, Y, ABCD, T), cascading,
mixed-mode conversion and passivity / causality / reciprocity checks."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class Network:
    """Frequency-domain N-port. ``s`` has shape (n_freq, n_ports, n_ports)."""

    f: np.ndarray
    s: np.ndarray
    z0: float = 50.0
    name: str = ""

    def __post_init__(self):
        self.f = np.asarray(self.f, dtype=float)
        self.s = np.asarray(self.s, dtype=complex)
        if self.s.ndim != 3 or self.s.shape[0] != self.f.size or self.s.shape[1] != self.s.shape[2]:
            raise ValueError(f"bad S shape {self.s.shape} for {self.f.size} frequencies")

    @property
    def n_ports(self) -> int:
        return self.s.shape[1]

    def sij(self, i: int, j: int) -> np.ndarray:
        """1-based port indices, as in S21."""
        return self.s[:, i - 1, j - 1]

    def db(self, i: int, j: int) -> np.ndarray:
        return 20 * np.log10(np.maximum(np.abs(self.sij(i, j)), 1e-300))

    def interpolate(self, f_new: np.ndarray) -> "Network":
        """Interpolate magnitude (linear) and unwrapped phase onto ``f_new``."""
        mag = np.abs(self.s)
        ph = np.unwrap(np.angle(self.s), axis=0)
        out = np.empty((f_new.size,) + self.s.shape[1:], complex)
        for i in range(self.n_ports):
            for j in range(self.n_ports):
                m = np.interp(f_new, self.f, mag[:, i, j])
                p = np.interp(f_new, self.f, ph[:, i, j])
                out[:, i, j] = m * np.exp(1j * p)
        return Network(f_new, out, self.z0, self.name)

    def reorder(self, order: list[int]) -> "Network":
        """Permute ports; ``order`` lists old 1-based port numbers in their new order."""
        idx = np.asarray(order) - 1
        return Network(self.f, self.s[:, idx][:, :, idx], self.z0, self.name)

    def renormalize(self, z_new: float) -> "Network":
        return Network(self.f, z_to_s(s_to_z(self.s, self.z0), z_new), z_new, self.name)


# ------------------------------------------------------------- conversions

def s_to_z(s: np.ndarray, z0: float) -> np.ndarray:
    n = s.shape[-1]
    i = np.eye(n)
    return z0 * np.linalg.solve(i - s, i + s)


def z_to_s(z: np.ndarray, z0: float) -> np.ndarray:
    n = z.shape[-1]
    i = np.eye(n)
    return (z - z0 * i) @ np.linalg.inv(z + z0 * i)


def s_to_y(s: np.ndarray, z0: float) -> np.ndarray:
    return np.linalg.inv(s_to_z(s, z0))


def s2_to_abcd(s: np.ndarray, z0: float) -> np.ndarray:
    s11, s12, s21, s22 = s[..., 0, 0], s[..., 0, 1], s[..., 1, 0], s[..., 1, 1]
    d = 2 * s21
    a = ((1 + s11) * (1 - s22) + s12 * s21) / d
    b = z0 * ((1 + s11) * (1 + s22) - s12 * s21) / d
    c = ((1 - s11) * (1 - s22) - s12 * s21) / (z0 * d)
    dd = ((1 - s11) * (1 + s22) + s12 * s21) / d
    return np.stack([np.stack([a, b], -1), np.stack([c, dd], -1)], -2)


def abcd_to_s2(m: np.ndarray, z0: float) -> np.ndarray:
    a, b, c, d = m[..., 0, 0], m[..., 0, 1], m[..., 1, 0], m[..., 1, 1]
    den = a + b / z0 + c * z0 + d
    s11 = (a + b / z0 - c * z0 - d) / den
    s12 = 2 * (a * d - b * c) / den
    s21 = 2 / den
    s22 = (-a + b / z0 - c * z0 + d) / den
    return np.stack([np.stack([s11, s12], -1), np.stack([s21, s22], -1)], -2)


def s_to_t(s: np.ndarray) -> np.ndarray:
    """2N-port S (ports 1..N left, N+1..2N right) -> transfer (T) matrix for cascading."""
    n = s.shape[-1] // 2
    s11, s12 = s[..., :n, :n], s[..., :n, n:]
    s21, s22 = s[..., n:, :n], s[..., n:, n:]
    s21i = np.linalg.inv(s21)
    t11 = s12 - s11 @ s21i @ s22
    t12 = s11 @ s21i
    t21 = -s21i @ s22
    t22 = s21i
    return np.concatenate([np.concatenate([t11, t12], -1), np.concatenate([t21, t22], -1)], -2)


def t_to_s(t: np.ndarray) -> np.ndarray:
    n = t.shape[-1] // 2
    t11, t12 = t[..., :n, :n], t[..., :n, n:]
    t21, t22 = t[..., n:, :n], t[..., n:, n:]
    t22i = np.linalg.inv(t22)
    s11 = t12 @ t22i
    s12 = t11 - t12 @ t22i @ t21
    s21 = t22i
    s22 = -t22i @ t21
    return np.concatenate([np.concatenate([s11, s12], -1), np.concatenate([s21, s22], -1)], -2)


def cascade(*nets: Network) -> Network:
    """Cascade 2N-port networks (left ports 1..N, right ports N+1..2N)."""
    base = nets[0]
    t = s_to_t(base.s)
    for n in nets[1:]:
        if n.n_ports != base.n_ports or not np.allclose(n.f, base.f):
            raise ValueError("networks must share port count and frequency grid")
        if not np.isclose(n.z0, base.z0):
            n = n.renormalize(base.z0)
        t = t @ s_to_t(n.s)
    return Network(base.f, t_to_s(t), base.z0, "+".join(x.name for x in nets))


# ------------------------------------------------------------- mixed mode

def mixed_mode(net: Network, pairs: tuple[tuple[int, int], ...] = ((1, 3), (2, 4))) -> dict[str, np.ndarray]:
    """Single-ended 4-port -> mixed-mode blocks Sdd, Sdc, Scd, Scc (each n_freq x 2 x 2).

    ``pairs`` lists the (positive, negative) single-ended ports of each differential
    port, 1-based. The default assumes ports 1/3 at the input and 2/4 at the output.
    """
    if net.n_ports != 4:
        raise ValueError("mixed-mode conversion needs a 4-port")
    m = np.zeros((4, 4))
    for k, (p, n) in enumerate(pairs):
        m[k, p - 1], m[k, n - 1] = 1, -1  # differential
        m[k + 2, p - 1], m[k + 2, n - 1] = 1, 1  # common
    m /= np.sqrt(2)
    smm = m @ net.s @ m.T
    return {"dd": smm[:, :2, :2], "dc": smm[:, :2, 2:], "cd": smm[:, 2:, :2], "cc": smm[:, 2:, 2:]}


# ------------------------------------------------------------- checks

def passivity_violation(net: Network) -> tuple[bool, float]:
    """Passive iff the largest singular value of S <= 1 at every frequency."""
    sv = np.linalg.svd(net.s, compute_uv=False)[:, 0]
    worst = float(sv.max())
    return worst <= 1 + 1e-9, worst


def enforce_passivity(net: Network) -> Network:
    """Simple singular-value clipping (keeps singular vectors)."""
    u, sv, vh = np.linalg.svd(net.s)
    sv = np.minimum(sv, 1.0 - 1e-9)
    return Network(net.f, u @ (sv[..., None] * vh), net.z0, net.name)


def reciprocity_error(net: Network) -> float:
    return float(np.max(np.abs(net.s - np.swapaxes(net.s, 1, 2))))


def causality_check(h: np.ndarray, f: np.ndarray, n_pad: int = 8, tol: float = 0.02):
    """Check a transfer function sampled on a uniform grid starting near DC.

    The response is extended to DC, mirrored to negative frequencies, windowed and
    transformed to the time domain. For a causal system with a pure delay ``td``
    there should be negligible energy before t = 0 (we inspect the last quarter of
    the periodic impulse response, which corresponds to negative time).

    Returns (is_causal, ratio of non-causal to total energy).
    """
    df = f[1] - f[0]
    if not np.allclose(np.diff(f), df, rtol=1e-6):
        raise ValueError("uniform frequency grid required")
    n0 = int(round(f[0] / df))
    hh = np.r_[np.full(n0, h[0]), h] if n0 > 0 else h.copy()
    if f[0] > 0 and n0 > 0:
        hh[0] = np.real(hh[0])
    w = np.hanning(2 * hh.size)[hh.size:]
    spec = hh * w
    m = n_pad * spec.size
    imp = np.fft.irfft(spec, 2 * m)
    e = imp ** 2
    neg = e[int(1.5 * m):].sum()
    ratio = float(neg / e.sum())
    return ratio < tol, ratio


# ------------------------------------------------------------- model generators

def rlgc_line(f: np.ndarray, length_m: float, r_dc=5.0, r_skin=2e-4, l=3.5e-7, g_tan=0.02, c=1.4e-10,
              z0_ref: float = 50.0) -> Network:
    """2-port lossy transmission line from per-metre RLGC.

    R(f) = r_dc + r_skin * sqrt(f) (ohm/m), G(f) = 2 pi f C tan_delta (S/m).
    """
    w = 2 * np.pi * np.maximum(f, 1.0)
    r = r_dc + r_skin * np.sqrt(np.maximum(f, 0)) * (1 + 1j)  # skin effect also adds internal inductance
    g = w * c * g_tan
    zs = r + 1j * w * l
    yp = g + 1j * w * c
    gamma = np.sqrt(zs * yp)
    zc = np.sqrt(zs / yp)
    gl = gamma * length_m
    abcd = np.empty((f.size, 2, 2), complex)
    abcd[:, 0, 0] = np.cosh(gl)
    abcd[:, 0, 1] = zc * np.sinh(gl)
    abcd[:, 1, 0] = np.sinh(gl) / zc
    abcd[:, 1, 1] = np.cosh(gl)
    return Network(f, abcd_to_s2(abcd, z0_ref), z0_ref, f"line{length_m * 1e3:.1f}mm")


def coupled_line_4port(f: np.ndarray, length_m: float, z_even=55.0, z_odd=45.0, eps_even=4.0, eps_odd=3.6,
                       loss_db_per_m_per_sqrtghz=8.0, z0_ref=50.0) -> Network:
    """Symmetric lossy coupled pair from even/odd-mode lines (ports 1,3 input; 2,4 output)."""
    c0 = 299792458.0
    out = np.zeros((f.size, 4, 4), complex)
    modes = {}
    for name, zc, eps in (("e", z_even, eps_even), ("o", z_odd, eps_odd)):
        beta = 2 * np.pi * f * np.sqrt(eps) / c0
        alpha = loss_db_per_m_per_sqrtghz * np.sqrt(f / 1e9) / 8.686
        gl = (alpha + 1j * beta) * length_m
        abcd = np.empty((f.size, 2, 2), complex)
        abcd[:, 0, 0] = np.cosh(gl)
        abcd[:, 0, 1] = zc * np.sinh(gl)
        abcd[:, 1, 0] = np.sinh(gl) / zc
        abcd[:, 1, 1] = np.cosh(gl)
        modes[name] = abcd_to_s2(abcd, z0_ref)
    se, so = modes["e"], modes["o"]
    # single-ended ports: 1 = in+, 2 = out+, 3 = in-, 4 = out-
    for (a, b) in ((0, 0), (0, 1), (1, 0), (1, 1)):
        p, q = (0, 1)[a], (0, 1)[b]
        out[:, p, q] = 0.5 * (se[:, a, b] + so[:, a, b])
        out[:, p + 2, q + 2] = 0.5 * (se[:, a, b] + so[:, a, b])
        out[:, p, q + 2] = 0.5 * (se[:, a, b] - so[:, a, b])
        out[:, p + 2, q] = 0.5 * (se[:, a, b] - so[:, a, b])
    return Network(f, out, z0_ref, f"pair{length_m * 1e3:.1f}mm")
