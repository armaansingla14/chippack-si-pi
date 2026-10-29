"""2-D finite-element thermo-mechanical solver (bilinear Q4, plane stress or plane strain).

Models a package cross-section made of rectangular material blocks on a
structured mesh and computes the response to a uniform temperature change
(CTE mismatch): displacements, warpage of the bottom surface and stresses.
Symmetry is used: the left edge (x = 0) is the centre line (u_x = 0) and one
node is pinned vertically.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla


@dataclass(frozen=True)
class Material:
    name: str
    E: float  # Pa
    nu: float
    cte: float  # 1/K


SILICON = Material("Si", 130e9, 0.28, 2.6e-6)
COPPER = Material("Cu", 117e9, 0.34, 17e-6)
ORGANIC_SUB = Material("organic substrate", 25e9, 0.30, 15e-6)
UNDERFILL = Material("underfill", 8e9, 0.33, 30e-6)
MOLD = Material("mold compound", 20e9, 0.30, 9e-6)
GLASS = Material("glass core", 75e9, 0.22, 3.5e-6)


@dataclass
class Block:
    material: Material
    x0: float
    x1: float
    y0: float
    y1: float


def _d_matrix(m: Material, plane: str) -> tuple[np.ndarray, np.ndarray]:
    E, nu, a = m.E, m.nu, m.cte
    if plane == "stress":
        c = E / (1 - nu ** 2)
        D = c * np.array([[1, nu, 0], [nu, 1, 0], [0, 0, (1 - nu) / 2]])
        eps_t = a * np.array([1, 1, 0])
    else:
        c = E / ((1 + nu) * (1 - 2 * nu))
        D = c * np.array([[1 - nu, nu, 0], [nu, 1 - nu, 0], [0, 0, (1 - 2 * nu) / 2]])
        eps_t = (1 + nu) * a * np.array([1, 1, 0])
    return D, eps_t


_GP = np.array([-1, 1]) / np.sqrt(3)


def _q4_matrices(a: float, b: float, D: np.ndarray, eps_t: np.ndarray):
    """Stiffness and unit-dT thermal load for an a x b rectangle (thickness 1)."""
    Ke = np.zeros((8, 8))
    fe = np.zeros(8)
    Bs = []
    for xi in _GP:
        for eta in _GP:
            dN_dxi = 0.25 * np.array([-(1 - eta), (1 - eta), (1 + eta), -(1 + eta)])
            dN_deta = 0.25 * np.array([-(1 - xi), -(1 + xi), (1 + xi), (1 - xi)])
            dN_dx = dN_dxi * 2 / a
            dN_dy = dN_deta * 2 / b
            B = np.zeros((3, 8))
            B[0, 0::2] = dN_dx
            B[1, 1::2] = dN_dy
            B[2, 0::2] = dN_dy
            B[2, 1::2] = dN_dx
            w = a * b / 4
            Ke += B.T @ D @ B * w
            fe += B.T @ D @ eps_t * w
            Bs.append(B)
    return Ke, fe, Bs


@dataclass
class FEAResult:
    x: np.ndarray  # node x coords (nx+1)
    y: np.ndarray  # node y coords (ny+1)
    u: np.ndarray  # (ny+1, nx+1, 2)
    von_mises: np.ndarray  # (ny, nx) element-average
    sxx: np.ndarray
    syy: np.ndarray
    sxy: np.ndarray
    elem_material: np.ndarray

    def warpage(self, row: int = 0) -> float:
        """Out-of-plane deflection range along a horizontal node row (default: bottom surface)."""
        uy = self.u[row, :, 1]
        return float(uy.max() - uy.min())

    def curvature(self, row: int = 0) -> float:
        """Least-squares curvature of the deflected node row (u_y ~ k x^2 / 2)."""
        c = np.polyfit(self.x, self.u[row, :, 1], 2)
        return float(2 * c[0])


def solve(blocks: list[Block], dT: float, nx: int = 120, ny_per_mm: float | None = None,
          y_breaks_n: int = 4, plane: str = "strain") -> FEAResult:
    """Mesh the bounding box of ``blocks`` (cells without material are void) and solve.

    The y-mesh resolves every block interface with ``y_breaks_n`` elements per layer.
    """
    x_max = max(b.x1 for b in blocks)
    ys = sorted({b.y0 for b in blocks} | {b.y1 for b in blocks})
    y_nodes = [ys[0]]
    for y0, y1 in zip(ys[:-1], ys[1:]):
        y_nodes.extend(np.linspace(y0, y1, y_breaks_n + 1)[1:])
    y = np.array(y_nodes)
    xs = sorted({b.x0 for b in blocks} | {b.x1 for b in blocks})
    x = np.linspace(0, x_max, nx + 1)
    x = np.unique(np.r_[x, xs])
    nxe, nye = x.size - 1, y.size - 1
    xc = 0.5 * (x[1:] + x[:-1])
    yc = 0.5 * (y[1:] + y[:-1])
    mat = np.full((nye, nxe), -1)
    mats: list[Material] = []
    for b in blocks:
        if b.material not in mats:
            mats.append(b.material)
        mi = mats.index(b.material)
        sel = (yc[:, None] > b.y0) & (yc[:, None] < b.y1) & (xc[None, :] > b.x0) & (xc[None, :] < b.x1)
        mat[sel] = mi
    n_nodes = (nxe + 1) * (nye + 1)

    def nid(i, j):
        return j * (nxe + 1) + i

    rows, cols, vals = [], [], []
    f = np.zeros(2 * n_nodes)
    cache = {}
    used = np.zeros(n_nodes, bool)
    for j in range(nye):
        for i in range(nxe):
            mi = mat[j, i]
            if mi < 0:
                continue
            a, bb = x[i + 1] - x[i], y[j + 1] - y[j]
            key = (mi, round(a, 15), round(bb, 15))
            if key not in cache:
                D, et = _d_matrix(mats[mi], plane)
                cache[key] = _q4_matrices(a, bb, D, et)
            Ke, fe, _ = cache[key]
            nodes = [nid(i, j), nid(i + 1, j), nid(i + 1, j + 1), nid(i, j + 1)]
            used[nodes] = True
            dofs = np.array([[2 * n, 2 * n + 1] for n in nodes]).ravel()
            rows.append(np.repeat(dofs, 8))
            cols.append(np.tile(dofs, 8))
            vals.append(Ke.ravel())
            f[dofs] += fe * dT
    K = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                      shape=(2 * n_nodes, 2 * n_nodes))
    fixed = set()
    for j in range(nye + 1):
        fixed.add(2 * nid(0, j))  # symmetry: u_x = 0 on x = 0
    j_pin = next(j for j in range(nye + 1) if used[nid(0, j)])
    fixed.add(2 * nid(0, j_pin) + 1)  # pin vertical rigid motion
    for n in np.nonzero(~used)[0]:
        fixed.update((2 * n, 2 * n + 1))
    free = np.setdiff1d(np.arange(2 * n_nodes), np.fromiter(fixed, int))
    u = np.zeros(2 * n_nodes)
    u[free] = spla.spsolve(K[free][:, free].tocsc(), f[free])
    U = u.reshape(nye + 1, nxe + 1, 2)
    # stresses at element centre
    sxx = np.zeros((nye, nxe))
    syy = np.zeros_like(sxx)
    sxy = np.zeros_like(sxx)
    szz_extra = np.zeros_like(sxx)
    for j in range(nye):
        for i in range(nxe):
            mi = mat[j, i]
            if mi < 0:
                continue
            a, bb = x[i + 1] - x[i], y[j + 1] - y[j]
            D, et = _d_matrix(mats[mi], plane)
            B = np.zeros((3, 8))
            dN_dx = 0.25 * np.array([-1, 1, 1, -1]) * 2 / a
            dN_dy = 0.25 * np.array([-1, -1, 1, 1]) * 2 / bb
            B[0, 0::2] = dN_dx
            B[1, 1::2] = dN_dy
            B[2, 0::2] = dN_dy
            B[2, 1::2] = dN_dx
            nodes = [nid(i, j), nid(i + 1, j), nid(i + 1, j + 1), nid(i, j + 1)]
            ue = np.array([[u[2 * n], u[2 * n + 1]] for n in nodes]).ravel()
            s = D @ (B @ ue - et * dT)
            sxx[j, i], syy[j, i], sxy[j, i] = s
            if plane == "strain":
                m_ = mats[mi]
                szz_extra[j, i] = m_.nu * (s[0] + s[1]) - m_.E * m_.cte * dT
    szz = szz_extra
    vm = np.sqrt(0.5 * ((sxx - syy) ** 2 + (syy - szz) ** 2 + (szz - sxx) ** 2) + 3 * sxy ** 2)
    vm[mat < 0] = 0
    return FEAResult(x, y, U, vm, sxx, syy, sxy, mat)


def timoshenko_curvature(m1: Material, t1: float, m2: Material, t2: float, dT: float, plane: str = "stress"):
    """Bimetal strip curvature (Timoshenko 1925). Layer 1 is the bottom layer.

    For plane strain the biaxial-free (wide plate in one direction) correction
    E -> E/(1-nu^2), alpha -> (1+nu) alpha is applied.
    """
    def eff(m):
        if plane == "stress":
            return m.E, m.cte
        return m.E / (1 - m.nu ** 2), (1 + m.nu) * m.cte
    E1, a1 = eff(m1)
    E2, a2 = eff(m2)
    m = t1 / t2
    n = E1 / E2
    h = t1 + t2
    k = 6 * (a2 - a1) * dT * (1 + m) ** 2 / (h * (3 * (1 + m) ** 2 + (1 + m * n) * (m ** 2 + 1 / (m * n))))
    return -k  # sign: positive k bends toward layer 1 when a2 > a1 and dT > 0
