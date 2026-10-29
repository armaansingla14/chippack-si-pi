"""Steady-state finite-volume thermal solver for layered packages ("CFD-lite").

3-D conduction on a structured grid of stacked layers, with
- per-layer (optionally per-cell) isotropic or orthotropic conductivity,
- volumetric heat sources (power maps),
- convective top surface to either a fixed ambient or a 1-D coolant stream that
  flows in +x and heats up by energy balance (upwind finite volume), and
- convective bottom surface (board side). Lateral faces are adiabatic.

It is a reduced-order model: the coolant is a lumped advecting stream with a
given heat-transfer coefficient, not a Navier-Stokes solution.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla


@dataclass
class Layer:
    name: str
    thickness: float  # m
    k: float | np.ndarray  # W/mK, scalar or (nx, ny) map
    nz: int = 2
    k_z: float | np.ndarray | None = None  # through-plane conductivity if different


@dataclass
class Coolant:
    """Fluid stream over the top surface (flowing in +x)."""

    mass_flow: float = 0.02  # kg/s across the full width
    cp: float = 4180.0  # J/kgK (water)
    t_in: float = 25.0  # degC


@dataclass
class ThermalModel:
    lx: float
    ly: float
    nx: int
    ny: int
    layers: list[Layer]  # bottom -> top
    h_top: float = 2e4
    h_bottom: float = 10.0
    t_ambient: float = 25.0
    coolant: Coolant | None = None
    power: dict[int, np.ndarray] = field(default_factory=dict)  # layer index -> (nx, ny) W map

    def _grid(self):
        z_edges = [0.0]
        cell_layer = []
        for li, L in enumerate(self.layers):
            for _ in range(L.nz):
                z_edges.append(z_edges[-1] + L.thickness / L.nz)
                cell_layer.append(li)
        return np.array(z_edges), np.array(cell_layer)

    def solve(self):
        nx, ny = self.nx, self.ny
        dx, dy = self.lx / nx, self.ly / ny
        z_edges, cl = self._grid()
        nz = cl.size
        dz = np.diff(z_edges)
        n_solid = nx * ny * nz
        use_fluid = self.coolant is not None
        n = n_solid + (nx * ny if use_fluid else 0)

        def idx(i, j, k):
            return (k * ny + j) * nx + i

        kxy = np.empty((nx, ny, nz))
        kz = np.empty((nx, ny, nz))
        for k, li in enumerate(cl):
            L = self.layers[li]
            kxy[:, :, k] = np.broadcast_to(L.k, (nx, ny))
            kz[:, :, k] = np.broadcast_to(L.k_z if L.k_z is not None else L.k, (nx, ny))
        rows, cols, vals = [], [], []
        b = np.zeros(n)

        def add(r, c, v):
            rows.append(r)
            cols.append(c)
            vals.append(v)

        I, J, K = np.meshgrid(np.arange(nx), np.arange(ny), np.arange(nz), indexing="ij")
        # x-direction conductances
        for axis, (d, area_fn) in enumerate(((dx, lambda k: dy * dz[k]), (dy, lambda k: dx * dz[k]))):
            if axis == 0:
                a, bq = (slice(0, nx - 1), slice(None)), (slice(1, nx), slice(None))
            else:
                a, bq = (slice(None), slice(0, ny - 1)), (slice(None), slice(1, ny))
            k1, k2 = kxy[a], kxy[bq]
            area = np.broadcast_to(dz[None, None, :], k1.shape) * (dy if axis == 0 else dx)
            g = area / (0.5 * d / k1 + 0.5 * d / k2)
            p = idx(I[a], J[a], K[a]).ravel()
            q = idx(I[bq], J[bq], K[bq]).ravel()
            gv = g.ravel()
            rows.extend(np.r_[p, q, p, q])
            cols.extend(np.r_[p, q, q, p])
            vals.extend(np.r_[gv, gv, -gv, -gv])
        # z-direction
        k1, k2 = kz[:, :, :-1], kz[:, :, 1:]
        g = dx * dy / (0.5 * dz[None, None, :-1] / k1 + 0.5 * dz[None, None, 1:] / k2)
        p = idx(I[:, :, :-1], J[:, :, :-1], K[:, :, :-1]).ravel()
        q = idx(I[:, :, 1:], J[:, :, 1:], K[:, :, 1:]).ravel()
        gv = g.ravel()
        rows.extend(np.r_[p, q, p, q])
        cols.extend(np.r_[p, q, q, p])
        vals.extend(np.r_[gv, gv, -gv, -gv])
        # bottom convection to ambient
        gb = 1.0 / (0.5 * dz[0] / kz[:, :, 0] + 1.0 / self.h_bottom) * dx * dy
        pb = idx(I[:, :, 0], J[:, :, 0], 0).ravel()
        rows.extend(pb)
        cols.extend(pb)
        vals.extend(gb.ravel())
        np.add.at(b, pb, gb.ravel() * self.t_ambient)
        # top convection
        gt = (1.0 / (0.5 * dz[-1] / kz[:, :, -1] + 1.0 / self.h_top) * dx * dy).ravel()
        pt = idx(I[:, :, -1], J[:, :, -1], nz - 1).ravel()
        rows.extend(pt)
        cols.extend(pt)
        vals.extend(gt)
        if not use_fluid:
            np.add.at(b, pt, gt * self.t_ambient)
        else:
            c = self.coolant
            mcp = c.mass_flow * c.cp / ny  # per y-strip
            fi = n_solid + (J[:, :, 0] * nx + I[:, :, 0]).ravel()
            # solid <-> fluid coupling
            rows.extend(np.r_[pt, fi, fi])
            cols.extend(np.r_[fi, pt, fi])
            vals.extend(np.r_[-gt, -gt, gt])
            # fluid advection (upwind): mcp * (Tf_i - Tf_{i-1})
            ii = I[:, :, 0].ravel()
            rows.extend(fi)
            cols.extend(fi)
            vals.extend(np.full(fi.size, mcp))
            up = ii > 0
            rows.extend(fi[up])
            cols.extend(fi[up] - 1)
            vals.extend(np.full(up.sum(), -mcp))
            np.add.at(b, fi[~up], mcp * c.t_in)
        # heat sources: spread over the cells of the layer
        for li, pmap in self.power.items():
            ks = np.nonzero(cl == li)[0]
            for k in ks:
                b[idx(I[:, :, 0], J[:, :, 0], k).ravel()] += (np.asarray(pmap) / ks.size).ravel()
        A = sp.csr_matrix((np.asarray(vals, float), (np.asarray(rows), np.asarray(cols))), shape=(n, n))
        t = spla.spsolve(A, b)
        ts = t[:n_solid].reshape(nz, ny, nx).transpose(2, 1, 0)
        tf = t[n_solid:].reshape(ny, nx).T if use_fluid else None
        q_top = float(np.sum(gt * (t[pt] - (t[n_solid + (J[:, :, 0] * nx + I[:, :, 0]).ravel()] if use_fluid
                                            else self.t_ambient))))
        q_bot = float(np.sum(gb.ravel() * (t[pb] - self.t_ambient)))
        return ThermalResult(ts, tf, z_edges, cl, q_top, q_bot, sum(float(np.sum(p)) for p in self.power.values()))


@dataclass
class ThermalResult:
    t: np.ndarray  # (nx, ny, nz) cell temperatures
    t_fluid: np.ndarray | None
    z_edges: np.ndarray
    cell_layer: np.ndarray
    q_top: float
    q_bottom: float
    p_total: float

    def layer_max(self, li: int) -> float:
        return float(self.t[:, :, self.cell_layer == li].max())

    @property
    def t_max(self) -> float:
        return float(self.t.max())

    @property
    def energy_balance_error(self) -> float:
        return abs(self.q_top + self.q_bottom - self.p_total) / max(self.p_total, 1e-12)


def block_map(nx: int, ny: int, lx: float, ly: float, x0: float, y0: float, w: float, h: float, power: float):
    """Uniform power over a rectangle [x0, x0+w] x [y0, y0+h] (area-weighted on the grid)."""
    xe = np.linspace(0, lx, nx + 1)
    ye = np.linspace(0, ly, ny + 1)
    fx = np.clip(np.minimum(xe[1:], x0 + w) - np.maximum(xe[:-1], x0), 0, None)
    fy = np.clip(np.minimum(ye[1:], y0 + h) - np.maximum(ye[:-1], y0), 0, None)
    a = np.outer(fx, fy)
    return power * a / a.sum()
