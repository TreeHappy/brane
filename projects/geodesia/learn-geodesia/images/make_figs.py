#!/usr/bin/env python3
"""Figures for learn-geodesia.

Run from this directory:  python3 make_figs.py

The surface is the same paraboloid used in the Lean demo:
    Phi(u, v) = (u, v, 0.3 (u^2 + v^2)).
Everything is computed with small explicit RK4 / gradient-descent loops so the
figures can be reproduced without SciPy.
"""

from __future__ import annotations

import math
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm

HERE = os.path.dirname(os.path.abspath(__file__))
C = 0.3  # z = C (u^2 + v^2)


# ---------------------------------------------------------------------------
# Geometry
# ---------------------------------------------------------------------------
def phi(q):
    u, v = q
    return np.array([u, v, C * (u * u + v * v)])


def ju(q):
    u, _ = q
    return np.array([1.0, 0.0, 2 * C * u])


def jv(q):
    _, v = q
    return np.array([0.0, 1.0, 2 * C * v])


def jac(q):
    return np.column_stack([ju(q), jv(q)])


def metric(q):
    u, v = q
    return np.array([
        [1 + 4 * C * C * u * u, 4 * C * C * u * v],
        [4 * C * C * u * v, 1 + 4 * C * C * v * v],
    ])


def christoffel(q):
    u, v = q
    gu = np.array([[8 * C * C * u, 4 * C * C * v],
                   [4 * C * C * v, 0.0]])
    gv = np.array([[0.0, 4 * C * C * u],
                   [4 * C * C * u, 8 * C * C * v]])
    gi = np.linalg.inv(metric(q))
    g = np.zeros((2, 2, 2))
    for k in range(2):
        for i in range(2):
            for j in range(2):
                total = 0.0
                for m in range(2):
                    dg_jm = (gu if i == 0 else gv)[j, m]
                    dg_im = (gu if j == 0 else gv)[i, m]
                    dg_ij = (gu if m == 0 else gv)[i, j]
                    total += gi[k, m] * (dg_jm + dg_im - dg_ij)
                g[k, i, j] = 0.5 * total
    return g


def geodesic_rhs(q, p):
    g = christoffel(q)
    acc = np.zeros(2)
    for k in range(2):
        acc[k] = -sum(g[k, i, j] * p[i] * p[j]
                      for i in range(2) for j in range(2))
    return acc


def rk4_step(q, p, h):
    k1q, k1p = p, geodesic_rhs(q, p)
    k2q = p + 0.5 * h * k1p
    k2p = geodesic_rhs(q + 0.5 * h * k1q, p + 0.5 * h * k1p)
    k3q = p + 0.5 * h * k2p
    k3p = geodesic_rhs(q + 0.5 * h * k2q, p + 0.5 * h * k2p)
    k4q = p + h * k3p
    k4p = geodesic_rhs(q + h * k3q, p + h * k3p)
    qn = q + h / 6 * (k1q + 2 * k2q + 2 * k3q + k4q)
    pn = p + h / 6 * (k1p + 2 * k2p + 2 * k3p + k4p)
    return qn, pn


def integrate(q, p, length, n=400):
    h = length / n
    qs = [q.copy()]
    ps = [p.copy()]
    for _ in range(n):
        q, p = rk4_step(q, p, h)
        qs.append(q.copy())
        ps.append(p.copy())
    return np.array(qs), np.array(ps)


def tangent_frame(q):
    e1 = ju(q) / np.linalg.norm(ju(q))
    raw = jv(q) - np.dot(jv(q), e1) * e1
    e2 = raw / np.linalg.norm(raw)
    return e1, e2


def initial_velocity(q, alpha):
    e1, e2 = tangent_frame(q)
    w = math.cos(alpha) * e1 + math.sin(alpha) * e2
    return np.linalg.inv(metric(q)) @ jac(q).T @ w


def normal(q):
    n = np.cross(ju(q), jv(q))
    return n / np.linalg.norm(n)


# ---------------------------------------------------------------------------
# Shooting
# ---------------------------------------------------------------------------
def shoot(q0, alpha, length, n=400):
    p0 = initial_velocity(q0, alpha)
    qs, ps = integrate(q0, p0, length, n)
    return qs


def shooting_solve(q0, target, n=400, iters=8):
    chord = np.linalg.norm(phi(target) - phi(q0))
    d = phi(target) - phi(q0)
    xi = np.linalg.inv(metric(q0)) @ jac(q0).T @ d
    w = jac(q0) @ xi
    e1, e2 = tangent_frame(q0)
    alpha0 = math.atan2(np.dot(w, e2), np.dot(w, e1))

    x = np.array([alpha0, 1.15 * chord])
    history = []
    paths = []

    def residual(z):
        qs = shoot(q0, z[0], z[1], n)
        return qs[-1] - target

    for _ in range(iters):
        r = residual(x)
        history.append(np.linalg.norm(r))
        paths.append(shoot(q0, x[0], x[1], n))
        if np.linalg.norm(r) < 1e-12:
            break
        jac_num = np.zeros((2, 2))
        eps = 1e-6
        for j in range(2):
            z = x.copy()
            z[j] += eps
            jac_num[:, j] = (residual(z) - r) / eps
        try:
            dx = np.linalg.solve(jac_num, -r)
        except np.linalg.LinAlgError:
            break
        step = 1.0
        for _ in range(20):
            if np.linalg.norm(residual(x + step * dx)) < np.linalg.norm(r):
                break
            step *= 0.5
        x = x + step * dx
    history.append(np.linalg.norm(residual(x)))
    paths.append(shoot(q0, x[0], x[1], n))
    return x, history, paths


# ---------------------------------------------------------------------------
# Relaxation
# ---------------------------------------------------------------------------
def poly_length(qs):
    return sum(np.linalg.norm(phi(qs[i + 1]) - phi(qs[i]))
               for i in range(len(qs) - 1))


def poly_grad(qs):
    grad = np.zeros_like(qs)
    for i in range(1, len(qs) - 1):
        d0 = phi(qs[i]) - phi(qs[i - 1])
        d1 = phi(qs[i + 1]) - phi(qs[i])
        t0 = d0 / np.linalg.norm(d0)
        t1 = d1 / np.linalg.norm(d1)
        grad[i] = jac(qs[i]).T @ (t0 - t1)
    return grad


def relax(q0, q1, n=24, iters=3000):
    t = np.linspace(0, 1, n + 1)
    qs = np.linspace(q0, q1, n + 1)
    d = q1 - q0
    perp = np.array([-d[1], d[0]])
    perp /= np.linalg.norm(perp)
    qs = qs + 0.45 * np.sin(np.pi * t)[:, None] * perp

    snaps = {}
    lengths = []
    for it in range(iters + 1):
        if it in (0, 50, 250, 500, 1000, 2000, iters):
            snaps[it] = qs.copy()
        lengths.append(poly_length(qs))
        if it == iters:
            break
        grad = poly_grad(qs)
        if np.linalg.norm(grad) < 1e-12:
            break
        step = 0.5
        l0 = poly_length(qs)
        for _ in range(30):
            cand = qs - step * grad
            l1 = poly_length(cand)
            if l1 < l0 - 1e-4 * step * np.sum(grad * grad):
                break
            step *= 0.5
        qs = cand
    return snaps, lengths


# ---------------------------------------------------------------------------
# Plot helpers
# ---------------------------------------------------------------------------
def draw_surface(ax, ulim=(-0.15, 1.65), vlim=(-0.35, 1.35), alpha=0.30):
    u = np.linspace(*ulim, 70)
    v = np.linspace(*vlim, 70)
    uu, vv = np.meshgrid(u, v)
    zz = C * (uu * uu + vv * vv)
    ax.plot_surface(uu, vv, zz, cmap="Greys", alpha=alpha,
                    linewidth=0, antialiased=True, rcount=40, ccount=40)
    ax.set_xlabel("u")
    ax.set_ylabel("v")
    ax.set_zlabel("z")
    ax.view_init(elev=28, azim=-125)


def plot_path(ax, qs, **kwargs):
    pts = np.array([phi(q) for q in qs])
    ax.plot(pts[:, 0], pts[:, 1], pts[:, 2], **kwargs)


def arrow(ax, p, d, color, label=None, scale=0.55):
    ax.quiver(p[0], p[1], p[2], d[0], d[1], d[2],
              color=color, linewidth=2.0, arrow_length_ratio=0.18,
              length=scale)
    if label:
        tip = p + scale * d
        ax.text(tip[0], tip[1], tip[2], " " + label, color=color,
                fontsize=9, weight="bold")


# ---------------------------------------------------------------------------
# Figure 1: fundamental forms
# ---------------------------------------------------------------------------
def fig_fundamental_forms():
    fig = plt.figure(figsize=(8.2, 6.4))
    ax = fig.add_subplot(111, projection="3d")
    draw_surface(ax, (-0.05, 1.95), (-1.0, 1.0), alpha=0.28)

    p = np.array([1.0, 0.0])
    P = phi(p)
    eu, ev = ju(p), jv(p)
    n = normal(p)

    # tangent plane patch
    s = np.linspace(-0.55, 0.55, 2)
    t = np.linspace(-0.55, 0.55, 2)
    ss, tt = np.meshgrid(s, t)
    patch = (P[:, None, None]
             + ss[None, :, :] * eu[:, None, None] / np.linalg.norm(eu)
             + tt[None, :, :] * ev[:, None, None])
    ax.plot_surface(patch[0], patch[1], patch[2], color="#8ecae6",
                    alpha=0.35, linewidth=0)

    arrow(ax, P, eu / np.linalg.norm(eu), "#1d3557", r"$\Phi_u$", 0.75)
    arrow(ax, P, ev / np.linalg.norm(ev), "#1d3557", r"$\Phi_v$", 0.75)
    arrow(ax, P, n, "#e63946", r"$n$", 0.85)

    # A non-geodesic curve through P, with its acceleration split into
    # normal (bending) and tangential (geodesic curvature) parts.
    ts = np.linspace(-0.9, 0.9, 160)
    qs = np.column_stack([1.0 + 0.4 * (np.cos(ts) - 1.0), 0.4 * np.sin(ts)])
    plot_path(ax, qs, color="#ff006e", linewidth=2.8, label="curve")

    q = np.array([1.0, 0.0])
    qdot = np.array([0.0, 0.4])
    qddot = np.array([-0.4, 0.0])
    gamma_ddot = (qddot[0] * ju(q) + qddot[1] * jv(q)
                  + qdot[0] ** 2 * np.array([0, 0, 2 * C])
                  + 2 * qdot[0] * qdot[1] * np.zeros(3)
                  + qdot[1] ** 2 * np.array([0, 0, 2 * C]))
    normal_part = np.dot(gamma_ddot, n) * n
    tangential_part = gamma_ddot - normal_part

    arrow(ax, P, gamma_ddot, "black", r"$\gamma''$", 0.9)
    arrow(ax, P, normal_part, "#0077b6", r"$II(\xi,\xi)n$", 2.5)
    arrow(ax, P, tangential_part, "#f77f00", r"$\kappa_g$ part", 1.1)
    ax.scatter(*P, color="black", s=35, depthshade=False)

    ax.set_title("First and second fundamental forms\n"
                 r"$I=J^{\mathrm{T}}J$ (ruler),  $II_{ij}=\langle n,\Phi_{ij}\rangle$ (bending)",
                 fontsize=11)
    ax.legend(loc="upper left", fontsize=9)
    ax.set_box_aspect((1, 1, 0.65))
    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(HERE, f"fundamental-forms.{ext}"),
                    dpi=170, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Figure 2: shooting
# ---------------------------------------------------------------------------
def fig_shooting():
    q0 = np.array([0.4, 0.2])
    target = np.array([1.2, 0.9])
    x, hist, paths = shooting_solve(q0, target)
    chord = np.linalg.norm(phi(target) - phi(q0))
    L0 = 1.15 * chord
    alpha0 = x[0]

    fig = plt.figure(figsize=(12.2, 5.2))
    ax = fig.add_subplot(121, projection="3d")
    draw_surface(ax)

    # Several shots with the *initial* length: visible misses.
    for k, da in enumerate(np.linspace(-0.28, 0.28, 5)):
        qs = shoot(q0, alpha0 + da, L0)
        plot_path(ax, qs, color="#adb5bd", linewidth=1.0, alpha=0.9,
                  linestyle="--")
        ax.scatter(*phi(qs[-1]), color="#6c757d", s=12, depthshade=False)

    plot_path(ax, paths[-1], color="#d00000", linewidth=3.0,
              label="converged geodesic")
    ax.scatter(*phi(q0), color="#2a9d8f", s=70, depthshade=False,
               label="source $a$")
    ax.scatter(*phi(target), color="#e63946", s=90, marker="*",
               depthshade=False, label="target $b$")
    ax.set_title("Shooting: grey = misses at the initial $L$, red = converged",
                 fontsize=10)
    ax.legend(loc="upper left", fontsize=8)
    ax.set_box_aspect((1, 1, 0.62))

    ax2 = fig.add_subplot(122)
    ax2.semilogy(range(len(hist)), hist, "o-", color="#d00000")
    ax2.set_xlabel("Newton iteration")
    ax2.set_ylabel(r"$\|R\|=\|q(L)-b\|$")
    ax2.set_title(f"Converged: $\\alpha={x[0]:.5f}$, $L={x[1]:.5f}$",
                  fontsize=10)
    ax2.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(HERE, f"shooting.{ext}"),
                    dpi=170, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Figure 3: relaxation
# ---------------------------------------------------------------------------
def fig_relaxation():
    q0 = np.array([0.4, 0.2])
    q1 = np.array([1.2, 0.9])
    snaps, lengths = relax(q0, q1)

    fig = plt.figure(figsize=(12.2, 5.2))
    ax = fig.add_subplot(121, projection="3d")
    draw_surface(ax)

    keys = sorted(snaps)
    colors = cm.viridis(np.linspace(0.05, 0.9, len(keys)))
    for color, it in zip(colors, keys):
        plot_path(ax, snaps[it], color=color, linewidth=2.0,
                  label=f"iter {it}" if it in (0, 250, 1000, keys[-1]) else None)
    ax.scatter(*phi(q0), color="#2a9d8f", s=70, depthshade=False)
    ax.scatter(*phi(q1), color="#e63946", s=90, marker="*", depthshade=False)
    ax.set_title("Relaxation: a bent path tightens into a geodesic", fontsize=10)
    ax.legend(loc="upper left", fontsize=8)
    ax.set_box_aspect((1, 1, 0.62))

    ax2 = fig.add_subplot(122)
    ax2.semilogy(lengths, color="#0077b6", linewidth=2.0)
    ax2.set_xlabel("iteration")
    ax2.set_ylabel("polyline length")
    ax2.set_title("Length decreases monotonically", fontsize=10)
    ax2.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(HERE, f"relaxation.{ext}"),
                    dpi=170, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    fig_fundamental_forms()
    fig_shooting()
    fig_relaxation()
    print("wrote fundamental-forms, shooting, relaxation (png + svg)")
