#!/usr/bin/env python3
"""Figures for learn-geodesia.

Run from this directory:  python3 make_figs.py

The surface is the same paraboloid used in the Lean demo:
    Phi(u, v) = (u, v, 0.3 (u^2 + v^2)).
Everything is computed with small explicit RK4 / gradient-descent loops so the
figures can be reproduced without SciPy.

Figures for chapter 1 are rendered in a book-illustration style using the
Catppuccin Latte palette.
"""

from __future__ import annotations

import math
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# Catppuccin Latte palette (https://github.com/catppuccin/catppuccin).
LATTE = {
    "rosewater": "#dc8a78", "flamingo": "#dd7878", "pink": "#ea76cb",
    "mauve": "#8839ef",     "red": "#d20f39",      "maroon": "#e64553",
    "peach": "#fe640b",     "yellow": "#df8e1d",   "green": "#40a02b",
    "teal": "#179299",      "sky": "#04a5e5",      "sapphire": "#209fb5",
    "blue": "#1e66f5",      "lavender": "#7287fd",
    "text": "#4c4f69",      "subtext1": "#5c5f77", "subtext0": "#6c6f85",
    "overlay2": "#7c7f93",  "overlay1": "#8c8fa1", "overlay0": "#9ca0b0",
    "surface2": "#acb0be",  "surface1": "#bcc0cc", "surface0": "#ccd0da",
    "base": "#eff1f5",      "mantle": "#e6e9ef",   "crust": "#dce0e8",
}

# Book-style rendering: LaTeX-flavoured serif labels on Latte canvases.
matplotlib.rcParams.update({
    "mathtext.fontset": "cm",
    "font.family": "serif",
    "axes.unicode_minus": False,
    "figure.facecolor": LATTE["base"],
    "axes.facecolor": LATTE["base"],
    "savefig.facecolor": LATTE["base"],
})

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
# Catppuccin Latte + book-style helpers: off-white canvas, shaded translucent
# surface with parameter grid-lines, round-shaft arrows with cone heads,
# matte point markers, dashed construction lines, serif math labels.
# ---------------------------------------------------------------------------
SURF = lambda U, V: np.stack([U, V, C * (U * U + V * V)])
DOMU = (-0.08, 2.02)
DOMV = (-1.06, 1.06)


def style3d(ax, xlim, ylim, zlim, elev=22, azim=-125, aspect=(1, 1, 0.62)):
    ax.set_proj_type("ortho")
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_zlim(*zlim)
    ax.view_init(elev=elev, azim=azim)
    ax.set_box_aspect(aspect)
    ax.set_facecolor(LATTE["base"])
    ax.set_axis_off()
    ax.computed_zorder = False   # draw strictly in our zorder handoff


def shaded_surface(ax, S, ulim, vlim, color=LATTE["crust"], alpha=0.96,
                   n=90):
    U = np.linspace(*ulim, n)
    V = np.linspace(*vlim, n)
    uu, vv = np.meshgrid(U, V)
    X, Y, Z = S(uu, vv)
    pc = ax.plot_surface(X, Y, Z, color=color, alpha=alpha, shade=True,
                         linewidth=0, antialiased=True)
    pc.set_zorder(1)      # sit behind everything drawn after it
    return pc


def iso_lines(ax, S, ulim, vlim, nu=10, nv=9, color=LATTE["overlay0"],
              lw=0.8, alpha=0.35):
    for u0 in np.linspace(*ulim, nu):
        V = np.linspace(*vlim, 80)
        U = np.full_like(V, u0)
        X, Y, Z = S(U, V)
        ax.plot(X, Y, Z, color=color, lw=lw, alpha=alpha, zorder=2)
    for v0 in np.linspace(*vlim, nv):
        U = np.linspace(*ulim, 80)
        V = np.full_like(U, v0)
        X, Y, Z = S(U, V)
        ax.plot(X, Y, Z, color=color, lw=lw, alpha=alpha, zorder=2)


def rim_lines(ax, S, ulim, vlim, color=LATTE["surface2"], lw=1.1, alpha=0.9):
    """Crisp outline along the four boundary curves of the patch."""
    U = np.linspace(*ulim, 90)
    V = np.linspace(*vlim, 90)
    for v0 in (vlim[0], vlim[1]):
        X, Y, Z = S(U, np.full_like(U, v0))
        ax.plot(X, Y, Z, color=color, lw=lw, alpha=alpha, zorder=3)
    for u0 in (ulim[0], ulim[1]):
        X, Y, Z = S(np.full_like(V, u0), V)
        ax.plot(X, Y, Z, color=color, lw=lw, alpha=alpha, zorder=3)


def _perp_basis(u):
    h = np.array([0.0, 0.0, 1.0])
    if abs(u @ h) > 0.9:
        h = np.array([0.0, 1.0, 0.0])
    e1 = np.cross(u, h)
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(u, e1)
    return e1, e2


def fancy_arrow(ax, base, vec, color, lw=2.2, head=None, alpha=1.0):
    """A book-style arrow: round-capped shaft + filled cone head."""
    base = np.asarray(base, float)
    vec = np.asarray(vec, float)
    L = np.linalg.norm(vec)
    if L < 1e-9:
        return base
    u = vec / L
    if head is None:
        head = min(0.15, 0.28 * max(L, 1e-9))
    shaft_end = base + vec - head * u
    ax.plot([base[0], shaft_end[0]], [base[1], shaft_end[1]],
            [base[2], shaft_end[2]], color=color, lw=lw, alpha=alpha,
            solid_capstyle="round", zorder=6)
    e1, e2 = _perp_basis(u)
    apex = base + vec
    ring = (shaft_end[None, :]
            + (0.5 * head) * (np.outer(np.cos(np.linspace(0, 2 * np.pi, 5)), e1)
                              + np.outer(np.sin(np.linspace(0, 2 * np.pi, 5)), e2)))
    tris = [[apex, ring[i], ring[i + 1]] for i in range(4)]
    ax.add_collection3d(Poly3DCollection(
        tris, facecolors=color, edgecolors="none", alpha=alpha, zorder=7))
    return apex


def ball(ax, p, s=36, color=LATTE["text"]):
    ax.scatter([p[0]], [p[1]], [p[2]], s=s, color=color, depthshade=False,
               edgecolors=LATTE["base"], linewidths=0.9, zorder=9)


def mlabel(ax, tip, text, color=None, offset=(0.10, -0.06, 0.08),
           fontsize=11.5):
    ax.text(tip[0] + offset[0], tip[1] + offset[1], tip[2] + offset[2],
            text, color=color or LATTE["text"], fontsize=fontsize, zorder=10)


# ---------------------------------------------------------------------------
# Figure 1: fundamental forms
# ---------------------------------------------------------------------------
def fig_fundamental_forms():
    """The chapter's anchor picture: surface with parameter grid, tangent
    plane, the tangent/normal vectors, and a curve whose acceleration is
    split into a normal (bending) and tangent (turning) part."""
    fig = plt.figure(figsize=(8.8, 6.8))
    ax = fig.add_subplot(111, projection="3d")
    style3d(ax, (-0.20, 2.12), (-1.24, 1.24), (-0.04, 1.75),
            elev=22, azim=-125, aspect=(1, 1, 0.60))
    shaded_surface(ax, SURF, DOMU, DOMV)
    iso_lines(ax, SURF, DOMU, DOMV)
    rim_lines(ax, SURF, DOMU, DOMV)

    p = np.array([1.0, 0.0])
    P = phi(p)
    eu, ev = ju(p), jv(p)
    n = normal(p)
    wu, wv = eu / np.linalg.norm(eu), ev / np.linalg.norm(ev)

    # tangent plane: translucent patch + dashed outline
    combos = [(-0.55, -0.55), (0.55, -0.55), (0.55, 0.55), (-0.55, 0.55)]
    pts = [P + a * wu + b * wv for a, b in combos]
    ax.add_collection3d(Poly3DCollection([pts], facecolors=LATTE["lavender"],
                                         edgecolors="none", alpha=0.16,
                                         zorder=5))
    ring = pts + [pts[0]]
    ax.plot([q[0] for q in ring], [q[1] for q in ring],
            [q[2] for q in ring], color=LATTE["blue"], lw=1.0,
            ls=(0, (3, 2)), alpha=0.75, zorder=5)

    # non-geodesic curve through P, and its acceleration split
    ts = np.linspace(-0.85, 0.85, 200)
    r = 0.75
    qs = np.column_stack([1.0 + r * (np.cos(ts) - 1.0), r * np.sin(ts)])
    plot_path(ax, qs, color=LATTE["mauve"], linewidth=2.8, zorder=5)

    qdot = np.array([0.0, r])
    qddot = np.array([-r, 0.0])
    gamma2 = (qddot[0] * eu + qddot[1] * ev
              + qdot[0] ** 2 * np.array([0, 0, 2 * C])
              + qdot[1] ** 2 * np.array([0, 0, 2 * C]))
    pn = (gamma2 @ n) * n
    tang = gamma2 - pn

    split = 1.35
    tip2 = fancy_arrow(ax, P, split * gamma2, LATTE["text"], lw=2.2)
    tipn = fancy_arrow(ax, P, split * pn, LATTE["teal"], lw=2.2)
    tipt = fancy_arrow(ax, P, split * tang, LATTE["peach"], lw=2.2)
    for a_, b_ in ((tip2, tipn), (tip2, tipt), (tipn, tipt)):
        ax.plot([a_[0], b_[0]], [a_[1], b_[1]], [a_[2], b_[2]],
                color=LATTE["overlay1"], lw=0.9, ls=(0, (2, 2.2)),
                alpha=0.85, zorder=4)

    fancy_arrow(ax, P, 0.78 * wu, LATTE["blue"], lw=2.2)
    fancy_arrow(ax, P, 0.78 * wv, LATTE["blue"], lw=2.2)
    fancy_arrow(ax, P, 0.92 * n, LATTE["red"], lw=2.2)
    ball(ax, P, s=70)

    # the vector J d of a parameter direction d = (1, −1) — green
    jd = jac(p) @ np.array([1.0, -1.0])
    jd = jd / np.linalg.norm(jd)
    tip_jd = fancy_arrow(ax, P, 0.62 * jd, LATTE["green"], lw=2.2)

    mlabel(ax, P + 0.78 * wu, r"$\Phi_u$", LATTE["blue"],
           offset=(0.10, 0.05, 0.06))
    mlabel(ax, P + 0.78 * wv, r"$\Phi_v$", LATTE["blue"],
           offset=(-0.30, -0.06, 0.02))
    mlabel(ax, P + 0.92 * n, r"$n$", LATTE["red"], offset=(0.08, 0.00, 0.10))
    mlabel(ax, tip_jd, r"$J\,d$", LATTE["green"], offset=(0.06, -0.14, 0.08))
    mlabel(ax, P - 0.55 * wu - 0.55 * wv, r"$T_pS$", LATTE["blue"],
           offset=(-0.16, -0.10, -0.14), fontsize=11)
    mlabel(ax, tip2, r"$\gamma''$", LATTE["text"], offset=(0.05, 0.12, 0.04))
    mlabel(ax, tipn, r"$II(d,d)\,n$", LATTE["teal"],
           offset=(0.16, 0.06, 0.16))
    mlabel(ax, tipt, r"$T\,\gamma''$", LATTE["peach"],
           offset=(-0.06, -0.22, -0.02))
    i_mid = int(np.argmin(np.abs(ts - 0.75)))     # upper tail of the arc
    mlabel(ax, phi(qs[i_mid]), r"$\gamma(s)$", LATTE["mauve"],
           offset=(-0.26, 0.02, -0.02), fontsize=10.5)

    ax.set_title("The tangent space $T_pS$ and the two measurements\n"
                 r"$I=J^{\mathrm{T}}J$ measures inside $T_pS$,"
                 r"  $II=\langle n,\Phi_{ij}\rangle$ up along $n$",
                 fontsize=11, color=LATTE["text"])
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
        plot_path(ax, qs, color=LATTE["overlay0"], linewidth=1.0, alpha=0.9,
                  linestyle="--")
        ax.scatter(*phi(qs[-1]), color=LATTE["overlay1"], s=12,
                   depthshade=False)

    plot_path(ax, paths[-1], color=LATTE["red"], linewidth=3.0,
              label="converged geodesic")
    ax.scatter(*phi(q0), color=LATTE["teal"], s=70, depthshade=False,
               label="source $a$")
    ax.scatter(*phi(target), color=LATTE["maroon"], s=90, marker="*",
               depthshade=False, label="target $b$")
    ax.set_title("Shooting: misses at the initial $L$ (grey),"
                 " converged geodesic (red)",
                 fontsize=10)
    ax.legend(loc="upper left", fontsize=8)
    ax.set_box_aspect((1, 1, 0.62))

    ax2 = fig.add_subplot(122)
    ax2.semilogy(range(len(hist)), hist, "o-", color=LATTE["red"])
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
    latte_fade = LinearSegmentedColormap.from_list(
        "latte_fade", [LATTE["blue"], LATTE["lavender"], LATTE["mauve"],
                       LATTE["red"]])
    colors = latte_fade(np.linspace(0.05, 0.9, len(keys)))
    for color, it in zip(colors, keys):
        plot_path(ax, snaps[it], color=color, linewidth=2.0,
                  label=f"iter {it}" if it in (0, 250, 1000, keys[-1]) else None)
    ax.scatter(*phi(q0), color=LATTE["teal"], s=70, depthshade=False)
    ax.scatter(*phi(q1), color=LATTE["maroon"], s=90, marker="*",
               depthshade=False)
    ax.set_title("Relaxation: a bent path tightens into a geodesic", fontsize=10)
    ax.legend(loc="upper left", fontsize=8)
    ax.set_box_aspect((1, 1, 0.62))

    ax2 = fig.add_subplot(122)
    ax2.semilogy(lengths, color=LATTE["blue"], linewidth=2.0)
    ax2.set_xlabel("iteration")
    ax2.set_ylabel("polyline length")
    ax2.set_title("Length decreases monotonically", fontsize=10)
    ax2.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(HERE, f"relaxation.{ext}"),
                    dpi=170, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Figure 4: normal curvature — κ(d) depends on the direction you walk
# ---------------------------------------------------------------------------
def fig_normal_curvature():
    """Left: the two principal sections of the paraboloid through one point,
    with their normal-bending vectors (lengths proportional to κ).
    Right: κ(θ) = II(d,d)/|Jd|² swept over all parameter directions; the
    extremes of the curve are the principal curvatures."""
    fig = plt.figure(figsize=(12.4, 5.4))
    ax = fig.add_subplot(121, projection="3d")
    style3d(ax, (-0.20, 2.12), (-1.24, 1.24), (-0.04, 1.75),
            elev=22, azim=-125, aspect=(1, 1, 0.60))
    shaded_surface(ax, SURF, DOMU, DOMV)
    iso_lines(ax, SURF, DOMU, DOMV)
    rim_lines(ax, SURF, DOMU, DOMV)

    p = np.array([1.0, 0.0])
    P = phi(p)
    eu, ev = ju(p), jv(p)
    n = normal(p)
    wu = eu / np.linalg.norm(eu)
    wv = ev / np.linalg.norm(ev)

    # I, II and the direction-resolved normal curvature at P
    I_ = metric(p)
    b = n @ np.array([0.0, 0.0, 2 * C])   # n·Φuu = n·Φvv at this symmetric point
    II_ = np.array([[b, 0.0], [0.0, b]])
    k_u = b / I_[0, 0]                    # κ for d ∥ e_u  (along Φ_u)
    k_v = b / I_[1, 1]                    # κ for d ∥ e_v  (along Φ_v)

    # the two principal sections of the surface through P
    ts_u = np.linspace(0.25, 1.72, 160)
    plot_path(ax, np.column_stack([ts_u, np.zeros_like(ts_u)]),
              color=LATTE["blue"], linewidth=1.8)
    ts_v = np.linspace(-0.88, 0.88, 160)
    plot_path(ax, np.column_stack([np.ones_like(ts_v), ts_v]),
              color=LATTE["peach"], linewidth=1.8)

    # normal-bending vectors II(d,d) n — same direction, different lengths;
    # each κ vector is fanned slightly off P along its own section so the
    # three parallel shafts stay legible, while n itself stays anchored at P
    base_u = P + 0.20 * wu
    base_v = P + 0.20 * wv
    tipu = fancy_arrow(ax, base_u, k_u * 2.8 * n, LATTE["blue"], lw=2.2)
    tipv = fancy_arrow(ax, base_v, k_v * 2.8 * n, LATTE["peach"], lw=2.2)
    tipn = fancy_arrow(ax, P, 1.05 * n, LATTE["red"], lw=2.2)
    ball(ax, P, s=70)
    mlabel(ax, tipn, r"$n$", LATTE["red"], offset=(0.18, 0.08, 0.02))
    mlabel(ax, tipu, r"$\kappa_u$", LATTE["blue"],
           offset=(0.16, -0.10, -0.02), fontsize=10.5)
    mlabel(ax, tipv, r"$\kappa_v$", LATTE["peach"],
           offset=(0.02, 0.16, 0.06), fontsize=10.5)

    ax.text2D(0.02, 0.98,
              "vectors: normal bending $II(d,d)\\,n$ for the two sections\n"
              rf"$\kappa_u={k_u:.4f}$ along $\Phi_u$ (blue),"
              rf"  $\kappa_v={k_v:.4f}$ along $\Phi_v$ (peach)",
              transform=ax.transAxes, va="top", fontsize=8.5,
              color=LATTE["subtext1"])
    ax.set_title("One tangent space, two directions,\n"
                 "two different bending rates",
                 fontsize=10.5, color=LATTE["text"])

    # polar sweep of κ(θ): rotate d=(cosθ, sinθ) through all directions
    ax2 = fig.add_subplot(122, projection="polar")
    th = np.linspace(0.0, 2.0 * np.pi, 480)
    kap = np.array([
        (np.array([c, s]) @ II_ @ np.array([c, s]))
        / (np.array([c, s]) @ I_ @ np.array([c, s]))
        for c, s in zip(np.cos(th), np.sin(th))
    ])
    H = np.mean(np.linalg.eigvalsh(np.linalg.solve(I_, II_)))
    ax2.plot(th, kap, color=LATTE["text"], linewidth=2.6,
             label=r"$\kappa(\theta)$", zorder=5)
    ax2.plot(th, np.full_like(th, H), "--", color=LATTE["overlay0"],
             linewidth=1.3, label=rf"mean  $H={H:.4f}$")
    ax2.scatter([0.0, np.pi], [k_u, k_u], s=42, color=LATTE["blue"], zorder=6)
    ax2.scatter([np.pi / 2, 3 * np.pi / 2], [k_v, k_v], s=42,
                color=LATTE["peach"], zorder=6)
    ax2.annotate(rf"$\kappa_u={k_u:.3f}$", xy=(0.0, k_u), xytext=(9, 15),
                 textcoords="offset points", fontsize=9, color=LATTE["blue"])
    ax2.annotate(rf"$\kappa_v={k_v:.3f}$", xy=(np.pi / 2, k_v), xytext=(9, 15),
                 textcoords="offset points", fontsize=9, color=LATTE["peach"])
    ax2.set_xticks([0.0, np.pi / 2, np.pi, 3 * np.pi / 2])
    ax2.set_xticklabels([r"along $\Phi_u$", r"along $\Phi_v$",
                         r"along $-\Phi_u$", r"along $-\Phi_v$"],
                        fontsize=9, color=LATTE["subtext1"])
    ax2.set_ylim(0.30, 0.58)
    ax2.set_yticks([0.35, 0.45, 0.55])
    ax2.set_yticklabels(["0.35", "0.45", "0.55"])
    ax2.tick_params(colors=LATTE["subtext0"], labelsize=7)
    ax2.grid(color=LATTE["surface0"], lw=0.8)
    ax2.spines["polar"].set_color(LATTE["surface1"])
    ax2.set_facecolor(LATTE["base"])
    ax2.set_title(r"$\kappa(\theta)=\dfrac{II(d,d)}{|J\,d|^2}$ swept over"
                  " all directions —\nextremes are the principal curvatures",
                  fontsize=10.5, pad=22, color=LATTE["text"])
    ax2.legend(loc="lower left", bbox_to_anchor=(-0.10, -0.18), fontsize=8,
               frameon=False, labelcolor=LATTE["subtext1"])
    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(HERE, f"normal-curvature.{ext}"),
                    dpi=170, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Figure 5: the shape operator — the normal tips as you slide
# ---------------------------------------------------------------------------
def fig_shape_operator():
    """Slide along the two principal directions through a point and sample
    the unit normal: how fast the arrow tips turn is exactly what −Sξ encodes."""
    fig = plt.figure(figsize=(8.8, 6.8))
    ax = fig.add_subplot(111, projection="3d")
    style3d(ax, (-0.20, 2.12), (-1.24, 1.24), (-0.04, 1.75),
            elev=22, azim=-125, aspect=(1, 1, 0.60))
    shaded_surface(ax, SURF, DOMU, DOMV)
    iso_lines(ax, SURF, DOMU, DOMV)
    rim_lines(ax, SURF, DOMU, DOMV)

    p = np.array([1.0, 0.0])
    P = phi(p)
    eu, ev = ju(p), jv(p)
    n = normal(p)
    I_ = metric(p)

    # unit tangent vectors for the two principal directions at this point,
    # converted to parameter directions by the chapter's own formula
    # ξ = I⁻¹ Jᵀ w   (the §1.2 "normal equations" box)
    w1 = eu / np.linalg.norm(eu)
    w2 = ev / np.linalg.norm(ev)
    xi1 = np.linalg.solve(I_, jac(p).T @ w1)
    xi2 = np.linalg.solve(I_, jac(p).T @ w2)

    fancy_arrow(ax, P, 0.78 * w1, LATTE["blue"], lw=2.2)
    fancy_arrow(ax, P, 0.78 * w2, LATTE["blue"], lw=2.2)
    fancy_arrow(ax, P, 1.00 * n, LATTE["red"], lw=2.2)

    scale = 0.34
    for xi, color in ((xi1, LATTE["red"]), (xi2, LATTE["teal"])):
        ts = np.linspace(-0.9, 0.9, 7)
        ts = ts[np.abs(ts) > 0.18]           # keep the center arrows unique
        plot_path(ax, np.column_stack([p[0] + ts * xi[0],
                                       p[1] + ts * xi[1]]),
                  color=LATTE["blue"], linewidth=1.0, alpha=0.45)
        track = np.linspace(-0.9, 0.9, 120)
        tips = np.array([phi(p + t * xi) + scale * normal(p + t * xi)
                         for t in track])
        ax.plot(tips[:, 0], tips[:, 1], tips[:, 2], color=color,
                linewidth=1.0, alpha=0.6, ls=(0, (2, 2.2)))
        for t in ts:
            Q = phi(p + t * xi)
            fancy_arrow(ax, Q, scale * normal(p + t * xi), color,
                        lw=1.5, alpha=0.88)
    ball(ax, P, s=70)

    mlabel(ax, P + 0.78 * w1, r"$\xi_1$", LATTE["blue"],
           offset=(0.10, 0.05, 0.06))
    mlabel(ax, P + 0.78 * w2, r"$\xi_2$", LATTE["blue"],
           offset=(-0.30, -0.06, 0.02))
    mlabel(ax, P + 0.92 * n, r"$n$", LATTE["red"], offset=(0.08, 0.00, 0.10))
    ax.text2D(0.02, 0.98,
              r"sampled unit normals $n(q)$ while sliding along $\xi_1$ (red)"
              "\n" r"and $\xi_2$ (teal); the dashed tip track is the turn"
              "\n" r"that $-S\xi$ measures",
              transform=ax.transAxes, va="top", fontsize=8.5,
              color=LATTE["subtext1"])
    ax.set_title("The shape operator, one map on $T_pS$:\n"
                 r"slide along $\xi$, the normal tips at rate $-S\xi$",
                 fontsize=10.5, color=LATTE["text"])
    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(HERE, f"shape-operator.{ext}"),
                    dpi=170, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Figure 6: the sign of K on three surfaces
# ---------------------------------------------------------------------------
def _signature_panel(ax, surf, q0, domain, title, zlim, elev, azim, aspect,
                     legend=False):
    """One panel in book style: shaded surface with parameter grid-lines,
    point P, unit normal, the two principal directions, plus a κ/K readout.
    surf is a vectorized (u, v) -> stack of arrays function; derivatives
    are central differences of surf itself (h = 1e-5)."""
    ulim, vlim = domain
    h = 1e-5
    q1, q2 = q0

    def PH(u, v):
        return np.asarray(surf(np.asarray(u, dtype=float),
                               np.asarray(v, dtype=float)))

    P = PH(q1, q2)
    du = (PH(q1 + h, q2) - PH(q1 - h, q2)) / (2 * h)
    dv = (PH(q1, q2 + h) - PH(q1, q2 - h)) / (2 * h)
    duu = (PH(q1 + h, q2) - 2 * P + PH(q1 - h, q2)) / h**2
    dvv = (PH(q1, q2 + h) - 2 * P + PH(q1, q2 - h)) / h**2
    duv = (PH(q1 + h, q2 + h) - PH(q1 + h, q2 - h)
           - PH(q1 - h, q2 + h) + PH(q1 - h, q2 - h)) / (4 * h**2)

    n = np.cross(du, dv)
    n /= np.linalg.norm(n)
    I_ = np.array([[du @ du, du @ dv], [dv @ du, dv @ dv]])
    II_ = np.array([[n @ duu, n @ duv], [n @ duv, n @ dvv]])
    S = np.linalg.solve(I_, II_)
    vals, vecs = np.linalg.eigh(S)
    K = np.linalg.det(II_) / np.linalg.det(I_)

    pad = 0.14
    style3d(ax, (ulim[0] - pad, ulim[1] + pad),
            (vlim[0] - pad, vlim[1] + pad), zlim, elev=elev, azim=azim,
            aspect=aspect)
    shaded_surface(ax, surf, ulim, vlim)
    iso_lines(ax, surf, ulim, vlim, nu=6, nv=5, alpha=0.28)
    rim_lines(ax, surf, ulim, vlim, lw=0.9)

    # normal + principal-direction arrows at P
    fancy_arrow(ax, P, 0.62 * n, LATTE["red"], lw=2.2)
    for idx in (0, 1):
        w = du * vecs[0, idx] + dv * vecs[1, idx]
        fancy_arrow(ax, P, 0.56 * w / np.linalg.norm(w), LATTE["blue"],
                    lw=2.0)
    ball(ax, P, s=70)

    mlabel(ax, P + 0.62 * n, r"$n$", LATTE["red"], offset=(0.06, 0.02, 0.07))
    w_tip1 = du * vecs[0, 0] + dv * vecs[1, 0]
    w_tip2 = du * vecs[0, 1] + dv * vecs[1, 1]
    mlabel(ax, P + 0.56 * w_tip1 / np.linalg.norm(w_tip1), r"$\xi_1$",
           LATTE["blue"], offset=(0.05, -0.02, 0.07), fontsize=10.5)
    mlabel(ax, P + 0.56 * w_tip2 / np.linalg.norm(w_tip2), r"$\xi_2$",
           LATTE["blue"], offset=(0.05, 0.04, 0.07), fontsize=10.5)

    k1, k2 = vals
    if abs(K) < 1e-9:
        stats = rf"$\kappa_1={k1:+.2f}$  $\kappa_2={k2:+.2f}$" "\n" r"$K=0$"
    else:
        stats = rf"$\kappa_1={k1:+.2f}$  $\kappa_2={k2:+.2f}$" f"\n$K={K:+.2f}$"
    ax.text2D(0.03, 0.96, stats, transform=ax.transAxes, va="top",
              fontsize=9, color=LATTE["subtext1"])
    ax.set_title(title, fontsize=10.5, color=LATTE["text"])

    if legend:
        from matplotlib.lines import Line2D
        handles = [
            Line2D([0], [0], color=LATTE["red"], lw=2,
                   label=r"unit normal $n$"),
            Line2D([0], [0], color=LATTE["blue"], lw=2,
                   label=r"principal directions $\xi_1,\xi_2$"),
        ]
        ax.legend(handles=handles, loc="lower left", fontsize=8,
                  frameon=False, labelcolor=LATTE["subtext1"])


def fig_signatures():
    fig = plt.figure(figsize=(13.8, 5.4))
    R = 0.7

    # bowl: z = 0.35 (u² + v²)
    ax1 = fig.add_subplot(131, projection="3d")
    _signature_panel(
        ax1,
        lambda U, V: np.stack([U, V, 0.35 * (U * U + V * V)]),
        (0.62, 0.28), ((-1.05, 1.05), (-1.05, 1.05)),
        "bowl — $K>0$: same-way bending",
        zlim=(-0.02, 0.80), elev=25, azim=-125, aspect=(1, 1, 0.6),
        legend=True)

    # saddle: z = 0.5 (u² − v²)
    ax2 = fig.add_subplot(132, projection="3d")
    _signature_panel(
        ax2,
        lambda U, V: np.stack([U, V, 0.5 * (U * U - V * V)]),
        (0.6, 0.3), ((-1.0, 1.0), (-1.0, 1.0)),
        "saddle — $K<0$: opposite bending",
        zlim=(-0.62, 0.62), elev=25, azim=-125, aspect=(1, 1, 0.6))

    # cylinder of radius R: (θ, y) ↦ (R cosθ, y, −R sinθ); the minus sign
    # makes Φ_θ × Φ_y point outward, so κ₁ = −1/R and κ₂ = 0
    ax3 = fig.add_subplot(133, projection="3d")
    _signature_panel(
        ax3,
        lambda U, V: np.stack([R * np.cos(U), V, -R * np.sin(U)]),
        (2.4, 0.0), ((0.25, 4.35), (-0.5, 0.5)),
        "cylinder — $K=0$: one straight direction",
        zlim=(-0.80, 0.80), elev=16, azim=-75, aspect=(1, 1, 0.7))

    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(HERE, f"curvature-signatures.{ext}"),
                    dpi=170, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Figure 7: shooting a fan of geodesics around a Gaussian hill
#
# The Jacobi metric of a unit-mass particle in a potential V is
#     g_ij = 2 (E - V) delta_ij.
# Its geodesics are the spatial trajectories of the mechanical system, so we
# can shoot them by integrating q'' = -grad V with the initial speed fixed by
# the energy.  A fan of equal-angle shots is the forward map of Chapter 5.
# ---------------------------------------------------------------------------
HILL_A = 1.0
HILL_SIGMA = 0.70
HILL_E = 1.10
HILL_Q0 = np.array([1.50, 1.50])
HILL_DOM = 2.80
HILL_NRAY = 40


def hill_potential(q):
    return HILL_A * np.exp(-np.sum(q * q, axis=-1) / (2 * HILL_SIGMA ** 2))


def hill_grad(q):
    return -hill_potential(q)[..., None] * q / (HILL_SIGMA ** 2)


def _hill_rhs(q, p):
    return p, -hill_grad(q)


def _hill_rk4(q, p, h):
    k1q, k1p = _hill_rhs(q, p)
    k2q, k2p = _hill_rhs(q + 0.5 * h * k1q, p + 0.5 * h * k1p)
    k3q, k3p = _hill_rhs(q + 0.5 * h * k2q, p + 0.5 * h * k2p)
    k4q, k4p = _hill_rhs(q + h * k3q, p + h * k3p)
    return (q + h / 6 * (k1q + 2 * k2q + 2 * k3q + k4q),
            p + h / 6 * (k1p + 2 * k2p + 2 * k3p + k4p))


def shoot_hill_fan(alphas, h=0.0025, tmax=40.0, rmax=HILL_DOM, lmax=10.0):
    """One geodesic shot per launch angle, integrated in the time parameter.

    The spatial curve is the same as the Jacobi geodesic; only its
    parameterisation differs.  Each path is stopped when it leaves the
    plotting square or exceeds a generous arc-length cap.
    """
    alphas = np.asarray(alphas, float)
    n = len(alphas)
    q = np.tile(HILL_Q0, (n, 1))
    speed = np.sqrt(np.maximum(2 * (HILL_E - hill_potential(q)), 0.0))
    p = speed[:, None] * np.stack([np.cos(alphas), np.sin(alphas)], axis=1)
    nsteps = int(tmax / h)
    traj = np.empty((nsteps + 1, n, 2))
    traj[0] = q
    arc = np.zeros(n)
    alive = np.ones(n, bool)
    death = np.full(n, nsteps + 1)
    for step in range(nsteps):
        qn, pn = _hill_rk4(q, p, h)
        qn[~alive] = q[~alive]
        pn[~alive] = p[~alive]
        arc += np.linalg.norm(qn - q, axis=1)
        q, p = qn, pn
        traj[step + 1] = q
        done = alive & ((np.linalg.norm(q, axis=1) > rmax) | (arc > lmax))
        death[done] = step + 1
        alive &= ~done
        if not alive.any():
            break
    return [traj[:death[i], i, :] for i in range(n)]


def fig_jacobi_fan():
    """A fan of shortest paths near a Gaussian potential hill.

    The shots are geodesics of g = 2(E-V) delta.  With E above the summit,
    every shot can cross the hill, but the hill deflects most of them; the
    shot aimed at the summit is the only one that reaches the shadow behind
    the hill.  The figure is the forward map whose Jacobian Chapter 5
    differentiates.
    """
    cmap = LinearSegmentedColormap.from_list(
        "latte_hill", [LATTE["base"], LATTE["crust"], LATTE["surface0"],
                       LATTE["surface1"], LATTE["yellow"], LATTE["peach"]])
    fig, ax = plt.subplots(figsize=(8.6, 7.8))

    g = np.linspace(-HILL_DOM, HILL_DOM, 620)
    X, Y = np.meshgrid(g, g)
    Z = HILL_A * np.exp(-(X * X + Y * Y) / (2 * HILL_SIGMA ** 2))
    ax.imshow(np.ma.masked_less(Z, 0.08), origin="lower",
              extent=[-HILL_DOM, HILL_DOM, -HILL_DOM, HILL_DOM],
              cmap=cmap, vmin=0.08, vmax=HILL_A, alpha=0.90,
              interpolation="bilinear", zorder=0)
    ax.contour(X, Y, Z, levels=[0.08], colors=[LATTE["surface2"]],
               linewidths=0.8, alpha=0.70, zorder=2)
    ax.contour(X, Y, Z, levels=np.linspace(0.12, 0.92, 9),
               colors=LATTE["overlay0"], linewidths=0.45, alpha=0.30,
               zorder=2)

    alphas = 2 * np.pi * np.arange(HILL_NRAY) / HILL_NRAY
    paths = shoot_hill_fan(alphas)
    alpha_c = math.atan2(-HILL_Q0[1], -HILL_Q0[0]) % (2 * np.pi)
    ic = int(np.argmin(np.abs((alphas - alpha_c + np.pi)
                             % (2 * np.pi) - np.pi)))
    segs = [p for j, p in enumerate(paths) if j != ic]
    ax.add_collection(LineCollection(
        segs, colors=LATTE["blue"], linewidths=0.9, alpha=0.55,
        zorder=4, capstyle="round"))

    p_c = paths[ic]
    ax.plot(p_c[:, 0], p_c[:, 1], color=LATTE["mauve"], lw=2.4,
            alpha=0.98, zorder=6, solid_capstyle="round")

    ax.scatter(0, 0, s=30, color=LATTE["text"], edgecolors=LATTE["base"],
               linewidths=0.8, zorder=8)
    ax.annotate(r"$V=A$", (0, 0), xytext=(10, 8),
                textcoords="offset points", fontsize=10,
                color=LATTE["text"], zorder=9)
    ax.scatter(*HILL_Q0, s=80, color=LATTE["red"], edgecolors=LATTE["base"],
               linewidths=0.9, zorder=8)
    ax.annotate(r"$q_0$", HILL_Q0, xytext=(8, 6),
                textcoords="offset points", fontsize=14,
                color=LATTE["red"], zorder=9)

    ax.set_xlim(-HILL_DOM, HILL_DOM)
    ax.set_ylim(-HILL_DOM, HILL_DOM)
    ax.set_aspect("equal")
    ax.set_xticks([-2, 0, 2])
    ax.set_yticks([-2, 0, 2])
    ax.set_xlabel(r"$x$", fontsize=12)
    ax.set_ylabel(r"$y$", fontsize=12)
    ax.tick_params(colors=LATTE["subtext0"], labelsize=8.5)
    for s in ax.spines.values():
        s.set_color(LATTE["surface2"])
        s.set_linewidth(0.8)
    ax.set_title(r"Shooting a fan of shortest paths near a Gaussian hill",
                 fontsize=13.5, color=LATTE["text"], pad=8)

    handles = [Line2D([], [], color=LATTE["blue"], lw=1.2, alpha=0.75,
                      label=r"geodesic shot"),
               Line2D([], [], color=LATTE["mauve"], lw=2.0,
                      label=r"shot through the summit")]
    ax.legend(handles=handles, loc="upper center",
              bbox_to_anchor=(0.5, -0.095), ncol=2, frameon=False,
              fontsize=8.5, handlelength=2.0, columnspacing=1.6,
              labelcolor=LATTE["subtext0"])
    ax.text(0.025, 0.975,
            rf"$A={HILL_A:.2f}$, $\sigma={HILL_SIGMA:.2f}$" "\n"
            rf"$E={HILL_E:.2f}>V_{{\max}}={HILL_A:.2f}$" "\n"
            rf"{HILL_NRAY} shots, equal angles",
            transform=ax.transAxes, va="top", ha="left", fontsize=8.5,
            color=LATTE["subtext0"], zorder=10,
            bbox=dict(boxstyle="round,pad=0.34", fc=LATTE["base"],
                      ec=LATTE["surface2"], alpha=0.9))
    ax.text(0.5, -0.155,
            r"$g_{ij}=2(E-V)\,\delta_{ij}$,   "
            r"$V(x,y)=A\,e^{-(x^2+y^2)/2\sigma^2}$",
            transform=ax.transAxes, ha="center", va="top", fontsize=10.5,
            color=LATTE["subtext0"])

    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(HERE, f"jacobi-fan.{ext}"),
                    dpi=180, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Figure 8: walking the embedded hill — the tangent-plane step
#
# The same Gaussian hill as fig_jacobi_fan, now as an embedded surface
#     Phi(x, y) = (x, y, h(x, y)),   h = A exp(-(x^2 + y^2) / 2 sigma^2),
# walked in ambient R^3 without ever writing the geodesic ODE: straight
# step, project the chord tip back onto the surface, re-tangent, repeat.
# The comparison curve is integrated with RK4 on the vector-route identity
# I (u'', v'') = -g specialised to a height field,
#     I_ij = delta_ij + h_i h_j,
#     g_k  = h_k (h_xx u'^2 + 2 h_xy u' v' + h_yy v'^2),
# the height-field twin of Geodesic.lean's gVector / accelVector.
# ---------------------------------------------------------------------------
HILL_START = np.array([1.40, 0.90])
HILL_LAUNCH_ANGLE = math.radians(215.0)
HILL_WALK_LENGTH = 1.6


def hill_height(x, y):
    """The hill's height, elementwise (scalars or arrays)."""
    return HILL_A * np.exp(-(x * x + y * y) / (2 * HILL_SIGMA ** 2))


def hill_dx(x, y):
    return -x / HILL_SIGMA ** 2 * hill_height(x, y)


def hill_dy(x, y):
    return -y / HILL_SIGMA ** 2 * hill_height(x, y)


def hill_dxx(x, y):
    return (x * x / HILL_SIGMA ** 4 - 1.0 / HILL_SIGMA ** 2) * hill_height(x, y)


def hill_dxy(x, y):
    return x * y / HILL_SIGMA ** 4 * hill_height(x, y)


def hill_dyy(x, y):
    return (y * y / HILL_SIGMA ** 4 - 1.0 / HILL_SIGMA ** 2) * hill_height(x, y)


def HILLSURF(U, V):
    return np.stack([U, V, hill_height(U, V)])


def hill_normal(x, y):
    n = np.array([-hill_dx(x, y), -hill_dy(x, y), 1.0])
    return n / np.linalg.norm(n)


def hill_frame(x, y):
    """Orthonormal tangent frame (e1 along +x, e2 = Φ_y orthogonalized):
    the height-field twin of `Geodesic.tangentFrame`."""
    e1 = np.array([1.0, 0.0, hill_dx(x, y)])
    e1 /= np.linalg.norm(e1)
    raw = np.array([0.0, 1.0, hill_dy(x, y)])
    e2 = raw - (raw @ e1) * e1
    return e1, e2 / np.linalg.norm(e2)


def hill_metric(x, y):
    gx, gy = hill_dx(x, y), hill_dy(x, y)
    return np.array([[1.0 + gx * gx, gx * gy], [gx * gy, 1.0 + gy * gy]])


def hill_gvector(x, y, p):
    """g_k = <Φ_k, Φ_ij> q'^i q'^j for the graph: <Φ_k, Φ_ij> = h_k h_ij,
    since Φ_k = (e_k, h_k) and Φ_ij = (0, 0, h_ij)."""
    d = (hill_dxx(x, y) * p[0] ** 2 + 2.0 * hill_dxy(x, y) * p[0] * p[1]
         + hill_dyy(x, y) * p[1] ** 2)
    return np.array([hill_dx(x, y) * d, hill_dy(x, y) * d])


def hill_accel(x, y, p):
    """(u'', v'') = -I⁻¹ g — the height-field twin of `Geodesic.accelVector`."""
    return -np.linalg.solve(hill_metric(x, y), hill_gvector(x, y, p))


def _hill_geo_rhs(q, p):
    return p, hill_accel(q[0], q[1], p)


def _hill_geo_rk4(q, p, h):
    k1q, k1p = _hill_geo_rhs(q, p)
    k2q, k2p = _hill_geo_rhs(q + 0.5 * h * k1q, p + 0.5 * h * k1p)
    k3q, k3p = _hill_geo_rhs(q + 0.5 * h * k2q, p + 0.5 * h * k2p)
    k4q, k4p = _hill_geo_rhs(q + h * k3q, p + h * k3p)
    return (q + h / 6 * (k1q + 2 * k2q + 2 * k3q + k4q),
            p + h / 6 * (k1p + 2 * k2p + 2 * k3p + k4p))


def hill_launch():
    """Start point, unit tangent direction w and parameter velocity ξ of the
    §5.9 walk: ξ = I⁻¹ Jᵀ w (`Geodesic.tangentVecToParams`'s formula)."""
    x, y = HILL_START
    e1, e2 = hill_frame(x, y)
    w = math.cos(HILL_LAUNCH_ANGLE) * e1 + math.sin(HILL_LAUNCH_ANGLE) * e2
    JT = np.array([[1.0, 0.0, hill_dx(x, y)], [0.0, 1.0, hill_dy(x, y)]])
    xi = np.linalg.solve(hill_metric(x, y), JT @ w)
    return np.array([x, y, hill_height(x, y)]), w, xi


def hill_geodesic(n=400):
    """RK4 on I (u'', v'') = -g for the hill surface, from the §5.9 launch
    data; returns the parameter path (u(s), v(s))."""
    _, _, xi = hill_launch()
    h = HILL_WALK_LENGTH / n
    q, p = HILL_START.copy(), xi
    pts = [q.copy()]
    for _ in range(n):
        q, p = _hill_geo_rk4(q, p, h)
        pts.append(q.copy())
    return np.array(pts)


def hill_walk(hstep, n_geo=400):
    """Walk the embedded hill without any geodesic ODE.

    Each step: move straight p' = p + hstep·w in ambient R^3 (the chord
    leaves the surface), drop the tip vertically back onto the graph, then
    re-tangent — kill the normal part of w in the new tangent plane (the
    `Surface.tangentialPart` split) and renormalize.  Returns the ambient
    polyline, the per-step chord misses, the per-step direction turns, and
    the walker's endpoint error against the RK4 geodesic."""
    P0, w0, _ = hill_launch()
    steps = int(round(HILL_WALK_LENGTH / hstep))
    p, w = P0.copy(), w0.copy()
    pts = [p.copy()]
    misses, turns = [], []
    for _ in range(steps):
        tip = p + hstep * w                       # straight chord
        x, y = tip[0], tip[1]                     # vertical projection
        q = np.array([x, y, hill_height(x, y)])
        misses.append(float(np.linalg.norm(tip - q)))
        n = hill_normal(x, y)
        wt = w - (w @ n) * n                      # tangential part
        w_new = wt / np.linalg.norm(wt)
        turns.append(math.degrees(math.acos(np.clip(w @ w_new, -1.0, 1.0))))
        p, w = q, w_new
        pts.append(q)
    geo = hill_geodesic(n_geo)
    geo_amb = np.column_stack([geo[:, 0], geo[:, 1],
                               hill_height(geo[:, 0], geo[:, 1])])
    err = float(np.linalg.norm(pts[-1] - geo_amb[-1]))
    return np.array(pts), geo_amb, misses, turns, err


def fig_hill_walk():
    """§5.9's picture: the embedded Gaussian hill, the start point with its
    tangent plane, the first straight chord floating its h²-sized gap above
    the surface, and the piecewise-linear walker (peach) against the RK4
    geodesic (red) of the same surface — near-identical on purpose: the
    walk and the ODE meet at the same curve."""
    walk, geo, misses, turns, err = hill_walk(0.1)
    P0, w0, xi0 = hill_launch()
    n0 = hill_normal(*HILL_START)
    e1, e2 = hill_frame(*HILL_START)

    fig = plt.figure(figsize=(9.6, 7.4))
    ax = fig.add_subplot(111, projection="3d")
    ulim, vlim = (-0.42, 1.86), (-0.42, 1.36)
    style3d(ax, (ulim[0] - 0.06, ulim[1] + 0.06),
            (vlim[0] - 0.06, vlim[1] + 0.06), (-0.03, 1.08),
            elev=24, azim=-122, aspect=(1, 1, 0.62))
    shaded_surface(ax, HILLSURF, ulim, vlim, color=LATTE["crust"], alpha=0.96)
    iso_lines(ax, HILLSURF, ulim, vlim, nu=7, nv=6)
    rim_lines(ax, HILLSURF, ulim, vlim)

    # tangent plane at the start: translucent patch + dashed outline
    r = 0.30
    patch = [P0 + a * e1 + b * e2 for a, b in
             [(-r, -r), (r, -r), (r, r), (-r, r)]]
    ax.add_collection3d(Poly3DCollection([patch], facecolors=LATTE["lavender"],
                                         edgecolors="none", alpha=0.18,
                                         zorder=5))
    ring = patch + [patch[0]]
    ax.plot([q[0] for q in ring], [q[1] for q in ring], [q[2] for q in ring],
            color=LATTE["blue"], lw=1.0, ls=(0, (3, 2)), alpha=0.8, zorder=5)

    # the RK4 geodesic (red) and the walker (peach, with chord vertices)
    ax.plot(geo[:, 0], geo[:, 1], geo[:, 2], color=LATTE["red"], lw=3.0,
            zorder=4, solid_capstyle="round")
    ax.plot(walk[:, 0], walk[:, 1], walk[:, 2], color=LATTE["peach"], lw=1.5,
            zorder=5, solid_capstyle="round")
    ax.scatter(walk[:, 0], walk[:, 1], walk[:, 2], s=9, color=LATTE["peach"],
               depthshade=False, zorder=6)

    # the first chord: straight, floating its miss above the surface
    tip1 = P0 + 0.1 * w0
    q1 = np.array([tip1[0], tip1[1], hill_height(tip1[0], tip1[1])])
    ax.plot([P0[0], tip1[0]], [P0[1], tip1[1]], [P0[2], tip1[2]],
            color=LATTE["text"], lw=1.6, ls=(0, (4, 2.2)), zorder=6)
    ax.plot([tip1[0], q1[0]], [tip1[1], q1[1]], [tip1[2], q1[2]],
            color=LATTE["teal"], lw=2.4, zorder=7)
    ball(ax, tip1, s=26, color=LATTE["teal"])

    fancy_arrow(ax, P0, 0.34 * w0, LATTE["blue"], lw=2.0, head=0.06)
    fancy_arrow(ax, P0, 0.80 * n0, LATTE["red"], lw=2.2)
    ball(ax, P0, s=70)
    ball(ax, walk[-1], s=70, color=LATTE["maroon"])

    mlabel(ax, P0 + 0.34 * w0, r"$w$", LATTE["blue"],
           offset=(0.02, 0.10, 0.12))
    mlabel(ax, P0 + 0.80 * n0, r"$n$", LATTE["red"],
           offset=(0.08, 0.02, 0.08))
    mlabel(ax, P0 - r * e1 - r * e2, r"$T_pS$", LATTE["blue"],
           offset=(-0.16, -0.10, -0.10), fontsize=10.5)
    mlabel(ax, walk[-1], "shared\nendpoint", LATTE["maroon"],
           offset=(-0.55, 0.10, -0.08), fontsize=9.5)

    ax.text2D(0.02, 0.97,
              "16 straight chords of $h=0.1$, projected and re-tangented\n"
              rf"first chord misses the surface by {misses[0]:.4f}"
              rf" — turn {turns[0]:.2f}$^\circ$"
              "\n" rf"steepest turn {max(turns):.2f}$^\circ$, total bend"
              f" {sum(turns):.1f}$^\\circ$"
              "\n" rf"walker vs ODE endpoint: {err:.4f} apart"
              r" (first order in $h$)",
              transform=ax.transAxes, va="top", fontsize=8.5,
              color=LATTE["subtext1"])
    handles = [
        Line2D([0], [0], color=LATTE["red"], lw=2.4,
               label="RK4 geodesic of the surface (ODE route)"),
        Line2D([0], [0], color=LATTE["peach"], lw=1.8,
               label="tangent-plane walk, 16 chords"),
        Line2D([0], [0], color=LATTE["text"], lw=1.6, ls=(0, (4, 2.2)),
               label=r"first chord: straight in $\mathbb{R}^3$"),
    ]
    ax.legend(handles=handles, loc="lower left", fontsize=8, frameon=False,
              labelcolor=LATTE["subtext1"])
    ax.set_title("Walking a hill: straight step, project, re-tangent\n"
                 r"every chord is straight in $\mathbb{R}^3$ — the bend"
                 " comes from the tilted tangent planes",
                 fontsize=10.5, color=LATTE["text"])
    fig.tight_layout()

    # side-view inset: the first steps along w — straight chords, curving
    # surface.  The vertical gap at each chord tip is the O(h²) miss of §5.9;
    # plotted to scale (no exaggeration), which is the point: it is tiny but
    # systematic, and it is what the projection + re-tangenting repair.
    ins = fig.add_axes((0.635, 0.53, 0.31, 0.21))
    smax = 0.15
    ss = np.linspace(0.0, smax, 200)
    xy0 = P0[:2] + ss[:, None] * w0[:2]
    chord_z = P0[2] + ss * w0[2]
    surf_z = hill_height(xy0[:, 0], xy0[:, 1])
    ins.plot(ss, surf_z, color=LATTE["overlay0"], lw=2.6, zorder=2,
             label="surface profile")
    ins.plot(ss, chord_z, color=LATTE["text"], lw=1.4, ls=(0, (4, 2.2)),
             zorder=3, label="straight chord")
    i1 = np.searchsorted(ss, 0.1)
    ins.scatter([ss[i1]], [chord_z[i1]], s=26, color=LATTE["teal"], zorder=5)
    ins.scatter([0.0], [P0[2]], s=30, color=LATTE["text"], zorder=5)
    ins.annotate("", xy=(ss[i1], surf_z[i1]), xytext=(ss[i1], chord_z[i1]),
                 arrowprops=dict(arrowstyle="-|>", color=LATTE["teal"],
                                 lw=1.6, mutation_scale=9))
    ins.text(ss[i1] + 0.005, (chord_z[i1] + surf_z[i1]) / 2,
             f"gap {misses[0]:.4f}", color=LATTE["teal"], fontsize=8.5,
             va="center")
    ins.legend(loc="lower right", fontsize=7.5, frameon=False,
               labelcolor=LATTE["subtext0"])
    ins.set_xlim(0.0, smax)
    ins.set_ylim(P0[2] - 0.002, float(surf_z[-1]) + 0.004)
    ins.set_xticks([0.0, 0.05, 0.10, 0.15])
    ins.set_xlabel("distance along $w$", fontsize=8, color=LATTE["subtext0"])
    ins.set_ylabel("$z$", fontsize=8, color=LATTE["subtext0"])
    ins.tick_params(labelsize=7, colors=LATTE["subtext0"])
    ins.set_facecolor(LATTE["mantle"])
    for sp in ins.spines.values():
        sp.set_color(LATTE["surface2"])
        sp.set_linewidth(0.8)
    ins.set_title("the chord floats above the curving surface",
                  fontsize=8.5, color=LATTE["text"], pad=4)

    for ext in ("png", "svg"):
        fig.savefig(os.path.join(HERE, f"hill-walk.{ext}"),
                    dpi=170, bbox_inches="tight")
    plt.close(fig)
    nrm = np.linalg.norm([hill_dx(*HILL_START), hill_dy(*HILL_START), 1.0])
    kap = float((hill_dxx(*HILL_START) * xi0[0] ** 2
                 + 2 * hill_dxy(*HILL_START) * xi0[0] * xi0[1]
                 + hill_dyy(*HILL_START) * xi0[1] ** 2) / nrm)
    print(f"hill-walk: endpoint ({geo[-1, 0]:.6f}, {geo[-1, 1]:.6f}), "
          f"walker error at h=0.1 {err:.6f}, first miss {misses[0]:.4f}, "
          f"first turn {turns[0]:.3f} deg, launch kappa_n {kap:.4f} rad/unit")


# ---------------------------------------------------------------------------
# Figure: tangent space schematic (2-D book style)
# ---------------------------------------------------------------------------
def fancy2(ax, base, vec, color, lw=2.2, head=0.055, zorder=5):
    """A 2-D arrow with a filled triangular head."""
    b = np.asarray(base, float)
    v = np.asarray(vec, float)
    n_ = v / np.linalg.norm(v)
    perp = np.array([-n_[1], n_[0]])
    tip = b + v
    left = tip - head * n_ + 0.6 * head * perp
    right = tip - head * n_ - 0.6 * head * perp
    ax.plot([b[0], tip[0]], [b[1], tip[1]], color=color, lw=lw,
            solid_capstyle="round", zorder=zorder)
    ax.fill([tip[0], left[0], right[0]], [tip[1], left[1], right[1]],
            color=color, zorder=zorder)
    return tip


def ball2(ax, p, s=42, color=LATTE["text"]):
    ax.scatter(*p, s=s, color=color, zorder=6)


def lab2(ax, p, text, color, dx=0.0, dy=0.0, fontsize=11):
    ax.text(p[0] + dx, p[1] + dy, text, color=color, fontsize=fontsize,
            ha="left", va="center", zorder=7)


def fig_tangent_space():
    """Jänich-style schematic: a bent sheet seen from the side, its flat
    tangent spaces attached at three points, and a zoom bubble showing the
    unique split of an ambient vector into tangent + normal parts."""
    fig, ax = plt.subplots(figsize=(9.2, 4.6))

    # the sheet, seen from the side: z = C u² profile (parabolic arc)
    us = np.linspace(-1.9, 1.9, 300)
    zs = C * us * us
    ax.plot(us, zs + 0.55, color=LATTE["subtext1"], lw=3.0,
            solid_capstyle="round", zorder=3)
    ax.fill_between(us, zs + 0.55, zs + 0.53, color=LATTE["lavender"],
                    alpha=0.35, zorder=1)

    # three points with their tangent spaces (flat world = line segment),
    # the tangent vector Φu (blue) and the normal n (red) at each
    probe = [-1.35, 0.0, 1.35]
    for u0 in probe:
        P = np.array([u0, C * u0 * u0 + 0.55])
        slope = 2 * C * u0
        tdir = np.array([1.0, slope]) / np.hypot(1.0, slope)
        ndir = np.array([-slope, 1.0]) / np.hypot(1.0, slope)
        r = 0.75
        ax.plot([P[0] - r * tdir[0], P[0] + r * tdir[0]],
                [P[1] - r * tdir[1], P[1] + r * tdir[1]],
                color=LATTE["blue"], lw=1.4, ls=(0, (4, 2.2)),
                alpha=0.9, zorder=4)
        fancy2(ax, P, 0.55 * tdir, LATTE["blue"])
        fancy2(ax, P, 0.55 * ndir, LATTE["red"])
        ball2(ax, P, s=36)
    lab2(ax, np.array([probe[0], C * probe[0] ** 2 + 0.55]) + np.array([-0.62, 0.28]),
         r"$\Phi_u$", LATTE["blue"], fontsize=11)
    lab2(ax, np.array([probe[0], C * probe[0] ** 2 + 0.55]) + np.array([0.02, 0.30]),
         r"$n$", LATTE["red"], fontsize=11)
    lab2(ax, np.array([probe[1], C + 0.55]) + np.array([0.72, -0.44]),
         r"$T_pS$ — flat", LATTE["blue"], fontsize=11.5)
    lab2(ax, np.array([probe[2], C * probe[2] ** 2 + 0.55]) + np.array([-0.95, 0.26]),
         r"$T_qS$ — just as flat", LATTE["blue"], fontsize=10)

    # zoom bubble at the middle point: one ambient vector w and its split
    zx, zy = 1.35, 2.62
    zw = 2.6
    circle = plt.Circle((zx, zy), 0.92, fill=True,
                        facecolor=LATTE["mantle"], edgecolor=LATTE["surface2"],
                        lw=1.2, zorder=2)
    ax.add_patch(circle)
    # connector from the sheet point to the bubble
    ax.plot([probe[1], zx - 0.30], [C + 0.55, zy - 0.88],
            color=LATTE["overlay1"], lw=1.0, ls=(0, (2, 2)), zorder=1)

    Z = np.array([zx, zy])            # bubble center = the point p
    tdir = np.array([1.0, 0.0])
    ndir = np.array([0.0, 1.0])
    fancy2(ax, Z - 0.80 * tdir, 1.6 * tdir, LATTE["blue"])
    fancy2(ax, Z, 0.66 * ndir, LATTE["red"])
    w = 0.86 * np.array([0.42, 0.91])             # ambient vector, tilted
    wn = (w @ ndir) * ndir                        # (w·n) n  — teal
    wt = w - wn                                   # T w — peach
    fancy2(ax, Z, w, LATTE["text"], lw=1.9, zorder=6)
    fancy2(ax, Z, 0.72 * wn, LATTE["teal"])
    fancy2(ax, Z, 0.72 * wt, LATTE["peach"])
    lab2(ax, Z + w + np.array([-0.02, 0.10]), r"$w$", LATTE["text"], fontsize=11)
    lab2(ax, Z + 0.72 * wn + np.array([0.12, 0.05]), r"$(w\cdot n)\,n$",
         LATTE["teal"], fontsize=10.5)
    lab2(ax, Z + 0.72 * wt + np.array([0.02, -0.16]), r"$T\,w$",
         LATTE["peach"], fontsize=11)
    lab2(ax, Z + np.array([-0.02, 0.74]), r"$n$", LATTE["red"], fontsize=11)

    ax.set_xlim(-2.6, 3.15)
    ax.set_ylim(0.35, 3.7)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title("The sheet is curved; its tangent spaces are flat\n"
                 "and every ambient vector splits uniquely",
                 fontsize=11.5, color=LATTE["text"])
    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(HERE, f"tangent-space.{ext}"),
                    dpi=170, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Figure: pullback metric + biproduct split (CT coda, 2-D book style)
# ---------------------------------------------------------------------------
def _ct_box(ax, x, y, w, h, lines, face, edge):
    ax.add_patch(plt.Rectangle((x, y), w, h, facecolor=face,
                               edgecolor=edge, lw=1.3, zorder=2))
    for k, (txt, col, fs) in enumerate(lines):
        ax.text(x + w / 2, y + h - (k + 1) * h / (len(lines) + 1), txt,
                color=col, fontsize=fs, ha="center", va="center", zorder=3)


def fig_pullback():
    """The coda picture: the metric as the pullback I = Φ*δ of the ambient
    inner product along the morphism J, and the split of an ambient vector as
    the biproduct ℝ³ = T_pS ⊕ ℝn."""
    fig, ax = plt.subplots(figsize=(9.6, 5.8))

    # top diagram: pullback — δ on ℝ³, pulled back to I on T_pS
    _ct_box(ax, 0.35, 2.75, 3.05, 1.5,
            [(r"$T\,\mathbb{R}^2 \to T_pS$", LATTE["text"], 12),
             (r"$d$: a parameter direction", LATTE["subtext0"], 10.5)],
            LATTE["mantle"], LATTE["blue"])
    _ct_box(ax, 5.35, 2.75, 3.05, 1.5,
            [(r"$\mathbb{R}^3$", LATTE["text"], 12),
             (r"$\delta$: the plain dot product", LATTE["subtext0"], 10.5)],
            LATTE["mantle"], LATTE["subtext1"])

    ax.annotate("", xy=(5.25, 3.85), xytext=(3.55, 3.85),
                arrowprops=dict(arrowstyle="-|>", color=LATTE["green"],
                                lw=2.2, mutation_scale=16))
    ax.text(4.4, 4.06, r"$J$", color=LATTE["green"], fontsize=13,
            ha="center")
    ax.text(4.4, 3.52, r"$J\,d = d_x\,\Phi_u + d_y\,\Phi_v$",
            color=LATTE["subtext0"], fontsize=10, ha="center")

    # the pullback: I = Jᵀ J read backwards along J
    ax.annotate("", xy=(3.55, 2.45), xytext=(5.25, 2.45),
                arrowprops=dict(arrowstyle="-|>", color=LATTE["peach"],
                                lw=2.2, mutation_scale=16,
                                linestyle=(0, (5, 2))))
    ax.text(4.4, 2.08, r"$I = \Phi^{*}\delta = J^{\mathrm{T}}J$ — the ruler,"
            " pulled back", color=LATTE["peach"], fontsize=11, ha="center")

    # bottom diagram: the biproduct split of an ambient vector w
    _ct_box(ax, 2.45, 0.85, 2.1, 0.75,
            [(r"$w \in \mathbb{R}^3$", LATTE["text"], 11.5)],
            LATTE["base"], LATTE["subtext1"])
    _ct_box(ax, 0.55, -0.05, 1.75, 0.85,
            [(r"$T_pS$", LATTE["text"], 12.5),
             (r"$T\,w$ (peach)", LATTE["peach"], 10)],
            LATTE["mantle"], LATTE["peach"])
    _ct_box(ax, 4.95, -0.05, 1.75, 0.85,
            [(r"$\mathbb{R}n$", LATTE["text"], 12.5),
             (r"$(w\cdot n)\,n$ (teal)", LATTE["teal"], 10)],
            LATTE["mantle"], LATTE["teal"])

    # projections (solid) and re-injections (dashed, slightly offset)
    ax.annotate("", xy=(1.45, 0.88), xytext=(3.05, 0.80),
                arrowprops=dict(arrowstyle="-|>", color=LATTE["peach"],
                                lw=2.0, mutation_scale=14))
    ax.text(1.95, 1.10, r"$T$", color=LATTE["peach"], fontsize=12)
    ax.annotate("", xy=(5.80, 0.88), xytext=(3.95, 0.80),
                arrowprops=dict(arrowstyle="-|>", color=LATTE["teal"],
                                lw=2.0, mutation_scale=14))
    ax.text(4.85, 1.10, r"$(w\cdot n)$", color=LATTE["teal"], fontsize=11)
    ax.annotate("", xy=(2.85, 0.82), xytext=(1.65, 0.92),
                arrowprops=dict(arrowstyle="-|>", color=LATTE["overlay1"],
                                lw=1.1, mutation_scale=10,
                                linestyle=(0, (2, 2))))
    ax.annotate("", xy=(4.15, 0.82), xytext=(5.60, 0.92),
                arrowprops=dict(arrowstyle="-|>", color=LATTE["overlay1"],
                                lw=1.1, mutation_scale=10,
                                linestyle=(0, (2, 2))))

    ax.text(4.4, -0.72, r"$\mathbb{R}^3 \cong T_pS \oplus \mathbb{R}n$:"
            " every vector splits uniquely, and the two pieces reassemble"
            " it losslessly",
            color=LATTE["subtext0"], fontsize=10.5, ha="center")

    ax.set_xlim(0.1, 8.9)
    ax.set_ylim(-1.0, 4.45)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title("The metric is a pullback; the split is a biproduct\n"
                 "(the category-theory reading of §1.2 and §1.5)",
                 fontsize=11.5, color=LATTE["text"])
    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(HERE, f"pullback.{ext}"),
                    dpi=170, bbox_inches="tight")
    plt.close(fig)


def fig_height_map():
    """Ch. 3 anchor: the paraboloid hill as a topographic height map.

    Left: contours of equal height h = C(u²+v²); three paths leave (1,0)
    with the same initial parameter velocity (1,0) — the constant-radius
    ring rides one contour, the parameter-straight line and the geodesic
    cross every contour head-on, but only the geodesic ends where the ODE
    says.  Right: the geodesic's ruler-speed |J q̇| stays flat (E held at
    0.680) while its parameter speed sinks; the straight line's ruler-cost
    rises — same walk, growing price."""
    cmap = LinearSegmentedColormap.from_list(
        "latte_topo", [LATTE["base"], LATTE["crust"], LATTE["surface0"],
                       LATTE["surface1"], LATTE["yellow"], LATTE["peach"]])
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(12.4, 5.6))

    # ---- left: the height map ------------------------------------------------
    g = np.linspace(0.0, 2.3, 420)
    U, V = np.meshgrid(g, g)
    Z = C * (U * U + V * V)
    cs = ax.contourf(U, V, Z, levels=14, cmap=cmap, alpha=0.95)
    lines = ax.contour(U, V, Z, levels=[0.3, 0.9, 1.5, 2.1, 2.7],
                       colors=LATTE["overlay0"], linewidths=0.8)
    ax.clabel(lines, inline=True, fontsize=8, fmt="%.1f",
              colors=[LATTE["subtext0"]])

    # the three paths of §2.9 from (1,0), initial parameter velocity (1,0)
    t = np.linspace(0.0, 1.0, 200)
    ring = np.column_stack([np.cos(t), np.sin(t)])            # circle r = 1
    straight = np.column_stack([1.0 + t, np.zeros_like(t)])   # straight line
    qs, ps = integrate(np.array([1.0, 0.0]), np.array([1.0, 0.0]), 1.0, n=400)
    ax.plot(ring[:, 0], ring[:, 1], color=LATTE["overlay1"], lw=2.0,
            linestyle=(0, (2, 2)))
    ax.plot(straight[:, 0], straight[:, 1], color=LATTE["overlay1"], lw=2.0,
            linestyle=(0, (5, 2)))
    ax.plot(qs[:, 0], qs[:, 1], color=LATTE["green"], lw=2.6)
    ax.plot([1.0], [0.0], marker="o", ms=7, color=LATTE["surface1"], zorder=5)
    ax.plot([2.0], [0.0], marker="o", ms=6, color=LATTE["overlay1"], zorder=5)
    ax.plot([qs[-1, 0]], [qs[-1, 1]], marker="o", ms=7, color=LATTE["green"],
            zorder=5)

    ax.text(0.28, 0.13, "start $(1,0)$,\nvelocity $(1,0)$",
            color=LATTE["subtext0"], fontsize=9)
    ax.text(2.24, 0.20, "straight\n$u=2.0$", color=LATTE["overlay1"],
            fontsize=9, ha="right")
    ax.text(0.30, 0.62, "ring rides\none contour", color=LATTE["overlay1"],
            fontsize=9, ha="center")
    ax.text(1.66, -0.30, "geodesic: $T\\gamma''=0$, ends at\n"
            "$u=1.879182$ — speed hoarded\nwhere steps cost most",
            color=LATTE["green"], fontsize=9)
    ax.text(0.06, 2.10, "contours: equal height $h = 0.3(u^2{+}v^2)$",
            color=LATTE["subtext0"], fontsize=8.5)
    ax.set_title("The height map is the ruler", fontsize=13,
                 color=LATTE["text"])
    ax.set_xlabel("$u$", color=LATTE["subtext0"])
    ax.set_ylabel("$v$", color=LATTE["subtext0"])
    ax.set_xlim(0.0, 2.3)
    ax.set_ylim(-0.45, 2.3)
    ax.set_aspect("equal")

    # ---- right: what is conserved -------------------------------------------
    ss = np.linspace(0.0, 1.0, 201)
    h = 1.0 / 200
    q, p = np.array([1.0, 0.0]), np.array([1.0, 0.0])
    geo_ruler, geo_param = [], []
    for i, s in enumerate(ss):
        I = metric(q)
        geo_ruler.append(np.sqrt(p @ I @ p))
        geo_param.append(np.linalg.norm(p))
        if i < len(ss) - 1:
            q, p = rk4_step(q, p, h)
    straight_ruler = np.sqrt(1.0 + 4 * C * C * (1.0 + ss) ** 2)

    ax2.plot(ss, geo_ruler, color=LATTE["green"], lw=2.6,
             label=r"geodesic: $|J\,\dot q|$ held — $E=0.680000$")
    ax2.plot(ss, straight_ruler, color=LATTE["overlay1"], lw=2.0,
             linestyle=(0, (5, 2)),
             label="straight line: same walk, rising price")
    ax2.plot(ss, geo_param, color=LATTE["peach"], lw=2.4,
             label=r"geodesic: parameter speed $|\dot q|$ sinks")
    ax2.set_xlabel("parameter time $s$", color=LATTE["subtext0"])
    ax2.set_title("What the geodesic conserves", fontsize=13,
                  color=LATTE["text"])
    ax2.legend(fontsize=9, loc="upper left", framealpha=0.0,
               labelcolor=LATTE["subtext0"])
    ax2.tick_params(colors=LATTE["subtext0"])
    for spine in ax2.spines.values():
        spine.set_color(LATTE["overlay0"])

    fig.text(0.5, 0.015,
             "Left: contours of the hill $h=0.3(u^2{+}v^2)$ — three paths "
             "from the same start and velocity.  Right: the geodesic keeps "
             "its ruler-cost per parameter second constant; the straight\n"
             "line pays more and more as the contours spread apart.",
             ha="center", fontsize=9.5, color=LATTE["subtext0"])
    fig.tight_layout(rect=(0.0, 0.055, 1.0, 1.0))
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(HERE, f"height-map.{ext}"),
                    dpi=170, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    fig_tangent_space()
    fig_fundamental_forms()
    fig_shooting()
    fig_relaxation()
    fig_normal_curvature()
    fig_shape_operator()
    fig_signatures()
    fig_jacobi_fan()
    fig_hill_walk()
    fig_pullback()
    fig_height_map()
    print("wrote tangent-space, fundamental-forms, shooting, relaxation, "
          "normal-curvature, shape-operator, curvature-signatures, "
          "jacobi-fan, hill-walk, pullback, height-map (png + svg)")
