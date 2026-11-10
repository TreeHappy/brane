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
    plane, the tangent/normal arrows, and a curve whose acceleration is
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

    mlabel(ax, P + 0.78 * wu, r"$\Phi_u$", LATTE["blue"],
           offset=(0.10, 0.05, 0.06))
    mlabel(ax, P + 0.78 * wv, r"$\Phi_v$", LATTE["blue"],
           offset=(-0.30, -0.06, 0.02))
    mlabel(ax, P + 0.92 * n, r"$n$", LATTE["red"], offset=(0.08, 0.00, 0.10))
    mlabel(ax, tip2, r"$\gamma''$", LATTE["text"], offset=(0.05, 0.12, 0.04))
    mlabel(ax, tipn, r"$II(d,d)\,n$", LATTE["teal"],
           offset=(0.16, 0.06, 0.16))
    mlabel(ax, tipt, r"$T\,\gamma''$", LATTE["peach"],
           offset=(-0.06, -0.22, -0.02))
    i_mid = int(np.argmin(np.abs(ts - 0.75)))     # upper tail of the arc
    mlabel(ax, phi(qs[i_mid]), r"$\gamma(s)$", LATTE["mauve"],
           offset=(-0.26, 0.02, -0.02), fontsize=10.5)

    ax.set_title("The two fundamental forms at one point\n"
                 r"$I=J^{\mathrm{T}}J$ measures in the tangent plane,"
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


# ---------------------------------------------------------------------------
# Figure 4: normal curvature — κ(d) depends on the direction you walk
# ---------------------------------------------------------------------------
def fig_normal_curvature():
    """Left: the two principal sections of the paraboloid through one point,
    with their normal-bending arrows (lengths proportional to κ).
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

    # normal-bending arrows II(d,d) n — same direction, different lengths;
    # each κ arrow is fanned slightly off P along its own section so the
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
              "arrows: normal bending $II(d,d)\\,n$ for the two sections\n"
              rf"$\kappa_u={k_u:.4f}$ along $\Phi_u$ (blue),"
              rf"  $\kappa_v={k_v:.4f}$ along $\Phi_v$ (peach)",
              transform=ax.transAxes, va="top", fontsize=8.5,
              color=LATTE["subtext1"])
    ax.set_title("One point, two directions,\ntwo different bending rates",
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
    the unit normal: how fast the arrow tips is exactly what −Sξ encodes."""
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

    # unit tangent arrows for the two principal directions at this point,
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
    ax.set_title("Shape operator $S = I^{-1}II$: slide along $\\xi$,\n"
                 r"and the normal arrow tips at rate $-S\xi$",
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


if __name__ == "__main__":
    fig_fundamental_forms()
    fig_shooting()
    fig_relaxation()
    fig_normal_curvature()
    fig_shape_operator()
    fig_signatures()
    print("wrote fundamental-forms, shooting, relaxation, normal-curvature, "
          "shape-operator, curvature-signatures (png + svg)")
