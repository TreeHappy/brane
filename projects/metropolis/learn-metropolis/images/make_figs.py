#!/usr/bin/env python3
"""
Generates the figures used in TUTORIAL.md section 4 ("Why only differences?").

Run:  python3 make_figs.py        (needs matplotlib)
Each figure is written as .png and .svg next to this script.
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
plt.rcParams.update({
    "font.size": 11,
    "axes.spines.top": False,
    "axes.spines.right": False,
})

BLUE, ORANGE, GREEN, GRAY, RED = "#2563eb", "#ea580c", "#16a34a", "#6b7280", "#dc2626"


def save(fig, name):
    fig.savefig(os.path.join(HERE, name + ".png"), dpi=160, bbox_inches="tight")
    fig.savefig(os.path.join(HERE, name + ".svg"), bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------- fig 1
# Why we work in logs: long products of per-observation probabilities
# collapse towards 0, while the log turns them into plain sums.
def fig1():
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.6))
    n = np.arange(1, 51)

    # pretend each observation contributes a factor of about 0.3
    p = 0.3 ** n

    ax = axes[0]
    ax.plot(n, p, color=BLUE, lw=2)
    ax.set_title("probability space: products")
    ax.set_xlabel("number of observations $n$")
    ax.set_ylabel(r"$0.3^{\,n}$")
    ax.annotate("crashes towards 0\n(precision loss, then underflow)",
                xy=(30, 0.3 ** 30), xytext=(18, 0.10),
                arrowprops=dict(arrowstyle="->", color=GRAY))

    ax = axes[1]
    ax.plot(n, np.log(p), color=GREEN, lw=2)
    ax.set_title("log space: sums")
    ax.set_xlabel("number of observations $n$")
    ax.set_ylabel(r"$\log(0.3^{\,n}) = n \cdot \log 0.3$")
    ax.annotate("stays in an ordinary,\ncomfortable range",
                xy=(30, 30 * np.log(0.3)), xytext=(20, -30),
                arrowprops=dict(arrowstyle="->", color=GRAY))

    fig.suptitle(r"$\log$ turns long products into plain sums:  $\log(ab) = \log a + \log b$",
                 y=1.06, fontsize=12)
    save(fig, "fig1_why_logs")


# ---------------------------------------------------------------- fig 2
# Adding a constant (here: subtracting log Z) slides the curve; the
# difference between any two points stays the same.
def fig2():
    fig, ax = plt.subplots(figsize=(8.8, 4.6))
    mu = np.linspace(0, 6, 400)
    f = -(mu - 3.0) ** 2          # a hump-shaped log-density (toy)
    c = -6.0                      # the constant being subtracted (= -log Z)

    ax.plot(mu, f, color=BLUE, lw=2, label=r"$\log U(\mu)$   (what the code computes)")
    ax.plot(mu, f + c, color=ORANGE, lw=2, ls="--",
            label=r"$\log U(\mu) - \log Z$   (true log posterior)")

    # two points and their difference, shown on BOTH curves
    x1, x2 = 2.0, 4.0
    y1, y2 = -(x1 - 3) ** 2, -(x2 - 3) ** 2     # y2 > y1 (x2 is nearer the peak)

    for (xx, yy, col) in ((x1, y1, BLUE), (x2, y2, BLUE),
                          (x1, y1 + c, ORANGE), (x2, y2 + c, ORANGE)):
        ax.plot([xx], [yy], "o", color=col, ms=6)
        ax.hlines(yy, xx, 7.6, color=col, ls=":", lw=1)

    # difference arrows at the right edge — same length in both columns
    ax.annotate("", xy=(7.0, y2), xytext=(7.0, y1),
                arrowprops=dict(arrowstyle="<->", color=BLUE, lw=2))
    ax.annotate("", xy=(8.1, y2 + c), xytext=(8.1, y1 + c),
                arrowprops=dict(arrowstyle="<->", color=ORANGE, lw=2))
    ax.text(6.9, (y1 + y2) / 2 + 0.35, r"$\Delta$", color=BLUE, fontsize=15, ha="right")
    ax.text(8.0, (y1 + y2) / 2 + 0.35, r"$\Delta$", color=ORANGE, fontsize=15, ha="right")
    ax.text(4.0, -11.4,
            r"same length!    $(f\,y + c) - (f\,x + c) \;=\; f\,y - f\,x$",
            fontsize=11.5, color=GRAY, ha="center")

    ax.set_xlabel(r"$\mu$")
    ax.set_ylim(-13.5, 1.2)
    ax.set_xlim(0, 8.6)
    ax.legend(loc="lower left", frameon=False)
    ax.set_title("Subtracting the normalizing constant just slides the whole curve down —\n"
                 "the difference between any two points is untouched", fontsize=12)
    save(fig, "fig2_shift_differences")


# ---------------------------------------------------------------- fig 3
# A fully worked numeric example: probability space vs log space.
def fig3():
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.9))
    for ax in axes:
        ax.axis("off")

    def panel(ax, title, lines, color):
        ax.set_title(title, fontsize=13, color=color)
        body = "\n".join(lines)
        ax.text(0.5, 0.5, body, ha="center", va="center", fontsize=12.5,
                family="monospace",
                bbox=dict(boxstyle="round,pad=0.9", fc="white", ec=color, lw=1.5))

    panel(axes[0],
          "probability space",
          [r"posterior$(\theta) = $ unnorm$(\theta) / Z$",
           "",
           r"unnorm$(\theta)  = 100$",
           r"unnorm$(\theta') = 10$",
           r"$Z = 1000$",
           "",
           r"post$(\theta)  = 100/1000 = 0.1$",
           r"post$(\theta') = \ 10/1000 = 0.01$",
           "",
           r"ratio $= 0.1 / 0.01 = \mathbf{10}$"],
          BLUE)

    panel(axes[1],
          "log space   (base 10)",
          [r"$\log$ post$(\theta) = \log$ unnorm$(\theta) - \log Z$",
           "",
           r"$\log$ unnorm$(\theta)  = \log 100 = 2$",
           r"$\log$ unnorm$(\theta') = \log 10\  = 1$",
           r"$\log Z = \log 1000 = 3$",
           "",
           r"$\log$ post$(\theta)  = 2 - 3 = -1$",
           r"$\log$ post$(\theta') = 1 - 3 = -2$",
           "",
           r"$\Delta = (-1) - (-2) = \mathbf{1} = \log \mathbf{10}$"],
          GREEN)

    axes[0].annotate("take the log", xy=(1.02, 0.5), xycoords="axes fraction",
                     xytext=(0.0, 0.5), textcoords="axes fraction",
                     arrowprops=dict(arrowstyle="->", color=GRAY, lw=2),
                     ha="center", va="center", fontsize=12, color=GRAY)

    fig.suptitle("Same model, same answer — in logs every division becomes a subtraction",
                 fontsize=13, y=1.02)
    save(fig, "fig3_worked_example")


# ---------------------------------------------------------------- fig 4
# The coin flip: the test  u < alpha  becomes  log u < Delta  on a log scale.
def fig4():
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.0))

    # ---- left: acceptance probability as a function of Delta -------------
    ax = axes[0]
    d = np.linspace(-5, 2.5, 400)
    alpha = np.minimum(1.0, np.exp(d))
    ax.plot(d, alpha, color=BLUE, lw=2)
    ax.fill_between(d, 0, alpha, color=GREEN, alpha=0.10)
    ax.axvline(0, color=GRAY, lw=1, ls=":")
    ax.set_xlabel(r"proposed  $\Delta = \log U(\theta') - \log U(\theta)$")
    ax.set_ylabel(r"accept probability  $\alpha$")
    ax.text(-4.8, 0.70, r"$\Delta \ll 0$:   worse proposal" "\n" r"→ almost surely reject",
            fontsize=10, color=RED)
    ax.text(0.15, 0.28, r"$\Delta \geq 0$:   better or equal" "\n" r"→ always accept  ($\alpha = 1$)",
            fontsize=10, color=GREEN)
    ax.text(0.05, 1.05, r"$\Delta = 0$:  50/50", fontsize=10, color=GRAY, ha="left")
    ax.set_title(r"$\alpha = \min(1,\ e^{\Delta})$", fontsize=12)

    # ---- right: the SAME test on two number lines ------------------------
    ax = axes[1]
    a = 0.70                       # one concrete acceptance probability
    delta = np.log10(a)            # its log (base 10): log10(0.7) ≈ -0.155
    ax.set_xlim(-1.30, 0.30)
    ax.set_ylim(-1.4, 2.3)
    ax.axis("off")
    ax.set_title("same coin flip, log coordinates", fontsize=12)

    def hline(y, x0, x1, label, color):
        ax.annotate("", xy=(x1, y), xytext=(x0, y),
                    arrowprops=dict(arrowstyle="->", color=color))
        ax.text(x1 + 0.03, y, label, fontsize=12, color=color, va="center")

    # top number line: u in [0, 1]
    hline(1.7, 0, 1.0, r"$u$", BLUE)
    ax.add_patch(plt.Rectangle((0, 1.62), a, 0.16, color=GREEN, alpha=0.5))
    ax.add_patch(plt.Rectangle((a, 1.62), 1 - a, 0.16, color=RED, alpha=0.35))
    ax.text(a / 2, 1.94, r"accept:  $u < \alpha$", fontsize=11, color=GREEN, ha="center")
    ax.text((1 + a) / 2, 1.94, "reject", fontsize=11, color=RED, ha="center")
    ax.vlines(a, 1.50, 1.84, color=GRAY, lw=1)
    ax.text(0.0, 1.13, rf"$\alpha = {a}$   (i.e. some coin we must beat)", fontsize=10.5, color=BLUE)

    # bottom number line: log u in (-inf, 0]
    hline(-0.2, -1.2, 0.0, r"$\log_{10} u$", ORANGE)
    ax.add_patch(plt.Rectangle((-1.2, -0.28), -delta + 1.2 - 0.0, 0.16, color=RED, alpha=0.35))
    ax.add_patch(plt.Rectangle((-delta, -0.28), 0 + delta, 0.16, color=GREEN, alpha=0.5))
    ax.text(-1.25, -0.40, r"reject:  $\log u > \Delta$", fontsize=11, color=RED, ha="left")
    ax.text(delta - 0.04/2, -0.40, "accept:  $\\log u < \\Delta$", fontsize=11,
            color=GREEN, ha="right")
    ax.vlines(delta, -0.36, -0.02, color=GRAY, lw=1)
    ax.text(delta, -0.60, rf"$\Delta = \log_{{10}} \alpha = {delta:+.3f}$",
            fontsize=10.5, ha="center", color=GRAY)

    # vertical connector: alpha on top line ↔ Delta on bottom line
    ax.plot([delta, delta], [-0.36, 0.90], color=GRAY, lw=1, ls=":")
    ax.text(-0.60, 0.45, r"log is increasing:" "\n" r"$u < \alpha \iff \log u < \log \alpha = \Delta$",
            fontsize=11, color=GRAY, ha="center")

    fig.canvas.draw()  # ensure layout before save

    save(fig, "fig4_coin_flip")


fig1()
fig2()
fig3()
fig4()
print("done")
