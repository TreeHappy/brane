# Brane

A playground for exploring mathematics — one self-contained *topic project* at
a time, each with code, docs, and (where possible) machine-checked proofs.

## Philosophy

- **Topics are independent.** Each `projects/<topic>/` has its own build
  (`lakefile`), docs, and CI entry; topics never depend on each other.
- **Learn by building.** Every topic ships a `learn-*/` folder: tutorials with
  exercises, worked proofs, and figures — written for someone meeting the
  material for the first time.
- **Prove what matters.** Lean isn't just here for running code: numerical
  claims get formal statements when the algorithms' correctness hinges on
  them (see the sampler's normalization-cancellation proof).
- **Deterministic by default.** Seeded RNGs and fixed datasets, so exercises
  can be verified by exact output match.

## Projects

| Project | What's explored | Docs entry point |
|---|---|---|
| [`metropolis`](projects/metropolis) | Probabilistic programming & MCMC: random-walk Metropolis a generic `Solver` architecture, log-space math, and Lean proofs of why the normalizing constant cancels | [`projects/metropolis/learn-metropolis/index.md`](projects/metropolis/learn-metropolis/index.md) |

*(next topics* — e.g. complex analysis, group theory, information theory —
*join here as new `projects/<topic>` folders.)*

## Repo layout

```
projects/<topic>/        ← self-contained topic project
  lakefile.toml          ← Lean build config (Lean 4 toolchain)
  *.lean                 ← the exploration code
  learn-*/               ← the tutorial folder for this topic
.github/workflows/       ← CI: builds each Lean project on push
```

## CI

GitHub Actions validates every Lean project on push via
[`leanprover/lean-action`](https://github.com/leanprover/lean-action) — see
the project table above for build commands you can run locally.
