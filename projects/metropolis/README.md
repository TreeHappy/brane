# 🎲 metropolis — probabilistic programming in Lean 4

Metropolis–Hastings MCMC with a generic solver architecture, plus
machine-checked proofs of the core invariants. ✨ First project in the
[Brane](../../README.md) math exploration repo.

## 🗂️ Layout

```
📄 Solver.lean           Model (unnormalized log-density), generic Solver
                         architecture, RNG helpers, runChain, metropolis
📄 Main.lean             A concrete model (μ ~ N(0,5), 50 obs from N(3,1)) + executable
📁 learn-metropolis/     The learning folder — start here: index.md
 ├ 📖 TUTORIAL.md        Ch. 1: Metropolis core, log-space, chain anatomy, proofs
 ├ 🦓 SAMPLERS.md        Ch. 2: the sampler zoo (Gibbs, slice, MALA, HMC, NUTS, …)
 ├ 🧾 Proofs.lean        Lean proof that the normalizing constant cancels
 └ 🖼️ images/            Figures + make_figs.py
```

## 🏗️ Build & run

```bash
lake build
lake exe metropolis                        # 👉 run the sampler (deterministic)
lake env lean learn-metropolis/Proofs.lean # ✅ check the proofs (no `sorry`)
```

Expected (seeded, so exact):

```
true μ = 3.0
posterior mean = ~3.09    posterior std = ~0.14
acceptance rate = ~0.47   num samples = 5000
```

## 🎓 Learning track

Read [`learn-metropolis/index.md`](learn-metropolis/index.md) first — it maps
Chapters 1–5 (🧱 core → 🦓 sampler zoo → 🔧 second solver → 🔬 chain
diagnostics → 🚀 beyond MCMC). Chapter 1 lives in
[`TUTORIAL.md`](learn-metropolis/TUTORIAL.md), Chapter 2 in
[`SAMPLERS.md`](learn-metropolis/SAMPLERS.md).
