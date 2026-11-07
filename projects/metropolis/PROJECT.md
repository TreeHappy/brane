# metropolis

A small Lean 4 project that makes Metropolis–Hastings MCMC learnable: the simplest
sampler, a generic solver architecture, and machine-checked proofs of why it works.

## Layout

```
Solver.lean           Model (unnormalized log-density), the generic Solver
                      architecture, randFloat/randNormal, runChain, metropolis
Main.lean             A concrete model (μ ~ N(0,5), 50 obs from N(3,1)) + executable
learn-metropolis/     The learning folder:
├── TUTORIAL.md       Hands-on walkthrough of Metropolis in Lean (start here)
├── Proofs.lean       Lean proof that the normalizing constant cancels
└── images/           Figures (fig1 why logs, fig2 shift differences,
                      fig3 worked example) + make_figs.py
```

## Build & run

```bash
lake build
lake exe metropolis
```

Output (deterministic — the RNG is seeded):

```
true μ = 3.0
posterior mean = ~3.09
posterior std  = ~0.14
acceptance rate = ~0.49
num samples = 5000
```

Check the proofs:

```bash
lake env lean learn-metropolis/Proofs.lean
```

Every theorem should depend only on Lean's standard kernel axioms
(`propext`, `Classical.choice`, `Quot.sound`) — no `sorry`.

## Key ideas

- A model is any unnormalized log-density: `abbrev Model := Array Float → Float`.
- Metropolis only needs *differences* of log-densities; the normalizing constant
  cancels — see `learn-metropolis/Proofs.lean`.
- A solver is anything with a `step : Model → Array Float → StdGen → StepResult`;
  `runChain` is solver-agnostic, so Gibbs/HMC can be added without touching the loop.

## Learning path

Read [`learn-metropolis/TUTORIAL.md`](learn-metropolis/TUTORIAL.md), then try its
exercises (tune σ, add a second parameter, thinning, swap the model, write a new solver).
Next steps on the roadmap: HMC, multiple chains + R-hat, reparameterization.
