# 🧮 Brane

> A playground for exploring mathematics — 🎲 probabilistic programming,
> 🤖 machine-checked proofs, 🔬 general math deep dives — one self-contained
> *topic project* at a time.

## ✨ Philosophy

- 🧱 **Topics are independent.** Each `projects/<topic>/` brings its own build
  (`lakefile`), docs, and CI entry; topics never depend on each other.
- 🏗️ **Learn by building.** Every topic ships a `learn-*/` folder: tutorials
  with exercises, worked proofs, and figures — written for someone meeting
  the material for the first time.
- ✅ **Prove what matters.** Lean isn't just here for running code: numerical
  claims get formal statements when the algorithms' correctness hinges on
  them (see the sampler's normalization-cancellation proof).
- 🎯 **Deterministic by default.** Seeded RNGs and fixed datasets, so
  exercises can be verified by exact output match.

## 📚 Projects

| Project | What's explored | Docs entry point |
|---|---|---|
| 🎲 [`metropolis`](projects/metropolis) | Probabilistic programming & MCMC: random-walk Metropolis, a generic `Solver` architecture, log-space math, and Lean proofs of why the normalizing constant cancels | 📖 [`learn-metropolis/index.md`](projects/metropolis/learn-metropolis/index.md) |

*(🚧 next topics* — e.g. complex analysis, group theory, information theory —
*join here as new `projects/<topic>` folders.)*

## 🗺️ Repo layout

```
📁 projects/<topic>/      ← self-contained topic project
 ├ 📄 lakefile.toml       ← Lean build config (Lean 4 toolchain)
 ├ 📄 *.lean              ← the exploration code
 └ 📁 learn-*/            ← the tutorial folder for this topic
📁 scripts/               ← host convenience scripts (start-brane.ps1)
📄 mise.toml              ← pinned tools + OCI dev image + tasks
📄 .github/workflows/     ← CI: builds each Lean project on push
```

## 💪 CI

GitHub Actions validates every Lean project on push via
[`leanprover/lean-action`](https://github.com/leanprover/lean-action) — see
the project table above for build commands you can run locally.

## 🧰 Tools & containers

Assumption: [`mise`](https://mise.jdx.dev) (and, for containers, Podman) are
already installed on the host — nothing here installs them for you.

Tool versions are pinned in the root [`mise.toml`](mise.toml) (Lean 4.35.0 +
lake, Python). The usual round-trip:

```bash
mise run build        # 🏗️ lake build every projects/*
mise run test         # ✅ lean-check every learn-*/Proofs.lean
mise run solve        # 🎲 build + run the metropolis sampler
```

### 🪟 Using the Lean compiler on Windows (Podman)

`mise.toml` also declares an OCI dev image (`[oci]` section, built with
[`mise oci`](https://mise.jdx.dev/dev-tools/mise-oci.html)) containing the
pinned Lean 4 toolchain and Python — ☕ one layer per tool, so a tool-version
bump only rebuilds that tool's layer.

**One-command route (Windows PowerShell):**

```powershell
.\scripts\start-brane.ps1
# ➕ options:
.\scripts\start-brane.ps1 -Rebuild
.\scripts\start-brane.ps1 -Command "lake build && lake exe metropolis"
```

Assumes `mise`, `wsl`, and `podman` are installed on the host; the script
builds the image inside WSL (Linux binaries), packs it into an OCI archive,
`podman load`s it, drops you into a ✨ bash shell with the repo mounted at
`/workspace`, and reports every step.

**Manual route (any Linux, e.g. WSL):**

```bash
mise run image-build             # 📦 writes the OCI layout to ./mise-oci/
mise run image-run -- bash       # 🐚 interactive shell; Podman preferred,
                                 #    repo mounted at /workspace
```

Because `mise oci build` copies host-installed tools into the image, build it
from a **Linux** environment (WSL is fine) — see the
[architecture caveat](https://mise.jdx.dev/dev-tools/mise-oci.html#requirements).

Inside the container:

```bash
cd projects/metropolis && lake build && lake exe metropolis
```

🏷️ The image is tagged `ghcr.io/TreeHappy/brane:latest`; publish it with
`mise run image-push` (after `podman login ghcr.io`).
