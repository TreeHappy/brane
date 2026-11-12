# Brane house conventions

Brane is a playground for exploring mathematics one self-contained *topic
project* at a time: probabilistic programming, machine-checked proofs,
general math deep dives. Everything lives under `projects/<topic>/`, and
every topic ships a `learn-<topic>/` tutorial written for someone meeting
the material for the first time.

## Repo layout

| Path | What it is |
|---|---|
| `projects/<topic>/` | self-contained topic project (own `lakefile.toml`, own CI entry) |
| `projects/<topic>/*.lean` | the exploration code (Lean 4 + Mathlib) |
| `projects/<topic>/learn-<topic>/` | the tutorial folder for this topic |
| `projects/<topic>/.github/workflows/lean_action_ci.yml` | CI for the topic (copy `metropolis`'s for new topics) |
| `scripts/` | host convenience scripts (`start-brane.ps1`) |
| `mise.toml` | pinned tools (Lean 4.35.0, Python) + tasks |

Topics are independent and never depend on each other.

## Topic project anatomy

- `lakefile.toml` and `lean-toolchain` pinned to the same version as
  `mise.toml` (currently Lean 4.35.0).
- Lean code as `Topic/*.lean` modules plus a `Main.lean` demo executable.
- `learn-<topic>/` tutorial folder (below).

## learn-<topic>/ anatomy

- **`index.md`** — entry point: emoji title; plain-language intro (explain
  *why* the formulas look the way they do); the central object and one
  boxed key equation; 🧭 reading map (mermaid flowchart) plus suggested
  paths by reader goal ("I never understood X → read Ch. a, b"); 📚 sources
  (primary links + license + "these notes are an independent explanation");
  🗺️ code map table mapping files to concepts; 🔤 notation section
  including a dictionary table mapping older textbook names → names used
  here; "the N pictures to keep in mind".
- **Numbered chapters `01-*.md`** — one idea per chapter: intuition →
  picture → formal statement → where it lives in the code → "try it"
  prompt.
- **Side supplements in CAPS filenames** (`GA_VIEW.md`, `SAMPLERS.md`,
  `TUTORIAL.md`) — alternate viewpoints, deliberately not numbered
  chapters.
- **`Proofs.lean`** — machine-checked proofs for claims the algorithms'
  correctness hinges on (e.g. why the normalizing constant cancels in
  Metropolis).
- **`images/make_figs.py`** — matplotlib script that regenerates every
  figure as both `.svg` and `.png`; never paste stale images by hand.

## Prose style

- Everyday language first, formalism second. No hidden vocabulary: define
  every symbol before use and keep the notation table current.
- Why-first: say what a formula *measures*, not just what it says.
  Blockquote callouts for asides like "why the literature leads with X".
- ```math fenced blocks for display math, inline $…$ otherwise. GitHub's
  renderer eats a backslash in inline `$…$`, so never put `\\` row
  separators, `\begin{…}` matrices or literal `|` inside inline math —
  move them into a ```math block (or reword: pairs, `\mathrm{diag}`,
  component lists). Official reference: docs.github.com "Writing
  mathematical expressions" (MathJax; the `` $`…`$ `` backtick form for
  markdown-overlapping characters).
- Tables for dictionaries and code maps; emoji section headers; mermaid
  for algorithm/reading diagrams.
- Deterministic examples: seeded RNG, fixed datasets, exact output match.
  Quote *real* program output, never invented numbers.
- Reference code with file links (and line ranges where helpful); snippets
  come from the actual files, not retyped from memory.

## Verification

```bash
mise run build                     # lake build every projects/*
mise run test                      # lean-check every learn-*/Proofs.lean
cd projects/<topic> && lake build && lake exe <topic>   # run a topic demo
cd projects/<topic>/learn-<topic>/images && python3 make_figs.py   # figures
```

## Lean code style

- Derivatives are data, not finite differences (see `Diff`/`Smooth2`).
- Thread `StdGen` explicitly; no global random state.
- Prefer Mathlib over hand-rolled; state theorems exactly, prove only the
  generality the demo needs.
- Fix compile errors at the root cause, not by weakening the statement.

## Housekeeping

- A new topic adds: the project folder, its CI workflow (copy
  `metropolis`'s), and a row in the `README.md` project table.
- Keep chapters self-contained; update the reading map, notation table and
  code map whenever chapters change.
