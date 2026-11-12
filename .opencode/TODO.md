# Brane working TODO

Cross-session work list. Editor maintains it; mark items as they land.
Position marker at the bottom says where we resume.

## Now — rework of `learn-geodesia/01-02-*` (Ch. 1–2)

- [x] **Prover: vector decls.** `Surface.tangentVec`, `Surface.tangentialPart`,
  `Surface.normalPart` added to `Geodesia/Surface.lean`; `lake build` ✅,
  `mise run test` ✅; verified outputs captured (J(1,0), J(1,1), split of
  w=(0,0,0.6) incl. reconstruction; regression lines unchanged).
- [x] **Chapter rewrite.** New text on disk: chunk 1 (header, dictionary,
  notation bridge, §1.1). Remaining appends at the `<!-- §APPEND -->` marker:
  — completed in chunks 1–6 (commits 04d3147…7bad859); §1.1–§1.8, §2.1–2.9,
  GA & CT coda all landed; all old verified fences preserved byte-identical
  (mechanically checked with `awk`+`comm`).
  - [x] chunk 2 — §1.2 (inner product on the tangent space; classical names
    demoted; all four verified example blocks preserved)
  - [x] chunk 3 — §1.3 (acceleration arrows) + §1.4 (shape operator as
    endomorphism of T_pS)
  - [x] chunk 4 — §1.5 (named split + new checked example), §1.6, §1.7
    (toolbox + new decls), §1.8 (+ new exercise 8), Part 2 header, §2.1–2.2
    (+ embedding Caution box)
  - [x] chunk 5 — §2.3 (executable covariant-derivative note) … §2.9
    (+ closing note extended to mention the new vector decls)
  - [x] chunk 6 — **GA & CT coda** (Cl₃: unit 2-blade B̂, sandwich split,
    wedge test, 1-forms d²=0; CT: Vect, morphism J, pullback I = Φ*δ,
    biproduct split, duality Jᵀ; three-coats mapping table) + pullback
    figure reference + closing "By the way"
- [x] **Intrinsic-vs-embedding aside.** Caution box in §1.1 (Hilbert, flat
  torus, Nash pointer, ADM/brane counterpoint) + Caution box in §2.2
  (vector route needs the embedding, symbol route does not) + pointers
  to Ch. 8.
- [ ] **Extend Ch. 8 §8.2** with the sharp examples it currently lacks:
  Hilbert's theorem (H² has *no* isometric embedding in ℝ³), the flat torus
  (intrinsically flat, needs ℝ⁴), Nash (embedding exists but in high dim
  and is not observable), and the counterpoint: where physics *does* assume
  an embedding (ADM 3+1 split, brane worldvolumes) the extrinsic machinery
  (II, S, H) comes right back.
- [x] **Figures** (`learn-geodesia/images/make_figs.py`) — all done:
  `fig_tangent_space()` (2-D schematic: sheet profile, three tangent spaces,
  zoom bubble with the split) and `fig_pullback()` (metric as pullback,
  split as biproduct) added; `fig_fundamental_forms` retuned (green J d
  vector + T_pS label + tangent-space-first title); titles of
  `fig_normal_curvature` + `fig_shape_operator` retitled; `fig_shooting` +
  `fig_relaxation` recolored to Latte (grey misses → overlay0, red
  converged → red, viridis → latte fade; alt-text semantics preserved);
  `__main__` updated; all 9 figures regenerated (png + svg), visually
  inspected for label collisions.
- [x] **NEXT SLOT — Terminology pass: retire the generic "arrow".** — done:
  every "arrow" in `01-02-*.md` re-evaluated by type and renamed
  (Φu,Φv → tangent vectors; γ″ → acceleration vector; Φij → curvature
  vectors; n → normal vector; w ∈ ℝ³ → ambient vector; basis arrows →
  basis vectors; arrow field → vector field; leftovers → vector); the
  two protected verified artifacts kept byte-identical (`cos θ (arrows)`
  snippet label + Output line, and the snippet comment
  `-- the arrows: J ξ = Φ_x`); fence-internal changes are comments and
  one mermaid label only (verified via fence-diff); awkward mechanical
  results hand-fixed (tautologies "vector math: vectors", "vector world
  … its own vectors", "four-vector table" → "table of four vectors",
  "surface vector γ″" → "acceleration vector γ″", dictionary rows →
  tangent/curvature vector). Also: `make_figs.py` rendered text
  ("arrows: normal bending…" → "vectors: …") + math-referring
  comments/docstrings updated (glyph-drawing "arrow"s kept);
  `fig_normal_curvature` regenerated and inspected; the three
  `/-- -/` docstrings in `Geodesia/Surface.lean` updated (comment-only,
  `lake build` passes); one-line "A note on words" added to the chapter
  header. Final sweep: only protected artifacts + LaTeX
  `\Longrightarrow` false positives remain.
- [x] **Rework Ch. 3** (`03-geodesics-and-lagrangians.md`) — done: full house
  header (epigraph picking up Ch. 2 §2.10 hand-off, what-this-chapter-is,
  reproducing block, dictionary + notation bridge incl. `GeoState`/`accel`/
  `energy`/`speed`/`tangentFrame`/`shoot`, note on words); §3.1–3.7 numbering
  preserved (external cites: `08` §3.4, `GA_VIEW` §3.3 — both still resolve);
  §3.2 gains Signature+Implementation of `Geodesic.accel` + verified eq-check
  ($g_u = 0.360000$, $u'' = -0.264706 = -\Gamma^u_{uu}$ of §2.7, residual
  $I\cdot\text{accel}+g = (0,0)$); §3.4 gains the setup fence (three surfaces
  + `pConserved`) and three verified conservation tables (paraboloid
  $E = 0.680000$/speed $1.166190$ at 5 checkpoints; polar $p_\theta =
  1.350000 = r^2\dot\theta$; sphere $p_\varphi = 0.424844$; sphere drift
  $\theta: 1.0 \to 1.834986$ woven into §3.5); new figure `height-map.png`
  ("height map is the topic": contours + three §2.9 paths | conserved
  ruler-speed vs straight line's rising price) embedded in §3.4; §3.6
  garbled $g_k$ pseudo-formula replaced by pointer to §3.2's box; new
  GA & CT coda (mirror wedge test $\gamma''\wedge n = 0$; Noether pairing /
  momentum of the $S^1$ symmetry; three-coats table; skippable note).
  All 23 retired-"arrow" occurrences retyped (2 remaining are the
  note-on-words meta-mentions). Certification: the chapter's exact embedded
  Lean (setup + eq-check + 3 tables) run as one scratch — all four quoted
  Output blocks byte-exact; render sweep: 0 odd-$ lines, 0 fragile inline
  (incl. two line-spanning inline spans fixed pre-commit). Ride-alongs in
  `02-tangent-parts.md`: coda GA gloss corrected ("carries no normal blade"
  → "the acceleration is purely normal") + ex. 9 inline-| hazard fixed;
  index.md node C mentions the height map.
- [ ] **Verify chapters.** Run every checked snippet via `lake env lean`
  scratch: old outputs must print unchanged (prover already confirmed I, II,
  n, κ, Γ lines), new §1.1/§1.5 outputs as quoted; confirm §1.x/§2.x
  cross-references resolve — now across the split files `01-fundamental-forms.md`
  and `02-tangent-parts.md`. (Ch. 3's embedded snippets certified at rewrite
  time.)
- [x] **index.md** — done during the split: reading-map node B retitled
  ("Turning tangent vectors", + hill walk), all retired-"arrow" terms fixed
  across index.md (incl. the pictures list), dictionary/code-map rows still
  accurate for the two new files (no file-links pointed at the old merged
  name).
- [ ] **House instructions** — add GA & CT coda convention to
  `.opencode/instructions/brane.md` (optional coda per chapter; deep dives
  live in CAPS supplements).
- [ ] **Commit-point report** — summarize diff, verification results, figure
  changes.

## Then — rollout after Ch. 1–2 lands

- [ ] **GA & CT coda for chapters 3–8** (one by one; scribe drafts, prover
  checks claims, librarian gathers GA/CT sources). Ch. 3 first: wedge test
  γ″∧n = 0 meets the geodesic equation.
- [ ] **GA_VIEW.md** — add pointer "previewed in the Ch. 1–2 coda"; keep as
  the deep-dive home. Decide later whether CT material ever needs its own
  `CT_VIEW.md` or lives in codas.
- [ ] **Filename decision** — keep `01-02-fundamental-forms-and-tangent-parts.md`
  (title inside is reworked); renaming would churn links from index/GA_VIEW/
  Ch. 3+. Default: keep.
- [ ] **CI gap** — `geodesia` has no `.github/workflows/lean_action_ci.yml`
  (README claims CI validates every project). Copy metropolis's workflow.

## Future topic ideas (from README + this session)

- [ ] New math topics flagged in README: complex analysis, group theory,
  information theory.
- [ ] metropolis: GA/CT coda (probability vs. measure; sampler as
  Markov-chain functor) — after the coda convention is proven on geodesia.
- [ ] Possibly a small "sources" pass: librarian digests for each topic's
  📚 section (Jänich *Vector Analysis* as a style/structure reference for
  the vector-first track).

## Resume marker

Next up: **verify chapters** (run all checked snippets across the split
files), then house instructions (GA & CT coda convention), 08 §8.2 extension,
GA_VIEW pointer, commit-point report. Split + §2.9/§2.10 + **Ch. 3 rework**
done (Ch. 3's snippets certified byte-exact at rewrite time).
