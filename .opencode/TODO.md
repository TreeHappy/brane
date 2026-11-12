# Brane working TODO

Cross-session work list. Editor maintains it; mark items as they land.
Position marker at the bottom says where we resume.

## Now — rework of `learn-geodesia/01-02-*` (Ch. 1–2)

- [x] **Prover: vector decls.** `Surface.tangentVec`, `Surface.tangentialPart`,
  `Surface.normalPart` added to `Geodesia/Surface.lean`; `lake build` ✅,
  `mise run test` ✅; verified outputs captured (J(1,0), J(1,1), split of
  w=(0,0,0.6) incl. reconstruction; regression lines unchanged).
- [ ] **Chapter rewrite.** New text on disk: chunk 1 (header, dictionary,
  notation bridge, §1.1). Remaining appends at the `<!-- §APPEND -->` marker:
  - [x] chunk 2 — §1.2 (inner product on the tangent space; classical names
    demoted; all four verified example blocks preserved)
  - [ ] chunk 3 — §1.3 (acceleration arrows) + §1.4 (shape operator as
    endomorphism of T_pS)
  - [ ] chunk 4 — §1.5 (named split + new checked example), §1.6, §1.7
    (toolbox + new decls), §1.8 (+ new exercise 8), Part 2 header, §2.1–2.2
  - [ ] chunk 5 — §2.3 (executable covariant-derivative note) … §2.9
  - [ ] chunk 6 — **GA & CT coda** (Cl₃: 2-blade B, sandwich split, wedge
    test, 1-forms as dual pairs; CT: Vect, morphism J, pullback I = Φ*δ,
    biproduct split, duality Jᵀ, d²=0; mapping table) + closing note
- [ ] **Intrinsic-vs-embedding aside.** Fold the "why not just assume the
  extra dimension?" exchange in: compact Caution box in §1.1/§2.2 + pointer
  to Ch. 8.
- [ ] **Extend Ch. 8 §8.2** with the sharp examples it currently lacks:
  Hilbert's theorem (H² has *no* isometric embedding in ℝ³), the flat torus
  (intrinsically flat, needs ℝ⁴), Nash (embedding exists but in high dim
  and is not observable), and the counterpoint: where physics *does* assume
  an embedding (ADM 3+1 split, brane worldvolumes) the extrinsic machinery
  (II, S, H) comes right back.
- [ ] **Figures** (`learn-geodesia/images/make_figs.py`):
  - [ ] add `fig_tangent_space()` — 2-D Jänich-style schematic (sheet profile,
    three tangent spaces, zoom bubble with the split)
  - [ ] add `fig_pullback()` — coda picture (metric as pullback I = Φ*δ;
    split as biproduct ℝ³ = T_pS ⊕ ℝn)
  - [ ] retune `fig_fundamental_forms` (add Jd arrow + T_pS label + title),
    titles of `fig_normal_curvature` + `fig_shape_operator`
  - [ ] recolor `fig_shooting` + `fig_relaxation` to Latte (semantic colors
    preserved: grey misses → overlay0, red converged → red; alt texts in
    05/06 stay valid)
  - [ ] run `python3 make_figs.py`, regenerate png+svg for all
- [ ] **Verify chapter.** Run every checked snippet via `lake env lean`
  scratch: old outputs must print unchanged (prover already confirmed I, II,
  n, κ, Γ lines), new §1.1/§1.5 outputs as quoted; confirm §1.x/§2.x
  cross-references resolve.
- [ ] **index.md** — reading-map node A ("The tangent space …"), suggested
  path 1 reword, GA_VIEW pointer ("preview in the coda").
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

Next up: chapter chunk 3 (§1.3 + §1.4), then continue down the Now list in order.
