import Mathlib

/-!
# Why the two geodesic accelerations agree

`Geodesia/Geodesic.lean` computes the parameter-space acceleration of a
geodesic in two ways:

* `accelChristoffel` reads the Christoffel symbols `Γ^k_ij` from the metric
  and evaluates the coordinate formula

  ```
  a^k = −Σ_{i,j} Γ^k_ij p^i p^j.
  ```

* `accelVector` differentiates the ambient curve `Φ(u(s), v(s))` twice and
  requires the tangential part of `Φ''` to vanish.  Writing
  `g_k = ⟨Φ_k, Φ_uu⟩ (u')² + 2⟨Φ_k, Φ_uv⟩ u'v' + ⟨Φ_k, Φ_vv⟩ (v')²`,
  that condition is `I (u'', v'') = −g`, so
  `(u'', v'') = −I⁻¹ g`.

This file proves, over `ℝ`, that the two formulas are the same:

```
accelChristoffel = accelVector.
```

A second half (for chapter 05, §5.9 "Walking a hill") checks the two
algebraic facts behind the discrete walker: the *graph identity* — on a
graph surface, the parameter velocity of an ambient tangent vector is
literally its `(x, y)`-components (`toParams_graphPhi_eq_xy`) — and the
*Gauss split* the re-tangent step runs on: `w = T w + (w·n) n` with
`T w` orthogonal to a unit normal (`gauss_split`, `tangentialPart_orthogonal`,
`normalVec_orthogonal`).

## Why the real numbers, not `Float`

The library itself runs on `Float`, because the demo wants fast numerical
output.  `Float` rounds at every operation, so two algebraically equal
formulas need not be *bitwise* equal, and no exact proof of equality could
exist.  The algebra below is the same algebra the program runs; the demo in
`Main.lean` checks that the `Float` implementations agree to displayed
precision, while this file checks the underlying identity exactly.

Using Mathlib's `EuclideanSpace ℝ (Fin 3)` and `Matrix` means the inner
products, matrix-vector products, and finite-sum manipulations come with
Mathlib's lemmas instead of hand-rolled coordinate expansions.  The metric
is the Gram matrix of the two tangent vectors, exactly as in
`Surface.metricI`.

The `#print axioms` commands at the bottom show that the main theorem uses
only Lean's standard kernel axioms — no `sorry`, no extra assumptions beyond
the named symmetry hypothesis on the second partials.
-/

noncomputable section

open scoped InnerProductSpace
open Matrix

namespace Geodesia
namespace GeodesicProof

/-- Ambient space `ℝ³`, as a Mathlib Euclidean space. -/
abbrev Vec3 := EuclideanSpace ℝ (Fin 3)

/-- First fundamental form: the Gram matrix of the two tangent vectors
`Φ_u`, `Φ_v`.  This is the real-valued counterpart of `Surface.metricI`. -/
def metric (Phi : Fin 2 → Vec3) : Matrix (Fin 2) (Fin 2) ℝ :=
  fun i j => ⟪Phi i, Phi j⟫_ℝ

/-- Derivative of the metric entry `I_ij` in the direction `l`, via the
product rule.  Counterpart of `Surface.dMetricEntry`. -/
def dMetric (Phi : Fin 2 → Vec3) (d2Phi : Fin 2 → Fin 2 → Vec3)
    (l i j : Fin 2) : ℝ :=
  ⟪d2Phi l i, Phi j⟫_ℝ + ⟪Phi i, d2Phi l j⟫_ℝ

/-- The first-kind symbol `⟨Φ_k, Φ_ij⟩`; every Christoffel formula in this
file ultimately reduces to it. -/
def gamma1 (Phi : Fin 2 → Vec3) (d2Phi : Fin 2 → Fin 2 → Vec3)
    (k i j : Fin 2) : ℝ :=
  ⟪Phi k, d2Phi i j⟫_ℝ

/-- Christoffel symbol of the second kind, in the same index convention as
`Surface.christoffel`:

`Γ^k_ij = ½ Σ_l I^{kl} (∂_i I_{jl} + ∂_j I_{il} − ∂_l I_{ij})`. -/
def christoffel (Phi : Fin 2 → Vec3) (d2Phi : Fin 2 → Fin 2 → Vec3)
    (Iinv : Matrix (Fin 2) (Fin 2) ℝ) (k i j : Fin 2) : ℝ :=
  (1 / 2 : ℝ) * ∑ l, Iinv k l *
    (dMetric Phi d2Phi i j l + dMetric Phi d2Phi j i l
      - dMetric Phi d2Phi l i j)

/-- **Lowering an index.**  The derivative of the metric is a sum of two
first-kind symbols.  This is the product rule written with inner products. -/
theorem dMetric_eq_gamma1 (Phi : Fin 2 → Vec3)
    (d2Phi : Fin 2 → Fin 2 → Vec3) (l i j : Fin 2) :
    dMetric Phi d2Phi l i j
      = gamma1 Phi d2Phi j l i + gamma1 Phi d2Phi i l j := by
  simp only [dMetric, gamma1]
  rw [real_inner_comm (d2Phi l i) (Phi j)]

/-- The first-kind symbol is symmetric in its lower indices. -/
theorem gamma1_symm (Phi : Fin 2 → Vec3) (d2Phi : Fin 2 → Fin 2 → Vec3)
    (hsymm : ∀ i j, d2Phi i j = d2Phi j i) (k i j : Fin 2) :
    gamma1 Phi d2Phi k i j = gamma1 Phi d2Phi k j i := by
  simp only [gamma1, hsymm i j]

/-- The cyclic combination in the Christoffel definition collapses to
`2 ⟨Φ_l, Φ_ij⟩`.  This is the algebraic heart of the coordinate formula. -/
theorem cyclic_term (Phi : Fin 2 → Vec3)
    (d2Phi : Fin 2 → Fin 2 → Vec3)
    (hsymm : ∀ i j, d2Phi i j = d2Phi j i) (l i j : Fin 2) :
    dMetric Phi d2Phi i j l + dMetric Phi d2Phi j i l
        - dMetric Phi d2Phi l i j
      = 2 * gamma1 Phi d2Phi l i j := by
  rw [dMetric_eq_gamma1 Phi d2Phi i j l,
      dMetric_eq_gamma1 Phi d2Phi j i l,
      dMetric_eq_gamma1 Phi d2Phi l i j]
  have h1 := gamma1_symm Phi d2Phi hsymm l i j
  have h2 := gamma1_symm Phi d2Phi hsymm j i l
  have h3 := gamma1_symm Phi d2Phi hsymm i j l
  rw [← h1, h2, h3]
  ring

/-- Contracting the metric inverse with the cyclic combination gives the
first-kind symbol.  In coordinates: `Γ^k_ij = Σ_l I^{kl} ⟨Φ_l, Φ_ij⟩`. -/
theorem christoffel_eq_gamma1 (Phi : Fin 2 → Vec3)
    (d2Phi : Fin 2 → Fin 2 → Vec3)
    (Iinv : Matrix (Fin 2) (Fin 2) ℝ)
    (hsymm : ∀ i j, d2Phi i j = d2Phi j i) (k i j : Fin 2) :
    christoffel Phi d2Phi Iinv k i j =
      ∑ l, Iinv k l * gamma1 Phi d2Phi l i j := by
  simp only [christoffel]
  rw [Finset.mul_sum]
  apply Finset.sum_congr rfl
  intro l _
  rw [cyclic_term Phi d2Phi hsymm l i j]
  ring

/-- The vector-calculus acceleration before inverting the metric:

`g_k = Σ_{i,j} ⟨Φ_k, Φ_ij⟩ p^i p^j`. -/
def gvec (Phi : Fin 2 → Vec3) (d2Phi : Fin 2 → Fin 2 → Vec3)
    (p : Fin 2 → ℝ) : Fin 2 → ℝ :=
  fun k => ∑ i, ∑ j, gamma1 Phi d2Phi k i j * p i * p j

/-- The Christoffel route: `a^k = −Σ_{i,j} Γ^k_ij p^i p^j`. -/
def accelChristoffel (Phi : Fin 2 → Vec3)
    (d2Phi : Fin 2 → Fin 2 → Vec3)
    (Iinv : Matrix (Fin 2) (Fin 2) ℝ) (p : Fin 2 → ℝ) : Fin 2 → ℝ :=
  fun k => -∑ i, ∑ j, christoffel Phi d2Phi Iinv k i j * p i * p j

/-- The vector-calculus route: `a = −I⁻¹ g`. -/
def accelVector (Phi : Fin 2 → Vec3) (d2Phi : Fin 2 → Fin 2 → Vec3)
    (Iinv : Matrix (Fin 2) (Fin 2) ℝ) (p : Fin 2 → ℝ) : Fin 2 → ℝ :=
  -(Iinv *ᵥ gvec Phi d2Phi p)

/-- Moving the `l`-sum to the outside in a triple sum. -/
theorem sum_comm3 (F : Fin 2 → Fin 2 → Fin 2 → ℝ) :
    (∑ i, ∑ j, ∑ l, F i j l) = ∑ l, ∑ i, ∑ j, F i j l := by
  calc
    (∑ i, ∑ j, ∑ l, F i j l) = ∑ j, ∑ i, ∑ l, F i j l := by
      rw [Finset.sum_comm]
    _ = ∑ j, ∑ l, ∑ i, F i j l := by
      apply Finset.sum_congr rfl
      intro j _
      rw [Finset.sum_comm]
    _ = ∑ l, ∑ j, ∑ i, F i j l := by
      rw [Finset.sum_comm]
    _ = ∑ l, ∑ i, ∑ j, F i j l := by
      apply Finset.sum_congr rfl
      intro l _
      rw [Finset.sum_comm]

/-- **The two accelerations agree.**  Expanding the Christoffel route and
using `christoffel_eq_gamma1`, both sides become the same double sum. -/
theorem accelChristoffel_eq_accelVector (Phi : Fin 2 → Vec3)
    (d2Phi : Fin 2 → Fin 2 → Vec3)
    (Iinv : Matrix (Fin 2) (Fin 2) ℝ)
    (hsymm : ∀ i j, d2Phi i j = d2Phi j i) (p : Fin 2 → ℝ)
    (k : Fin 2) :
    accelChristoffel Phi d2Phi Iinv p k
      = accelVector Phi d2Phi Iinv p k := by
  simp only [accelChristoffel, accelVector,
    christoffel_eq_gamma1 Phi d2Phi Iinv hsymm]
  congr 1
  rw [Matrix.mulVec_apply_eq_sum]
  simp only [gvec, Finset.sum_mul, Finset.mul_sum]
  rw [sum_comm3]
  apply Finset.sum_congr rfl
  intro l _
  apply Finset.sum_congr rfl
  intro i _
  apply Finset.sum_congr rfl
  intro j _
  ring

/-- **The vector route is exactly the geodesic condition.**  Applying the
metric to `accelVector` returns `−g`.  A geodesic is defined by
`I a = −g`, so this is the equation the coordinate formula solves. -/
theorem metric_mulVec_accelVector (Phi : Fin 2 → Vec3)
    (d2Phi : Fin 2 → Fin 2 → Vec3)
    (Iinv : Matrix (Fin 2) (Fin 2) ℝ)
    (hInv : metric Phi * Iinv = 1) (p : Fin 2 → ℝ) :
    metric Phi *ᵥ accelVector Phi d2Phi Iinv p =
      -gvec Phi d2Phi p := by
  simp only [accelVector, Matrix.mulVec_neg, Matrix.mulVec_mulVec, hInv,
    Matrix.one_mulVec]

/-! ## The ambient acceleration is purely normal

Parameterize a curve on the surface as `s ↦ Φ(u(s), v(s))` with velocity
`p = (u', v')` and parameter acceleration `a = (u'', v'')`.  The chain rule
gives the ambient second derivative

```
Φ'' = Σ_i Φ_i u''_i + Σ_{i,j} Φ_ij u'_i u'_j.
```

The claim behind chapter 05's "a straight walk on a hill bends only off the
surface" is that when `a` comes from the `accelVector` route, `Φ''` has no
tangential component: `⟨Φ_k, Φ''⟩ = 0` for both tangent directions `Φ_k`,
so the ambient acceleration is purely normal.  Note the two sums cancel
against each other exactly because `accelVector` was built as `−I⁻¹ g`:
pairing the first sum with `Φ_k` gives `(I a)_k = −g_k`, pairing the second
gives `g_k`. -/

/-- The ambient second derivative of the parameter curve
`s ↦ Φ(u(s), v(s))` with velocity `p = (u', v')` and parameter acceleration
`a = (u'', v'')`, by the chain rule:

`Φ'' = Σ_i Φ_i • u''_i + Σ_{i,j} Φ_ij • (u'_i u'_j)`.

The first sum is the chain rule through the first partials (the part the
curve steers); the second is the grid-bending part (how the parameter grid
itself curves).  Note that in Lean the scalar comes first in `•`:
`a i • Phi i` is the vector `Φ_i` scaled by `u''_i`. -/
def gamma2 (Phi : Fin 2 → Vec3) (d2Phi : Fin 2 → Fin 2 → Vec3)
    (a p : Fin 2 → ℝ) : Vec3 :=
  (∑ i, a i • Phi i) + (∑ i, ∑ j, (p i * p j) • d2Phi i j)

/-- **The ambient acceleration is purely normal.**  Pairing the ambient
second derivative `Φ''` of the curve with either tangent vector `Φ_k` gives
zero, i.e. the tangential part of `Φ''` vanishes: the curve is straight
inside the surface and bent only off it.  The proof pairs each of the two
sums in `gamma2` with `Φ_k`; the first gives `(I a)_k = −g_k` (by
`metric_mulVec_accelVector`), the second gives `g_k`, and they cancel. -/
theorem tangential_gamma2_eq_zero (Phi : Fin 2 → Vec3)
    (d2Phi : Fin 2 → Fin 2 → Vec3)
    (Iinv : Matrix (Fin 2) (Fin 2) ℝ)
    (hInv : metric Phi * Iinv = 1)
    (p : Fin 2 → ℝ) (k : Fin 2) :
    ⟪Phi k, gamma2 Phi d2Phi (accelVector Phi d2Phi Iinv p) p⟫_ℝ = 0 := by
  -- Unfold `gamma2` and split the inner product over the two sums.
  show ⟪Phi k, (∑ i, accelVector Phi d2Phi Iinv p i • Phi i)
      + ∑ i, ∑ j, (p i * p j) • d2Phi i j⟫_ℝ = 0
  rw [inner_add_right]
  -- Pair the `a`-sum with `Φ_k`: this is exactly `I *ᵥ a` at `k`.
  have h1 : ⟪Phi k, ∑ i, accelVector Phi d2Phi Iinv p i • Phi i⟫_ℝ
      = (metric Phi *ᵥ accelVector Phi d2Phi Iinv p) k := by
    rw [Matrix.mulVec_apply_eq_sum, inner_sum]
    refine Finset.sum_congr rfl fun i _ => ?_
    rw [inner_smul_right]
    -- `metric Phi k i` is definitionally `⟪Phi k, Phi i⟫_ℝ`.
    have hmet : metric Phi k i = ⟪Phi k, Phi i⟫_ℝ := rfl
    rw [hmet]
    ring
  -- Pair the grid-bending sum with `Φ_k`: this is `gvec` at `k`, since
  -- `gamma1 Phi d2Phi k i j = ⟪Phi k, d2Phi i j⟫_ℝ` by definition.
  have h2 : ⟪Phi k, ∑ i, ∑ j, (p i * p j) • d2Phi i j⟫_ℝ
      = gvec Phi d2Phi p k := by
    show ⟪Phi k, ∑ i, ∑ j, (p i * p j) • d2Phi i j⟫_ℝ
      = ∑ i, ∑ j, gamma1 Phi d2Phi k i j * p i * p j
    rw [inner_sum]
    refine Finset.sum_congr rfl fun i _ => ?_
    rw [inner_sum]
    refine Finset.sum_congr rfl fun j _ => ?_
    rw [inner_smul_right]
    -- `gamma1 Phi d2Phi k i j` is definitionally `⟪Phi k, d2Phi i j⟫_ℝ`.
    have hgam : gamma1 Phi d2Phi k i j = ⟪Phi k, d2Phi i j⟫_ℝ := rfl
    rw [hgam]
    ring
  -- The two sums cancel: `I a = −g` is exactly how `accelVector` was built.
  rw [h1, h2, metric_mulVec_accelVector Phi d2Phi Iinv hInv p]
  simp

/-! ## The walker's two facts (chapter 05, §5.9 "Walking a hill")

The discrete walker of §5.9 walks the graph surface `Φ(x, y) = (x, y, h(x, y))`
with three moves: a straight ambient chord, a vertical projection back onto
the surface, and a re-tangenting step `w ↦ w - (w·n) n`.  Two algebraic facts
make this scheme well defined, and both are checked here over `ℝ` — the
`Float` library rounds at every operation, so no exact statement about it
could be proved (same story as in the first half of this file):

1. **Graph identity** (`toParams_graphPhi_eq_xy`).  For a graph surface, the
   parameter velocity of an ambient tangent vector `w` is literally the
   `(x, y)`-components of `w`: the conversion `ξ = I⁻¹ Jᵀ w` (the code's
   `Geodesic.tangentVecToParams`) returns exactly `(w.x, w.y)`.

2. **Tangency preservation** (`gauss_split`, `tangentialPart_orthogonal`,
   `normalVec_orthogonal`, `tangentialPart_normalVec_orthogonal`).  The Gauss
   split behind the re-tangent step is a genuine decomposition
   `w = T w + (w·n) n`, and `T w` is orthogonal to `n` — so the walker's
   updated direction is tangent at the new point.

**Not formalized on purpose.**  That the re-tangent update on `ξ` agrees with
one Euler step of the geodesic ODE `ξ' = ξ − h·Γ(ξ, ξ)` up to `O(h²)` is a
first-order numerical claim, out of Lean's `Float` scope.  The tutorial
verifies it numerically: differences 2.45e-3, 5.88e-4, 1.44e-4 at
`h = 0.1, 0.05, 0.025` — exactly `O(h²)` scaling. -/

/-! ### The graph identity: parameter velocity = the `(x, y)`-components

Everything here is at one point of the surface, so the graph's tangent
vectors are just two ambient vectors: `Φ_u = (1, 0, h_x)` and
`Φ_v = (0, 1, h_y)` — the unit vectors `e₁, e₂` in the first two slots,
plus the slope in the third. -/

/-- The tangent vectors of the graph `Φ(x, y) = (x, y, h(x, y))` at a point,
with `hx = h_x` and `hy = h_y` there: `Φ_u = (1, 0, hx)` and
`Φ_v = (0, 1, hy)`. -/
def graphPhi (hx hy : ℝ) : Fin 2 → Vec3 :=
  fun i => WithLp.toLp 2 (match i with
    | 0 => ![1, 0, hx]
    | 1 => ![0, 1, hy])

@[simp] lemma graphPhi_ofLp_zero (hx hy : ℝ) : (graphPhi hx hy 0).ofLp = ![1, 0, hx] := rfl

@[simp] lemma graphPhi_ofLp_one (hx hy : ℝ) : (graphPhi hx hy 1).ofLp = ![0, 1, hy] := rfl

/-- Parameter-space velocity `ξ = I⁻¹ Jᵀ w` at a point: the abstract
counterpart of `Geodesic.tangentVecToParams`.  As in `Surface.tangentVec`,
`J d = d.x • Φ_u + d.y • Φ_v`, so `Jᵀ w` is the vector of dot products
`⟨Φ_u, w⟩, ⟨Φ_v, w⟩`, and `tangentVecToParams` is exactly `I⁻¹` applied to
it. -/
def toParams (Phi : Fin 2 → Vec3) (Iinv : Matrix (Fin 2) (Fin 2) ℝ) (w : Vec3) :
    Fin 2 → ℝ :=
  Iinv *ᵥ (fun l => ⟪Phi l, w⟫_ℝ)

/-- **The graph identity.**  On a graph surface `Φ(x, y) = (x, y, h(x, y))`,
the parameter velocity of an ambient tangent vector `w` is literally the
`(x, y)`-components of `w`: if `w = a • Φ_u + b • Φ_v` is tangent — which for
the graph forces `a = w.ofLp 0` and `b = w.ofLp 1`, since `Φ_u = (1, 0, h_x)`
and `Φ_v = (0, 1, h_y)` — then the walker's conversion `ξ = I⁻¹ Jᵀ w`
(`Geodesic.tangentVecToParams`) returns exactly `(w.ofLp 0, w.ofLp 1)`, i.e.
`(w.x, w.y)`.  Equivalently: `J ξ = w` with `ξ = (w.x, w.y)`.

The proof is the same one-line algebra the tutorial tells: `Jᵀ w = I ξ` for
`ξ = (a, b)` (the metric is *defined* as the Gram matrix of `Φ_u, Φ_v`), so
`I⁻¹ Jᵀ w = ξ`. -/
theorem toParams_graphPhi_eq_xy (hx hy a b : ℝ)
    (Iinv : Matrix (Fin 2) (Fin 2) ℝ)
    (hInv : metric (graphPhi hx hy) * Iinv = 1)
    (w : Vec3) (hw : w = a • graphPhi hx hy 0 + b • graphPhi hx hy 1) :
    toParams (graphPhi hx hy) Iinv w = ![w.ofLp 0, w.ofLp 1] := by
  -- `I⁻¹` is a two-sided inverse of `I` (both are square matrices over ℝ).
  have hRI : Iinv * metric (graphPhi hx hy) = 1 := mul_eq_one_symm hInv
  -- Chain rule along the graph: `Jᵀ w = I ![a, b]` — pairing `w` with `Φ_l`
  -- gives `a ⟨Φ_l, Φ_u⟩ + b ⟨Φ_l, Φ_v⟩`, the `l`-th entry of `I ![a, b]`.
  have hJT : (fun l => ⟪graphPhi hx hy l, w⟫_ℝ)
      = metric (graphPhi hx hy) *ᵥ ![a, b] := by
    funext l
    rw [hw, inner_add_right, real_inner_smul_right, real_inner_smul_right,
      Matrix.mulVec_apply_eq_sum, Fin.sum_univ_two]
    -- `metric Phi l j` is definitionally `⟪Phi l, Phi j⟫_ℝ`.
    have h0 : metric (graphPhi hx hy) l 0
        = ⟪graphPhi hx hy l, graphPhi hx hy 0⟫_ℝ := rfl
    have h1 : metric (graphPhi hx hy) l 1
        = ⟪graphPhi hx hy l, graphPhi hx hy 1⟫_ℝ := rfl
    rw [h0, h1]
    simp
    ring
  -- `ξ = I⁻¹ Jᵀ w = I⁻¹ (I ![a, b]) = ![a, b]`.
  unfold toParams
  rw [hJT, Matrix.mulVec_mulVec, hRI, Matrix.one_mulVec]
  -- The graph tangent vectors carry `e₁, e₂` in their first two slots, so
  -- the `(x, y)`-components of `w` are `a, b`.
  simp [hw]

/-! ### Tangency preservation: the re-tangent step stays in the tangent plane

The Gauss split the walker's step 3 runs on: for a unit normal `n`,

```
w  =  T w + (w·n) n,      T w = w - (w·n) n,      ⟨T w, n⟩ = 0.
```

`T w` is the code's `Surface.tangentialPart`, `(w·n) n` its
`Surface.normalPart`, and `n` its `Surface.normalAt` — which is the
normalized cross product `Φ_u × Φ_v / |Φ_u × Φ_v|`. -/

/-- The normal part `(w·n) n` of the Gauss split; abstract counterpart of
`Surface.normalPart`. -/
def normalPart (n w : Vec3) : Vec3 := ⟪w, n⟫_ℝ • n

/-- The tangential part `T w = w - (w·n) n` of the Gauss split; abstract
counterpart of `Surface.tangentialPart`. -/
def tangentialPart (n w : Vec3) : Vec3 := w - ⟪w, n⟫_ℝ • n

/-- **The Gauss split is a decomposition**: `w = T w + (w·n) n` — the
re-tangent step keeps exactly one summand and discards the other. -/
theorem gauss_split (n w : Vec3) :
    tangentialPart n w + normalPart n w = w := by
  simp [tangentialPart, normalPart]

/-- **The re-tangent step kills only the normal part.**  For a *unit* normal
`n`, the tangential part is orthogonal to it: `⟨T w, n⟩ = 0`.  Unit-ness is
where the nondegeneracy of the surface patch enters — `Surface.normalAt`
divides by `|Φ_u × Φ_v|`, which is what makes `n` a unit vector. -/
theorem tangentialPart_orthogonal (n : Vec3) (hn : ‖n‖ = 1) (w : Vec3) :
    ⟪tangentialPart n w, n⟫_ℝ = 0 := by
  simp only [tangentialPart, inner_sub_left, real_inner_smul_left,
    real_inner_self_eq_norm_sq, hn]
  ring

/-- The cross product `Φ_u × Φ_v`, as an ambient vector: the unnormalized
part of the normal `Surface.normalAt` builds. -/
def crossVec (Phi : Fin 2 → Vec3) : Vec3 :=
  WithLp.toLp 2 (crossProduct (Phi 0).ofLp (Phi 1).ofLp)

lemma crossVec_ofLp (Phi : Fin 2 → Vec3) :
    (crossVec Phi).ofLp = crossProduct (Phi 0).ofLp (Phi 1).ofLp := rfl

/-- The unit normal `n = Φ_u × Φ_v / |Φ_u × Φ_v|`; abstract counterpart of
`Surface.normalAt`. -/
def normalVec (Phi : Fin 2 → Vec3) : Vec3 :=
  (1 / ‖crossVec Phi‖) • crossVec Phi

/-- The normal is a **unit** vector whenever the patch is nondegenerate
(`crossVec Phi ≠ 0`): this is `normalAt`'s division by `|Φ_u × Φ_v|`
paying off. -/
theorem norm_normalVec (Phi : Fin 2 → Vec3) (hne : crossVec Phi ≠ 0) :
    ‖normalVec Phi‖ = 1 := by
  have hc : 0 < ‖crossVec Phi‖ := norm_pos_iff.mpr hne
  rw [normalVec, norm_smul, Real.norm_eq_abs,
    abs_of_pos (div_pos one_pos hc)]
  exact div_mul_cancel₀ 1 (ne_of_gt hc)

/-- The normal is **orthogonal to the surface**: `⟨Φ_u, n⟩ = ⟨Φ_v, n⟩ = 0`.
Note this needs no nondegeneracy hypothesis — the cross product is
perpendicular to both of its factors. -/
theorem normalVec_orthogonal (Phi : Fin 2 → Vec3) (i : Fin 2) :
    ⟪Phi i, normalVec Phi⟫_ℝ = 0 := by
  -- `normalVec` is a scalar multiple of the cross product, and inner
  -- products are linear in the second slot.
  rw [normalVec, real_inner_smul_right]
  have hzero : ⟪Phi i, crossVec Phi⟫_ℝ = 0 := by
    -- Bridge to the plain dot product on `Fin 3 → ℝ` (`star` is trivial
    -- over `ℝ`), then use that a cross product is perpendicular to both
    -- of its factors.
    rw [EuclideanSpace.inner_eq_star_dotProduct, star_trivial, crossVec_ofLp,
      dotProduct_comm]
    fin_cases i
    · exact dot_self_cross _ _
    · exact dot_cross_self _ _
  rw [hzero, mul_zero]

/-- **The walker's re-tangent step stays tangent** (§5.9, step 3).  With the
surface's own unit normal `n = Φ_u × Φ_v / |Φ_u × Φ_v|` — `Surface.normalAt`
in the code — the tangential part of any ambient vector `w` is orthogonal to
`n`; combined with `normalVec_orthogonal` (`n` is perpendicular to both
tangent directions), `T w` lies in the tangent plane at the new point. -/
theorem tangentialPart_normalVec_orthogonal (Phi : Fin 2 → Vec3)
    (hne : crossVec Phi ≠ 0) (w : Vec3) :
    ⟪tangentialPart (normalVec Phi) w, normalVec Phi⟫_ℝ = 0 :=
  tangentialPart_orthogonal (normalVec Phi) (norm_normalVec Phi hne) w

end GeodesicProof
end Geodesia

-- The theorems depend only on the standard kernel axioms.
#print axioms Geodesia.GeodesicProof.accelChristoffel_eq_accelVector
#print axioms Geodesia.GeodesicProof.metric_mulVec_accelVector
#print axioms Geodesia.GeodesicProof.tangential_gamma2_eq_zero
#print axioms Geodesia.GeodesicProof.toParams_graphPhi_eq_xy
#print axioms Geodesia.GeodesicProof.gauss_split
#print axioms Geodesia.GeodesicProof.tangentialPart_orthogonal
#print axioms Geodesia.GeodesicProof.norm_normalVec
#print axioms Geodesia.GeodesicProof.normalVec_orthogonal
#print axioms Geodesia.GeodesicProof.tangentialPart_normalVec_orthogonal
