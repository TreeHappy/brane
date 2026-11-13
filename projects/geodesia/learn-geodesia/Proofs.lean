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

end GeodesicProof
end Geodesia

-- The main theorem depends only on the standard kernel axioms.
#print axioms Geodesia.GeodesicProof.accelChristoffel_eq_accelVector
#print axioms Geodesia.GeodesicProof.metric_mulVec_accelVector
#print axioms Geodesia.GeodesicProof.tangential_gamma2_eq_zero
