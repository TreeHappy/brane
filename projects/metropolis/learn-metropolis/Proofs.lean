import Std

/-!
# Why Metropolis only needs *differences* of log-density

This file proves, in Lean 4, the claim made in `TUTORIAL.md`:

> the normalizing constant of the posterior cancels out of every
> accept/reject decision, so the sampler never needs to compute it.

## The Bayes story

For data `D` and parameter `θ`, Bayes' theorem says

    posterior(θ | D) = p(D | θ) p(θ) / Z

where the **normalizing constant** (also called the *evidence* or *marginal
likelihood*) is

    Z = ∫ p(D | θ) p(θ) dθ.

`Z` is a single number: it does not depend on `θ`. Writing

    logU(θ) = log p(D | θ) + log p(θ)      -- the unnormalized log-density

we get

    log posterior(θ) = logU(θ) - log Z.

## Metropolis acceptance

A random-walk Metropolis step from `θ` to a proposal `θ'` accepts with
probability

    α(θ → θ') = min(1, posterior(θ') / posterior(θ))
              = min(1, (unnorm(θ') / Z) / (unnorm(θ) / Z))
              = min(1, unnorm(θ') / unnorm(θ)).       -- Z cancels!

In log space the acceptance test is

    Δ = log posterior(θ') - log posterior(θ)
      = (logU(θ') - log Z) - (logU(θ) - log Z)
      = logU(θ') - logU(θ).                          -- log Z cancels!

So the solver can use `logU` directly. `Main.lean` does exactly this:
`model` returns `prior + ll` (that is `logU`) and `step` accepts
iff `log u < Δ`.

## What is proved below

The proofs are done over `Rat` (exact rational numbers), not `Float`.
Core Lean 4 has `Rat` but not `ℝ` (the reals live in Mathlib), and `Float`
is binary floating point, so equalities there are approximate. The algebra
is identical over `ℝ`; `Rat` keeps this file self-contained and the proofs
exact.

1. `diff_shift_rat` — adding the same constant to both values does not
   change their difference.
2. `posterior_diff_rat` — subtracting `log Z` from both log-densities does
   not change their difference. This is the Bayes/normalization step.
3. `accept_prob_shift_rat` — any acceptance probability that is a function
   of `Δ` (e.g. `min 1 (exp Δ)`) is unchanged.
4. `decision_shift_rat` — the Boolean accept/reject decision is unchanged
   for every coin flip.
5. `ratio_rat` — at the density level, the same `Z` cancels in a ratio.

The `#print axioms` commands at the bottom show that these theorems depend
only on Lean's standard kernel axioms (`propext`, `Classical.choice`,
`Quot.sound`), never on `sorry` or `admit`.
-/

namespace Metropolis

/-- **Shift invariance of differences.**

Adding the same constant `c` to both `f x` and `f y` leaves their
difference unchanged. This is the algebraic heart of "the normalizing
constant does not matter". -/
theorem diff_shift_rat (f : Rat → Rat) (c x y : Rat) :
    (f y + c) - (f x + c) = f y - f x := by
  rw [Rat.sub_eq_add_neg, Rat.sub_eq_add_neg, Rat.neg_add]
  have h : c + (-f x + -c) = (c + -c) + -f x := by ac_rfl
  rw [Rat.add_assoc, h, Rat.add_neg_cancel, Rat.zero_add]

/-- **The normalizing constant cancels in log space.**

If the true log posterior is `logU θ - logZ`, then the difference of log
posteriors at two points equals the difference of the *unnormalized*
log-densities. This is the theorem that lets the sampler ignore `Z`. -/
theorem posterior_diff_rat (logU : Rat → Rat) (logZ x y : Rat) :
    (logU y - logZ) - (logU x - logZ) = logU y - logU x := by
  rw [Rat.sub_eq_add_neg, Rat.sub_eq_add_neg, Rat.sub_eq_add_neg,
      Rat.sub_eq_add_neg, Rat.neg_add, Rat.neg_neg]
  have h : logU y + -logZ + (-logU x + logZ)
      = logU y + (-logU x + (-logZ + logZ)) := by ac_rfl
  rw [h, Rat.neg_add_cancel, Rat.add_zero]

/-- **The acceptance probability is unchanged.**

`acceptProb` stands for any function of the log-ratio `Δ`, for example
`fun Δ => min 1 (exp Δ)`. Since its argument is unchanged by a constant
shift, the acceptance probability is unchanged too. (We keep `acceptProb`
abstract because `exp` is not part of core Lean.) -/
theorem accept_prob_shift_rat (f : Rat → Rat) (c x y : Rat)
    (acceptProb : Rat → Rat) :
    acceptProb ((f y + c) - (f x + c)) = acceptProb (f y - f x) := by
  rw [diff_shift_rat]

/-- **The coin flip is unchanged.**

The actual accept/reject decision compares `Δ` with a random threshold
(`log u` in `Solver.lean`). Any Boolean decision procedure `decide` gives
the same answer for the shifted and unshifted log-densities. -/
theorem decision_shift_rat (f : Rat → Rat) (c x y : Rat)
    (decide : Rat → Bool) :
    decide ((f y + c) - (f x + c)) = decide (f y - f x) := by
  rw [diff_shift_rat]

/-- **The normalizing constant cancels in a density ratio.**

At the level of densities (not logs), the posterior ratio is

    (p y / Z) / (p x / Z) = p y / p x,

provided `Z ≠ 0`. This is the non-log version of the same cancellation;
it is the form usually shown in Metropolis–Hastings derivations. -/
theorem ratio_rat (p : Rat → Rat) (x y Z : Rat) (hZ : Z ≠ 0) :
    (p y / Z) / (p x / Z) = p y / p x := by
  simp only [Rat.div_def]
  rw [Rat.inv_mul_rev, Rat.inv_inv]
  rw [← Rat.mul_assoc (p y * Z⁻¹) Z (p x)⁻¹]
  rw [Rat.mul_assoc (p y) Z⁻¹ Z]
  rw [Rat.inv_mul_cancel Z hZ, Rat.mul_one]

end Metropolis

/-! ## Axiom check

`#print axioms` reports the axioms a theorem ultimately relies on. Lean's
kernel has three standard axioms (`propext`, `Classical.choice`,
`Quot.sound`); seeing only those means the proof is complete — there is no
hidden `sorry`. -/

#print axioms Metropolis.diff_shift_rat
#print axioms Metropolis.posterior_diff_rat
#print axioms Metropolis.accept_prob_shift_rat
#print axioms Metropolis.decision_shift_rat
#print axioms Metropolis.ratio_rat
