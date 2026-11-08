/-!
# Geodesia — the metric, Christoffel symbols, curvature

For a surface patch `Φ : (x, y) ↦ Φ(x, y) ∈ ℝ³`:

## First fundamental form (the metric on the surface)

```
g = [[ E  F ]      E = ⟨Φx, Φx⟩,  F = ⟨Φx, Φy⟩,  G = ⟨Φy, Φy⟩
     [ F  G ]]
```
This is *the geometry*: every intrinsic measurement (lengths, angles, areas,
even Gaussian curvature) is computable from `g` alone.

## Christoffel symbols

```
Γ^k_ij = ½ Σ_m g^{km} (∂_i g_{jm} + ∂_j g_{im} - ∂_m g_{ij})
```
with `g^{-1}` the inverse metric. They correct naive differentiation: the
covariant derivative of a surface vector field is
`(∇_i V)^k = ∂_i V^k + Σ_j Γ^k_ij V^j`; the extra term comes from the
*turning* of the coordinate grid. Geodesics are paths with zero covariant
acceleration (see `Geodesic.lean`).

## Second fundamental form and curvature

```
II = [[ e  f ]      e = ⟨n, Φxx⟩, f = ⟨n, Φxy⟩ = ⟨n, Φyx⟩, g₂ = ⟨n, Φyy⟩
     [ f  g₂ ]]     n = unit normal
```
Gaussian curvature `K = det II / det I`; mean curvature
`H = (e G − 2 f F + g₂ E) / 2(E G − F²)`.
-/

import Geodesia.Surface

namespace Geodesia

/-- 2×2 symmetric Float matrix `[[a, b], [b, c]]`. -/
structure Mat2 where
  a : Float
  b : Float
  c : Float

namespace Mat2

def det (M : Mat2) : Float := M.a * M.c - M.b * M.b

/-- Inverse via Cramer's rule (positive-definite metrics have `det > 0`). -/
def inv (M : Mat2) : Mat2 :=
  let d := M.det
  ⟨M.c / d, -M.b / d, M.a / d⟩

/-- `(k, m)` entry of the inverse matrix. -/
def invEntry (M : Mat2) : Fin 2 → Fin 2 → Float
  | 0, 0 => M.a
  | 0, 1 => M.b
  | 1, 0 => M.b
  | 1, 1 => M.c
  -- all 4 pairs of `Fin 2 × Fin 2` are covered

end Mat2

namespace Surface

/-! ## Derivatives of the patch Φ. -/

/-- `∂Φ/∂i` for `i = 0` (the x-direction) or `i = 1` (the y-direction). -/
def dPhi (s : Surface) (i : Fin 2) (x y : Float) : Vec3 :=
  match i with
  | 0 => ⟨dd s.sx x y, dd s.sy x y, dd s.sz x y⟩
  | 1 => ⟨dd (s.sx x) y, dd (s.sy x) y, dd (s.sz x) y⟩

/-- `(i, j)` entry of the first fundamental form: `⟨∂_iΦ, ∂_jΦ⟩`. -/
def metricEntry (s : Surface) (x y : Float) : Fin 2 → Fin 2 → Float :=
  fun i j => (dPhi s i x y).dot (dPhi s j x y)

/-- The metric as a matrix. -/
def metricI (s : Surface) (x y : Float) : Mat2 :=
  ⟨metricEntry s x y 0 0, metricEntry s x y 0 1, metricEntry s x y 1 1⟩

namespace MetricI

def E (s : Surface) (x y : Float) : Float := metricEntry s x y 0 0
def Ff (s : Surface) (x y : Float) : Float := metricEntry s x y 0 1
def G (s : Surface) (x y : Float) : Float := metricEntry s x y 1 1

end MetricI

/-- `∂_j g_{ij}` — one metric derivative entry, direction `u = 0 (∂/∂x) or 1 (∂/∂y)`. -/
def dMetricEntry (s : Surface) (u : Fin 2) (i j : Fin 2) (x y : Float) : Float :=
  match u with
  | 0 => dd (fun t => metricEntry s t y i j) x
  | 1 => dd (fun t => metricEntry s x t i j) y

/-- `Γ^k_ij` at `(x, y)`: the Christoffel symbol of the second kind. -/
def christoffel (s : Surface) (x y : Float) :
    Fin 2 → Fin 2 → Fin 2 → Float :=
  fun k i j =>
    let gInv := (metricI s x y).inv
    let acc0 : Float :=
      0.5 * gInv.invEntry k 0 *
        (dMetricEntry s i 0 j x y + dMetricEntry s j 0 i x y - dMetricEntry s 0 i j x y)
    let acc1 : Float :=
      0.5 * gInv.invEntry k 1 *
        (dMetricEntry s i 1 j x y + dMetricEntry s j 1 i x y - dMetricEntry s 1 i j x y)
    acc0 + acc1

/-- `(i, j)` entry of the second fundamental form: `⟨n, ∂²_ij Φ⟩`. -/
def secondFormEntry (s : Surface) (x y : Float) : Fin 2 → Fin 2 → Float :=
  fun i j =>
    let n := normal s x y
    let d2comp (f : Float → Float → Float) : Float :=
      match i, j with
      | 0, 0 => ddd (fun t => f t y) x
      | 0, 1 => dd (fun t => dd (f t) y) x
      | 1, 0 => dd (fun t => dd (f t) y) x
      | 1, 1 => ddd (f x) y
      -- all 4 pairs of `Fin 2 × Fin 2` are covered
    n.dot ⟨d2comp s.sx, d2comp s.sy, d2comp s.sz⟩

/-- The second fundamental form as a matrix. -/
def metricII (s : Surface) (x y : Float) : Mat2 :=
  ⟨secondFormEntry s x y 0 0, secondFormEntry s x y 0 1, secondFormEntry s x y 1 1⟩

/-- Gaussian curvature `K = det II / det I` — intrinsic, lives in the metric alone. -/
def gaussK (s : Surface) (x y : Float) : Float :=
  (metricII s x y).det / (metricI s x y).det

/-- Mean curvature `H = (e G − 2 f F + g₂ E) / 2 (E G − F²)` — extrinsic. -/
def meanH (s : Surface) (x y : Float) : Float :=
  let E := MetricI.E s x y
  let Ff := MetricI.Ff s x y
  let G := MetricI.G s x y
  let II := metricII s x y
  let denom := 2 * (E * G - Ff * Ff)
  (II.a * G - 2 * II.b * Ff + II.c * E) / denom

end Surface

end Geodesia
