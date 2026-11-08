import Geodesia.Surface

/-!
# Geodesia — the metric, Christoffel symbols, curvature

For a surface patch `Φ : (u, v) ↦ Φ(u, v) ∈ ℝ³`, write `Φ_u = ∂Φ/∂u`,
`Φ_v = ∂Φ/∂v`, and similarly for the second partials.

## First fundamental form (the metric on the surface)

```
I = [[ E  F ]      E = ⟨Φ_u, Φ_u⟩,  F = ⟨Φ_u, Φ_v⟩,  G = ⟨Φ_v, Φ_v⟩
     [ F  G ]]
```

Equivalently `I = Jᵀ J`, where `J = [Φ_u  Φ_v]` is the 3×2 Jacobian of the
parameterization.  This is *the geometry*: every intrinsic measurement
(lengths, angles, areas, even Gaussian curvature) is computable from `I`
alone.

## Second fundamental form (how the surface bends in space)

```
II = [[ e  f ]      e = ⟨n, Φ_uu⟩,  f = ⟨n, Φ_uv⟩,  g₂ = ⟨n, Φ_vv⟩
      [ f  g₂ ]]     n = unit normal
```

The pair `(I, II)` is the complete local shape data.  The Gauss formula

```
Φ_ij = Γ^k_ij Φ_k + II_ij n
```

splits the second derivative of the embedding into a *tangent* part (the
Christoffel symbols) and a *normal* part (the second fundamental form).  That
single equation is the geometric meaning of both objects.

## Christoffel symbols

```
Γ^k_ij = ½ Σ_m I^{km} (∂_i I_{jm} + ∂_j I_{im} − ∂_m I_{ij})
```

They correct naive differentiation: the covariant derivative of a surface
vector field is `(∇_i V)^k = ∂_i V^k + Σ_j Γ^k_ij V^j`; the extra term comes
from the *turning* of the coordinate basis.  Geodesics are paths with zero
covariant acceleration (see `Geodesic.lean`).

## Curvature

The shape operator is `S = I⁻¹ II`; its eigenvalues are the principal
curvatures `κ₁, κ₂`.  Then

* Gaussian curvature `K = det S = det II / det I` — intrinsic;
* mean curvature `H = ½ tr S`.

For a graph `Φ(x, y) = (x, y, f(x, y))` this reduces to

```
I  = [[1 + f_x², f_x f_y], [f_x f_y, 1 + f_y²]]
II = Hessian(f) / √(1 + |∇f|²)
```

which is exactly what `gaussK`/`meanH` evaluate below.
-/

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

/-- `(k, m)` entry of a symmetric matrix. -/
def entry (M : Mat2) : Fin 2 → Fin 2 → Float
  | 0, 0 => M.a
  | 0, 1 => M.b
  | 1, 0 => M.b
  | 1, 1 => M.c
  -- all 4 pairs of `Fin 2 × Fin 2` are covered

end Mat2

namespace Surface

/-! ## First fundamental form. -/

/-- `(i, j)` entry of the first fundamental form: `⟨∂_iΦ, ∂_jΦ⟩`. -/
def metricEntry (s : Surface) (u v : Float) : Fin 2 → Fin 2 → Float :=
  fun i j => (dPhi s i u v).dot (dPhi s j u v)

/-- The metric as a matrix. -/
def metricI (s : Surface) (u v : Float) : Mat2 :=
  ⟨metricEntry s u v 0 0, metricEntry s u v 0 1, metricEntry s u v 1 1⟩

namespace MetricI

def E (s : Surface) (u v : Float) : Float := metricEntry s u v 0 0
def Ff (s : Surface) (u v : Float) : Float := metricEntry s u v 0 1
def G (s : Surface) (u v : Float) : Float := metricEntry s u v 1 1

end MetricI

/-- `∂_u I_{ij}` from the exact second partials:
`∂_u⟨Φ_i, Φ_j⟩ = ⟨Φ_{ui}, Φ_j⟩ + ⟨Φ_i, Φ_{uj}⟩`. -/
def dMetricEntry (s : Surface) (u i j : Fin 2) (x y : Float) : Float :=
  (d2Phi s u i x y).dot (dPhi s j x y)
    + (dPhi s i x y).dot (d2Phi s u j x y)

/-- `Γ^k_ij` at `(x, y)`: the Christoffel symbol of the second kind. -/
def christoffel (s : Surface) (x y : Float) :
    Fin 2 → Fin 2 → Fin 2 → Float :=
  fun k i j =>
    let gInv := (metricI s x y).inv
    let term0 :=
      dMetricEntry s i j 0 x y + dMetricEntry s j i 0 x y
        - dMetricEntry s 0 i j x y
    let term1 :=
      dMetricEntry s i j 1 x y + dMetricEntry s j i 1 x y
        - dMetricEntry s 1 i j x y
    0.5 * gInv.entry k 0 * term0 + 0.5 * gInv.entry k 1 * term1

/-! ## Second fundamental form. -/

/-- `(i, j)` entry of the second fundamental form: `⟨n, ∂²_ij Φ⟩`. -/
def secondFormEntry (s : Surface) (x y : Float) : Fin 2 → Fin 2 → Float :=
  fun i j => (normalAt s x y).dot (d2Phi s i j x y)

/-- The second fundamental form as a matrix. -/
def metricII (s : Surface) (x y : Float) : Mat2 :=
  ⟨secondFormEntry s x y 0 0, secondFormEntry s x y 0 1,
    secondFormEntry s x y 1 1⟩

/-- `I(w, w)` for a parameter-space vector `w`. -/
def firstForm (s : Surface) (x y : Float) (w : Vec2) : Float :=
  let I := metricI s x y
  I.a * w.x * w.x + 2 * I.b * w.x * w.y + I.c * w.y * w.y

/-- `II(w, w)` for a parameter-space vector `w`. -/
def secondForm (s : Surface) (x y : Float) (w : Vec2) : Float :=
  let II := metricII s x y
  II.a * w.x * w.x + 2 * II.b * w.x * w.y + II.c * w.y * w.y

/-- Normal curvature of the surface in direction `w`:
`κ_n = II(w, w) / I(w, w)`. -/
def normalCurvature (s : Surface) (x y : Float) (w : Vec2) : Float :=
  secondForm s x y w / firstForm s x y w

/-- The shape operator `S = I⁻¹ II`, applied to a parameter-space vector.
Its eigenvalues are the principal curvatures. -/
def shapeOperator (s : Surface) (x y : Float) (w : Vec2) : Vec2 :=
  let II := metricII s x y
  let gInv := (metricI s x y).inv
  let b1 := II.a * w.x + II.b * w.y
  let b2 := II.b * w.x + II.c * w.y
  ⟨gInv.a * b1 + gInv.b * b2, gInv.b * b1 + gInv.c * b2⟩

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

/-- Principal curvatures `(κ₁, κ₂)` from `K` and `H`: `κ = H ± √(H² − K)`.
Returns `(min, max)`. -/
def principalK (s : Surface) (x y : Float) : Float × Float :=
  let H := meanH s x y
  let K := gaussK s x y
  let disc0 := H * H - K
  let disc := Float.sqrt (if disc0 < 0 then 0 else disc0)
  (H - disc, H + disc)

end Surface

end Geodesia
