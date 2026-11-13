import Geodesia.Curvature

/-!
# Geodesia — geodesics as an initial-value problem

A geodesic is a curve `s ↦ (u(s), v(s))` in the parameter plane whose
covariant acceleration vanishes.  In coordinates the geodesic equation is

```
u'' + Γ¹₁₁ (u')² + 2Γ¹₁₂ u'v' + Γ¹₂₂ (v')² = 0
v'' + Γ²₁₁ (u')² + 2Γ²₁₂ u'v' + Γ²₂₂ (v')² = 0
```

with the Christoffel symbols from `Curvature.lean`.  Writing

```
u' = du,  v' = dv,  u'' = a_u,  v'' = a_v
```

turns the second-order system into the first-order system

```
(u, v, du, dv)' = (du, dv, a_u, a_v).
```

`rhs` is exactly the map on the right: it leaves the two velocities alone and
asks an acceleration function for the last two slots.  `rhsWith` makes that
acceleration argument explicit, so the *same* integrator can drive either of
the two accelerations below.

## Two ways to compute the acceleration

1. **Christoffel route** (`accelChristoffel`).  Read `Γ^k_ij` from the
   metric and evaluate

   ```
   a^k = −Σ_{i,j} Γ^k_ij p^i p^j,     p = (du, dv).
   ```

   This is the textbook coordinate formula.

2. **Vector-calculus route** (`accelVector`).  Differentiate the ambient
   curve `Φ(u(s), v(s))` twice.  The chain rule gives

   ```
   Φ'' = Φ_u u'' + Φ_v v'' + Φ_uu (u')² + 2Φ_uv u'v' + Φ_vv (v')².
   ```

   A geodesic has no *tangential* acceleration, so `⟨Φ_k, Φ''⟩ = 0` for
   `k = u, v`.  That is the linear system `I (u'', v'') = −g`, where
   `I` is the first fundamental form and

   ```
   g_k = ⟨Φ_k, Φ_uu⟩ (u')² + 2⟨Φ_k, Φ_uv⟩ u'v' + ⟨Φ_k, Φ_vv⟩ (v')².
   ```

   `gVector` computes `g`, and `accelVector` inverts `I` and negates.
   No Christoffel symbol appears: this route uses only dot products of the
   embedding and its exact partials.

The two formulas are the same identity — substituting the definition of `Γ`
into the coordinate equation produces the vector equation.  The identity is
machine-checked over `ℝ` in `learn-geodesia/Proofs.lean`.

## Integration

`eulerStep` is one forward-Euler step: follow the tangent line of the state
for time `h`.  `rk4Step` takes four such slope samples and averages them, so
it is fourth-order accurate.  `eulerIntegrate` and `integrate` repeat their
single step a prescribed number of times.

`shoot` is the first half of a shooting method: pick a unit tangent direction
on the surface, convert it to parameter-space velocity with `G⁻¹ Jᵀ`, integrate
for a prescribed length, and return the ambient endpoint.  Root-finding on the
direction is the other half; see `learn-geodesia/05-shooting.md`.
-/

namespace Geodesia

/-- State of the first-order geodesic system: two parameters and their rates. -/
structure GeoState where
  u : Float
  v : Float
  du : Float
  dv : Float

namespace GeoState

def add (a b : GeoState) : GeoState :=
  ⟨a.u + b.u, a.v + b.v, a.du + b.du, a.dv + b.dv⟩

def smul (c : Float) (a : GeoState) : GeoState :=
  ⟨c * a.u, c * a.v, c * a.du, c * a.dv⟩

/-- Compact demo printing: `(u = …, v = …, du = …, dv = …)`. -/
instance : ToString GeoState :=
  ⟨fun st => s!"(u = {st.u}, v = {st.v}, du = {st.du}, dv = {st.dv})"⟩

end GeoState

namespace Geodesic

/-- Coordinate-space acceleration `(u'', v'')` from the Christoffel symbols:

`a^k = −Σ_{i,j} Γ^k_ij p^i p^j`, with `p = (du, dv)`.

The sum is written out using the symmetry `Γ^k_01 = Γ^k_10`, so the mixed
term appears twice.  This is the coordinate form of the geodesic equation. -/
def accelChristoffel (s : Surface) (st : GeoState) : Vec2 :=
  let Γ := Surface.christoffel s st.u st.v
  let g1uu := Γ 0 0 0
  let g1uv := Γ 0 0 1
  let g1vv := Γ 0 1 1
  let g2uu := Γ 1 0 0
  let g2uv := Γ 1 0 1
  let g2vv := Γ 1 1 1
  let du := st.du
  let dv := st.dv
  ⟨-(g1uu * du * du + 2 * g1uv * du * dv + g1vv * dv * dv),
   -(g2uu * du * du + 2 * g2uv * du * dv + g2vv * dv * dv)⟩

/-- The bending part of the ambient acceleration,

`g_k = ⟨Φ_k, Φ_uu⟩ (u')² + 2⟨Φ_k, Φ_uv⟩ u'v' + ⟨Φ_k, Φ_vv⟩ (v')²`

for `k = u, v`.  A geodesic satisfies `I (u'', v'') = −g`; see
`accelVector`. -/
def gVector (s : Surface) (st : GeoState) : Vec2 :=
  let Φu := Surface.dPhi s 0 st.u st.v
  let Φv := Surface.dPhi s 1 st.u st.v
  let Φuu := Surface.d2Phi s 0 0 st.u st.v
  let Φuv := Surface.d2Phi s 0 1 st.u st.v
  let Φvv := Surface.d2Phi s 1 1 st.u st.v
  let du := st.du
  let dv := st.dv
  ⟨Φu.dot Φuu * du * du + 2 * Φu.dot Φuv * du * dv + Φu.dot Φvv * dv * dv,
   Φv.dot Φuu * du * du + 2 * Φv.dot Φuv * du * dv + Φv.dot Φvv * dv * dv⟩

/-- Vector-calculus acceleration: require the tangential part of `Φ''` to
vanish.  Since `Φ'' = Φ_u u'' + Φ_v v'' + (grid-bending terms)`, the two
tangential conditions are `I (u'', v'') = −gVector`, so inverting the metric
gives `(u'', v'') = −I⁻¹ gVector`.

This is the same acceleration as `accelChristoffel`, but computed from dot
products instead of Christoffel symbols. -/
def accelVector (s : Surface) (st : GeoState) : Vec2 :=
  let I := Surface.metricI s st.u st.v
  let gInv := I.inv
  let g := gVector s st
  ⟨-(gInv.a * g.x + gInv.b * g.y), -(gInv.b * g.x + gInv.c * g.y)⟩

/-- Package a chosen acceleration function into the first-order system:

`(u, v, du, dv) ↦ (du, dv, a_u, a_v)`.

Keeping the acceleration abstract is what lets `eulerStep`, `rk4Step` and
`integrate` run either `accelChristoffel` or `accelVector`. -/
def rhsWith (accel : Surface → GeoState → Vec2) (s : Surface)
    (st : GeoState) : GeoState :=
  let a := accel s st
  ⟨st.du, st.dv, a.x, a.y⟩

/-- Right-hand side of the coordinate (Christoffel) geodesic system. -/
def rhs (s : Surface) (st : GeoState) : GeoState :=
  rhsWith accelChristoffel s st

/-- Right-hand side of the vector-calculus geodesic system. -/
def rhsVector (s : Surface) (st : GeoState) : GeoState :=
  rhsWith accelVector s st

/-- One forward-Euler step: follow the tangent line of the state for time `h`.

Euler is the simplest integrator; `rk4Step` below is the accurate one. -/
def eulerStep (s : Surface) (h : Float) (st : GeoState) : GeoState :=
  GeoState.add st (GeoState.smul h (rhs s st))

/-- Integrate `n` Euler steps of size `h`. -/
def eulerIntegrate (s : Surface) (h : Float) (n : Nat) (st : GeoState) : GeoState :=
  match n with
  | 0 => st
  | n + 1 => eulerIntegrate s h n (eulerStep s h st)

/-- One classical fourth-order Runge–Kutta step. -/
def rk4Step (s : Surface) (h : Float) (st : GeoState) : GeoState :=
  let k1 := rhs s st
  let k2 := rhs s (GeoState.add st (GeoState.smul (h / 2) k1))
  let k3 := rhs s (GeoState.add st (GeoState.smul (h / 2) k2))
  let k4 := rhs s (GeoState.add st (GeoState.smul h k3))
  let avg :=
    GeoState.add (GeoState.add k1 (GeoState.smul 2 k2))
      (GeoState.add (GeoState.smul 2 k3) k4)
  GeoState.add st (GeoState.smul (h / 6) avg)

/-- Integrate `n` RK4 steps of size `h`. -/
def integrate (s : Surface) (h : Float) (n : Nat) (st : GeoState) : GeoState :=
  match n with
  | 0 => st
  | n + 1 => integrate s h n (rk4Step s h st)

/-- Kinetic energy `½ I(ξ, ξ)`; constant along an affinely parametrized
geodesic (a first integral of the ODE). -/
def energy (s : Surface) (st : GeoState) : Float :=
  0.5 * Surface.firstForm s st.u st.v ⟨st.du, st.dv⟩

/-- Speed `|γ'| = √(2E)`. -/
def speed (s : Surface) (st : GeoState) : Float :=
  Float.sqrt (2 * energy s st)

/-- An orthonormal frame `(e₁, e₂)` for the tangent plane at `(u, v)`. -/
def tangentFrame (s : Surface) (u v : Float) : Vec3 × Vec3 :=
  let e1 := Vec3.normalize (Surface.dPhi s 0 u v)
  let raw := Surface.dPhi s 1 u v
  let e2 := Vec3.normalize (Vec3.sub raw (Vec3.smul (raw.dot e1) e1))
  (e1, e2)

/-- Convert an ambient tangent vector `w` into parameter-space velocity:
`ξ = I⁻¹ Jᵀ w`, exactly the formula in Eberly's note. -/
def tangentVecToParams (s : Surface) (u v : Float) (w : Vec3) : Vec2 :=
  let gInv := (Surface.metricI s u v).inv
  let a := (Surface.dPhi s 0 u v).dot w
  let b := (Surface.dPhi s 1 u v).dot w
  ⟨gInv.a * a + gInv.b * b, gInv.b * a + gInv.c * b⟩

/-- Integrate from `(u, v)` in the unit tangent direction `angle`, for arc
length `length`, and return the final state.  `steps` controls the RK4 grid. -/
def shootState (s : Surface) (u v : Float) (angle length : Float)
    (steps : Nat) : GeoState :=
  let frame := tangentFrame s u v
  let w := Vec3.add (Vec3.smul (Float.cos angle) frame.1)
    (Vec3.smul (Float.sin angle) frame.2)
  let xi := tangentVecToParams s u v w
  let h := length / Float.ofNat steps
  integrate s h steps ⟨u, v, xi.x, xi.y⟩

/-- The ambient endpoint of the shot. -/
def shoot (s : Surface) (u v : Float) (angle length : Float)
    (steps : Nat) : Vec3 :=
  let st := shootState s u v angle length steps
  Surface.point s st.u st.v

end Geodesic

end Geodesia
