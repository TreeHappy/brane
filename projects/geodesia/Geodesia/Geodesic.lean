import Geodesia.Curvature

/-!
# Geodesia — geodesics as an initial-value problem

The geodesic equation in parameter space is

```
u'' + Γ¹₁₁ (u')² + 2Γ¹₁₂ u'v' + Γ¹₂₂ (v')² = 0
v'' + Γ²₁₁ (u')² + 2Γ²₁₂ u'v' + Γ²₂₂ (v')² = 0
```

with the Christoffel symbols from `Curvature.lean`.  Writing `p = u'` and
`q = v'` turns this into the first-order system

```
u' = p,  v' = q,  p' = F(u,v,p,q),  q' = G(u,v,p,q)
```

which `rk4Step` integrates.  `shoot` is the first half of a shooting method:
pick a unit tangent direction on the surface, convert it to parameter-space
velocity with `G⁻¹ Jᵀ`, integrate for a prescribed length, and return the
ambient endpoint.  Root-finding on the direction is the other half; see
`learn-geodesia/05-shooting.md`.
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

end GeoState

namespace Geodesic

/-- Parameter-space acceleration `(u'', v'')` from the geodesic equation. -/
def accel (s : Surface) (st : GeoState) : Vec2 :=
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

/-- Right-hand side of the first-order system. -/
def rhs (s : Surface) (st : GeoState) : GeoState :=
  let a := accel s st
  ⟨st.du, st.dv, a.x, a.y⟩

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
