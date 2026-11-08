/-!
# Geodesia — differential operators

The planar operators on a height field and on 2-vector fields:

* gradient  `∇h = (h_x, h_y)`         — direction of steepest climb
* Hessian   ` Hess h = [[h_xx, h_xy], [h_yx, h_xy]]` — the 2×2 of all
  second partials (symmetric because the field is C²)
* Laplacian `Δh = h_xx + h_yy`        — sum of curvature in both directions
* the field operators for `F : (x, y) ↦ (F₁, F₂)`:
  divergence  `div F = ∂F₁/∂x + ∂F₂/∂y` (scalar, how much a point is a source)
  curl        `curl F = ∂F₂/∂x - ∂F₁/∂y` (scalar here; how much the field
  swirls — no z-component can exist in 2-D)
* Jacobian of a map `(x, y) ↦ (u, v)` — area stretch of the reparametrization

The star identity of the plane: `curl (∇h) = 0` for C² `h`, because
`h_xy = h_yx` (Clairaut / symmetry of mixed partials).
-/

import Geodesia.Surface

namespace Geodesia

/-- Gradient of a height field at a point. -/
def grad (h : Height) (x y : Float) : Vec2 := ⟨Surface.hx h x y, Surface.hy h x y⟩

/-- A planar (2-component) vector field. -/
abbrev PlanarField := Float → Float → Vec2

/-- Divergence: how "source-like" a field is at a point. -/
def divergence (F : PlanarField) (x y : Float) : Float :=
  let e := Surface.eps
  (F (x + e) y).x / (2 * e) - (F (x - e) y).x / (2 * e)
  + (F x (y + e)).y / (2 * e) - (F x (y - e)).y / (2 * e)

/-- Scalar curl: how "swirly" a planar field is at a point. -/
def curl2 (F : PlanarField) (x y : Float) : Float :=
  let e := Surface.eps
  (F (x + e) y).y / (2 * e) - (F (x - e) y).y / (2 * e)
  - (F x (y + e)).x / (2 * e) + (F x (y - e)).x / (2 * e)

/-- Hessian as a symmetric 2×2 Float matrix `[ [h_xx, h_xy], [h_xy, h_yy] ]`. -/
def hessian (h : Height) (x y : Float) : Float × Float × Float × Float :=
  (Surface.hxx h x y, Surface.hxy h x y, Surface.hxy h x y, Surface.hyy h x y)

/-- Laplacian: trace of the Hessian. -/
def laplacian (h : Height) (x y : Float) : Float :=
  Surface.hxx h x y + Surface.hyy h x y

/-- Jacobian determinant of a planar map `Φ : (x, y) ↦ (u(x,y), v(x,y))`. -/
def jacobianDet (u : Height) (v : Height) (x y : Float) : Float :=
  let ux := Surface.hx u x y
  let uy := Surface.hy u x y
  let vx := Surface.hx v x y
  let vy := Surface.hy v x y
  ux * vy - uy * vx

end Geodesia
