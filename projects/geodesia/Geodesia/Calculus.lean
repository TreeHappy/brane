import Geodesia.Surface

/-!
# Geodesia — differential operators

The planar operators on a height field and on 2-vector fields.  Because the
height field and the field components are `Smooth2` values, every operator
below is exact: it reads the stored derivatives instead of differencing
samples.

* gradient  `∇h = (h_x, h_y)`         — direction of steepest climb
* Hessian   ` Hess h = [[h_xx, h_xy], [h_xy, h_yy]]` — the 2×2 of all
  second partials (symmetric because the field is C²)
* Laplacian `Δh = h_xx + h_yy`        — sum of curvature in both directions
* divergence `div F` and scalar curl `curl F` of a planar field
* Jacobian of a map `(x, y) ↦ (u, v)` — area stretch of the reparametrization

The star identity of the plane: `curl (∇h) = 0` for C² `h`, because
`h_xy = h_yx` (Clairaut / symmetry of mixed partials).
-/

namespace Geodesia

/-- A planar (2-component) vector field with exact partial derivatives. -/
structure PlanarField where
  px : Smooth2
  py : Smooth2

/-- Gradient of a height field at a point. -/
def grad (h : Height) (x y : Float) : Vec2 := Surface.grad h x y

/-- Divergence: how "source-like" a field is at a point. -/
def divergence (F : PlanarField) (x y : Float) : Float :=
  F.px.fx x y + F.py.fy x y

/-- Scalar curl: how "swirly" a planar field is at a point. -/
def curl2 (F : PlanarField) (x y : Float) : Float :=
  F.py.fx x y - F.px.fy x y

/-- Hessian as a symmetric 2×2 Float matrix `[ [h_xx, h_xy], [h_xy, h_yy] ]`. -/
def hessian (h : Height) (x y : Float) : Float × Float × Float × Float :=
  (Surface.hxx h x y, Surface.hxy h x y, Surface.hxy h x y, Surface.hyy h x y)

/-- Laplacian: trace of the Hessian. -/
def laplacian (h : Height) (x y : Float) : Float :=
  Surface.hxx h x y + Surface.hyy h x y

/-- Jacobian determinant of a planar map `Φ : (x, y) ↦ (u(x,y), v(x,y))`. -/
def jacobianDet (u v : Smooth2) (x y : Float) : Float :=
  u.fx x y * v.fy x y - u.fy x y * v.fx x y

/-! ## Differentiating along a curve

A curve is two unary `Diff` values, so its position **and** velocity are
available as functions.  The chain rule then gives derivatives of a height
field along the curve without ever forming a difference quotient. -/

/-- A parametrized planar curve `s ↦ (x(s), y(s))` with derivatives. -/
structure Curve2 where
  x : Diff
  y : Diff

namespace Curve2

/-- Position of the curve. -/
def point (c : Curve2) (s : Float) : Vec2 := ⟨c.x.f s, c.y.f s⟩

/-- Velocity of the curve. -/
def velocity (c : Curve2) (s : Float) : Vec2 := ⟨c.x.d s, c.y.d s⟩

end Curve2

/-- Directional derivative of a height field along a curve:
`d/ds h(c(s)) = h_x(c) x'(s) + h_y(c) y'(s)`. -/
def dAlong (h : Height) (c : Curve2) (s : Float) : Float :=
  let p := c.point s
  let v := c.velocity s
  h.fx p.x p.y * v.x + h.fy p.x p.y * v.y

end Geodesia
