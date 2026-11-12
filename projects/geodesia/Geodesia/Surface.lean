import Geodesia.Diff
import Geodesia.Vec

/-!
# Geodesia — height fields and surfaces

Two levels of the tower:

* a *height field* `h : Smooth2` is the `ℝ² → ℝ` case: map the domain
  `(x, y)` to a scalar (height / potential / temperature).
* a *surface* is the `ℝ² → ℝ³` case: an embedding `Φ` from the parameter
  plane into ambient space. The canonical trick in this library is the
  *terrain embedding* of a height field:
      Φ(x, y) = (x, y, h(x, y))
  so every height map is a bent sheet in 3-space, ready for surface geometry.

`Smooth2` carries the derivatives as data, so `Φ_u`, `Φ_v`, `Φ_uu`, … below
are exact whenever the height/surface was built from the `Diff` combinators.
A finite-difference fallback is available via `Smooth2.ofFiniteDiff`, but the
rest of the library never needs an `ε`.
-/

namespace Geodesia

/-- A height field over the parameter plane, with exact partial derivatives. -/
abbrev Height := Smooth2

/-- A surface (parametrized patch): two parameters, three coordinate
functions, each carrying its own partial derivatives. -/
structure Surface where
  sx : Smooth2
  sy : Smooth2
  sz : Smooth2

/-- The canonical terrain embedding: a height field lifted to a surface. -/
def Surface.ofHeight (h : Height) : Surface :=
  ⟨Smooth2.coordX, Smooth2.coordY, h⟩

namespace Surface

/-! ## Point and partial derivatives of `Φ`. -/

/-- The ambient point `Φ(u, v)`. -/
def point (s : Surface) (u v : Float) : Vec3 :=
  ⟨s.sx.f u v, s.sy.f u v, s.sz.f u v⟩

/-- `∂Φ/∂u` (`i = 0`) or `∂Φ/∂v` (`i = 1`). -/
def dPhi (s : Surface) (i : Fin 2) (u v : Float) : Vec3 :=
  match i with
  | 0 => ⟨s.sx.fx u v, s.sy.fx u v, s.sz.fx u v⟩
  | 1 => ⟨s.sx.fy u v, s.sy.fy u v, s.sz.fy u v⟩

/-- The parameter-direction arrow `J d = d.x • Φu + d.y • Φv`: the linear map
from the parameter plane to the tangent space at (u, v). -/
def tangentVec (s : Surface) (u v : Float) (d : Vec2) : Vec3 :=
  Vec3.add (Vec3.smul d.x (dPhi s 0 u v)) (Vec3.smul d.y (dPhi s 1 u v))

/-- `∂²Φ/∂i∂j`, using the exact second partials stored in `Smooth2`. -/
def d2Phi (s : Surface) (i j : Fin 2) (u v : Float) : Vec3 :=
  match i, j with
  | 0, 0 => ⟨s.sx.fxx u v, s.sy.fxx u v, s.sz.fxx u v⟩
  | 0, 1 => ⟨s.sx.fxy u v, s.sy.fxy u v, s.sz.fxy u v⟩
  | 1, 0 => ⟨s.sx.fxy u v, s.sy.fxy u v, s.sz.fxy u v⟩
  | 1, 1 => ⟨s.sx.fyy u v, s.sy.fyy u v, s.sz.fyy u v⟩

/-! ## Height-field helpers. -/

/-- First partials of the height field. -/
def hx (h : Height) (x y : Float) : Float := h.fx x y
def hy (h : Height) (x y : Float) : Float := h.fy x y

/-- Second pure partials (height field version of the Hessian diagonal). -/
def hxx (h : Height) (x y : Float) : Float := h.fxx x y
def hyy (h : Height) (x y : Float) : Float := h.fyy x y

/-- The mixed partial, read off the height field. -/
def hxy (h : Height) (x y : Float) : Float := h.fxy x y

/-- Gradient `∇h = (h_x, h_y)` — direction of steepest climb. -/
def grad (h : Height) (x y : Float) : Vec2 := ⟨hx h x y, hy h x y⟩

/-! ## Tangents and normal of a height field. -/

/-- Tangent along increasing `x`, sitting in the surface: `∂Φ/∂x`. -/
def tangentX (h : Height) (x y : Float) : Vec3 := ⟨1, 0, hx h x y⟩

/-- Tangent along increasing `y` on the surface: `∂Φ/∂y`. -/
def tangentY (h : Height) (x y : Float) : Vec3 := ⟨0, 1, hy h x y⟩

/-- Unit normal of a general surface: `n = Φ_u × Φ_v / |Φ_u × Φ_v|`. -/
def normalAt (s : Surface) (u v : Float) : Vec3 :=
  let c := (dPhi s 0 u v).cross (dPhi s 1 u v)
  Vec3.smul (1 / c.norm) c

/-- The normal part `(w · n) n` of an ambient arrow `w` — the bending
component of the Gauss split. -/
def normalPart (s : Surface) (u v : Float) (w : Vec3) : Vec3 :=
  Vec3.smul (w.dot (normalAt s u v)) (normalAt s u v)

/-- The tangential part `T w = w - (w · n) n` of an ambient arrow `w` — the
surface part of the Gauss split; `w = tangentialPart w + normalPart w`. -/
def tangentialPart (s : Surface) (u v : Float) (w : Vec3) : Vec3 :=
  Vec3.sub w (normalPart s u v w)

/-- Upward unit normal of a height field. -/
def normal (h : Height) (x y : Float) : Vec3 :=
  normalAt (ofHeight h) x y

/-- Slope magnitude `|∇h|` at a point — how steep the terrain is. -/
def slope (h : Height) (x y : Float) : Float :=
  let g := grad h x y
  g.norm

end Surface

end Geodesia
