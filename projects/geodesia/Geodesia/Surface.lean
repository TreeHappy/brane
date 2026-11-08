/-!
# Geodesia — height fields and surfaces

Two levels of the tower:

* a *height field* `h : Float → Float → Float` is the `ℝ² → ℝ` case:
  map the domain `(x, y)` to a scalar (height / potential / temperature).
* a *surface* is the `ℝ² → ℝ³` case: an embedding `Φ` from the parameter plane
  into ambient space. The canonical trick in this library is the *terrain
  embedding* of a height field:
      Φ(x, y) = (x, y, h(x, y))
  so every height map is a bent sheet in 3-space, ready for surface geometry.

Everything numeric is exact-for-the-Float-example and Differentiable in spirit:
the surface functions are C² (indeed C^∞), so all derivatives below exist.
The library approximates them with central finite differences.
-/

import Geodesia.Vec

namespace Geodesia

/-- A height field over the parameter plane: `(x, y) ↦ h(x, y)`. -/
abbrev Height := Float → Float → Float

/-- A surface (parametrized patch): two parameters, three coordinates. -/
structure Surface where
  sx : Float → Float → Float
  sy : Float → Float → Float
  sz : Float → Float → Float

/-- The canonical terrain embedding: a height field lifted to a surface. -/
def Surface.ofHeight (h : Height) : Surface :=
  ⟨fun x _ => x, fun _ y => y, h⟩

namespace Surface

attribute [local simp] Vec3.add Vec3.sub Vec3.smul

/-! ## Central finite differences.

`d/du f(u) ≈ (f(u + ε) - f(u - ε)) / 2ε`, second derivative via the standard
3-point stencil. The step is tiny enough for smooth functions to look
"exact" in Float output.
-/

def eps : Float := 1e-5

/-- First partial with respect to the first parameter, for a unary profile. -/
def dd (f : Float → Float) (u : Float) : Float :=
  (f (u + eps) - f (u - eps)) / (2 * eps)

/-- Second partial for a unary profile (convexity of a cross-section). -/
def ddd (f : Float → Float) (u : Float) : Float :=
  (f (u + eps) - 2 * f u + f (u - eps)) / (eps * eps)

/-- First partials of the height field. -/
def hx (h : Height) (x y : Float) : Float := dd (fun u => h u y) x
def hy (h : Height) (x y : Float) : Float := dd (fun v => h x v) y
/-- Second pure partials (height field version of the Hessian diagonal). -/
def hxx (h : Height) (x y : Float) : Float := ddd (fun u => h u y) x
def hyy (h : Height) (x y : Float) : Float := ddd (fun v => h x v) y
/-- The mixed partial, read off the height field. -/
def hxy (h : Height) (x y : Float) : Float := dd (fun u => hy h u y) x

/-! ## Tangents and normal of a height field. -/

/-- Tangent along increasing `x`, sitting in the surface: `∂Φ/∂x`. -/
def tangentX (h : Height) (x y : Float) : Vec3 := ⟨1, 0, hx h x y⟩

/-- Tangent along increasing `y` on the surface: `∂Φ/∂y`. -/
def tangentY (h : Height) (x y : Float) : Vec3 := ⟨0, 1, hy h x y⟩

/-- Upward unit normal: `n = Φx × Φy / |Φx × Φy|`. -/
def normal (h : Height) (x y : Float) : Vec3 :=
  let c := (tangentX h x y).cross (tangentY h x y)
  Vec3.smul (1 / c.norm) c

/-- Slope magnitude `|∇h|` at a point — how steep the terrain is. -/
def slope (h : Height) (x y : Float) : Float :=
  Float.sqrt ((hx h x y) ^ 2 + (hy h x y) ^ 2)

end Surface

end Geodesia
