/-!
# Geodesia — differentiation as data

The finite-difference stencil is convenient for arbitrary user functions, but
it is only an *approximation*.  This module makes the derivative itself part
of the value:

* `Diff` is a unary function `f` packaged with its first and second
  derivatives `f'`, `f''`.
* `Smooth2` is a bivariate function `f(x, y)` packaged with all first and
  second partials `f_x`, `f_y`, `f_xx`, `f_xy`, `f_yy`.

Every combinator below (`const`, `add`, `mul`, `sin`, `exp`, …) differentiates
symbolically using the sum, product and chain rules.  If you build a surface
from these pieces, every metric entry, Christoffel symbol and curvature is
computed from *exact* derivatives — no `ε` anywhere.

`Smooth2.ofFiniteDiff` is the escape hatch for a bare `Float → Float → Float`
whose derivatives you do not know.  It is the only place a finite-difference
step appears.
-/

namespace Geodesia

/-! ## Unary functions -/

/-- A differentiable unary function, with its derivative functions supplied. -/
structure Diff where
  f : Float → Float
  d : Float → Float
  d2 : Float → Float

namespace Diff

def const (c : Float) : Diff :=
  { f := fun _ => c, d := fun _ => 0, d2 := fun _ => 0 }

def id : Diff :=
  { f := fun x => x, d := fun _ => 1, d2 := fun _ => 0 }

def add (a b : Diff) : Diff :=
  { f := fun x => a.f x + b.f x
    d := fun x => a.d x + b.d x
    d2 := fun x => a.d2 x + b.d2 x }

def neg (a : Diff) : Diff :=
  { f := fun x => -a.f x, d := fun x => -a.d x, d2 := fun x => -a.d2 x }

def sub (a b : Diff) : Diff := add a (neg b)

def smul (c : Float) (a : Diff) : Diff :=
  { f := fun x => c * a.f x
    d := fun x => c * a.d x
    d2 := fun x => c * a.d2 x }

/-- Product rule, twice: `(fg)'' = f''g + 2f'g' + fg''`. -/
def mul (a b : Diff) : Diff :=
  { f := fun x => a.f x * b.f x
    d := fun x => a.d x * b.f x + a.f x * b.d x
    d2 := fun x =>
      a.d2 x * b.f x + 2 * a.d x * b.d x + a.f x * b.d2 x }

/-- Chain rule, twice. -/
def comp (g a : Diff) : Diff :=
  { f := fun x => g.f (a.f x)
    d := fun x => g.d (a.f x) * a.d x
    d2 := fun x =>
      g.d2 (a.f x) * a.d x * a.d x + g.d (a.f x) * a.d2 x }

def sin (a : Diff) : Diff :=
  { f := fun x => Float.sin (a.f x)
    d := fun x => Float.cos (a.f x) * a.d x
    d2 := fun x =>
      -Float.sin (a.f x) * a.d x * a.d x + Float.cos (a.f x) * a.d2 x }

def cos (a : Diff) : Diff :=
  { f := fun x => Float.cos (a.f x)
    d := fun x => -Float.sin (a.f x) * a.d x
    d2 := fun x =>
      -Float.cos (a.f x) * a.d x * a.d x - Float.sin (a.f x) * a.d2 x }

def exp (a : Diff) : Diff :=
  { f := fun x => Float.exp (a.f x)
    d := fun x => Float.exp (a.f x) * a.d x
    d2 := fun x =>
      Float.exp (a.f x) * a.d x * a.d x + Float.exp (a.f x) * a.d2 x }

end Diff

/-! ## Bivariate functions -/

/-- A `C²` bivariate function, with every partial derivative as a function.

The fields are the data; there is no hidden approximation.  Construct values
with the combinators in `namespace Smooth2` (or fill in the fields yourself
when you know the derivatives). -/
structure Smooth2 where
  f : Float → Float → Float
  fx : Float → Float → Float
  fy : Float → Float → Float
  fxx : Float → Float → Float
  fxy : Float → Float → Float
  fyy : Float → Float → Float

namespace Smooth2

def const (c : Float) : Smooth2 :=
  { f := fun _ _ => c
    fx := fun _ _ => 0
    fy := fun _ _ => 0
    fxx := fun _ _ => 0
    fxy := fun _ _ => 0
    fyy := fun _ _ => 0 }

/-- The first coordinate function `(x, y) ↦ x`. -/
def coordX : Smooth2 :=
  { f := fun x _ => x
    fx := fun _ _ => 1
    fy := fun _ _ => 0
    fxx := fun _ _ => 0
    fxy := fun _ _ => 0
    fyy := fun _ _ => 0 }

/-- The second coordinate function `(x, y) ↦ y`. -/
def coordY : Smooth2 :=
  { f := fun _ y => y
    fx := fun _ _ => 0
    fy := fun _ _ => 1
    fxx := fun _ _ => 0
    fxy := fun _ _ => 0
    fyy := fun _ _ => 0 }

def add (a b : Smooth2) : Smooth2 :=
  { f := fun x y => a.f x y + b.f x y
    fx := fun x y => a.fx x y + b.fx x y
    fy := fun x y => a.fy x y + b.fy x y
    fxx := fun x y => a.fxx x y + b.fxx x y
    fxy := fun x y => a.fxy x y + b.fxy x y
    fyy := fun x y => a.fyy x y + b.fyy x y }

def neg (a : Smooth2) : Smooth2 :=
  { f := fun x y => -a.f x y
    fx := fun x y => -a.fx x y
    fy := fun x y => -a.fy x y
    fxx := fun x y => -a.fxx x y
    fxy := fun x y => -a.fxy x y
    fyy := fun x y => -a.fyy x y }

def sub (a b : Smooth2) : Smooth2 := add a (neg b)

def smul (c : Float) (a : Smooth2) : Smooth2 :=
  { f := fun x y => c * a.f x y
    fx := fun x y => c * a.fx x y
    fy := fun x y => c * a.fy x y
    fxx := fun x y => c * a.fxx x y
    fxy := fun x y => c * a.fxy x y
    fyy := fun x y => c * a.fyy x y }

/-- Product rule in two variables. -/
def mul (a b : Smooth2) : Smooth2 :=
  { f := fun x y => a.f x y * b.f x y
    fx := fun x y => a.fx x y * b.f x y + a.f x y * b.fx x y
    fy := fun x y => a.fy x y * b.f x y + a.f x y * b.fy x y
    fxx := fun x y =>
      a.fxx x y * b.f x y + 2 * a.fx x y * b.fx x y + a.f x y * b.fxx x y
    fxy := fun x y =>
      a.fxy x y * b.f x y + a.fx x y * b.fy x y
        + a.fy x y * b.fx x y + a.f x y * b.fxy x y
    fyy := fun x y =>
      a.fyy x y * b.f x y + 2 * a.fy x y * b.fy x y + a.f x y * b.fyy x y }

def sin (a : Smooth2) : Smooth2 :=
  { f := fun x y => Float.sin (a.f x y)
    fx := fun x y => Float.cos (a.f x y) * a.fx x y
    fy := fun x y => Float.cos (a.f x y) * a.fy x y
    fxx := fun x y =>
      -Float.sin (a.f x y) * a.fx x y * a.fx x y
        + Float.cos (a.f x y) * a.fxx x y
    fxy := fun x y =>
      -Float.sin (a.f x y) * a.fx x y * a.fy x y
        + Float.cos (a.f x y) * a.fxy x y
    fyy := fun x y =>
      -Float.sin (a.f x y) * a.fy x y * a.fy x y
        + Float.cos (a.f x y) * a.fyy x y }

def cos (a : Smooth2) : Smooth2 :=
  { f := fun x y => Float.cos (a.f x y)
    fx := fun x y => -Float.sin (a.f x y) * a.fx x y
    fy := fun x y => -Float.sin (a.f x y) * a.fy x y
    fxx := fun x y =>
      -Float.cos (a.f x y) * a.fx x y * a.fx x y
        - Float.sin (a.f x y) * a.fxx x y
    fxy := fun x y =>
      -Float.cos (a.f x y) * a.fx x y * a.fy x y
        - Float.sin (a.f x y) * a.fxy x y
    fyy := fun x y =>
      -Float.cos (a.f x y) * a.fy x y * a.fy x y
        - Float.sin (a.f x y) * a.fyy x y }

def exp (a : Smooth2) : Smooth2 :=
  { f := fun x y => Float.exp (a.f x y)
    fx := fun x y => Float.exp (a.f x y) * a.fx x y
    fy := fun x y => Float.exp (a.f x y) * a.fy x y
    fxx := fun x y =>
      Float.exp (a.f x y) * a.fx x y * a.fx x y
        + Float.exp (a.f x y) * a.fxx x y
    fxy := fun x y =>
      Float.exp (a.f x y) * a.fx x y * a.fy x y
        + Float.exp (a.f x y) * a.fxy x y
    fyy := fun x y =>
      Float.exp (a.f x y) * a.fy x y * a.fy x y
        + Float.exp (a.f x y) * a.fyy x y }

/-! ### Fallback: a bare function with finite differences -/

/-- A finite-difference step used only by `ofFiniteDiff`.

Prefer `Smooth2` combinators or a hand-filled record whenever the derivative
is known: the exact route is both faster and noise-free. -/
def eps : Float := 1e-5

def fd1 (g : Float → Float) (t : Float) : Float :=
  (g (t + eps) - g (t - eps)) / (2 * eps)

def fd2 (g : Float → Float) (t : Float) : Float :=
  (g (t + eps) - 2 * g t + g (t - eps)) / (eps * eps)

/-- Wrap an arbitrary `f : Float → Float → Float` using central differences.

This is the *only* approximation in the library.  It exists so exploratory
functions still work; use the exact constructors for anything you want to
differentiate repeatedly or trust quantitatively. -/
def ofFiniteDiff (f : Float → Float → Float) : Smooth2 :=
  { f := f
    fx := fun x y => fd1 (fun u => f u y) x
    fy := fun x y => fd1 (fun v => f x v) y
    fxx := fun x y => fd2 (fun u => f u y) x
    fxy := fun x y => fd1 (fun u => fd1 (fun v => f u v) y) x
    fyy := fun x y => fd2 (fun v => f x v) y }

end Smooth2

end Geodesia
