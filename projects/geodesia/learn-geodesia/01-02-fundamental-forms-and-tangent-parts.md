# Geodesia — the tangent space (Ch. 1–2)

> *The surface is curved; the tangent space is flat. This chapter lives in
> the flat world — and takes two measurements of how the curved world bends
> it.*

> **What this chapter is.** One object, seen from several sides: the
> **tangent space** of a parameterized surface — the flat 2-D vector world
> attached to each of its points. The classical "two fundamental forms" are
> two things you can do in that world: measure with its inner product
> ($I = J^{\mathsf T}J$), and record how acceleration leaves it
> ($II_{ij} = n\cdot\Phi_{ij}$). Everything is ordinary vector calculus:
> arrows, dot products, one subtraction. The chapter is written like the
> Lean documentation (as in the Verso-built Lean Manual and *Theorem
> Proving in Lean*): every boxed formula is followed by the **Lean
> declarations that compute it** — a **Signature** block — and every
> numeric claim is a **Checked example**: the Lean snippet together with
> the output it really prints. Nothing in an `Output:` block is
> hand-computed; each was produced by running the snippet against the
> library. The math notation is kept in full ($\Phi$, $I$, $II$, $T$,
> $\Gamma$, …); the Lean identifiers are the computational shadows of the
> same objects.
>
> It merges the former `01-fundamental-forms.md` and `02-tangent-parts.md`.
> The section numbers of both chapters are preserved (`§1.1–1.8`,
> `§2.1–2.9`), so the cross-references from Chapters 3–8, `index.md` and
> `GA_VIEW.md` still resolve. Three small declarations were added to the
> library for this chapter — `Surface.tangentVec`, `Surface.tangentialPart`,
> `Surface.normalPart` — so that the vector statements of §1.1 and §1.5 are
> executable, not just prose.

**Reproducing the checked examples.** All snippets run against this project:

```bash
cd projects/geodesia
lake build
lake exe geodesia                   # the demo printed in §1.6
lake env lean <your-scratch>.lean   # any snippet below (it begins with `import Geodesia`)
```

The dictionary that the whole chapter keeps re-deriving:

| Math object | Lean declaration | File |
|---|---|---|
| $\Phi(u,v)$ | `Surface.point` | `Geodesia/Surface.lean` |
| $\Phi_i = \partial\Phi/\partial i$ | `Surface.dPhi` | `Geodesia/Surface.lean` |
| $\Phi_{ij} = \partial^2\Phi/\partial i\,\partial j$ | `Surface.d2Phi` | `Geodesia/Surface.lean` |
| $n$ | `Surface.normalAt` | `Geodesia/Surface.lean` |
| $J\,d$ — the arrow of a parameter direction | `Surface.tangentVec` | `Geodesia/Surface.lean` |
| $I = J^{\mathsf T}J$ — the inner product of the tangent space | `Surface.metricI`, `Surface.firstForm` | `Geodesia/Curvature.lean` |
| $II_{ij} = n\cdot\Phi_{ij}$ — the acceleration table | `Surface.metricII`, `Surface.secondForm` | `Geodesia/Curvature.lean` |
| $T w = w - (w\cdot n)\,n$ — the tangent part of the split | `Surface.tangentialPart` | `Geodesia/Surface.lean` |
| $(w\cdot n)\,n$ — the normal part of the split | `Surface.normalPart` | `Geodesia/Surface.lean` |
| $\kappa(d) = II(d,d)/I(d,d)$ | `Surface.normalCurvature` | `Geodesia/Curvature.lean` |
| $S = I^{-1}II$ | `Surface.shapeOperator` | `Geodesia/Curvature.lean` |
| $K$, $H$, $\kappa_{1,2}$ | `Surface.gaussK`, `Surface.meanH`, `Surface.principalK` | `Geodesia/Curvature.lean` |
| $\Gamma^k_{\;ij}$ (coordinates of $T\Phi_{ij}$) | `Surface.christoffel` | `Geodesia/Curvature.lean` |
| $\partial_u\,(\Phi_j\cdot\Phi_k)$ | `Surface.dMetricEntry` | `Geodesia/Curvature.lean` |
| $\xi = I^{-1}J^{\mathsf T}w$ | `Geodesic.tangentVecToParams` | `Geodesia/Geodesic.lean` |
| exact derivative arrows | `Smooth2` + its combinators | `Geodesia/Diff.lean` |

**The notation bridge.** The Lean code uses short identifiers where the text
uses math symbols. Every occurrence is translatable by this table — and each
**Signature** block below is followed by the implementation and its
vector-calculus reading, so any formula can be checked against the code:

| Lean notation | Vector-calculus meaning |
|---|---|
| `s.sx`, `s.sy`, `s.sz` | the coordinate functions $x(u,v)$, $y(u,v)$, $z(u,v)$ of $\Phi$ |
| `i : Fin 2` with values `0`, `1` | the parameter index: $0 \leftrightarrow u$, $1 \leftrightarrow v$ |
| `u v : Float` | the point $(u, v)$ of the parameter plane |
| `Surface.dPhi s i u v` | the arrow $\Phi_i(u,v)$ |
| `Surface.d2Phi s i j u v` | the arrow $\Phi_{ij}(u,v)$ |
| `Surface.tangentVec s u v d` | the tangent arrow $J\,d = a\,\Phi_u + b\,\Phi_v$ of the parameter direction $d=(a,b)$ |
| `Surface.tangentialPart s u v w` | the tangent part $T\,w$ of the split of an ambient arrow $w$ |
| `Surface.normalPart s u v w` | the normal part $(w\cdot n)\,n$ of the split of $w$ |
| `v.dot w`, `v.cross w`, `v.norm` | $v\cdot w$, $v\times w$, $\lvert v\rvert$ |
| `Vec3.smul a v` | the scaled arrow $a\,v$ |
| `Mat2` with fields `a`, `b`, `c` | the symmetric matrix $\begin{bmatrix}a & b\\ b & c\end{bmatrix}$ |
| `metricI` / `metricII` | $I = J^{\mathsf T}J$ / $II = [\,n\cdot\Phi_{ij}\,]$ |
| `firstForm` / `secondForm` | the quadratic forms $I(w, w) = w^{\mathsf T}Iw$ / $II(w, w) = w^{\mathsf T}II\,w$ |
| `Surface.christoffel s x y k i j` | the coefficient $\Gamma^k_{\;ij}$ |
| `GeoState` fields `u v du dv` | $q(s) = (u,v)$ and $q'(s) = (u',v')$ along a curve |

(Implementations are shown verbatim from the sources, with the fully
qualified names they get from their `namespace`.)

![The anchor picture: the paraboloid with its parameter grid; at one point the tangent space T_pS (dashed outline) carrying the tangent arrows Φu, Φv, a parameter direction d and its arrow Jd (green), and the normal n; a curve γ(s) whose acceleration γ″ splits into the tangent part Tγ″ (peach) and the normal part II(d,d)n (teal).](images/fundamental-forms.png)

---

## Part 1 — The tangent space: one flat world, two measurements

> **Goal.** Understand exactly what `metricI`, `metricII`, `gaussK`,
> `meanH` and `principalK` in `Curvature.lean` are computing — as three
> facts about one flat world: the tangent space has an inner product (that
> is $I$), acceleration has a normal component (that is $II$), and a linear
> operator compares the two (that is $S$, hence $K$ and $H$). Everything
> here is ordinary vector math: arrows, dot products, one subtraction.

### 1.1 The surface and its tangent space

> **Picture.** Take a sheet and bend it smoothly — no creases. At every
> point, lay a small flat plane against it, the way a coin rests on a ball.
> That little plane is the **tangent space** of the point: a complete 2-D
> vector world, flat as $\mathbb{R}^2$, with its own arrows, its own
> lengths and angles — and, because it sits inside 3-space, it *inherits*
> the dot products of the ambient world. The sheet is curved; its tangent
> spaces never are. Nearly everything in this chapter happens in those
> flat worlds.

![The tangent space, schematically: the sheet seen from the side (curved), with its flat tangent spaces T_pS attached at three points, each carrying the tangent arrow Φu and the normal n; the zoom bubble at one point shows an ambient arrow w split uniquely into its tangent part T w (peach) and its normal part (w·n)n (teal) — the flat world of one point.](images/tangent-space.png)

A surface is a function

```math
\Phi : \mathbb{R}^2 \to \mathbb{R}^3,
\qquad \Phi(u,v) = \big(x(u,v),\ y(u,v),\ z(u,v)\big).
```

You feed it a point of the parameter plane and it returns a point of 3-space.
The two numbers $(u,v)$ are just coordinates — like row and column indices
painted on the surface.

Differentiate $\Phi$ with respect to each of its two inputs. Each derivative
is a genuine 3-D arrow:

```math
\Phi_u = \begin{bmatrix}\partial x/\partial u\\ \partial y/\partial u\\ \partial z/\partial u\end{bmatrix},
\qquad
\Phi_v = \begin{bmatrix}\partial x/\partial v\\ \partial y/\partial v\\ \partial z/\partial v\end{bmatrix}.
```

Put them side by side as the columns of a $3\times 2$ matrix, the
**Jacobian** $J$:

```math
J = \begin{bmatrix}\Phi_u & \Phi_v\end{bmatrix}.
```

**The tangent space.** Wherever you stand, the two arrows span a plane —
and this is the sentence the whole chapter stands on: *every linear
combination of them is again the velocity of a real curve on the surface.*
Lift the parameter-plane line $t \mapsto p + t\,(a,b)$ to
$\gamma(t) = \Phi\big(p + t\,(a,b)\big)$ and differentiate at $t = 0$: the
velocity is exactly $a\,\Phi_u + b\,\Phi_v$. So the span deserves the name

```math
T_pS \;=\; \{\, a\,\Phi_u + b\,\Phi_v \;:\; a, b \in \mathbb{R} \,\}
\;\subset\;\mathbb{R}^3,
```

the **tangent space** at $p$. It is a genuine vector space — sums and
scalings of tangent arrows stay tangent, because a span always is — and it
is *flat*: no curvature of the sheet leaks into it. Curvature will enter
only as how these flat worlds tilt relative to each other from point to
point (Part 2), never inside one of them.

> **Caution.** *Is it a problem that this flat world floats inside an extra
> dimension?* For embedded surfaces: none — the ambient $\mathbb{R}^3$ is a
> gift. The normal $n$, the projection $T$, the whole of §1.3 come free
> with it, and that is why this chapter can work in honest arrows. Physics
> is stricter, and for good reasons: spacetime comes with no surrounding
> space to project into, a hypothetical embedding would be unobservable
> scaffolding, and some perfectly honest surfaces have *no* embedding at
> all — Hilbert's theorem forbids the hyperbolic plane in $\mathbb{R}^3$,
> and the intrinsically flat torus needs $\mathbb{R}^4$. That is why
> Chapter 8 rebuilds everything from the ruler $I$ alone — not because the
> extra dimension is forbidden, but because it is not always there when you
> need it. And where physicists *do* assume one (brane worldvolumes, the
> ADM 3+1 split), the extrinsic objects of this chapter come right back.

**Signature** — `Geodesia/Surface.lean` stores the arrows, not approximations
of them (the `Fin 2` index `0, 1` plays the role of $u, v$):

```lean
structure Surface where
  sx : Smooth2          -- each coordinate function carries its own
  sy : Smooth2          -- exact partial derivatives
  sz : Smooth2

def Surface.ofHeight (h : Height) : Surface
    -- the terrain embedding Φ(x, y) = (x, y, h(x, y))
def Surface.point (s : Surface) (u v : Float) : Vec3
    -- Φ(u, v)
def Surface.dPhi (s : Surface) (i : Fin 2) (u v : Float) : Vec3
    -- Φ_u (i = 0) and Φ_v (i = 1)
def Surface.d2Phi (s : Surface) (i j : Fin 2) (u v : Float) : Vec3
    -- Φ_ij: the curvature arrows
```

**Implementation.** The arrows are built component-wise from the exact
partials stored in `Smooth2`:

```lean
def Surface.point (s : Surface) (u v : Float) : Vec3 :=
  ⟨s.sx.f u v, s.sy.f u v, s.sz.f u v⟩

def Surface.dPhi (s : Surface) (i : Fin 2) (u v : Float) : Vec3 :=
  match i with
  | 0 => ⟨s.sx.fx u v, s.sy.fx u v, s.sz.fx u v⟩
  | 1 => ⟨s.sx.fy u v, s.sy.fy u v, s.sz.fy u v⟩

def Surface.d2Phi (s : Surface) (i j : Fin 2) (u v : Float) : Vec3 :=
  match i, j with
  | 0, 0 => ⟨s.sx.fxx u v, s.sy.fxx u v, s.sz.fxx u v⟩
  | 0, 1 => ⟨s.sx.fxy u v, s.sy.fxy u v, s.sz.fxy u v⟩
  | 1, 0 => ⟨s.sx.fxy u v, s.sy.fxy u v, s.sz.fxy u v⟩
  | 1, 1 => ⟨s.sx.fyy u v, s.sy.fyy u v, s.sz.fyy u v⟩
```

**Vector-calculus reading.** Component $k$ of `dPhi s i` is $\partial_i x_k$:
e.g. $\Phi_u = (\partial_u x,\ \partial_u y,\ \partial_u z)$, the displayed
matrix column by column. Both mixed-partial rows read the same field `fxy`
because $\Phi_{uv} = \Phi_{vu}$ (Clairaut's theorem). Note what is *absent*:
the code never materializes the $3\times 2$ Jacobian $J$ itself — every
consumer below needs only its columns and their dot products. (Curves
$\gamma$ are carried by `GeoState` $= (q, q')$ and integrated in Chapter 4.)

#### What $J$ does to a direction

Introduce a real parameter $s \in \mathbb{R}$ and a differentiable
**parameter-plane curve**

```math
c(s) := \bigl(u(s),\, v(s)\bigr) \in \mathbb{R}^2,
```

where $u(s), v(s) \in \mathbb{R}$ are the parameter coordinates along the
curve. Composing with $\Phi$ gives the corresponding **surface curve**

```math
\gamma(s) := \Phi\!\bigl(c(s)\bigr) = \Phi\!\bigl(u(s),\, v(s)\bigr) \in \mathbb{R}^3.
```

Write derivatives with respect to $s$ using a prime,
$u'(s) := \frac{du}{ds}$, $v'(s) := \frac{dv}{ds}$, and
$\gamma'(s) := \frac{d\gamma}{ds}$. Differentiating $\gamma$ by the chain
rule gives

```math
\gamma'(s) = \frac{d\gamma}{ds}
= \Phi_u\, u'(s) + \Phi_v\, v'(s)
= J \begin{bmatrix}u'(s)\\ v'(s)\end{bmatrix}.
```

So $J$ converts a 2-D direction $(a,b)$ in the parameter plane into the
actual 3-D direction you move on the surface:

```math
\boxed{\;J\begin{bmatrix}a\\ b\end{bmatrix} = a\,\Phi_u + b\,\Phi_v
\;\in\; T_pS.\;}
```

That is all "parameter space" means here: the numbers $(a,b)$ live in the
flat plane, and $J$ carries them into the flat tangent world. No other
machinery.

**Signature** — the linear map from parameter directions to tangent arrows
now has a name in the code:

```lean
def Surface.tangentVec (s : Surface) (u v : Float) (d : Vec2) : Vec3
    -- the tangent arrow J d = d.x • Φu + d.y • Φv
```

**Checked example.** Handing it the coordinate directions reproduces the
basis arrows; a diagonal direction gives the diagonal arrow:

```lean
#eval show IO Unit from do
  IO.println s!"J(1,0) = {Surface.tangentVec surf 1 0 ⟨1, 0⟩}"
  IO.println s!"J(1,1) = {Surface.tangentVec surf 1 0 ⟨1, 1⟩}"
```

**Output:**
```text
J(1,0) = (1.000000, 0.000000, 0.600000)
J(1,1) = (1.000000, 1.000000, 0.600000)
```

`J(1,0)` is literally the $\Phi_x$ of the four-arrow table below — the same
arrow, now obtained by feeding the *direction* to the linear map. And
`J(1,1) = Φ_x + Φ_y`: the arrow of the diagonal direction, reused by the
angles example of §1.2.

#### A graph is a special case

For a height field $z = f(x,y)$ the surface is $\Phi(x,y) = (x, y, f(x,y))$,
and the two tangent vectors are

```math
\Phi_x = \begin{bmatrix}1\\ 0\\ f_x\end{bmatrix},
\qquad
\Phi_y = \begin{bmatrix}0\\ 1\\ f_y\end{bmatrix}.
```

This is the surface built by `Surface.ofHeight` in the code.

**Implementation.** The terrain embedding, and the two graph tangent arrows —
literally the displayed matrices:

```lean
def Surface.ofHeight (h : Height) : Surface :=
  ⟨Smooth2.coordX, Smooth2.coordY, h⟩       -- Φ(x, y) = (x, y, h(x, y))

def Surface.tangentX (h : Height) (x y : Float) : Vec3 := ⟨1, 0, hx h x y⟩
def Surface.tangentY (h : Height) (x y : Float) : Vec3 := ⟨0, 1, hy h x y⟩
```

i.e. $\Phi_x = (1, 0, f_x)$ and $\Phi_y = (0, 1, f_y)$, with $f_x$, $f_y$ the
exact partials of $h$.

**The running example.** Most checked examples in this chapter use the
paraboloid $z = 0.3\,(x^2+y^2)$, built from exact `Smooth2` combinators in
`Main.lean` (§1.7 explains the combinator machinery):

```lean
import Geodesia
open Geodesia

def paraboloid : Smooth2 :=
  Smooth2.smul 0.3
    (Smooth2.add (Smooth2.mul Smooth2.coordX Smooth2.coordX)
      (Smooth2.mul Smooth2.coordY Smooth2.coordY))

def surf : Surface := Surface.ofHeight paraboloid
```

**Checked example.** The four arrows at $(x,y) = (1,0)$ — compare
$f_x = 0.6$, $f_{xx} = 0.6$ below:

```lean
#eval show IO Unit from do
  IO.println s!"Φ    = {Surface.point surf 1 0}"
  IO.println s!"Φ_x  = {Surface.dPhi surf 0 1 0}"
  IO.println s!"Φ_y  = {Surface.dPhi surf 1 1 0}"
  IO.println s!"Φ_xx = {Surface.d2Phi surf 0 0 1 0}"
```

**Output:**
```text
Φ    = (1.000000, 0.000000, 0.300000)
Φ_x  = (1.000000, 0.000000, 0.600000)
Φ_y  = (0.000000, 1.000000, 0.000000)
Φ_xx = (0.000000, 0.000000, 0.600000)
```

### 1.2 The inner product of the tangent space — how to measure on a bent sheet

> **Picture.** Last section gave the tangent space its arrows. A vector
> space with arrows but no way to compare them is a dark room: you can hold
> things, not measure them. What turns $T_pS$ into a place where "length",
> "angle" and "area" mean something is one piece of equipment: an **inner
> product**. And here is the pleasant surprise — we don't have to invent
> one. The tangent space sits inside $\mathbb{R}^3$, and $\mathbb{R}^3$
> already owns the plain dot product $v \cdot w$. The tangent space simply
> *inherits* it.

Concretely: take the two basis arrows $\Phi_u, \Phi_v$ and record the three
dot products among them —

```math
\Phi_u\cdot\Phi_u, \qquad
\Phi_u\cdot\Phi_v, \qquad
\Phi_v\cdot\Phi_v.
```

Stacked into a symmetric $2\times 2$ matrix, this is the *matrix of the
inner product in the basis $\{\Phi_u,\Phi_v\}$* — the **Gram matrix** of two
vectors, the same object you meet as the normal-equation matrix
$J^{\mathsf T}J$ in least squares:

```math
\boxed{\;I = J^{\mathsf T}J
=
\begin{bmatrix}
\Phi_u\cdot\Phi_u & \Phi_u\cdot\Phi_v\\
\Phi_v\cdot\Phi_u & \Phi_v\cdot\Phi_v
\end{bmatrix}\;}
```

> **What the old books call it.** "The *first fundamental form*." The name
> sounds like a sacrament; the object is a Gram matrix. The three entries
> were historically named $E$, $F$, $G$ — letters which tell you nothing
> (the code keeps them as `MetricI.E`, `MetricI.Ff`, `MetricI.G` for
> cross-reading older texts). The content is one sentence: *to measure on
> the surface, dot the ambient arrows.*

**Signature** — `Geodesia/Curvature.lean`:

```lean
structure Mat2 where
  a : Float        -- the symmetric matrix [[a, b], [b, c]]
  b : Float
  c : Float

def Mat2.det (M : Mat2) : Float          -- det I = ac − b²
def Mat2.inv (M : Mat2) : Mat2           -- Cramer's rule (det > 0 here)

def Surface.metricEntry (s : Surface) (u v : Float) : Fin 2 → Fin 2 → Float
    -- I_ij = ⟨∂_iΦ, ∂_jΦ⟩, one dot product
def Surface.metricI (s : Surface) (u v : Float) : Mat2
    -- the whole Gram matrix I = JᵀJ
def Surface.firstForm (s : Surface) (x y : Float) (w : Vec2) : Float
    -- the quadratic form I(w, w)
```

**Implementation.** The symmetric matrix is stored as three `Float`s, and
every entry of the boxed matrix is one `metricEntry` dot product:

```lean
def Mat2.det (M : Mat2) : Float := M.a * M.c - M.b * M.b
    -- det [[a, b], [b, c]] = ac − b²

def Mat2.inv (M : Mat2) : Mat2 :=
  let d := M.det
  ⟨M.c / d, -M.b / d, M.a / d⟩
    -- Cramer's rule: I⁻¹ = (1/det I) [[c, −b], [−b, a]]

def Surface.metricEntry (s : Surface) (u v : Float) : Fin 2 → Fin 2 → Float :=
  fun i j => (dPhi s i u v).dot (dPhi s j u v)

def Surface.metricI (s : Surface) (u v : Float) : Mat2 :=
  ⟨metricEntry s u v 0 0, metricEntry s u v 0 1, metricEntry s u v 1 1⟩
```

**Vector-calculus reading.** `metricEntry s u v i j` is one entry of the
boxed matrix, $I_{ij} = \Phi_i\cdot\Phi_j$; `metricI` stacks the three
independent ones ($b = I_{12} = \Phi_u\cdot\Phi_v = I_{21}$, by symmetry of
the dot product). `Mat2.inv` is the closed-form inverse of a symmetric
$2\times 2$ matrix,

```math
I^{-1} = \frac{1}{\det I}\begin{bmatrix}c & -b\\ -b & a\end{bmatrix},
```

so every use of $I^{-1}$ below (the shape operator, the Christoffel symbols,
the $I^{-1}J^{\mathsf T}w$ solve) is four divisions — no linear-algebra
library anywhere.

#### What the inner product buys you

Three purchases, all of them plain expansions of $v\cdot w$.

**Lengths.** A parameter direction $(a,b)$ becomes the tangent arrow
$J\,(a,b) = a\Phi_u + b\Phi_v$, and its squared length is

```math
\bigl|J\begin{bmatrix}a\\ b\end{bmatrix}\bigr|^2
= \bigl(a\Phi_u + b\Phi_v\bigr)\cdot\bigl(a\Phi_u + b\Phi_v\bigr)
= \begin{bmatrix}a & b\end{bmatrix} I \begin{bmatrix}a\\ b\end{bmatrix}.
```

So $I$ is the **ruler** of the tangent space: it converts the flat
parameter-plane numbers into true lengths on the bent sheet.

**Implementation.** The quadratic form is that expansion, written out entry
by entry:

```lean
def Surface.firstForm (s : Surface) (x y : Float) (w : Vec2) : Float :=
  let I := metricI s x y
  I.a * w.x * w.x + 2 * I.b * w.x * w.y + I.c * w.y * w.y
```

i.e. for $w = (w_x, w_y)$:

```math
I(w, w) \;=\; a\,w_x^2 \,+\, 2\,b\,w_x w_y \,+\, c\,w_y^2
\;=\; w^{\mathsf T} I\, w.
```

**Angles.** Two tangent directions $\xi$, $\eta$ (2-vectors in the parameter
plane) meet at the angle of their arrows $J\xi$, $J\eta$ — and the inner
product computes it without ever forming the arrows:

```math
\cos\theta = \frac{(J\xi)\cdot(J\eta)}{|J\xi|\,|J\eta|}
= \frac{\xi^{\mathsf T} I\, \eta}{\sqrt{\xi^{\mathsf T} I\, \xi}\;\sqrt{\eta^{\mathsf T} I\, \eta}}.
```

**Checked example.** No dedicated function computes $\cos\theta$, but both
routes of the boxed formula are checkable with the library's pieces. Take
$\xi = (1,0)$, $\eta = (1,1)$ at $(1,0)$ on the paraboloid. The **arrow
route** forms $J\xi = \Phi_x$ and $J\eta = \Phi_x + \Phi_y$ and uses plain
dot products; the **ruler route** recovers the bilinear form by polarizing
the quadratic form, $\xi^{\mathsf T}I\eta = \tfrac12\bigl(I(\xi{+}\eta,\xi{+}\eta) - I(\xi,\xi) - I(\eta,\eta)\bigr)$,
using `firstForm` from above:

```lean
#eval show IO Unit from do
  -- the arrows: J ξ = Φ_x, J η = Φ_x + Φ_y for ξ = (1, 0), η = (1, 1)
  let Jxi : Vec3 := Surface.dPhi surf 0 1 0
  let Jeta : Vec3 := Vec3.add (Surface.dPhi surf 0 1 0) (Surface.dPhi surf 1 1 0)
  IO.println s!"cos θ (arrows) = {(Jxi.dot Jeta) / (Jxi.norm * Jeta.norm)}"
  -- the ruler route via polarization of I(w, w) = wᵀ I w
  let q : Vec2 → Float := fun w => Surface.firstForm surf 1 0 w
  let xi : Vec2 := ⟨1, 0⟩
  let eta : Vec2 := ⟨1, 1⟩
  let sum : Vec2 := Vec2.add xi eta
  let bil := (q sum - q xi - q eta) / 2
  IO.println s!"cos θ (ruler)   = {bil / (Float.sqrt (q xi) * Float.sqrt (q eta))}"
  -- the squared-length claim: I(w, w) = |J w|²
  IO.println s!"I((1,0),(1,0))  = {q xi}"
  IO.println s!"|Φ_x|²          = {Jxi.dot Jxi}"
```

**Output:**
```text
cos θ (arrows) = 0.759125
cos θ (ruler)   = 0.759125
I((1,0),(1,0))  = 1.360000
|Φ_x|²          = 1.360000
```

Two claims verified at once: the angle formula (the two routes agree) and
the squared-length formula ($I(w,w) = |Jw|^2$: here
$I((1,0),(1,0)) = 1.36 = |\Phi_x|^2 = 1 + f_x^2$).

**Area.** The two basis arrows span a parallelogram whose area is the norm
of their cross product — and that equals the determinant of the Gram matrix:

```math
|\Phi_u \times \Phi_v|^2 = \det I
\qquad\Longrightarrow\qquad
dA = \sqrt{\det I}\; du\, dv.
```

#### Recovering a 2-D direction from a 3-D arrow

> **By the way.** This subsection is really about *coordinates*: given an
> arrow $w$ that happens to live in $T_pS$, what are its components in the
> basis $\Phi_u, \Phi_v$? The answer needs the inner product — which is why
> it lives here and not in §1.1.

Expand $w$ in the basis:

```math
w = a\,\Phi_u + b\,\Phi_v = J\begin{bmatrix}a\\ b\end{bmatrix}.
```

$J$ is $3\times 2$, not square — but $w$ is known to lie in the span, so a
solution exists and the normal-equations trick finds it: multiply by
$J^{\mathsf T}$ and use $J^{\mathsf T}J = I$:

```math
\boxed{\;\begin{bmatrix}a\\ b\end{bmatrix} = I^{-1} J^{\mathsf T} w\;}
```

This is exactly what `Geodesic.tangentVecToParams` computes — give it a
tangent arrow $w$ and it returns the 2 coordinates whose arrows reproduce it:

```lean
def Geodesic.tangentVecToParams (s : Surface) (u v : Float) (w : Vec3) : Vec2
    -- ξ = I⁻¹ Jᵀ w
```

**Implementation.** Two dot products form $J^{\mathsf T}w$, then the stored
Cramer inverse applies $I^{-1}$:

```lean
def Geodesic.tangentVecToParams (s : Surface) (u v : Float) (w : Vec3) : Vec2 :=
  let gInv := (Surface.metricI s u v).inv
  let a := (Surface.dPhi s 0 u v).dot w        -- b₀ = Φ_u · w
  let b := (Surface.dPhi s 1 u v).dot w        -- b₁ = Φ_v · w
  ⟨gInv.a * a + gInv.b * b, gInv.b * a + gInv.c * b⟩
```

**Vector-calculus reading.** With $b := J^{\mathsf T}w =
\begin{bmatrix}\Phi_u\cdot w\\ \Phi_v\cdot w\end{bmatrix}$, the code returns
$I^{-1}b$ component-wise: $\xi_1 = (I^{-1})_{11}b_0 + (I^{-1})_{12}b_1$,
$\xi_2 = (I^{-1})_{21}b_0 + (I^{-1})_{22}b_1$ — the boxed solve
$\xi = I^{-1}J^{\mathsf T}w$, with no matrix product ever formed.

#### The graph metric again, by dot products

For $\Phi = (x, y, f)$ the three dot products are immediate:

```math
\Phi_x\cdot\Phi_x = 1 + f_x^2,\qquad
\Phi_x\cdot\Phi_y = f_x f_y,\qquad
\Phi_y\cdot\Phi_y = 1 + f_y^2,
```

so

```math
I=
\begin{bmatrix}
1+f_x^2 & f_x f_y\\
f_x f_y & 1+f_y^2
\end{bmatrix}.
```

Cross-check it the vector way: a plane step $(dx,dy)$ lifts to the arrow
$(dx, dy, f_x dx + f_y dy)$, whose squared norm is
$dx^2 + dy^2 + (f_x dx + f_y dy)^2$ — the same quadratic form.

**Checked example.** Three claims verified at once on `surf` at $(1,0)$:
the graph formula ($I_{11} = 1 + 0.6^2 = 1.36$), the determinant identity
$\det I = |\Phi_x\times\Phi_y|^2$ (Exercise 2), and the $I^{-1}J^{\mathsf T}w$
round-trip, fed the non-unit arrow $\Phi_x$ itself:

```lean
#eval show IO Unit from do
  let I := Surface.metricI surf 1 0
  IO.println s!"I            = [[{I.a}, {I.b}], [{I.b}, {I.c}]]"
  IO.println s!"det I        = {I.det}"
  let c := (Surface.dPhi surf 0 1 0).cross (Surface.dPhi surf 1 1 0)
  IO.println s!"|Φ_x × Φ_y|² = {c.dot c}"
  IO.println s!"I⁻¹ Jᵀ Φ_x   = {Geodesic.tangentVecToParams surf 1 0 (Surface.dPhi surf 0 1 0)}"
```

**Output:**
```text
I            = [[1.360000, 0.000000], [0.000000, 1.000000]]
det I        = 1.360000
|Φ_x × Φ_y|² = 1.360000
I⁻¹ Jᵀ Φ_x   = (1.000000, 0.000000)
```

The round-trip returns $(1, 0)$: $\Phi_x$ is indeed the first basis arrow.

### 1.3 The normal, and how acceleration leaves the tangent space

> **Picture.** §1.2 equipped the tangent space with a ruler. But the
> tangent space can't tell you *everything* about the bent sheet it sits
> on — a flat plane has no opinion about bending. The missing information
> lives in one extra arrow, the **unit normal** $n$: the direction that
> points *away* from the flat world, perpendicular to everything in $T_pS$.
> With $T_pS$ and $n$ together, 3-space splits cleanly in two: every ambient
> arrow is *tangent part plus normal part*. The normal is the gatekeeper of
> that split — and, one derivative later, the yardstick of bending: an
> acceleration can move you along the surface freely, but the piece of it
> that leaves the tangent plane must go in the $n$ direction, and $n$ alone
> measures how much.

```math
n = \frac{\Phi_u \times \Phi_v}{|\Phi_u \times \Phi_v|}.
```

Differentiate $\Phi$ once more. Besides the two tangent arrows $\Phi_u$,
$\Phi_v$, the surface now has three **curvature arrows**:

```math
\Phi_{uu} = \frac{\partial^2\Phi}{\partial u^2},\qquad
\Phi_{uv} = \frac{\partial^2\Phi}{\partial u\,\partial v},\qquad
\Phi_{vv} = \frac{\partial^2\Phi}{\partial v^2}.
```

Each is an ordinary 3-D arrow, and each can be tested the same way: *how
much of it points along $n$?* That is one dot product. Three arrows, three
dot products, one symmetric table:

```math
\boxed{\;II =
\begin{bmatrix}
n\cdot\Phi_{uu} & n\cdot\Phi_{uv}\\
n\cdot\Phi_{uv} & n\cdot\Phi_{vv}
\end{bmatrix}\;}
```

> **What the old books call it.** The "*second fundamental form*", entries
> $e$, $f$, $g$. Again the object is simpler than the ceremony: it is the
> table of *normal components of the curvature arrows* — dot each
> $\Phi_{ij}$ with $n$, write it down. (The code's other second-form entry
> names live in `MetricII`.) Note the deep symmetry with §1.2: there,
> tangent arrows dotted with *each other*; here, curvature arrows dotted
> with the *gatekeeper*.

**Signature** — `Geodesia/Surface.lean` and `Geodesia/Curvature.lean`:

```lean
def Surface.normalAt (s : Surface) (u v : Float) : Vec3
    -- n = Φ_u × Φ_v / |Φ_u × Φ_v|
def Surface.secondFormEntry (s : Surface) (x y : Float) : Fin 2 → Fin 2 → Float
    -- II_ij = ⟨n, ∂²_ij Φ⟩, one dot product
def Surface.metricII (s : Surface) (x y : Float) : Mat2
    -- the whole bending matrix II
def Surface.secondForm (s : Surface) (x y : Float) (w : Vec2) : Float
    -- the quadratic form II(w, w)
def Surface.normalCurvature (s : Surface) (x y : Float) (w : Vec2) : Float
    -- κ_n = II(w, w) / I(w, w)
```

**Implementation.** The gatekeeper and the bending table:

```lean
def Surface.normalAt (s : Surface) (u v : Float) : Vec3 :=
  let c := (dPhi s 0 u v).cross (dPhi s 1 u v)
  Vec3.smul (1 / c.norm) c

def Surface.secondFormEntry (s : Surface) (x y : Float) : Fin 2 → Fin 2 → Float :=
  fun i j => (normalAt s x y).dot (d2Phi s i j x y)

def Surface.metricII (s : Surface) (x y : Float) : Mat2 :=
  ⟨secondFormEntry s x y 0 0, secondFormEntry s x y 0 1,
    secondFormEntry s x y 1 1⟩
```

**Vector-calculus reading.** `cross` is the determinant expansion
$u\times v = (u_yv_z - u_zv_y,\ u_zv_x - u_xv_z,\ u_xv_y - u_yv_x)$, so
`normalAt` is the boxed $n = \dfrac{\Phi_u\times\Phi_v}{|\Phi_u\times\Phi_v|}$
(the `smul (1 / c.norm)` scales to unit length). `secondFormEntry s x y i j`
is one entry of the boxed bending table, $II_{ij} = n\cdot\Phi_{ij}$; the
quadratic form mirrors §1.2's `firstForm`:

```lean
def Surface.secondForm (s : Surface) (x y : Float) (w : Vec2) : Float :=
  let II := metricII s x y
  II.a * w.x * w.x + 2 * II.b * w.x * w.y + II.c * w.y * w.y

def Surface.normalCurvature (s : Surface) (x y : Float) (w : Vec2) : Float :=
  secondForm s x y w / firstForm s x y w
```

That is $II(w,w) = e\,w_x^2 + 2f\,w_xw_y + g\,w_y^2$, and the normal
curvature is the boxed ratio in one line,

```math
\kappa(d) \;=\; \frac{II(d,d)}{I(d,d)}
\;=\; \frac{n\cdot\gamma''(s)}{|J\,d|^2}\qquad (d = \gamma'(s)\ \text{in parameters}),
```

so `normalCurvature` **is** the $n\cdot\gamma''(s)$ formula of this section,
normalized — the code never differentiates $\gamma$; it just evaluates the
two quadratic forms.

#### The acceleration test, done by hand

Take the same surface curve $\gamma(s) = \Phi(u(s), v(s)) \in \mathbb{R}^3$
and differentiate it twice with respect to $s$, writing
$u''(s) := \frac{d^2u}{ds^2}$, $v''(s) := \frac{d^2v}{ds^2}$, and
$\gamma''(s) := \frac{d^2\gamma}{ds^2}$ — the **acceleration arrow**:

```math
\gamma''(s) = \frac{d^2\gamma}{ds^2}
= \Phi_u\, u''(s) + \Phi_v\, v''(s)
+ \Phi_{uu}\,u'(s)^2 + 2\,\Phi_{uv}\,u'(s)v'(s) + \Phi_{vv}\,v'(s)^2.
```

Now dot it with $n$. The first two terms die ($n \perp \Phi_u, \Phi_v$ —
the gatekeeper again), and what survives is exactly the bending table
applied to the velocity:

```math
\boxed{\;n\cdot\gamma''(s) = n\cdot\frac{d^2\gamma}{ds^2}
= \begin{bmatrix}u'(s) & v'(s)\end{bmatrix} II \begin{bmatrix}u'(s)\\ v'(s)\end{bmatrix}\;}
```

Let $d := \begin{bmatrix}u'(s)\\ v'(s)\end{bmatrix} \in \mathbb{R}^2$ be the
parameter velocity. So $II(d,d)$ *is* the normal component of the
acceleration when you move in surface direction $d$. Dividing by the
squared speed $|J\,d|^2$ normalizes it — that ratio is the **normal
curvature** `Surface.normalCurvature`:

```math
\kappa(d) = \frac{II(d,d)}{|J\, d|^2}.
```

In plain vector terms it answers: *if I walk across the surface in the
direction of arrow $J d$, how fast does my acceleration arrow bend away
from the tangent plane?* The sign depends on which way $n$ points, but the
product of the two principal curvatures below does not.

![One point of the paraboloid, two principal normal sections: the blue section along Φu and the peach section along Φv bend at different rates κu and κv (arrow lengths proportional to the bending); the polar plot sweeps κ(θ) over all parameter directions, its extremes being the principal curvatures.](images/normal-curvature.png)

**Local height above the tangent plane.** Put a point at the origin, make
the tangent plane the $xy$-plane so $n = (0,0,1)$. Then the third coordinate
of $\Phi_{uu}$ is exactly $n\cdot\Phi_{uu} = II_{11}$ etc., and Taylor's
theorem gives

```math
z(u,v) = \tfrac12\,\bigl(II_{11}u^2 + 2\,II_{12}uv + II_{22}v^2\bigr) + O(3).
```

So $II$ **is** the quadratic approximation of the surface above its tangent
plane: a bowl, a saddle, or flat in one direction — you can read all three
straight off the signs of $II$.

#### Second equivalent formula: derivative of the normal

Differentiate $n \cdot \Phi_u = 0$ with respect to $v$:

```math
\frac{\partial n}{\partial v}\cdot\Phi_u + n\cdot\Phi_{uv} = 0
\qquad\Longrightarrow\qquad
II_{12} = -\frac{\partial n}{\partial v}\cdot\Phi_u.
```

The same holds entry by entry. So $II$ also measures **how fast the normal
arrow turns as you move** — two ways of seeing the same bending. (This is
the identity the shape operator of §1.4 is built on.)

> **What stays paper-only here.** The Taylor development
> $z(u,v) = \tfrac12\,(II_{11}u^2 + 2\,II_{12}uv + II_{22}v^2) + O(3)$ and
> the identity $II_{12} = -(\partial_v n)\cdot\Phi_u$ have no code: the
> library differentiates the *coordinate functions* (`Smooth2`), never
> `normalAt` itself, so there is no stored $\partial_v n$ arrow to dot
> with. Both remain checkable from printed arrows: the Taylor model by
> comparing `Surface.point` against $\Phi$ plus the quadratic model near the
> tangency point, the derivative-of-normal sign by finite-differencing
> `Surface.normalAt`.

**Checked example.** The normal, $II$, and the two normal curvatures along
the coordinate directions — the two extremes of the polar plot above:

```lean
#eval show IO Unit from do
  let II := Surface.metricII surf 1 0
  IO.println s!"II     = [[{II.a}, {II.b}], [{II.b}, {II.c}]]"
  IO.println s!"n      = {Surface.normalAt surf 1 0}"
  IO.println s!"κ(1,0) = {Surface.normalCurvature surf 1 0 ⟨1, 0⟩}"
  IO.println s!"κ(0,1) = {Surface.normalCurvature surf 1 0 ⟨0, 1⟩}"
```

**Output:**
```text
II     = [[0.514496, 0.000000], [0.000000, 0.514496]]
n      = (-0.514496, 0.000000, 0.857493)
κ(1,0) = 0.378306
κ(0,1) = 0.514496
```

Exactly the §1.6 paraboloid arithmetic: $II_{11} = 0.6/\sqrt{1.36} = 0.5145$,
and $\kappa(1,0) = 0.5145/1.36 = 0.3783$.

### 1.4 The shape operator — one linear map on the tangent space

> **Picture.** Two tables now sit on the same flat world: $I$, the ruler
> (§1.2), and $II$, the bending table (§1.3). Comparing matrices is
> bookkeeping; comparing *maps* is geometry. There is one linear map from
> the tangent space to itself that does the comparison in one stroke:

```math
\boxed{\;II = I\,S
\qquad\Longleftrightarrow\qquad
S = I^{-1} II\;}
```

$S$ takes a tangent direction and returns a tangent direction — an
**endomorphism of $T_pS$** — and it encodes exactly one question: *if I
slide in surface direction $\xi$, how fast and in which way does the normal
arrow turn?* Concretely, $-S\xi$ is the rate of change of $n$ along $\xi$:

```math
-\,S\,\begin{bmatrix}a\\ b\end{bmatrix}
= a\,\frac{\partial n}{\partial u} + b\,\frac{\partial n}{\partial v},
\qquad \bigl(a,b\bigr)\ \text{a parameter direction}.
```

![Sliding along the two principal directions ξ1 (red arrows) and ξ2 (teal arrows), the sampled unit normals tip at different rates; the dashed tracks of the arrow tips are exactly the turns that −Sξ measures.](images/shape-operator.png)

$S$ is a $2\times 2$ map acting on the 2-D tangent space, so it has two
eigenvalues, the **principal curvatures** $\kappa_1, \kappa_2$, and two
eigendirections — the directions of greatest and least bending. Its two
basic invariants get the historical names:

```math
K = \det S = \frac{\det II}{\det I}
\qquad\text{(Gaussian curvature)},
```

```math
H = \tfrac12\,\mathrm{tr}\,S = \tfrac12\,\mathrm{tr}\!\bigl(I^{-1}II\bigr)
\qquad\text{(mean curvature)}.
```

The eigenvalues come back from the trace/determinant pair:

```math
\kappa_{1,2} = H \pm \sqrt{H^2 - K}.
```

- $K>0$: both directions bend the same way (bowl, ellipsoid).
- $K<0$: directions bend oppositely (saddle).
- $K=0$: at least one direction is straight (plane, cylinder, cone).
- $H=0$: the two bendings cancel (soap films).

![The three sign regimes of the Gaussian curvature K: a bowl bends the same way in both principal directions (K>0), a saddle bends oppositely (K<0), and a cylinder has one straight direction (K=0); each panel shows the unit normal n and the two principal directions ξ1, ξ2 at one point, with the κ1, κ2, K readout.](images/curvature-signatures.png)

> **What the old books call it.** The map $S$ is usually introduced as
> $-Dn$, "*the Weingarten map*", and the identity $II = IS$ is buried under
> the name. The code states it the honest way: $S$ is defined as the map
> that makes $II = IS$ true, and the normal-turning picture is a theorem
> about it, not its definition. (Which came first is a matter of taste;
> that they agree is Exercise 7.)

**Signature** — `Geodesia/Curvature.lean`:

```lean
def Surface.shapeOperator (s : Surface) (x y : Float) (w : Vec2) : Vec2
    -- S w = I⁻¹ (II w); its eigenvalues are the principal curvatures
def Surface.gaussK (s : Surface) (x y : Float) : Float
    -- K = det II / det I — intrinsic
def Surface.meanH (s : Surface) (x y : Float) : Float
    -- H = (eG − 2fF + gE) / 2(EG − F²)
def Surface.principalK (s : Surface) (x y : Float) : Float × Float
    -- (κ_min, κ_max) = H ± √(H² − K)
```

One wrinkle worth knowing: `Surface.principalK` computes
$\kappa_{1,2} = H \pm \sqrt{H^2-K}$ rather than calling an eigensolver —
with a $2\times 2$ map the eigenvalues come from trace and determinant,
which is all `Mat2` needs.

**Implementation.** All four curvature numbers are closed-form $2\times 2$
algebra — no eigensolver, no SVD:

```lean
def Surface.shapeOperator (s : Surface) (x y : Float) (w : Vec2) : Vec2 :=
  let II := metricII s x y
  let gInv := (metricI s x y).inv
  let b1 := II.a * w.x + II.b * w.y
  let b2 := II.b * w.x + II.c * w.y
  ⟨gInv.a * b1 + gInv.b * b2, gInv.b * b1 + gInv.c * b2⟩

def Surface.gaussK (s : Surface) (x y : Float) : Float :=
  (metricII s x y).det / (metricI s x y).det

def Surface.meanH (s : Surface) (x y : Float) : Float :=
  let E := MetricI.E s x y
  let Ff := MetricI.Ff s x y
  let G := MetricI.G s x y
  let II := metricII s x y
  let denom := 2 * (E * G - Ff * Ff)
  (II.a * G - 2 * II.b * Ff + II.c * E) / denom

def Surface.principalK (s : Surface) (x y : Float) : Float × Float :=
  let H := meanH s x y
  let K := gaussK s x y
  let disc0 := H * H - K
  let disc := Float.sqrt (if disc0 < 0 then 0 else disc0)
  (H - disc, H + disc)
```

**Vector-calculus reading.**

- `shapeOperator` computes $Sw = I^{-1}(II\,w)$ in the two natural steps:
  first the bending table acts on $w$,
  $b = II\,w = \begin{bmatrix}e\,w_x + f\,w_y\\ f\,w_x + g\,w_y\end{bmatrix}$
  (`b1`, `b2`), then the ruler inverse acts, $Sw = I^{-1}b$ (the same
  Cramer expansion as §1.2). The defining property $II = IS$ is then
  automatic: $I\,(Sw) = I\,I^{-1}(II\,w) = II\,w$, which is Exercise 7.
- `gaussK` is the boxed ratio $K = \det II/\det I$ — one line.
- `meanH` looks like a different formula because it *is* the trace,
  expanded. From
  $I^{-1} = \dfrac{1}{EG-F^2}\begin{bmatrix}G & -F\\ -F & E\end{bmatrix}$:

  ```math
  \mathrm{tr}\,\bigl(I^{-1}II\bigr)
  = \frac{G\,e \,-\, 2F\,f \,+\, E\,g}{EG - F^2}
  \qquad\Longrightarrow\qquad
  H = \frac{eG - 2fF + gE}{2\,(EG - F^2)},
  ```

  which is literally the code line — the code never builds the matrix
  $I^{-1}II$, only its trace.
- `principalK` is the boxed $\kappa_{1,2} = H \pm \sqrt{H^2 - K}$; the `if`
  clamps a tiny negative discriminant (a `Float` rounding artifact) to zero.

**Checked example.** $S$ applied to the two coordinate directions, then
$K$, $H$ and the eigenvalue pair — note $S(1,0)$ and $S(0,1)$ reproduce
`normalCurvature` from §1.3 (here $S$ is diagonal, so $\Phi_x$, $\Phi_y$
are the eigendirections):

```lean
#eval show IO Unit from do
  IO.println s!"S(1,0)       = {Surface.shapeOperator surf 1 0 ⟨1, 0⟩}"
  IO.println s!"S(0,1)       = {Surface.shapeOperator surf 1 0 ⟨0, 1⟩}"
  IO.println s!"K, H         = {Surface.gaussK surf 1 0}, {Surface.meanH surf 1 0}"
  let (kmin, kmax) := Surface.principalK surf 1 0
  IO.println s!"κ_min, κ_max = {kmin}, {kmax}"
```

**Output:**
```text
S(1,0)       = (0.378306, 0.000000)
S(0,1)       = (0.000000, 0.514496)
K, H         = 0.194637, 0.446401
κ_min, κ_max = 0.378306, 0.514496
```

### 1.5 The split of an arrow — tangent part and normal part

> **Picture.** Put §1.1 and §1.3 together: at a point, 3-space is exactly
> the flat world $T_pS$ plus one extra direction, $n$. Two flat pieces, one
> arrow to be placed — so every ambient arrow $w$ splits **uniquely** into a
> piece that stays in the flat world and a piece that leaves it:

```math
\boxed{\;
\begin{aligned}
T &: \mathbb{R}^3 \to \mathbb{R}^3, & T(w) &:= w - (w\cdot n)\,n, \qquad w \in \mathbb{R}^3,\\[6pt]
&& \Phi_{ij} &= T(\Phi_{ij}) + II_{ij}\,n, \qquad i,j \in \{u,v\}.
\end{aligned}\;}
```

$T$ — the **tangential projection** — subtracts exactly the normal
component: $T(w)$ is the **tangent part** of $w$, $(w\cdot n)\,n$ its
**normal part**. Applied to a curvature arrow $\Phi_{ij}$, the second line
is the classical **Gauss split**: the tangent part $T(\Phi_{ij})$ is the
turning data, the normal coefficient $II_{ij}$ the bending data (§1.3).
Expanding $T(\Phi_{ij})$ in the basis $\{\Phi_u,\Phi_v\}$ yields the
coefficients textbooks call Christoffel symbols; Part 2 works with the
arrows directly.

**Signature.** The split now has named declarations in the library —
this is the vector machinery of the chapter, made executable:

```lean
def Surface.tangentialPart (s : Surface) (u v : Float) (w : Vec3) : Vec3
    -- T w = w − (w·n) n, the tangent part of the split
def Surface.normalPart (s : Surface) (u v : Float) (w : Vec3) : Vec3
    -- (w·n) n, the normal part of the split

Vec3.sub w (Vec3.smul (w.dot n) n)   -- T w = w − (w·n) n, in bare primitives
```

**Implementation.** Each part is one `Vec3` line, verbatim from
`Geodesia/Surface.lean` (the `namespace Surface` supplies the prefix):

```lean
def normalPart (s : Surface) (u v : Float) (w : Vec3) : Vec3 :=
  Vec3.smul (w.dot (normalAt s u v)) (normalAt s u v)

def tangentialPart (s : Surface) (u v : Float) (w : Vec3) : Vec3 :=
  Vec3.sub w (normalPart s u v w)
```

**Vector-calculus reading.** The three `Vec3` primitives underneath
(`Geodesia/Vec.lean`) are the textbook definitions — `dot` is the
coordinate expansion of $u\cdot v = u_xv_x + u_yv_y + u_zv_z$, `smul` is
scalar multiplication $a\,u$, `sub` is the component-wise difference
$u - v$:

```lean
def Vec3.dot (u v : Vec3) : Float := u.x * v.x + u.y * v.y + u.z * v.z
def Vec3.smul (a : Float) (u : Vec3) : Vec3 := ⟨a * u.x, a * u.y, a * u.z⟩
def Vec3.sub (u v : Vec3) : Vec3 := ⟨u.x - v.x, u.y - v.y, u.z - v.z⟩
```

So `normalPart s u v w` is $(w\cdot n)\,n$ coordinate by coordinate, and
`tangentialPart s u v w` is the boxed $T(w) = w - (w\cdot n)\,n$:

```math
T(w) = w - (w\cdot n)\,n
= \begin{bmatrix}w_x\\ w_y\\ w_z\end{bmatrix}
- (w_xn_x + w_yn_y + w_zn_z)\begin{bmatrix}n_x\\ n_y\\ n_z\end{bmatrix},
```

It lands in the tangent plane by construction:

```math
T(w)\cdot n = w\cdot n - (w\cdot n)\,(n\cdot n) = 0
\qquad (\text{since } |n| = 1).
```

No test is needed that the result stays tangent — the algebra forbids
anything else. The same three primitives build `Surface.normalAt` (§1.3,
where `smul` rescales the cross product to unit length) and the orthonormal
shooting frame `Geodesic.tangentFrame` (a Gram–Schmidt step).

**Checked example.** Split the arrow $w = (0, 0, 0.6)$ — the paraboloid's
$\Phi_{vv}$ at $(1,0)$ — using the named parts, and check the reconstruction:

```lean
#eval show IO Unit from do
  let w : Vec3 := ⟨0, 0, 0.6⟩
  IO.println s!"w·n         = {w.dot (Surface.normalAt surf 1 0)}"
  IO.println s!"T(w)        = {Surface.tangentialPart surf 1 0 w}"
  IO.println s!"(w·n)n      = {Surface.normalPart surf 1 0 w}"
  IO.println s!"T(w) + (w·n)n = {Vec3.add (Surface.tangentialPart surf 1 0 w) (Surface.normalPart surf 1 0 w)}"
```

**Output:**
```text
w·n         = 0.514496
T(w)        = (0.264706, 0.000000, 0.158824)
(w·n)n      = (-0.264706, 0.000000, 0.441176)
T(w) + (w·n)n = (0.000000, 0.000000, 0.600000)
```

The last line returns $w$ exactly: tangent part plus normal part *is* the
arrow — the split is lossless, as it must be. The normal coefficient
$0.514496$ is the $II_{22}$ entry of §1.3, and the $\Phi_u$-coefficient
$0.264706$ of the tangent part is a Christoffel symbol — §2.7 reads this
same arrow out of the code.

**Checked example (in raw primitives).** The same split written with bare
`Vec3` lines — what you would compute before the named declarations existed:

```lean
#eval show IO Unit from do
  let n := Surface.normalAt surf 1 0
  let w : Vec3 := ⟨0, 0, 0.6⟩
  let Tw := Vec3.sub w (Vec3.smul (w.dot n) n)
  IO.println s!"w·n  = {w.dot n}"
  IO.println s!"T(w) = {Tw}"
```

**Output:**
```text
w·n  = 0.514496
T(w) = (0.264706, 0.000000, 0.158824)
```

> **One sentence summary.** $J=[\Phi_u\;\Phi_v]$, $I=J^{\mathsf T}J$ (ruler),
> $II_{ij}=n\cdot\Phi_{ij}$ (bending), and every arrow — curvature arrows
> included — splits as $w = T(w) + (w\cdot n)\,n$, so
> $\Phi_{ij}=T(\Phi_{ij})+II_{ij}n$ delivers the turning data $T(\Phi_{ij})$
> in the same sweep.

---

### 1.6 Worked examples as vector computations

Four classical surfaces, all measured with the same three moves: tangent
arrows → dot products ($I$), curvature arrows → dot products with $n$
($II$), and $S = I^{-1}II$.

#### Plane

```math
\Phi(u,v) = (u,v,0):\quad
\Phi_u = \begin{bmatrix}1\\0\\0\end{bmatrix},\quad
\Phi_v = \begin{bmatrix}0\\1\\0\end{bmatrix},\quad
\Phi_{uu}=\Phi_{uv}=\Phi_{vv}=0.
```

Then $I = J^{\mathsf T}J$ is the identity and every dot product with the
(curvature) vectors vanishes, so $II = 0$, $K=H=0$. Flatness is the
statement "all second derivatives are zero vectors".

#### Cylinder of radius $R$

```math
\Phi(\theta,z) = (R\cos\theta,\,R\sin\theta,\,z).
```

Tangent arrows:

```math
\Phi_\theta = \begin{bmatrix}-R\sin\theta\\ R\cos\theta\\ 0\end{bmatrix},
\qquad
\Phi_z = \begin{bmatrix}0\\ 0\\ 1\end{bmatrix},
```

so $I = \begin{bmatrix}R^2 & 0\\ 0 & 1\end{bmatrix}$. The single nonzero
curvature arrow is

```math
\Phi_{\theta\theta} = \begin{bmatrix}-R\cos\theta\\ -R\sin\theta\\ 0\end{bmatrix}
= -R\,\frac{\Phi_\theta\times\Phi_z}{|\Phi_\theta\times\Phi_z|},
```

which is parallel to the normal — so $n\cdot\Phi_{\theta\theta} = -R$ and
$n\cdot\Phi_{zz} = 0$, giving

```math
II = \begin{bmatrix}-R & 0\\ 0 & 0\end{bmatrix}.
```

Hence one direction bends ($-1/R$) and one is straight ($0$): $K = 0$. That
is why a cylinder can be cut and unrolled flat — nothing forces it to bend.

#### Sphere of radius $R$

```math
\Phi(\theta,\varphi)
= R\,(\sin\theta\cos\varphi,\ \sin\theta\sin\varphi,\ \cos\theta).
```

Tangent arrows have lengths $|R\Phi_\theta| = R$ and
$|\Phi_\varphi| = R\sin\theta$, perpendicular to each other, so

```math
I = \begin{bmatrix}R^2 & 0\\ 0 & R^2\sin^2\theta\end{bmatrix}.
```

The normal here is simply $n = \Phi/R$ (the arrow from the center to the
point). The curvature arrows all point along $-n$:

```math
\Phi_{\theta\theta} = -R\,n,
\qquad
\Phi_{\varphi\varphi} = -R\sin^2\theta\,n,
```

so

```math
II = \begin{bmatrix}-R & 0\\ 0 & -R\sin^2\theta\end{bmatrix}.
```

Then $S = I^{-1}II = -(1/R)$ times identity, so
$\kappa_1 = \kappa_2 = -1/R$, $K = 1/R^2$, $H = -1/R$. The sphere bends
equally in every direction; only the *sign* of the first order bending
depends on choosing the outward $n$. (The code's unit-sphere demo checks
exactly this: $K = 1$.)

#### Paraboloid graph $z = 0.3(x^2+y^2)$ at $(1,0)$

Here $f_x = 0.6$, $f_y = 0$, $f_{xx} = 0.6$, $f_{yy}=0.6$, $f_{xy}=0$, so

```math
\Phi_x = \begin{bmatrix}1\\0\\0.6\end{bmatrix},
\qquad
\Phi_y = \begin{bmatrix}0\\1\\0\end{bmatrix},
\qquad
n = \frac{1}{\sqrt{1.36}}\begin{bmatrix}-0.6\\0\\1\end{bmatrix}.
```

Dot products:

```math
I
= \begin{bmatrix}1.36 & 0\\ 0 & 1\end{bmatrix},
\qquad
II
= \frac{1}{\sqrt{1.36}}\begin{bmatrix}0.6 & 0\\ 0 & 0.6\end{bmatrix}
= \begin{bmatrix}0.5145&0\\ 0&0.5145\end{bmatrix}.
```

Notice $II$ is a multiple of the identity while $I$ is not, so

```math
S = I^{-1}II
=
\begin{bmatrix}
0.3783 & 0\\
0 & 0.5145
\end{bmatrix}.
```

The two bending rates differ — the **metric** decides which direction is
more curved, not $II$ alone. Then

```math
K = \frac{\det II}{\det I} = \frac{0.5145^2}{1.36} = 0.1946,
\qquad
H = \tfrac12(0.3783 + 0.5145) = 0.4464.
```

Every number in this example is a dot product of four 3-D arrows — exactly
what `lake exe geodesia` prints:

```text
Geodesia — exact differential geometry (no finite differences)
paraboloid z = 0.3(x² + y²) at (x, y) = (1.000000, 0.000000)
  Φ            = (1.000000, 0.000000, 0.300000)
  I            = [[1.360000, 0.000000], [0.000000, 1.000000]]
  Γ¹₁₁, Γ¹₁₂, Γ¹₂₂ = 0.264706, 0.000000, 0.264706
  Γ²₁₁, Γ²₁₂, Γ²₂₂ = 0.000000, 0.000000, 0.000000
  K, H          = 0.194637, 0.446401
  κ_min, κ_max  = 0.378306, 0.514496
geodesic from (1, 0) with parameter velocity (1, 0), s ∈ [0, 1]
  endpoint      = (1.879182, 0.000000, 1.059398)
  E(0), E(1)    = 0.680000, 0.680000
  shoot(angle=0.6, length=1) = (1.606884, 0.550855, 0.865655)
unit sphere at (θ, φ) = (1.000000, 0.200000)
  K (should be 1) = 1.000000
```

<!-- §APPEND -->
