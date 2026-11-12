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

<!-- §APPEND -->
