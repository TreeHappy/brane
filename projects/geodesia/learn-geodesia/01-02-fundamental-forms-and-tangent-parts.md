# Geodesia — the two fundamental forms and the tangent parts (Ch. 1–2)

> **What this chapter is.** One merged chapter, written in the style of the Lean
> documentation (as in the Verso-built Lean Manual and *Theorem Proving in Lean*):
> every boxed formula is followed by the **Lean declarations that compute it** — a
> **Signature** block — and every numeric claim is a **Checked example**: the Lean
> snippet together with the output it really prints. Nothing in an `Output:` block
> is hand-computed; each was produced by running the snippet against the library.
> The math notation is kept in full ($\Phi$, $I$, $II$, $T$, $\Gamma$, …); the Lean
> identifiers are the computational shadows of the same objects.
>
> It merges the former `01-fundamental-forms.md` and `02-tangent-parts.md`. The
> section numbers of both chapters are preserved (`§1.1–1.8`, `§2.1–2.9`), so the
> cross-references from Chapters 3–8, `index.md` and `GA_VIEW.md` still resolve.

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
| $I = J^{\mathsf T}J$, the ruler | `Surface.metricI`, `Surface.firstForm` | `Geodesia/Curvature.lean` |
| $II_{ij} = n\cdot\Phi_{ij}$, the bending | `Surface.metricII`, `Surface.secondForm` | `Geodesia/Curvature.lean` |
| $\kappa(d) = II(d,d)/I(d,d)$ | `Surface.normalCurvature` | `Geodesia/Curvature.lean` |
| $S = I^{-1}II$ | `Surface.shapeOperator` | `Geodesia/Curvature.lean` |
| $K$, $H$, $\kappa_{1,2}$ | `Surface.gaussK`, `Surface.meanH`, `Surface.principalK` | `Geodesia/Curvature.lean` |
| $T w = w - (w\cdot n)\,n$ | `Vec3.sub w (Vec3.smul (w.dot n) n)` | `Geodesia/Vec.lean` |
| $\Gamma^k_{\;ij}$ (coordinates of $T\Phi_{ij}$) | `Surface.christoffel` | `Geodesia/Curvature.lean` |
| $\partial_u\,(\Phi_j\cdot\Phi_k)$ | `Surface.dMetricEntry` | `Geodesia/Curvature.lean` |
| $\xi = I^{-1}J^{\mathsf T}w$ | `Geodesic.tangentVecToParams` | `Geodesia/Geodesic.lean` |
| exact derivative arrows | `Smooth2` + its combinators | `Geodesia/Diff.lean` |

![The paraboloid with its parameter grid; at one point the tangent plane (dashed outline) carrying the tangent arrows Φu, Φv and the normal n, plus a curve γ(s) whose second derivative γ″ splits into a tangent-plane part Tγ″ and a normal part II(d,d)n.](images/fundamental-forms.png)

---

## Part 1 — The two fundamental forms: dot products of tangent vectors

> **Goal.** Understand exactly what `metricI`, `metricII`, `gaussK`, `meanH`
> and `principalK` in `Curvature.lean` are computing. Everything here is
> ordinary vector math: a surface is a function, its derivatives are vectors,
> and the two forms are just dot products between those vectors.

### 1.1 The surface and its derivative matrix

A surface is a function

```math
\Phi : \mathbb{R}^2 \to \mathbb{R}^3,
\qquad \Phi(u,v) = \big(x(u,v),\ y(u,v),\ z(u,v)\big).
```

You feed it a point of the parameter plane and it returns a point of 3-space.
The two numbers $(u,v)$ are just coordinates — like row and column indices
painted on the surface.

Differentiate $\Phi$ with respect to each of its two inputs. Each derivative
is a genuine 3-D arrow (a "world-space" vector):

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

These two columns are the whole story of Part 1. Wherever you stand on
the surface, they are the two arrows that span the **tangent plane** at that
point. Everything else is dot products between them, between their
derivatives, and between their derivatives and a normal.

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
J\begin{bmatrix}a\\ b\end{bmatrix} = a\,\Phi_u + b\,\Phi_v.
```

That is all "parameter space" means here: the numbers $(a,b)$ live in the
flat plane, and $J$ carries them onto the surface. No other machinery.

#### A graph is a special case

For a height field $z = f(x,y)$ the surface is $\Phi(x,y) = (x, y, f(x,y))$,
and the two tangent vectors are

```math
\Phi_x = \begin{bmatrix}1\\ 0\\ f_x\end{bmatrix},
\qquad
\Phi_y = \begin{bmatrix}0\\ 1\\ f_y\end{bmatrix}.
```

This is the surface built by `Surface.ofHeight` in the code.

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

### 1.2 The first fundamental form $I$ — the Gram matrix of the tangent vectors

Given two tangent vectors $\Phi_u, \Phi_v$, the most basic vector quantity
you can form from them is their collection of dot products:

```math
\Phi_u\cdot\Phi_u, \qquad
\Phi_u\cdot\Phi_v, \qquad
\Phi_v\cdot\Phi_v.
```

Stack them into a symmetric $2\times 2$ matrix and you have the **first
fundamental form**:

```math
\boxed{\;I = J^{\mathsf T}J
=
\begin{bmatrix}
\Phi_u\cdot\Phi_u & \Phi_u\cdot\Phi_v\\
\Phi_v\cdot\Phi_u & \Phi_v\cdot\Phi_v
\end{bmatrix}\;}
```

It is simply the **Gram matrix** of the two tangent vectors — the same thing
you compute when you build a least-squares normal-equation matrix $J^{\mathsf T}J$.
(Older books name the three entries $E$, $F$, $G$; the code keeps them as
`MetricI.E`, `MetricI.Ff`, `MetricI.G`, but the entries are just dot products.)

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

#### What the dot products buy you

**Squared length on the surface.** A direction $(a,b)$ in the parameter
plane becomes the 3-D arrow $a\Phi_u + b\Phi_v$, and its squared length is

```math
\bigl|J\begin{bmatrix}a\\ b\end{bmatrix}\bigr|^2
= \bigl(a\Phi_u + b\Phi_v\bigr)\cdot\bigl(a\Phi_u + b\Phi_v\bigr)
= \begin{bmatrix}a & b\end{bmatrix} I \begin{bmatrix}a\\ b\end{bmatrix}.
```

so $I$ is the **ruler**: it turns parameter-plane directions into true
lengths on the bent surface. Notice how the formula is pure expansion of a
dot product.

**Angles.** Two surface directions $\xi$, $\eta$ (both 2-vectors) meet at
the angle of their arrows $J\xi$, $J\eta$:

```math
\cos\theta = \frac{(J\xi)\cdot(J\eta)}{|J\xi|\,|J\eta|}
= \frac{\xi^{\mathsf T} I\, \eta}{\sqrt{\xi^{\mathsf T} I\, \xi}\;\sqrt{\eta^{\mathsf T} I\, \eta}}.
```

**Area.** The two arrows $\Phi_u$, $\Phi_v$ span a parallelogram whose area
is the norm of their cross product, and that equals the determinant of $I$:

```math
|\Phi_u \times \Phi_v|^2 = \det I
\qquad\Longrightarrow\qquad
dA = \sqrt{\det I}\; du\, dv.
```

#### Recovering a 2-D direction from a 3-D arrow

Suppose you are handed a world-space tangent arrow $w$ (3 components) and
want its coordinates $(a,b)$ in the tangent basis. Expand:

```math
w = a\,\Phi_u + b\,\Phi_v = J\begin{bmatrix}a\\ b\end{bmatrix}.
```

$J$ is not square, so solve it by multiplying both sides by $J^{\mathsf T}$ —
the normal-equations trick. Using $J^{\mathsf T}J = I$:

```math
\boxed{\;\begin{bmatrix}a\\ b\end{bmatrix} = I^{-1} J^{\mathsf T} w\;}
```

This is exactly what `Geodesic.tangentVecToParams` computes — give it a
tangent arrow $w$ and it returns the 2 numbers whose arrows reproduce it:

```lean
def Geodesic.tangentVecToParams (s : Surface) (u v : Float) (w : Vec3) : Vec2
    -- ξ = I⁻¹ Jᵀ w
```

#### The graph metric again, by dot products

For $\Phi = (x, y, f)$ the dot products are immediate:

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

### 1.3 The second fundamental form $II$ — dot products with the normal

The tangent vectors describe the surface **locally flat**. To describe how
it bends we need, besides the tangent plane, the **unit normal** arrow:

```math
n = \frac{\Phi_u \times \Phi_v}{|\Phi_u \times \Phi_v|}.
```

Now differentiate once more. Besides the two tangent vectors $\Phi_u$,
$\Phi_v$, the surface has three **curvature vectors**:

```math
\Phi_{uu} = \frac{\partial^2\Phi}{\partial u^2},\qquad
\Phi_{uv} = \frac{\partial^2\Phi}{\partial u\,\partial v},\qquad
\Phi_{vv} = \frac{\partial^2\Phi}{\partial v^2}.
```

Each is a 3-D arrow. The part of such an arrow that bends away from the
tangent plane is its component along the normal, that is, its dot product
with $n$. Collect those three dot products into a matrix — the **second
fundamental form**:

```math
\boxed{\;II =
\begin{bmatrix}
n\cdot\Phi_{uu} & n\cdot\Phi_{uv}\\
n\cdot\Phi_{uv} & n\cdot\Phi_{vv}
\end{bmatrix}\;}
```

(Again, old books call the entries $e$, $f$, $g$; the entries are just
"curvature vector dotted with the normal", which is easier to remember.)

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

#### What the normal dot products buy you

**The bending of a curve.** Take the same surface curve
$\gamma(s) = \Phi(u(s), v(s)) \in \mathbb{R}^3$ and differentiate it twice
with respect to $s$, writing $u''(s) := \frac{d^2u}{ds^2}$,
$v''(s) := \frac{d^2v}{ds^2}$, and
$\gamma''(s) := \frac{d^2\gamma}{ds^2}$:

```math
\gamma''(s) = \frac{d^2\gamma}{ds^2}
= \Phi_u\, u''(s) + \Phi_v\, v''(s)
+ \Phi_{uu}\,u'(s)^2 + 2\,\Phi_{uv}\,u'(s)v'(s) + \Phi_{vv}\,v'(s)^2.
```

Dot this with the normal $n$. The $\Phi_u$, $\Phi_v$ terms die because
$n\perp\Phi_u,\Phi_v$, and what survives is

```math
\boxed{\;n\cdot\gamma''(s) = n\cdot\frac{d^2\gamma}{ds^2}
= \begin{bmatrix}u'(s) & v'(s)\end{bmatrix} II \begin{bmatrix}u'(s)\\ v'(s)\end{bmatrix}\;}
```

a quadratic form again. Let $d := \begin{bmatrix}u'(s)\\ v'(s)\end{bmatrix} \in \mathbb{R}^2$
be the parameter-plane velocity (a surface direction). So $II(d,d)$ is the
normal component of the second derivative when you move in surface direction
$d$. Dividing by the squared speed $|J\, d|^2$ normalizes it: that ratio is
the **normal curvature** `Surface.normalCurvature`:

```math
\kappa(d) = \frac{II(d,d)}{|J\, d|^2}.
```

It answers, in plain vector terms: *if I walk across the surface in the
direction of arrow $J d$, how fast does my arrow bend away from the tangent
plane?* The sign depends on which way $n$ points, but the product of the two
principal curvatures below does not.

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
arrow turns as you move** — two ways of seeing the same bending.

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

### 1.4 The shape operator, principal curvatures, $K$ and $H$

Both forms are matrices on the same 2-D tangent space: $I$ is the metric
(Gram) matrix, $II$ is the bending matrix. To compare them directly, relate
both to one linear map $S$:

```math
\boxed{\;II = I\,S
\qquad\Longleftrightarrow\qquad
S = I^{-1} II\;}
```

$S$ is the **shape operator**. Geometrically, $S$ encodes: "if I move in
surface direction $\xi$, how fast and in which way does the normal arrow
turn?" Concretely, negative $S\xi$ is the rate of change of the normal as
you slide in direction $\xi$:

```math
-\,S\,\begin{bmatrix}a\\ b\end{bmatrix}
= a\,\frac{\partial n}{\partial u} + b\,\frac{\partial n}{\partial v},
\qquad \bigl(a,b\bigr)\ \text{a parameter direction}.
```

![Sliding along the two principal directions ξ1 (red arrows) and ξ2 (teal arrows), the sampled unit normals tip at different rates; the dashed tracks of the arrow tips are exactly the turns that −Sξ measures.](images/shape-operator.png)

$S$ is a $2\times 2$ matrix acting on tangent directions, so it has two
eigenvalues, the **principal curvatures** $\kappa_1, \kappa_2$, and two
eigendirections — the directions of greatest and least bending. And its two
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
with a $2\times 2$ matrix the eigenvalues come from trace and determinant,
which is all `Mat2` needs.

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

### 1.5 The Gauss split

Define the tangential projection $T$ and recall $II_{ij} := n\cdot \Phi_{ij}$:

```math
\boxed{\;
\begin{aligned}
T &: \mathbb{R}^3 \to \mathbb{R}^3, & T(w) &:= w - (w\cdot n)\,n, \qquad w \in \mathbb{R}^3,\\[6pt]
&& \Phi_{ij} &= T(\Phi_{ij}) + II_{ij}\,n, \qquad i,j \in \{u,v\}.
\end{aligned}\;}
```

The second line is the **Gauss split**: $T(\Phi_{ij})$ is the tangent
(turning) part and $II_{ij}\,n$ the normal (bending) part. Expanding
$T(\Phi_{ij})$ in the basis $\{\Phi_u,\Phi_v\}$ yields the coefficients
textbooks call Christoffel symbols; Part 2 works with the arrows $T(\Phi_{ij})$ directly.

**Signature.** There is no dedicated `T` function in the library — the
projection is one subtraction with `Geodesia/Vec.lean` primitives:

```lean
Vec3.sub w (Vec3.smul (w.dot n) n)   -- T w = w − (w·n) n
```

**Checked example.** Take the paraboloid's second-derivative arrow in the
$v$-direction, $\Phi_{vv} = (0, 0, 0.6)$, and project it:

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

So $\Phi_{vv} = (0.264706,\,0,\,0.158824) + 0.514496\,n$: one arrow, two
projections. The normal coefficient $0.514496$ is the $II_{22}$ entry of
§1.3, and the $\Phi_u$-coefficient $0.264706$ of the tangent part is a
Christoffel symbol — §2.7 reads this same arrow out of the code.

> **One sentence summary.** $J=[\Phi_u\;\Phi_v]$, $I=J^{\mathsf T}J$ (ruler),
> $II_{ij}=n\cdot\Phi_{ij}$ (bending), and $\Phi_{ij}=T(\Phi_{ij})+II_{ij}n$
> delivers the turning data $T(\Phi_{ij})$ in the same sweep.

---

### 1.6 Worked examples as vector computations

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

### 1.7 How the code represents this (exact derivatives)

`Geodesia/Diff.lean` defines

```lean
structure Smooth2 where
  f   : Float → Float → Float
  fx  : Float → Float → Float
  fy  : Float → Float → Float
  fxx : Float → Float → Float
  fxy : Float → Float → Float
  fyy : Float → Float → Float
```

Each derivative is a **function**, so there is no finite-difference step in
the main path — the derivative arrows above are exact. You build a surface
from exact combinators:

```lean
def paraboloid : Smooth2 :=
  Smooth2.smul 0.3
    (Smooth2.add (Smooth2.mul Smooth2.coordX Smooth2.coordX)
                 (Smooth2.mul Smooth2.coordY Smooth2.coordY))
```

`Smooth2.mul` implements the product rule in both variables, `Smooth2.sin`
implements the chain rule twice, etc. If you only have a bare
`Float → Float → Float`, use `Smooth2.ofFiniteDiff` as a fallback; that is
the only place an `ε` appears.

A general (non-graph) surface is just three `Smooth2` coordinate functions —
that is how the sphere and the polar plane of Part 2 are built:

```lean
def sphere : Surface where
  sx := Smooth2.mul (Smooth2.sin Smooth2.coordX) (Smooth2.cos Smooth2.coordY)
  sy := Smooth2.mul (Smooth2.sin Smooth2.coordX) (Smooth2.sin Smooth2.coordY)
  sz := Smooth2.cos Smooth2.coordX
```

Then, in `Surface.lean` / `Curvature.lean`, the whole toolbox of this part:

```lean
Surface.dPhi s 0 x y         -- the tangent arrow Φ_u
Surface.dPhi s 1 x y         -- the tangent arrow Φ_v
Surface.d2Phi s 0 0 x y      -- the curvature arrow Φ_uu
Surface.normalAt s x y       -- n
Surface.metricI s x y        -- I  = Gram matrix of the tangent arrows
Surface.metricII s x y       -- II = normal dot products of curvature arrows
Surface.gaussK s x y         -- K = det II / det I
Surface.meanH s x y          -- H = ½ tr(I⁻¹ II)
Surface.principalK s x y     -- (κ_min, κ_max)
Surface.shapeOperator s x y w
```

### 1.8 Exercises

1. **Cone.** Take $\Phi(r,\theta) = (r\cos\theta, r\sin\theta, r)$.
   Compute the tangent arrows $\Phi_r$, $\Phi_\theta$ and $I = J^{\mathsf T}J$.
   Show that one bending direction is zero and hence $K = 0$.
2. **Determinant identity.** Prove $\det I = |\Phi_u\times\Phi_v|^2$ using
   only the cross-product identity
   $(a\times b)\cdot(c\times d) = (a\cdot c)(b\cdot d) - (a\cdot d)(b\cdot c)$.
   *(The §1.2 checked example verifies it numerically on the paraboloid.)*
3. **Graph curvature.** From the dot-product formulas in §1.2–§1.3 for
   $\Phi=(x,y,f)$, derive

   ```math
   K=\frac{f_{xx}f_{yy}-f_{xy}^2}{\bigl(1+f_x^2+f_y^2\bigr)^2}.
   ```

4. **Weingarten sign.** Differentiate $n\cdot\Phi_u = 0$ with respect to $v$
   and show $n\cdot\Phi_{uv} = -(\partial n/\partial v)\cdot\Phi_u$. What
   changes if you choose $-n$ instead of $n$?
5. **Cylinder vs sphere.** In §1.6, the cylinder's non-straight direction
   had $n\cdot\Phi_{\theta\theta}=-R$ and the straight one had $0$. Point
   out which matrix entry of $II$ encodes each, and why $K=0$ follows.
6. **Code experiment.** In `Main.lean`, move the paraboloid point from
   $(1,0)$ to $(2,0)$ and predict which of $I$, $II$, $K$, $\kappa_1$,
   $\kappa_2$ change and in which direction before running.
7. **Shape operator.** Verify numerically on the paraboloid that
   $II(d,d) = I(Sd, d)$ holds for several parameter directions $d$ — the
   defining property of $S = I^{-1}II$. *(Starter: `Surface.secondForm s 1 0 d`
   vs `Surface.firstForm s 1 0 (Surface.shapeOperator s 1 0 d)`.)*

---

## Part 2 — How tangent arrows turn: the tangential projection

> **Goal.** Express everything in **pure vector calculus**: the derivative of
> any surface-bound arrow is an ordinary 3-D arrow, and the *surface-visible*
> part is obtained by the tangential projection $T w = w - (w\cdot n)\,n$
> (defined in §1.5, computable in one `Vec3` line). Coefficients — the
> Christoffel symbols — appear only as the *coordinate bookkeeping* of that
> projection (§2.8), because that is what the Lean code stores.

### 2.1 The problem: basis arrows turn as you walk

In the $xy$-plane the basis arrows $e_1=(1,0)$, $e_2=(0,1)$ are the same
everywhere, so differentiating a vector field is component-wise. On a
parameterized surface, the tangent basis arrows are the *columns of the
Jacobian*, $\Phi_u$ and $\Phi_v$, and these genuinely change from point to
point. Differentiate them: each derivative $\Phi_{ui},\Phi_{vi},\dots$ (more
generally $\Phi_{ij}$) is again a 3-D arrow in $\mathbb{R}^3$.

Here is the central **vector** move. At any point, 3-space decomposes into the
tangent plane plus the normal arrow $n$, so every 3-D arrow splits uniquely:

```math
\boxed{\;
\Phi_{ij} \;=\; \underbrace{T\,\Phi_{ij}}_{\text{stays in the surface}}
\;+\; \underbrace{(n\cdot\Phi_{ij})\, n}_{\text{bending out of the surface}},
\qquad
T w \;=\; w - (w\cdot n)\,n
\;}
```

$T$ — *tangential projection* — is a linear map from 3-D arrows back onto the
tangent plane. Everything in this part is one of exactly two statements about
the arrows $\Phi_{ij}$:

- the **normal part** $(n\cdot\Phi_{ij})$ is the second fundamental form $II$
  (§1.3 — surface bending);
- the **tangent part** $T\Phi_{ij}$ is how much the surface-bound derivative
  of an arrow *unexpectedly* turns. That is the subject of this part.

No coefficients yet — these are honest arrows. Keep it that way as long as
possible. (§1.5's checked example already did one split:
$\Phi_{vv} = (0.264706, 0, 0.158824) + 0.514496\,n$.)

### 2.2 Why not Christoffels? (and why this track uses vectors)

You will meet $\Gamma^k_{\;ij}$ everywhere in textbooks, and it is fair to ask
why this track relegates them to a sub-role.

**Why the literature is symbol-heavy.** Classical surface theory and general
relativity grew up with *intrinsic* geometry: work with a quantity that lives
on a manifold with no ambient space at all. Inside a manifold with only
coordinates, there is nothing to project onto — no ambient derivative exists
— so you *must* track how the coordinate basis itself twists, and that means
carrying a 3-index array of correction numbers, $\Gamma^k_{\;ij}$, into every
derivative. They are the price of computing without an embedding. Index
gymnastics ($\nabla_i V^k = \partial_i V^k + \Gamma^k_{\;ij}V^j$ and the
higher contractions) is then a full telescope for that correction.

**Why vectors are preferable here.** This course always has the embedding
$\Phi : \mathbb{R}^2 \to \mathbb{R}^3$ in hand, and for embedded surfaces the
covariant derivative has a **one-line vector formula**:

```math
\boxed{\;\nabla_v \, w \;=\; D_v w - \bigl((D_v w)\cdot n\bigr)\,n \;=\; T(D_v w)\;}
```

i.e. differentiate the ambient arrows with the ordinary $\mathbb{R}^3$
derivative, then drop the normal piece. "Tangent + projection" replaces every
$\Gamma$ contraction. Concretely:

| | symbol route | vector route (this track) |
|---|---|---|
| what you differentiate | coordinate *coefficients* $V^k$ | the arrows $V$ themselves |
| correction mechanism | the $\Gamma$ terms (3-index bookkeeping) | subtract the normal part (`T(w) = w − (w·n)n`) |
| geometry read out | after summing repeated indices | directly, by looking at arrows |
| works without embedding | ✅ (its super-power) | needs $\Phi$ (fine — we have it) |
| code cost | build $\Gamma$ array, then contract | project + dot with basis arrows |

Where the symbols *are* the right tool, this track says so honestly: when you
only have the ruler $I$ (an abstract surface, a higher-dimensional manifold,
relativity), or when the textbook source you are following (Eberly's notes)
exploits the pure-$I$ formulas. And — since the Lean code does implement
`Surface.christoffel` — its job in the code is precisely to *store the
coordinates of $T\Phi_{ij}$* (§2.8, and see `Geodesia/Curvature.lean`). The
vector statement stays primary; the symbol array is a derived container.

### 2.3 The covariant derivative is the correction-free rate of turn

Take a surface-bound arrow field $V = V^u\,\Phi_u + V^v\,\Phi_v$. Its
ambient derivative in direction $i$ is ordinary calculus:

```math
\frac{\partial V}{\partial u_i}
= \frac{\partial V^u}{\partial u_i}\Phi_u
+ \frac{\partial V^v}{\partial u_i}\Phi_v
+ V^u\,\Phi_{ui} + V^v\,\Phi_{vi}.
```

The first two terms are expected changes of the arrow; the last two are the
surprise — *the basis arrows themselves turn* as the field's coefficients
move with the point. The surface-visible rate of turn is: project the whole
thing back.

```math
\boxed{\;\nabla_i V \;=\; T\!\left(\frac{\partial V}{\partial u_i}\right)\;}
```

Since $V^u\Phi_u+V^v\Phi_v$ is already tangent ($T$ passes through it
unchanged), the correction is exactly "differentiate the arrows, project".
In pencil-and-paper coordinates the same statement becomes the familiar
textbook pair

```math
\nabla_i V^k \;=\; \frac{\partial V^k}{\partial u_i}
+ \Gamma^k_{\;ij}\,V^j,
```

with the $\Gamma$-correction terms coming straight from
$T(\Phi_{ki}) = \Gamma^\cdot_{ki}\Phi_u + \cdots$ (i.e., from expanding
$T\Phi_{ij}$ in $\{\Phi_u,\Phi_v\}$). On a flat plane in standard coordinates
the basis arrows never turn: $\Gamma^k_{\;ij}=0$ and every $\nabla$ is an
ordinary partial. The same plane *in polar coordinates* does need
corrections; on a sphere every direction needs them (§2.5–2.6). This is the
whole narrative of the Christoffel story — but written in arrow-language it
is simply $T(\partial\text{-derivative})$, with the $\Gamma$ terms surfacing
only when a pencil demands them.

**Signature.** The library's consumer of the $\Gamma$ array is the geodesic
equation: a geodesic is a curve whose *tangent part* of $\gamma''$ vanishes,
$T\gamma'' = 0$ (Chapter 3). `Geodesia/Geodesic.lean` implements the
coordinate form of that single projection:

```lean
def Geodesic.accel (s : Surface) (st : GeoState) : Vec2
    -- (u″, v″) = −Γ^k_ij q′ⁱ q′ʲ: the coefficient form of T γ″ = 0
```

### 2.4 The Gauss formula, in arrow form

The one-sentence motto of the whole part:
*the second-derivative arrows split into tangent part $+$ normal part, and
every formula is one of those two projections.*

```mermaid
flowchart LR
    A["second-derivative arrows<br/>Φᵤᵤ, Φᵤᵥ, Φᵥᵥ (3-D)"] -->|"project out bending:<br/>(n·Φᵢⱼ) n"| B["II — second fundamental form"]
    A -->|"remove bending:<br/>T = w − (w·n)n"| C["tangent part: turning<br/>(covariant derivative)"]
    B -.->|"basis-coefficients<br/>(optional!)"| D["Γᵏᵢⱼ — the symbols"]
    C -.->|"basis-coefficients<br/>(optional!)"| D
```

Both routes — the "project-out" route and the "differentiate-the-metric"
route of §2.8 — measure the same arrows and agree numerically; they are the
same computation with and without coordinates. (§2.8's checked example makes
"agree numerically" exact: the gap is `0.000000`.)

### 2.5 Example: the flat plane in polar coordinates

Flat plane, wrapped geometry — the perfect test that tangent parts are about
the *basis arrows*, not about bending:

```math
\Phi(r,\theta) = (r\cos\theta,\ r\sin\theta,\ 0).
```

Basis arrows and Gram matrix:

```math
\Phi_r = (\cos\theta,\ \sin\theta,\ 0),\quad
\Phi_\theta = (-r\sin\theta,\ r\cos\theta,\ 0),\quad
I = \begin{bmatrix}1&0\\0&r^2\end{bmatrix}.
```

Differentiate the arrows (all ordinary $\mathbb{R}^3$ calculus):

```math
\frac{\partial\Phi_r}{\partial r}=0,\qquad
\frac{\partial\Phi_r}{\partial\theta}=\frac{1}{r}\,\Phi_\theta,\qquad
\frac{\partial\Phi_\theta}{\partial r}=\frac{1}{r}\,\Phi_\theta,\qquad
\frac{\partial\Phi_\theta}{\partial\theta}=-r\,\Phi_r.
```

The plane's normal is a constant vector $n=(0,0,1)$, so these derivatives are
*entirely tangent*: $T(\text{each row}) = $ each row. Two consequences worth
pausing over:

- the *basis arrows turn* — $\Phi_r$ turns into a $\Phi_\theta$ arrow when
  pushed around $\theta$ — despite the *plane being flat*;
- and the *arrow turning is exactly what cancels a curve's own turning*:
  for a trajectory that stays geometrically straight, the coordinate parts
  $\ddot r$, $\ddot\theta$ must push against these turning arrows so that
  the *surface arrow* $\gamma''$ ends up purely normal (Chapter 3 makes
  this precise).

For comparison, in coordinates: these arrow equations *are*
$\Gamma^\theta_{\;\theta r} = \Gamma^\theta_{\;r\theta} = 1/r$ and
$\Gamma^r_{\;\theta\theta}=-r$ — read off the same displayed rows, with
coefficients recorded in the basis. The vector statement comes first; the
coefficients are nothing more than these rows written as slips of paper.
The cross-check with the pure-$I$ Gram-derivative formula (§2.8) gives the
same numbers, as it must.

**Signature + Checked example.** The polar plane is three exact `Smooth2`
combinators — then the code's pure-$I$ route reproduces every coefficient of
the hand-derived rows above (at $r = 1.5$: $-r = -1.5$, $1/r = 0.666667$),
while $II = 0$ confirms the plane is flat:

```lean
def polarPlane : Surface where
  sx := Smooth2.mul Smooth2.coordX (Smooth2.cos Smooth2.coordY)
  sy := Smooth2.mul Smooth2.coordX (Smooth2.sin Smooth2.coordY)
  sz := Smooth2.const 0

#eval show IO Unit from do
  let I := Surface.metricI polarPlane 1.5 0.7
  IO.println s!"I (r = 1.5) = [[{I.a}, {I.b}], [{I.b}, {I.c}]]"
  let II := Surface.metricII polarPlane 1.5 0.7
  IO.println s!"II          = [[{II.a}, {II.b}], [{II.b}, {II.c}]]"
  let Γ := Surface.christoffel polarPlane 1.5 0.7
  IO.println s!"Γ^r_θθ      = {Γ 0 1 1}"
  IO.println s!"Γ^θ_rθ, Γ^θ_θr = {Γ 1 0 1}, {Γ 1 1 0}"
  IO.println s!"Γ^r_rr, Γ^r_rθ, Γ^θ_θθ = {Γ 0 0 0}, {Γ 0 0 1}, {Γ 1 1 1}"
```

**Output:**
```text
I (r = 1.5) = [[1.000000, 0.000000], [0.000000, 2.250000]]
II          = [[0.000000, 0.000000], [0.000000, 0.000000]]
Γ^r_θθ      = -1.500000
Γ^θ_rθ, Γ^θ_θr = 0.666667, 0.666667
Γ^r_rr, Γ^r_rθ, Γ^θ_θθ = 0.000000, 0.000000, 0.000000
```

Flat plane, nonzero $\Gamma$: the *arrows* turn, the *plane* does not bend.

### 2.6 Example: the unit sphere

Sphere $\Phi(\theta,\varphi)=
(\sin\theta\cos\varphi, \sin\theta\sin\varphi,\cos\theta)$,
normal $n = \Phi$ (the radial arrow), and Gram matrix
$I = \operatorname{diag}(1, \sin^2\theta)$.

Now differentiate the basis arrows — keeping the split intact:

```math
\frac{\partial\Phi_\theta}{\partial\theta}
= -\Phi \;=\; -\,n\quad (\text{entirely normal — } T\text{-part is 0}),
```

```math
\frac{\partial\Phi_\varphi}{\partial\theta}
= \cot\theta\,\Phi_\varphi\quad (\text{purely tangent}),
\qquad
\frac{\partial\Phi_\varphi}{\partial\varphi}
= -\sin\theta\cos\theta\,\Phi_\theta - \sin^2\theta\,n.
```

(Run the dot-product check yourself: $n\cdot\Phi_\theta = 0$ shows the first
split is clean, and $n\cdot \partial_\varphi\Phi_\varphi = -\sin^2\theta$
agrees with $-\sin^2\theta\,n$.)

This single row-set explains all the geometry the $\Gamma$ table used to
mumble about:

- a **meridian** ($\varphi$ fixed): its second-derivative arrow is
  $\Phi_{\theta\theta} = -n$ — normal only, tangent part zero → geodesic.
- around the **equator** ($\theta = \pi/2$):
  $\partial_\varphi\Phi_\varphi = -n$ too, so latitude circles at $\pi/2$
  are geodesics but no other latitudes are: at $\theta\ne\pi/2$ the tangent
  part $-\sin\theta\cos\theta\,\Phi_\theta$ is nonzero (the circle must
  accelerate inward, and the surface can only supply the normal fraction).
- general great circles are combinations, in the ordinary sense.

**Checked example.** The sphere's $I$ and its two nontrivial
$\Gamma$-coefficients at $(\theta,\varphi) = (1.0, 0.2)$ — exactly the two
arrow rows above ($\sin^2 1.0 = 0.708073$,
$-\sin\theta\cos\theta = -0.454649$, $\cot\theta = 0.642093$):

```lean
#eval show IO Unit from do
  let I := Surface.metricI sphere 1.0 0.2
  IO.println s!"I_sphere = [[{I.a}, {I.b}], [{I.b}, {I.c}]]"
  let Γ := Surface.christoffel sphere 1.0 0.2
  IO.println s!"Γ^θ_φφ   = {Γ 0 1 1}"  -- should be −sin θ cos θ
  IO.println s!"Γ^φ_θφ   = {Γ 1 0 1}"  -- should be cot θ
```

**Output:**
```text
I_sphere = [[1.000000, 0.000000], [0.000000, 0.708073]]
Γ^θ_φφ   = -0.454649
Γ^φ_θφ   = 0.642093
```

(`Γ^θ_θθ`, `Γ^θ_θφ`, `Γ^φ_θθ`, `Γ^φ_φφ` all come out $0$ — the $\Gamma$
table here is as sparse as the arrow rows.)

### 2.7 Example: the paraboloid

For $z = 0.3(x^2+y^2)$ at $(x,y)=(1,0)$ the code prints

```math
I=\begin{bmatrix}1.36&0\\0&1\end{bmatrix},
\qquad
\Gamma^u_{\;uu}=0.2647,\quad\Gamma^u_{\;uv}=0,\quad\Gamma^u_{\;vv}=0.2647.
```

Vector reading: the $v$-coordinate curve $x=\text{const}$ has second-derivative
arrow $(0,0,0.6)$, while the normal there is
$n=(-0.6,0,1)/\sqrt{1.36}$. Since $(0,0,0.6)$ is **not** parallel to $n$,
exactly the leftover — the tangent part $T(0,0,0.6)$ — has a nonzero
$\Phi_u$ component, and $0.2647$ is that arrow's coefficient in the user
direction (the coefficient of $\Phi_u$ inside $T\Phi_{vv}$, which in
coordinates is labelled $\Gamma^u_{\;vv}$). The same arrow's *normal*
coefficient is the $II_{22}$ entry from §1.3. One arrow, two
projections — that split pattern is the *entire subject*.

**Checked example.** Both halves of that story, printed by the code:
the $\Gamma^u$ row (§1.6 demo) and the explicit split of the same arrow
(§1.5):

```text
Γ¹₁₁, Γ¹₁₂, Γ¹₂₂ = 0.264706, 0.000000, 0.264706     -- coordinates of T Φ_vv (and T Φ_uu)
T(w) = (0.264706, 0.000000, 0.158824)               -- the arrow T Φ_vv itself, w = (0, 0, 0.6)
w·n  = 0.514496                                     -- the normal coefficient = II₂₂
```

### 2.8 The coordinate bookkeeping (what the code stores)

Everything above was arrow language. When you (or the Lean code) want actual
**numbers to feed an ODE**, the tangent parts have to be expressed in the
basis $\{\Phi_u,\Phi_v\}$ and the *coefficients* recorded. That's what the
following code computes — the coefficients are what tradition calls the
Christoffel symbols:

```lean
Surface.dPhi  s i x y   -- the basis arrows Φ_i
Surface.d2Phi s i j x y -- the second-derivative arrows Φ_ij
Surface.dMetricEntry s u i j x y
                        -- exact Gram-matrix derivative ∂_u (Φ_i·Φ_j)
Surface.metricI s x y   -- the Gram matrix I
Surface.christoffel s x y
                        -- coordinates of T(Φ_ij) in {Φ_u, Φ_v}:
                        -- Fin 2 → Fin 2 → Fin 2 → Float
```

Two equivalent coordinate formulas, both implemented exactly:

**(a) Direct projection.** For an unknown expansion
$T\Phi_{ij} = a\,\Phi_u + b\,\Phi_v$, dot each row with the basis arrows;
the normal-part drops, and the Gram matrix $I=J^{\mathsf T}J$ remains:

```math
I\begin{bmatrix}a\\b\end{bmatrix} =
\begin{bmatrix}\Phi_u\cdot\Phi_{ij}\\ \Phi_v\cdot\Phi_{ij}\end{bmatrix}
\;\;\Longrightarrow\;\;
\Gamma^k_{\;ij}=\bigl(I^{-1}\bigr)^{k\ell}
\bigl(\Phi_\ell\cdot\Phi_{ij}\bigr).
```

**(b) From the ruler alone.** Differentiate $I_{ij} = \Phi_i\cdot\Phi_j$,
cancel shared terms cyclically, and get the pure-$I$ formula

```math
\Gamma^k_{\;ij}
= \tfrac12\bigl(I^{-1}\bigr)^{k\ell}
\bigl(\partial_i I_{j\ell}+\partial_j I_{i\ell}-\partial_\ell I_{ij}\bigr),
```

which `Curvature.lean` computes with **exact** second partials
(`Surface.dMetricEntry`):

```math
\partial_u\,(\Phi_j\cdot\Phi_k)
= \Phi_{uj}\cdot\Phi_k + \Phi_j\cdot\Phi_{uk}.
```

This bare-$I$ route is the one the historical notes favor, since it works for
absolutely-intrinsic data (a metric alone, no embedding). The code implements
(b); route (a) is three dot products and a $2\times 2$ solve, so it makes a
portable independent check:

```lean
/-- Route (a): Γ^k_ij = (I⁻¹)^kℓ (Φ_ℓ · Φ_ij), straight from the projection. -/
def gammaDirect (s : Surface) (x y : Float) : Fin 2 → Fin 2 → Fin 2 → Float :=
  fun k i j =>
    let Iinv := (Surface.metricI s x y).inv
    let b0 := (Surface.dPhi s 0 x y).dot (Surface.d2Phi s i j x y)
    let b1 := (Surface.dPhi s 1 x y).dot (Surface.d2Phi s i j x y)
    Iinv.entry k 0 * b0 + Iinv.entry k 1 * b1

/-- Squared gap between the two routes, over all eight (k, i, j). -/
def gammaGap (s : Surface) (x y : Float) : Float :=
  let Γ := Surface.christoffel s x y
  let d := gammaDirect s x y
  (Γ 0 0 0 - d 0 0 0) * (Γ 0 0 0 - d 0 0 0)
    + (Γ 0 0 1 - d 0 0 1) * (Γ 0 0 1 - d 0 0 1)
    + (Γ 0 1 0 - d 0 1 0) * (Γ 0 1 0 - d 0 1 0)
    + (Γ 0 1 1 - d 0 1 1) * (Γ 0 1 1 - d 0 1 1)
    + (Γ 1 0 0 - d 1 0 0) * (Γ 1 0 0 - d 1 0 0)
    + (Γ 1 0 1 - d 1 0 1) * (Γ 1 0 1 - d 1 0 1)
    + (Γ 1 1 0 - d 1 1 0) * (Γ 1 1 0 - d 1 1 0)
    + (Γ 1 1 1 - d 1 1 1) * (Γ 1 1 1 - d 1 1 1)
```

**Checked example.** The two routes agree exactly on every surface of this
chapter — the arrow statement (a) and the ruler statement (b) are the same
computation:

```lean
#eval show IO Unit from do
  IO.println s!"paraboloid:  gap = {gammaGap surf 1 0}"
  IO.println s!"polar plane: gap = {gammaGap polarPlane 1.5 0.7}"
  IO.println s!"sphere:      gap = {gammaGap sphere 1.0 0.2}"
```

**Output:**
```text
paraboloid:  gap = 0.000000
polar plane: gap = 0.000000
sphere:      gap = 0.000000
```

> One sentence to keep both views in place: *the Christoffel symbols are the
> coordinates of the tangent part of the second-derivative arrows: numbers
> read off by whichever formula is convenient, never required by the vector
> story itself.*

### 2.9 Exercises

1. **Polar, both ways.** From §2.5's arrow rows read
   $T(\partial_\theta\Phi_r) = (1/r)\Phi_\theta$ and
   $T(\partial_\theta\Phi_\theta) = -r\,\Phi_r$; then check the same numbers
   via the pure-$I$ formula of §2.8 with $I=\mathrm{diag}(1,r^2)$.
   *(The §2.5 checked example already ran the second half for you.)*
2. **Standard coordinates.** Show that every tangent-part coefficient
   vanishes for the plane in standard $(x,y)$ coordinates — and explain why
   that matches §2.5's nonzero polar numbers: the *plane* is flat in both;
   the *arrows* turn only after polar wrapping.
3. **Symmetry.** Use $\Phi_{uv}=\Phi_{vu}$ to show the tangent-part
   coordinates are symmetric in $i,j$: $\Gamma^k_{\;ij}=\Gamma^k_{\;ji}$.
   *(Check: `Surface.christoffel s x y k 0 1 == Surface.christoffel s x y k 1 0`.)*
4. **First-kind identity.** Prove
   $\partial_k\,(\Phi_i\cdot\Phi_j)
    = \Phi_{ik}\cdot\Phi_j + \Phi_i\cdot\Phi_{jk}$
   and use it to derive the pure-$I$ formula of §2.8(b) in full.
5. **Cylinder.** For $\Phi(\theta,z)=(R\cos\theta,R\sin\theta,z)$, compute
   each of $\partial_\theta\Phi_\theta=(-R\,n)$,
   $\partial_\theta\Phi_z = 0$-type rows and deduce that *every* tangent
   part vanishes. Is that consistent with the cylinder being flat?
6. **Sphere check.** Verify the three arrow splits of §2.6 from
   $I=\mathrm{diag}(1,\sin^2\theta)$ and the pure-$I$ route.
   *(The §2.6 checked example did the pure-$I$ half.)*
7. **Conserved pairing.** Show that on the sphere
   $\sin^2\theta\,\dot\varphi$ stays constant along any curve whose tangent
   part of $\gamma''$ vanishes (Chapter 3). Hint: differentiate
   $\sin^2\theta\,\dot\varphi$ in $s$ and substitute the two component
   equations.
8. **Code experiment.** Evaluate `Surface.christoffel` on the paraboloid at
   $(0,0)$ and $(2,0)$; translate the printed numbers back into
   *which arrow turned how much*, using §2.7's vector reading.

---

> **Where the outputs come from.** Every `Output:` block in this chapter was
> produced by running the adjacent snippet with `lake env lean` against the
> library, in the spirit of the Lean documentation (Verso): the examples are
> part of the build, not decoration. If the project ever adopts Verso, these
> ```lean blocks are exactly the pieces that would become checked examples.
