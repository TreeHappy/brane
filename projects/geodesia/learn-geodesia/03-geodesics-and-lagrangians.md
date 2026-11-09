# 3️⃣ Geodesics — curves whose bending is purely normal

> **Goal.** Define geodesics in pure vector language: the *tangential
> projection* $T = w - (w\cdot n)n$ of the curve's second derivative must
> vanish. Derive the geodesic equation twice — once from that vector
> condition, once from a length functional — in a form where **no
> Christoffel symbols appear**; then state the conserved quadratic form
> that makes the numerical integration stable.

---

## 3.1 Three equivalent definitions

A geodesic is characterised three ways, all the same equation:

1. **Geometric (vector).** $\;T\gamma'' = 0$: the second-derivative arrow of
   the surface curve points *only* along the normal — zero component inside
   the tangent plane.
2. **Variational.** A critical point of the length functional
   $\ell[\gamma] = \int |\gamma'|\, ds$ (or of the quadratic functional
   $\int \tfrac12\,|J\,\dot q|^2 ds$).
3. **Intrinsic straightness.** Cut the surface locally and lay it flat: the
   geodesic becomes a straight line — the straightest possible curve on the
   surface.

> **Caveat.** A geodesic is *locally* length-minimising, not globally
> shortest: on a sphere the long arc of a great circle is a geodesic too,
> and between antipodal points every great circle ties for the same length.

---

## 3.2 The vector condition, and the equation it produces

Let $\gamma(s) = \Phi(q(s))$ with direction numbers $\dot q = (\dot u,
\dot v)$. Differentiating twice with the chain rule:

```math
\gamma' = \dot u\,\Phi_u + \dot v\,\Phi_v,
```

```math
\gamma'' = \ddot u\,\Phi_u + \ddot v\,\Phi_v
+ \Phi_{uu}\,\dot u^2 + 2\,\Phi_{uv}\,\dot u\dot v + \Phi_{vv}\,\dot v^2.
```

Everything on the right is an ordinary 3-D arrow combination — no coordinate
machinery yet. The geometric definition says its **tangential projection
must vanish**, $T\gamma'' = 0$, i.e. $\gamma'' = \bigl(n\cdot\gamma''\bigr)n$.
To get a *differential equation for the two scalar functions* $(u,v)$, take
it the Chapter-2 way: dot each tangent basis arrow into $\gamma''$
(the normal part drops out of these dot products):

```math
\Phi_k \cdot \gamma'' \;=\; 0, \qquad k \in \{u, v\}.
```

Expand (orthogonality of $\Phi_u,\Phi_v$ to $n$ kills the second-half terms):

```math
\boxed{\;
\underbrace{I\,\ddot q \;+\; g\bigl(q, \dot q\bigr)}_{\text{a } 2\text{-vector}}
\;=\; 0
\qquad\text{with}\qquad
g_k \;=\;
\Phi_k\cdot\Bigl(
\Phi_{uu}\,\dot u^2 + 2\,\Phi_{uv}\,\dot u \dot v + \Phi_{vv}\,\dot v^2
\Bigr),\;}
```

where $I = J^{\mathsf T}J$ is the first fundamental form and $g$ is a
2-component vector assembled purely by *dotting basis arrows with
second-derivative arrows* — vectors, twice, no symbols. Solved for the
acceleration:

```math
\ddot q = -\,I^{-1}\,g.
```

That is the **geodesic equation** in this course's native currency: compare
by eye against the boxed form the $\Gamma$'s usually occupy,

```math
\underbrace{\Bigl(\ddot q^k + \Gamma^k_{\;ij}\dot q^i\dot q^j\Bigr)}_{tangent-plane part} = 0
\quad\Longleftrightarrow\quad
I\,\ddot q + g = 0,
```

— the two right-hand sides are identical because $g_k = \Gamma_{kij}\,d^id^j$,
where $\Gamma_{kij} = \Phi_k\cdot\Phi_{ij}$ are the "first-kind" coordinates of
the same tangent parts (Chapter 2 §2.8). The vector form is nicer for two
reasons: it does not need the $\Gamma$ array at all — just arrows and dots —
and it *exposes* what the equation is: *make the parameter acceleration
exactly cancel the turning of the basis arrows, so the surface arrow
$\gamma''$ is left purely normal.*

Note what does *not* appear: the normal part of $\gamma''$ never vanishes on
a curved surface. A geodesic is the curve that lets the surface's own
bending supply *all* of its acceleration.

---

## 3.3 The variational route, in vector terms

The same equation falls out of minimising length, again without symbols.
Use the quadratic integrand to avoid the square root:

```math
L\bigl(q,\dot q\bigr) = \tfrac12\,\dot q^{\mathsf T} I\bigl(q\bigr)\,\dot q
= \tfrac12\,\bigl|J(q)\,\dot q\bigr|^2.
```

Euler–Lagrange for $\int L\,ds$:

```math
\frac{d}{ds}\frac{\partial L}{\partial \dot q^k}
- \frac{\partial L}{\partial q^k} = 0.
```

Compute by pure calculus:

```math
\frac{\partial L}{\partial \dot q^k}
= \bigl(I\,\dot q\bigr)_k
= \Phi_k\cdot\gamma',
\qquad
\frac{\partial L}{\partial q^k}
= \dot q^i\dot q^j\, \Phi_{ik}\cdot\Phi_j.
```

Expanding the $s$-derivative of $\Phi_k \cdot \gamma'$ with the product rule
and collecting the $\ddot q$ terms produces exactly the same boxed equation
of §3.2 — the $\Gamma$ symbols that this derivation is traditionally dressed
in are, once more, just the coordinates of the tangent parts that the vector
arithmetic already handles directly. No hidden assumptions involved: it is
the standard fact that *critical points of length* satisfy the same
equation as *vanishing tangent projection* — and the split
"a tangent part and a normal part" of Chapter 2 is all it takes to go
between them.

Like any Lagrangian system, this one automatically transports its own
"energy":

```math
\frac{d}{ds}\Bigl(\tfrac12\,\dot q^{\mathsf T} I\, \dot q\Bigr) = 0,
```

i.e. the surface arrow $\gamma'$ keeps constant length along the curve — the
conserved quadratic form used as a numerical health check in the code.

---

## 3.4 The conserved quadratic form and reparameterization

Direct check straight from the boxed equation — notice it needs *no*
$\Gamma$'s, only $I\ddot q + g = 0$ and the dot-product expansion of $g$:

```math
\frac{d}{ds}\Bigl(\tfrac12\,\dot q^{\mathsf T} I\,\dot q\Bigr)
= \dot q^{\mathsf T}\bigl(I\,\ddot q + g\bigr) = 0.
```

Consequences:

- $|J\dot q| = \sqrt{2 \times (\text{that value})}$ is **constant**: a
  geodesic parameterized this way is automatically constant-speed.
- If $s \mapsto as+b$ is an affine reparameterization, the equation is
  unchanged; arc length is a conventional choice.
- For the raw length functional
  $\ell = \int \sqrt{\dot q^{\mathsf T} I\,\dot q}\, ds$, the equations only
  determine the *image* of the curve; picking constant speed reduces them
  to the same geodesic equation. The quadratic functional is the
  computationally convenient version of the same geometry.

**Conserved quantities from symmetry.** Whenever the Gram matrix $I$ does
not depend on a coordinate $q^\ell$, the pairing

```math
p_\ell = \bigl(I\,\dot q\bigr)_\ell = \Phi_\ell\cdot\gamma'
```

is conserved — i.e. the projection of the surface velocity arrow along a
symmetry direction of the surface cannot change:

- **Polar plane:** $\Phi_\theta\cdot\gamma' = r^2\dot\theta$ is constant —
  the straightness of the geodesic, seen in rotating coordinates.
- **Sphere:** $\sin^2\theta\,\dot\varphi$ is constant.
- **Surfaces of revolution** $\Phi(r,\theta) = (r\cos\theta, r\sin\theta, f(r))$:
  the same $\Phi_\theta$-pairing is constant (historically "Clairaut's
  relation").

This is a practical numerical check: a correct integrator must conserve
these pairings, and any measurable drift is a red flag.

Vector framing once more: for the coordinate function $q^\ell$, the pairing
$\Phi_\ell\cdot\gamma'$ is "the surface arrow's $\ell$-component against the
symmetric-direction basis arrow" — a coordinate-free statement of *the
projection of the velocity along a symmetry cannot change*.

---

## 3.5 Worked examples

### Plane in Cartesian coordinates

$I$ is the identity; $g = 0$ because every second-derivative arrow
$\Phi_{ij}$ vanishes; so $\ddot q = 0$: geodesics are straight lines at
constant parameter speed.

### Plane in polar coordinates

With the polar basis arrows of Chapter 2 (§2.5), the two boxed component
equations $I\ddot q = -g$ become

```math
\ddot r - r\,\dot\theta^2 = 0,
\qquad
\ddot\theta + \frac{2}{r}\,\dot r\,\dot\theta = 0.
```

These are straight lines written in rotating coordinates. The conserved
pairing $r^2\dot\theta$ is the projection of the (constant) direction arrow
of the straight line onto the rotating $\Phi_\theta$ arrow — the $\sim 1/r$
in the second equation comes from that arrow's changing length $r$.

### Sphere

With the sphere arrows of Chapter 2 (§2.6):

```math
\ddot\theta - \sin\theta\cos\theta\,\dot\varphi^2 = 0,
\qquad
\ddot\varphi + 2\cot\theta\,\dot\theta\dot\varphi = 0.
```

Solutions are great circles. Meridians ($\varphi = \text{const}$) give
$\ddot\theta = 0$; the equator ($\theta = \pi/2$) gives $\ddot\varphi = 0$.
The conserved pairing is $\sin^2\theta\,\dot\varphi$ — the projection of the
velocity against the rotation-symmetric arrow $\Phi_\varphi$.

### Paraboloid

For $z = f(x,y)$ the Gram matrix is
$I = \begin{bmatrix}1+f_x^2 & f_xf_y\\ f_xf_y & 1+f_y^2\end{bmatrix}$; the
$g$-vector is assembled from the second-derivative arrows
$\Phi_{xx}, \Phi_{xy}, \Phi_{yy}$ exactly as in §3.2's box; and no
closed-form solution exists in general — which is precisely why the
numerical methods of Chapters 4–6 exist. The code prints a conserved
quadratic form $\tfrac12\dot q^{\mathsf T}I\,\dot q$ that stays constant
along the integrated geodesic — a first check that the ODE and the
integrator agree.

### Cone

A cone is developable: cut it along a generator and unroll it to a flat
sector, where geodesics become straight lines. Beautiful check of the
intrinsic nature of geodesics: $K = 0$, and the "shortest path" is a
straight line after unrolling.

---

## 3.6 How the code implements it

```lean
Geodesic.accel s st      -- (u'', v'') from the tangent-part coordinates
Geodesic.rhs s st        -- first-order system (u, v, u', v')
Geodesic.rk4Step s h st  -- one RK4 step
Geodesic.integrate s h n st
Geodesic.energy s st     -- ½ · dot qᵀ I dot q  (the conserved quadratic form)
Geodesic.speed s st      -- |J·dotq|  =  √(2 · that value)
```

`accel` computes exactly the boxed equation of §3.2, up to coordinates. In
the code the tangent-part coefficients come pre-computed as the
`Surface.christoffel` array (those are, again, purely the coordinates of
$T\Phi_{ij}$ in the basis — the routine *is* the vector formula of Chapter
2 §2.8, in coefficient clothing), and `accel` contracts them with the
direction numbers. In vector form:

```math
\begin{bmatrix}\ddot u\\ \ddot v\end{bmatrix} = -\,I^{-1}\,g,
\qquad
g_k \;=\; \bigl(\text{the }k\text{-th tangent part, dotted with basis
arrows}\bigr)\cdot \dot q\otimes\dot q,
```

which is the line-for-line Lean recipe for "form the arrow sum
$\dot u^2\Phi_{uu} + 2\dot u\dot v\Phi_{uv} + \dot v^2\Phi_{vv}$, project
tangentially, read off coordinates" — the coordinate-rotations in §3.2's
first boxed equation and $\Gamma$-contractions tie up to the same result.

The demo integrates a geodesic on the paraboloid and prints the conserved
quadratic form holding to the printed precision — the discrete analogue of
§3.4's identity.

---

## 3.7 Exercises

1. **Normal acceleration.** Show that a geodesic satisfies
   $n\cdot\gamma'' = II(\dot q,\dot q)$ and that its *tangential projection*
   vanishes — then verify both numerically on the sphere.
2. **Variational algebra.** Fill in the symmetrization step in §3.3 and
   explain why contracting with $\dot q^i \dot q^j$ lets you replace
   $\partial_i g_{kj}$ by
   $\tfrac12(\partial_i g_{kj} + \partial_j g_{ki})$ inside $I$-derivative
   terms.
3. **Length vs quadratic functional.** Derive Euler–Lagrange for
   $\ell = \int\sqrt{\dot q^{\mathsf T}I\,\dot q}\,ds$ and show it reduces
   to the geodesic equation when $|\gamma'|$ is constant.
4. **Polar straight lines.** Show that
   $r\cos(\theta - \theta_0) = d$ satisfies both polar component equations.
5. **Sphere great circles.** Verify that
   $\gamma(s) = \cos s\,U + \sin s\,V$ with two orthonormal arrows $U,V$
   satisfies the sphere component equations — *first* by checking that its
   surface second derivative points purely radially (the vector test!),
   *then* from the components.
6. **Symmetry pairing.** Prove that $\partial_\ell\,I = 0$ implies
   $p_\ell = \Phi_\ell\cdot\gamma'$ is constant along a geodesic.
7. **Surface of revolution.** Derive the Gram matrix for
   $\Phi(r,\theta) = (r\cos\theta, r\sin\theta, f(r))$ and the corresponding
   conserved pairing.
8. **Affine parameter.** Show that if $q(s)$ solves the geodesic equation,
   then $q(as+b)$ also solves it. What happens to the conserved quadratic
   form?
9. **Code check.** Run `lake exe geodesia` and verify that the conserved
   quadratic form is unchanged after 100 RK4 steps. Then change the step
   size and see how the error behaves.
