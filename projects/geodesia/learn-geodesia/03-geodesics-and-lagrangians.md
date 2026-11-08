# 3️⃣ Geodesics — curves whose bending is purely normal

> **Goal.** Define geodesics purely with vectors, derive the geodesic
> equation twice (once by "the tangent part of $\gamma''$ must vanish",
> once from a length functional), and state the conserved quadratic form
> that makes the numerical integration stable.

---

## 3.1 Three equivalent definitions

A geodesic on a surface can be characterised in three ways, all leading to
the same equation:

1. **Geometric.** A curve whose second-derivative arrow $d^2\gamma/ds^2$
   has *no component inside the tangent plane* — it is entirely along the
   normal $n$.
2. **Variational.** A critical point of the length functional
   $\ell[\gamma] = \int |\gamma'|\, ds$ (or of the quadratic functional
   $\int \tfrac12\,|J\,\dot q|^2 ds$).
3. **Intrinsic straightness.** If you cut the surface locally and lay it
   flat, the geodesic becomes a straight line: its image is the straightest
   possible curve on the surface.

> **Caveat.** A geodesic is *locally* length-minimising. It need not be the
> global shortest path: on a sphere, the long arc of a great circle is also
> a geodesic, and between antipodal points every great circle is a geodesic
> of the same length.

---

## 3.2 The tangent part of the second derivative must vanish

Let $\gamma(s) = \Phi\big(u(s), v(s)\big)$ with direction numbers
$\dot q = (\dot u, \dot v)$. Differentiate twice (the chain rule again):

```math
\gamma' = \dot u\,\Phi_u + \dot v\,\Phi_v,
```

```math
\gamma'' = \ddot u\,\Phi_u + \ddot v\,\Phi_v
+ \Phi_{uu}\,\dot u^2 + 2\,\Phi_{uv}\,\dot u\dot v + \Phi_{vv}\,\dot v^2.
```

Now insert the Chapter 1/2 decomposition of each curvature arrow into its
tangent part and normal part. Collect terms:

```math
\gamma'' =
\underbrace{\Bigl(\ddot q^k + \Gamma^k_{\;ij}\dot q^i\dot q^j\Bigr)\Phi_k
}_{\text{tangent-plane part}}
+
\underbrace{II_{ij}\,\dot q^i\dot q^j\; n}_{\text{normal part}}.
```

Definition 1 says: the tangent-plane part is zero. Therefore

```math
\boxed{\;\ddot u
+ \Gamma^u_{\;uu}\dot u^2 + 2\Gamma^u_{\;uv}\dot u\dot v
+ \Gamma^u_{\;vv}\dot v^2 = 0,\;}
```

```math
\boxed{\;\ddot v
+ \Gamma^v_{\;uu}\dot u^2 + 2\Gamma^v_{\;uv}\dot u\dot v
+ \Gamma^v_{\;vv}\dot v^2 = 0.\;}
```

This is the **geodesic equation** — two coupled second-order equations in
the two parameter functions. Everything on the left is assembled from
vectors you already know how to build: the second derivative arrows, their
projections, and the rates of change of the coordinates.

Note what does *not* appear: the normal part never vanishes (on a curved
surface, $\gamma''$ must point along $n$). A geodesic is the curve that
lets the surface's own bending supply all of its acceleration.

---

## 3.3 The variational route, in vector terms

The same equation falls out of minimizing length. Work with the square of
the length integrand to avoid the square root:

```math
L\bigl(q,\dot q\bigr) = \frac12\,\dot q^{\mathsf T} I\!\bigl(q\bigr)\,\dot q
= \frac12\,\bigl|J(q)\,\dot q\bigr|^2.
```

The **Euler–Lagrange equations** for the functional $\int L\,ds$ are

```math
\frac{d}{ds}\frac{\partial L}{\partial \dot q^k}
- \frac{\partial L}{\partial q^k} = 0.
```

Compute the two pieces by pure calculus:

```math
\frac{\partial L}{\partial \dot q^k} = \bigl(I\,\dot q\bigr)_k
= \bigl(\Phi_k\cdot\gamma'\bigr),
```

```math
\frac{\partial L}{\partial q^k}
= \tfrac12\, \dot q^{i}\dot q^{j}\,
\frac{\partial I_{ij}}{\partial q^k}
= \dot q^{i}\dot q^{j}\,\Phi_{ik}\cdot\Phi_j.
```

Substituting, expanding the $s$-derivative with the product rule and
multiplying through by $I^{-1}$ yields precisely the geodesic equation
above. No hidden assumptions are involved: it is the standard fact that
"critical points of length" satisfy the same equation as "vanishing tangent
component of $\gamma''$", and the $\Gamma$ coefficients are exactly the
number-crunching device that expresses both.

Like any Lagrangian system, this one automatically transports the length
functional: solving it also shows

```math
\frac{d}{ds}\Bigl(\tfrac12\,\dot q^{\mathsf T} I\, \dot q\Bigr) = 0,
```

i.e. the arrow $\gamma'$ keeps constant length along the curve — the
conserved quadratic form used as a numerical health check in the code.

---

## 3.4 The conserved quadratic form and reparameterization

Direct check straight from the geodesic equation:

```math
\frac{d}{ds}\Bigl(\tfrac12\,\dot q^{\mathsf T} I\,\dot q\Bigr)
= \dot q^{\mathsf T} I\Bigl(\ddot q + \Gamma[\dot q,\dot q]\Bigr) = 0.
```

Consequences:

- $|J\dot q| = \sqrt{2 \times (\text{that value})}$ is **constant**. A
  geodesic parameterized this way is automatically constant-speed.
- If $s \mapsto as+b$ is an affine reparameterization, the equation is
  unchanged; arc length is simply the conventional choice.
- For the raw length functional
  $\ell = \int \sqrt{\dot q^{\mathsf T} I\,\dot q}\, ds$, the equations
  only determine the *image* of the curve, and picking constant speed
  reduces them to the same geodesic equation. The quadratic functional is
  the computationally convenient version of the same geometry.

**Conserved quantities from symmetry.** Whenever the Gram matrix $I$ does
not depend on a coordinate $q^\ell$
($\partial_\ell\, (\Phi_i\cdot\Phi_j) = 0$), the momentum-like pair

```math
p_\ell = \bigl(I\,\dot q\bigr)_\ell = \Phi_\ell\cdot\gamma'
```

is conserved:

- **Polar plane:** $\Phi_\theta\cdot\gamma' = r^2\dot\theta$ is constant —
  this is exactly the straightness of the geodesic written in rotating
  coordinates.
- **Sphere:** $\sin^2\theta\,\dot\varphi$ is constant.
- **Surfaces of revolution** $\Phi(r,\theta) = (r\cos\theta, r\sin\theta, f(r))$:
  the same $\Phi_\theta$-pairing is constant (historically
  "Clairaut's relation").

This is a practical numerical check: a correct integrator must conserve
these pairings, and any measurable drift is a red flag.

Note that for a free (uncoupled) parameter function, the pairing
$\Phi_\ell\cdot\gamma'$ is just "the $\ell$-component of the surface
arrow in the $\Phi_\ell$ basis" — a coordinate-free way of saying the
projection of the arrow along a symmetry direction cannot change.

---

## 3.5 Worked examples

### Plane in Cartesian coordinates

$I$ is the identity, all $\Gamma = 0$, so $\ddot q = 0$: geodesics are
straight lines at constant parameter speed.

### Plane in polar coordinates

```math
\ddot r - r\,\dot\theta^2 = 0,
\qquad
\ddot\theta + \frac{2}{r}\,\dot r\,\dot\theta = 0.
```

These are straight lines written in rotating coordinates. The conserved
pairing $r^2\dot\theta$ is the projection of the (constant) direction arrow
of the straight line onto the rotating $\Phi_\theta$ arrow — the ~$1/r$
dependence comes from that arrow's changing length $r$.

### Sphere

```math
\ddot\theta - \sin\theta\cos\theta\,\dot\varphi^2 = 0,
\qquad
\ddot\varphi + 2\cot\theta\,\dot\theta\dot\varphi = 0.
```

Solutions are great circles. Meridians ($\varphi = \text{const}$) give
$\ddot\theta = 0$; the equator ($\theta = \pi/2$) gives $\ddot\varphi = 0$.
The conserved pairing is $\sin^2\theta\,\dot\varphi$.

### Paraboloid

For $z = f(x,y)$ the Gram matrix is
$I = \begin{bmatrix}1+f_x^2 & f_xf_y\\ f_xf_y & 1+f_y^2\end{bmatrix}$, the
$\Gamma$'s follow from it, and no closed-form solution exists in general —
which is precisely why the numerical methods of Chapters 4–6 exist. The
code prints a conserved quadratic form that stays constant along the
integrated geodesic — a first check that the ODE and the integrator agree.

### Cone

A cone is developable: cut it along a generator and unroll it to a flat
sector. Geodesics become straight lines in the sector. This is a beautiful
check of the intrinsic nature of geodesics: $K = 0$, and the "shortest
path" is a straight line after unrolling.

---

## 3.6 How the code implements it

```lean
Geodesic.accel s st      -- (u'', v'') from the Γ terms
Geodesic.rhs s st        -- first-order system (u, v, u', v')
Geodesic.rk4Step s h st  -- one RK4 step
Geodesic.integrate s h n st
Geodesic.energy s st     -- ½ · dot qᵀ I dot q  (the conserved quadratic form)
Geodesic.speed s st      -- |J·dotq|  =  √(2 · that value)
```

`accel` is exactly, in vector form,

```math
\begin{bmatrix}\ddot u\\ \ddot v\end{bmatrix}
= -\begin{bmatrix}
\Gamma^u_{\;uu}\dot u^2 + 2\Gamma^u_{\;uv}\dot u\dot v + \Gamma^u_{\;vv}\dot v^2\\
\Gamma^v_{\;uu}\dot u^2 + 2\Gamma^v_{\;uv}\dot u\dot v + \Gamma^v_{\;vv}\dot v^2
\end{bmatrix}.
```

The demo integrates a geodesic on the paraboloid and prints a conserved
quadratic form $\tfrac12\dot q^{\mathsf T}I\,\dot q$ that stays constant to
the printed precision — the discrete analogue of the identity in §3.4.

---

## 3.7 Exercises

1. **Normal acceleration.** Show that a geodesic satisfies
   $n\cdot\gamma'' = II(d,d)$ where $d = \dot q$, and that its tangent-plane
   acceleration vanishes — then verify both numerically on the sphere.
2. **Variational algebra.** Fill in the symmetrization step in §3.3 and
   explain why contracting with $\dot q^i \dot q^j$ lets you replace
   $\partial_i g_{kj}$ by
   $\tfrac12(\partial_i g_{kj} + \partial_j g_{ki})$ inside $I$-derivative
   terms.
3. **Length vs quadratic functional.** Derive the Euler–Lagrange equation
   for $\ell = \int\sqrt{\dot q^{\mathsf T}I\,\dot q}\,ds$ and show it
   reduces to the geodesic equation when $|\gamma'|$ is constant.
4. **Polar straight lines.** Show that
   $r\cos(\theta - \theta_0) = d$ satisfies both polar geodesic equations.
5. **Sphere great circles.** Verify that
   $\gamma(s) = \cos s\,U + \sin s\,V$ with two orthonormal arrows $U,V$
   satisfies the sphere geodesic equations.
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
