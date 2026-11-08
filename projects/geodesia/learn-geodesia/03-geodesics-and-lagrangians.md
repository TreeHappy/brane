# 3️⃣ Geodesics, Lagrangians and the physics of free motion

> **Goal.** See the geodesic equation as **zero covariant acceleration**,
> derive it a second time from a Lagrangian, and recognise the Christoffel
> terms as the fictitious forces of curvilinear coordinates.

---

## 3.1 Three equivalent definitions

A geodesic on a surface can be characterised in three ways:

1. **Geometric.** A curve with zero *covariant* acceleration:
   $\nabla_{\dot\gamma}\dot\gamma=0$.
2. **Variational.** A critical point of the length functional
   $\ell[\gamma]=\int |\gamma'|\,ds$ (or of the energy functional
   $\int \tfrac12|\gamma'|^2ds$).
3. **Physical.** The trajectory of a free particle constrained to the
   surface (no tangential force; the surface may exert a normal force).

They all lead to the same ODE. The first is the cleanest conceptually; the
second is the most useful computationally; the third is the most intuitive.

> **Caveat.** A geodesic is *locally* length-minimising. It need not be the
> global shortest path: on a sphere, the long arc of a great circle is also a
> geodesic, and between antipodal points every great circle is a geodesic of
> the same length.

---

## 3.2 Zero covariant acceleration

Let $\gamma(s)=\Phi(q(s))$ with $q(s)=(u(s),v(s))$. Differentiate twice:

```math
\gamma'=\dot q^i\Phi_i,
```

```math
\gamma''=\ddot q^k\Phi_k+\dot q^i\dot q^j\Phi_{ij}.
```

Insert the Gauss formula $\Phi_{ij}=\Gamma^k_{\;ij}\Phi_k+II_{ij}n$:

```math
\gamma''
=\bigl(\ddot q^k+\Gamma^k_{\;ij}\dot q^i\dot q^j\bigr)\Phi_k
+II_{ij}\dot q^i\dot q^j\,n.
```

The **covariant acceleration** is the tangent part:

```math
\frac{D\dot\gamma}{ds}
=\bigl(\ddot q^k+\Gamma^k_{\;ij}\dot q^i\dot q^j\bigr)\Phi_k.
```

A geodesic has zero covariant acceleration, so

```math
\boxed{\;\ddot q^k+\Gamma^k_{\;ij}\dot q^i\dot q^j=0\;}
\qquad k=u,v.
\tag{geodesic equation}
```

The normal part is still there:

```math
n\cdot\gamma''=II_{ij}\dot q^i\dot q^j.
```

So a geodesic’s acceleration is **purely normal to the surface**. That is
exactly what “the surface pushes back” means: the surface supplies whatever
normal force is needed to keep the curve on it, but there is no tangential
force. On a sphere, geodesics are great circles and their acceleration points
to the centre — normal to the sphere.

If $|\gamma'|=1$, the normal curvature of the geodesic is
$|n\cdot\gamma''|=|II(\xi,\xi)|$, while its **geodesic curvature**
(tangential bending) is zero.

---

## 3.3 The Lagrangian derivation

This is the derivation that makes the Christoffel formula unsurprising.

The **kinetic energy** of a particle moving in parameter space is

```math
T=\tfrac12\,g_{ij}(q)\,\dot q^i\dot q^j,
\qquad
L=T
```

(no potential: free motion). The Euler–Lagrange equations are

```math
\frac{d}{ds}\frac{\partial L}{\partial\dot q^k}
-\frac{\partial L}{\partial q^k}=0.
```

Compute the pieces:

```math
\frac{\partial L}{\partial\dot q^k}=g_{kj}\dot q^j,
```

```math
\frac{\partial L}{\partial q^k}
=\tfrac12\,\partial_k g_{ij}\,\dot q^i\dot q^j.
```

Therefore

```math
\frac{d}{ds}\bigl(g_{kj}\dot q^j\bigr)
-\tfrac12\,\partial_k g_{ij}\,\dot q^i\dot q^j=0,
```

i.e.

```math
g_{kj}\ddot q^j
+\partial_i g_{kj}\,\dot q^i\dot q^j
-\tfrac12\,\partial_k g_{ij}\,\dot q^i\dot q^j=0.
```

The middle term is symmetric in $i,j$ after contraction with
$\dot q^i\dot q^j$, so replace it by its symmetrisation:

```math
\partial_i g_{kj}\,\dot q^i\dot q^j
=\tfrac12\bigl(\partial_i g_{kj}+\partial_j g_{ki}\bigr)\dot q^i\dot q^j.
```

Then

```math
g_{kj}\ddot q^j
+\tfrac12\bigl(\partial_i g_{kj}+\partial_j g_{ki}
-\partial_k g_{ij}\bigr)\dot q^i\dot q^j=0.
```

Contract with $g^{mk}$ and use the definition of the Christoffel symbol:

```math
\boxed{\;\ddot q^m+\Gamma^m_{\;ij}\dot q^i\dot q^j=0.\;}
```

The “mysterious” metric formula for $\Gamma$ is simply what you get when you
move the metric derivatives in the Euler–Lagrange equation to the other side.

---

## 3.4 Energy, arc length and affine parameters

Because $L$ does not depend explicitly on $s$, the **energy**

```math
E=\tfrac12\,g_{ij}\dot q^i\dot q^j
```

is conserved. Check directly from the geodesic equation:

```math
\frac{dE}{ds}
=g_{ij}\dot q^i\Bigl(\ddot q^j+\Gamma^j_{\;k\ell}\dot q^k\dot q^\ell\Bigr)=0.
```

(The bracket is zero by the geodesic equation; the factor is just the
equation contracted with $\dot q$.)

Consequences:

- $|\gamma'|=\sqrt{2E}$ is constant. A geodesic is automatically
  **constant-speed** in an affine parameter.
- If you reparameterise $s\mapsto as+b$, the equation is unchanged. Any
  affine parameter is equally valid; arc length is the conventional choice.
- For the **length** functional
  $\ell=\int\sqrt{g_{ij}\dot q^i\dot q^j}\,ds$, the Euler–Lagrange equations
  are reparameterisation-invariant and only determine the *image* curve. If
  you then choose arc length, they reduce to the same geodesic equation.
  The energy functional is the computationally convenient “gauge-fixed”
  version.

This distinction matters in the numerics: shooting with a fixed “time”
$L$ is really using the energy action, and the value of $E$ fixes the speed.

---

## 3.5 Physics dictionary

| Geometry | Mechanics |
|---|---|
| $q^i$ | generalised coordinates |
| $g_{ij}(q)$ | mass / inertia matrix |
| $\tfrac12 g_{ij}\dot q^i\dot q^j$ | kinetic energy |
| $\Gamma^k_{\;ij}$ | fictitious force coefficients |
| $\ddot q^k=-\Gamma^k_{\;ij}\dot q^i\dot q^j$ | “coordinate force” needed for free motion |
| geodesic | free-particle trajectory |
| $II$ | normal force supplied by the constraint surface |
| conserved $p_i=g_{ij}\dot q^j$ | momentum from a symmetry |

### Fictitious forces

In a general coordinate system, a free particle obeys

```math
\ddot q^k=-\Gamma^k_{\;ij}\dot q^i\dot q^j.
```

The right-hand side looks like a force but is entirely due to the
coordinates. Polar coordinates make this vivid:

```math
\ddot r=r\dot\theta^2,
\qquad
\ddot\theta=-\frac{2}{r}\dot r\dot\theta.
```

The first is centrifugal, the second is the polar-coordinate analogue of the
Coriolis term. The particle is moving in a straight line; the coordinates are
rotating.

On a genuinely curved surface the same interpretation holds: the
$\Gamma$-terms are the coordinate description of inertia, while the normal
force is the physical constraint that keeps the particle on the surface.

### Noether momentum

If the metric does not depend on a coordinate $q^\ell$, then

```math
p_\ell=g_{\ell j}\dot q^j
```

is conserved. Examples:

- **Polar plane:** $\partial_\theta g=0$, so $r^2\dot\theta$ is conserved
  (angular momentum).
- **Sphere:** $\partial_\varphi g=0$, so $\sin^2\theta\,\dot\varphi$ is
  conserved.
- **Surface of revolution:** $\Phi(r,\theta)=(r\cos\theta,r\sin\theta,f(r))$
  has $\partial_\theta g=0$, so $r^2\dot\theta$ is conserved. This is
  **Clairaut’s relation** for geodesics on a surface of revolution.

This is also a practical numerical check: a good integrator should conserve
these momenta to high accuracy.

---

## 3.6 Hamiltonian form (optional, for physics readers)

The Legendre transform gives

```math
p_i=\frac{\partial L}{\partial\dot q^i}=g_{ij}\dot q^j,
\qquad
H=p_i\dot q^i-L=\tfrac12\,g^{ij}p_ip_j.
```

Hamilton’s equations are

```math
\dot q^i=\frac{\partial H}{\partial p_i}=g^{ij}p_j,
```

```math
\dot p_k=-\frac{\partial H}{\partial q^k}
=-\tfrac12\,\partial_k g^{ij}\,p_ip_j.
```

They are equivalent to the geodesic equation. The geodesic flow lives on the
cotangent bundle $T^*S$; this is the starting point for symplectic
integrators. For long integrations, a symplectic method often conserves
energy and momenta better than a generic RK4, at the cost of a more
complicated implementation.

---

## 3.7 Worked examples

### Plane in Cartesian coordinates

$g=I$, $\Gamma=0$, so $\ddot q=0$: geodesics are straight lines.

### Plane in polar coordinates

```math
\ddot r-r\dot\theta^2=0,\qquad
\ddot\theta+\frac{2}{r}\dot r\dot\theta=0.
```

These are straight lines written in polar coordinates. The conserved
$r^2\dot\theta$ is the angular momentum of the straight motion.

### Sphere

```math
\ddot\theta-\sin\theta\cos\theta\,\dot\varphi^2=0,
\qquad
\ddot\varphi+2\cot\theta\,\dot\theta\dot\varphi=0.
```

Solutions are great circles. Meridians ($\varphi=\text{const}$) give
$\ddot\theta=0$; the equator ($\theta=\pi/2$) gives $\ddot\varphi=0$.
The conserved momentum is $\sin^2\theta\,\dot\varphi$.

### Paraboloid

For $z=f(x,y)$ the metric is
$I=\begin{bmatrix}1+f_x^2&f_xf_y\\f_xf_y&1+f_y^2\end{bmatrix}$, and the
geodesic equations are messy but mechanical. There is no closed-form solution
in general, which is precisely why the numerical methods of Chapters 5–6
exist. The code prints an energy that stays constant along the integrated
geodesic — a first check that the ODE and the integrator agree.

### Cone

A cone is developable: cut it along a generator and unroll it to a flat
sector. Geodesics become straight lines in the sector. This is a beautiful
check of the intrinsic nature of geodesics: $K=0$, and the “shortest path”
is a straight line after unrolling.

---

## 3.8 How the code implements it

```lean
Geodesic.accel s st      -- (u'', v'') from Γ
Geodesic.rhs s st        -- first-order system (u, v, u', v')
Geodesic.rk4Step s h st  -- one RK4 step
Geodesic.integrate s h n st
Geodesic.energy s st     -- ½ I(ξ, ξ)
Geodesic.speed s st      -- √(2E)
```

`accel` is exactly

```math
u''=-(\Gamma^u_{\;uu}\dot u^2+2\Gamma^u_{\;uv}\dot u\dot v
+\Gamma^u_{\;vv}\dot v^2),
```

```math
v''=-(\Gamma^v_{\;uu}\dot u^2+2\Gamma^v_{\;uv}\dot u\dot v
+\Gamma^v_{\;vv}\dot v^2).
```

The demo integrates a geodesic on the paraboloid and prints
`E(0) = E(1) = 0.68`, which is the discrete analogue of energy conservation.

---

## 3.9 Exercises

1. **Normal acceleration.** Show that a geodesic satisfies
   $n\cdot\gamma''=II(\dot q,\dot q)$ and that its tangential acceleration is
   zero.
2. **Lagrangian algebra.** Fill in the symmetrisation step in §3.3 and
   explain why $\dot q^i\dot q^j$ lets you replace
   $\partial_i g_{kj}$ by
   $\tfrac12(\partial_i g_{kj}+\partial_j g_{ki})$.
3. **Length vs energy.** Derive the Euler–Lagrange equation for
   $\ell=\int\sqrt{g_{ij}\dot q^i\dot q^j}\,ds$ and show that it reduces to
   the geodesic equation when $|\dot q|$ is constant.
4. **Polar straight lines.** Show that
   $r\cos(\theta-\theta_0)=d$ (a straight line) satisfies both polar
   geodesic equations.
5. **Sphere great circles.** Verify that
   $\gamma(s)=\cos s\,U+\sin s\,V$ with orthonormal $U,V$ satisfies the
   sphere geodesic equations.
6. **Noether.** Prove that $\partial_\ell g=0$ implies
   $p_\ell=g_{\ell j}\dot q^j$ is constant along a geodesic.
7. **Surface of revolution.** Derive the metric for
   $\Phi(r,\theta)=(r\cos\theta,r\sin\theta,f(r))$ and Clairaut’s relation.
8. **Affine parameter.** Show that if $q(s)$ solves the geodesic equation,
   then $q(as+b)$ also solves it. What happens to the energy?
9. **Code check.** Run `lake exe geodesia` and verify that the energy is
   unchanged after 100 RK4 steps. Then change the step size and see how the
   error behaves.
