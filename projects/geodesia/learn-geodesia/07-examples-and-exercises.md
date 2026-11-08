# 7️⃣ Worked examples, cheat sheet and exercises

This chapter ties the pieces together with concrete calculations. Every
formula here is implemented somewhere in `Geodesia/`, and every numerical
example can be reproduced with `lake exe geodesia` or the Python figure
script in `images/make_figs.py`.

---

## 7.1 Example A — the flat plane in polar coordinates

$$
\Phi(r,\theta)=(r\cos\theta,r\sin\theta,0).
$$

### Metric

$$
I=\begin{bmatrix}1&0\\0&r^2\end{bmatrix},
\qquad
I^{-1}=\begin{bmatrix}1&0\\0&1/r^2\end{bmatrix}.
$$

### Christoffel symbols

$$
\Gamma^r_{\;\theta\theta}=-r,
\qquad
\Gamma^\theta_{\;r\theta}=\Gamma^\theta_{\;\theta r}=\frac1r,
\qquad
\text{all others }0.
$$

### Geodesic equations

$$
\ddot r-r\dot\theta^2=0,
\qquad
\ddot\theta+\frac{2}{r}\dot r\dot\theta=0.
$$

### Verification: a straight line

A straight line at distance $d$ from the origin is

$$
r\cos(\theta-\theta_0)=d.
$$

Differentiate twice with respect to $s$ and use the two geodesic equations;
both are satisfied. The conserved quantity is

$$
r^2\dot\theta=\text{const},
$$

the angular momentum of the straight motion. This is the cleanest example of
**nonzero Christoffel symbols on a flat surface**: the symbols encode the
rotation of the polar basis, not curvature.

---

## 7.2 Example B — great circles on the unit sphere

$$
\Phi(\theta,\varphi)
=(\sin\theta\cos\varphi,\ \sin\theta\sin\varphi,\ \cos\theta).
$$

### Metric and symbols

$$
I=\begin{bmatrix}1&0\\0&\sin^2\theta\end{bmatrix},
$$

$$
\Gamma^\theta_{\;\varphi\varphi}=-\sin\theta\cos\theta,
\qquad
\Gamma^\varphi_{\;\theta\varphi}=\Gamma^\varphi_{\;\varphi\theta}
=\cot\theta.
$$

### Geodesic equations

$$
\ddot\theta-\sin\theta\cos\theta\,\dot\varphi^2=0,
$$

$$
\ddot\varphi+2\cot\theta\,\dot\theta\dot\varphi=0.
$$

### Verification: a great circle

Let $U,V$ be orthonormal vectors and set

$$
\gamma(s)=\cos s\,U+\sin s\,V.
$$

Then $\gamma''=-\gamma$, which is normal to the sphere, so the tangential
acceleration is zero: $\gamma$ is a geodesic. The conserved momentum is

$$
\sin^2\theta\,\dot\varphi=\text{const}.
$$

Meridians ($\varphi=\text{const}$) and the equator ($\theta=\pi/2$) are the
easiest special cases.

---

## 7.3 Example C — the paraboloid, shooting vs relaxation

$$
\Phi(u,v)=\bigl(u,v,0.3(u^2+v^2)\bigr).
$$

### Metric

$$
I=
\begin{bmatrix}
1+0.36u^2 & 0.36uv\\
0.36uv & 1+0.36v^2
\end{bmatrix}.
$$

### Curvature at $(1,0)$

From Chapter 1:

$$
I=\begin{bmatrix}1.36&0\\0&1\end{bmatrix},
\qquad
II=\begin{bmatrix}0.5145&0\\0&0.5145\end{bmatrix},
$$

$$
\kappa_1=0.3783,\quad
\kappa_2=0.5145,\quad
K=0.1946,\quad
H=0.4464.
$$

### Shooting

Source $a=(0.4,0.2)$, target $b=(1.2,0.9)$. The chord length is $1.2281$.
Damped Newton on $(\alpha,L)$ with unit speed:

| iteration | $\alpha$ | $L$ | $\|R\|$ |
|---:|---:|---:|---:|
| 0 | 0.67857208 | 1.41231284 | $1.26\times10^{-1}$ |
| 1 | 0.67976725 | 1.23649445 | $2.44\times10^{-3}$ |
| 2 | 0.67974028 | 1.23976961 | $8.87\times10^{-7}$ |
| 3 | 0.67974028 | 1.23977080 | $2.02\times10^{-14}$ |

Converged geodesic length: $L=1.23977$, tangent angle
$\alpha=0.67974$ rad.

### Relaxation

Start from a bent path (chord plus a perpendicular bump of amplitude $0.45$),
$N=24$ segments, and run length-gradient descent:

| iteration | polyline length |
|---:|---:|
| 0 | 1.5741 |
| 250 | 1.2751 |
| 500 | 1.2454 |
| 1000 | 1.23987 |
| 2000 | 1.23974 |
| final | 1.23974 |

The two methods agree on the geodesic length to about $3\times10^{-5}$
(polyline discretisation error). The final relaxed path deviates from the
parameter-space chord by up to $0.0915$; the geodesic is not the chord on a
curved surface.

### Code

```bash
cd projects/geodesia
lake build
lake exe geodesia
```

Expected exact-derivative demo output (abridged):

```
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
unit sphere at (θ, φ) = (1.000000, 0.200000)
  K (should be 1) = 1.000000
```

---

## 7.4 Example D — height fields and graph formulas

For $\Phi(x,y)=(x,y,f(x,y))$:

$$
I=
\begin{bmatrix}
1+f_x^2 & f_xf_y\\
f_xf_y & 1+f_y^2
\end{bmatrix},
$$

$$
II=\frac{1}{\sqrt{1+|\nabla f|^2}}
\begin{bmatrix}
f_{xx}&f_{xy}\\
f_{xy}&f_{yy}
\end{bmatrix},
$$

$$
K=\frac{f_{xx}f_{yy}-f_{xy}^2}{(1+|\nabla f|^2)^2},
$$

$$
H=\frac{(1+f_y^2)f_{xx}-2f_xf_yf_{xy}+(1+f_x^2)f_{yy}}
{2(1+|\nabla f|^2)^{3/2}}.
$$

These are the formulas behind `Surface.gaussK`, `Surface.meanH`,
`Surface.metricII`, and the exact-derivative `Smooth2` fields. For a height
field, the normal is

$$
n=\frac{(-f_x,-f_y,1)}{\sqrt{1+|\nabla f|^2}}.
$$

A useful sanity check: a plane $f=ax+by+c$ has zero Hessian, hence
$II=0$, $K=H=0$; a paraboloid $f=\tfrac12 k(x^2+y^2)$ has constant Hessian
$kI$, so $II=kI/\sqrt{1+k^2r^2}$.

---

## 7.5 Cheat sheet

### Fundamental forms

$$
I=J^{\mathsf T}J,\qquad
J=[\Phi_u\ \Phi_v],\qquad
I_{ij}=\langle\Phi_i,\Phi_j\rangle,
$$

$$
II_{ij}=\langle n,\Phi_{ij}\rangle
=- \langle n_i,\Phi_j\rangle,
\qquad
n=\frac{\Phi_u\times\Phi_v}{|\Phi_u\times\Phi_v|}.
$$

### Christoffel

$$
\Gamma_{kij}
=\tfrac12(\partial_i g_{jk}+\partial_j g_{ik}-\partial_k g_{ij}),
\qquad
\Gamma^k_{\;ij}=g^{k\ell}\Gamma_{\ell ij}.
$$

### Gauss formula

$$
\Phi_{ij}=\Gamma^k_{\;ij}\Phi_k+II_{ij}n.
$$

### Covariant derivative and geodesic

$$
\nabla_i V^k=\partial_i V^k+\Gamma^k_{\;ij}V^j,
$$

$$
\ddot q^k+\Gamma^k_{\;ij}\dot q^i\dot q^j=0.
$$

### Lagrangian

$$
L=\tfrac12 g_{ij}\dot q^i\dot q^j,
\qquad
\frac{d}{ds}\frac{\partial L}{\partial\dot q^k}
-\frac{\partial L}{\partial q^k}=0.
$$

### Curvature

$$
S=I^{-1}II,\qquad
K=\det S=\frac{\det II}{\det I},\qquad
H=\tfrac12\operatorname{tr}S,
$$

$$
\kappa_{1,2}=H\pm\sqrt{H^2-K}.
$$

### Numerical

$$
v_{k+1}=v_k-[DR(v_k)]^{-1}R(v_k),\qquad
R(v)=q(L;a,v)-b.
$$

$$
\frac{\partial q^k}{\partial t}
=q''^k+\Gamma^k_{\;ij}q'^i q'^j
\quad\text{(relaxation PDE)}.
$$

---

## 7.6 Exercises with hints

### Exercise 1 — cone

$\Phi(r,\theta)=(r\cos\theta,r\sin\theta,r)$. Compute $I$, $II$, $K$ and the
Christoffel symbols. Show that unrolling the cone turns geodesics into
straight lines.

*Hint:* $I=\mathrm{diag}(2,r^2)$ if the cone has slope $1$. The cone is
developable, so $K=0$.

### Exercise 2 — cylinder geodesics

$\Phi(\theta,z)=(R\cos\theta,R\sin\theta,z)$. Show that geodesics are helices
(including circles and vertical lines as degenerate cases). Which
Christoffel symbols are nonzero?

*Hint:* $I=\mathrm{diag}(R^2,1)$, and the geodesic equations are
$\ddot\theta=0$, $\ddot z=0$.

### Exercise 3 — Clairaut

For a surface of revolution
$\Phi(r,\theta)=(r\cos\theta,r\sin\theta,f(r))$, prove that
$r^2\dot\theta=\text{const}$ along a geodesic.

*Hint:* the metric does not depend on $\theta$; use Noether’s theorem.

### Exercise 4 — sphere antipodes

Show that all unit-speed geodesics from the north pole reach the south pole
at $s=\pi$. What does this imply for the shooting Jacobian?

*Hint:* great circles through the poles.

### Exercise 5 — variational Jacobian

Implement the variational equation for the paraboloid and compare
$DR=Y(L)$ with a finite-difference Jacobian. How close are they?

### Exercise 6 — relaxation Hessian

Derive the block tridiagonal Hessian of the discrete Dirichlet energy
$E=\sum\frac12|\Delta q_i|^2$. Write the linear system for a Newton step.

*Hint:* the Hessian is the discrete Laplacian with Dirichlet boundary
conditions.

### Exercise 7 — resampling

Implement arc-length resampling for a relaxed polyline. Why does it improve
convergence?

### Exercise 8 — curvature sign

For the unit sphere, compute $II$ with the inward normal. How do $II$, $H$,
$\kappa_i$ change? Which quantity is invariant?

### Exercise 9 — graph curvature

Verify the graph formula for $K$ by direct computation for
$f(x,y)=x^2-y^2$ (a saddle). Where is $K$ negative?

### Exercise 10 — code

Add a `geodesicLength` function to `Geodesic.lean` that sums
$\sqrt{g_{ij}\Delta q^i\Delta q^j}$ along an integrated geodesic, and check
that it matches `length` in `shootState`.

---

## 7.7 Extensions

- **Higher-dimensional manifolds.** Everything generalises: $I$ is $n\times n$,
  $J$ is $m\times n$, $\Gamma^k_{\;ij}$ has three indices, and the geodesic
  state has $2n$ components. Eberly’s Riemannian note is written for this
  case.
- **Anisotropic metrics.** If the metric is not $J^{\mathsf T}J$ but an
  arbitrary positive-definite matrix (e.g. a Finsler/Riemannian model), the
  same formulas apply with $g_{ij}$ given directly.
- **Meshes.** On a triangle mesh, replace derivatives by per-vertex estimates;
  see Eberly’s *Mesh Differential Geometry* and *Geodesics On Triangle
  Meshes*.
- **General relativity.** Geodesics in spacetime use the same equation with a
  Lorentzian metric; timelike, spacelike and null geodesics correspond to
  $E>0$, $E<0$ and $E=0$.
- **Optimal control.** The geodesic equation is the Euler–Lagrange equation
  of a free particle; adding a potential gives forced motion, and the same
  variational machinery produces the control equations.

---

## 7.8 Where to go next

1. Re-derive the Gauss formula and the metric formula for $\Gamma$ on paper.
2. Implement the variational Jacobian and a Newton shooting solver.
3. Implement hierarchical relaxation and compare it with gradient descent.
4. Try a surface with negative curvature (a saddle or a catenoid) and watch
   geodesics diverge.
5. Read Eberly’s *Computing Geodesics on a Riemannian Manifold* §2–§5 with
   these notes as a companion.
