# 2️⃣ Christoffel symbols — how a coordinate grid twists

> **Goal.** See $\Gamma^k_{\;ij}$ as the **tangent part of the derivative of
> the basis vectors**, derive the metric formula without magic, and recognise
> the same symbols as the “fictitious forces” of polar coordinates.

---

## 2.1 The problem: vectors live in different tangent planes

On a flat plane in Cartesian coordinates, differentiating a vector field is
easy: the basis $e_x,e_y$ is the same everywhere, so you differentiate the
components.

On a surface, or even on a plane in curvilinear coordinates, the basis
**changes from point to point**. The polar basis

$$
e_r=(\cos\theta,\sin\theta,0),\qquad
e_\theta=(-r\sin\theta,r\cos\theta,0)
$$

rotates as you move. If you differentiate a vector field
$V=V^r e_r+V^\theta e_\theta$ naively, you get

$$
\partial_\theta V
=(\partial_\theta V^r)e_r+(\partial_\theta V^\theta)e_\theta
+V^r\,\partial_\theta e_r+V^\theta\,\partial_\theta e_\theta.
$$

The last two terms are the **change of the basis itself**. Christoffel
symbols are exactly the coefficients of that change.

---

## 2.2 The moving-frame definition

Write the coordinate basis vectors as

$$
e_i=\Phi_i=\frac{\partial\Phi}{\partial q^i}.
$$

Their derivative in the $j$-direction is the second derivative of the
embedding:

$$
D_j e_i = \partial_j\Phi_i = \Phi_{ij}.
$$

Now decompose $\Phi_{ij}$ into a **tangent part** and a **normal part**:

$$
\boxed{\;\Phi_{ij}=\Gamma^k_{\;ij}\,e_k+II_{ij}\,n\;}
\tag{Gauss formula}
$$

- The normal coefficient is forced: dot with $n$ gives
  $II_{ij}=\langle n,\Phi_{ij}\rangle$.
- The tangent coefficients are the **Christoffel symbols of the second kind**:

$$
\Gamma^k_{\;ij}
=\langle \Phi_{ij},e_\ell\rangle\,g^{\ell k}.
$$

This is the geometric definition used by the code:
`Surface.d2Phi` supplies $\Phi_{ij}$, and `Surface.christoffel` projects it
onto the tangent basis using the inverse metric.

In words:

> $\Gamma^k_{\;ij}$ is **how much of $e_k$ appears when you move $e_i$ in the
> $j$-direction**. It measures the twisting of the coordinate grid, not the
> bending of the surface. The bending is $II$.

This is why the Gauss formula is the best single picture of the subject:

$$
\underbrace{\Phi_{ij}}_{\text{second derivative}}
=
\underbrace{\Gamma^k_{\;ij}e_k}_{\text{grid twist}}
+
\underbrace{II_{ij}n}_{\text{bending}}.
$$

---

## 2.3 First kind and the metric formula

Define the **Christoffel symbols of the first kind**

$$
\Gamma_{kij}=\langle \Phi_{ij},\Phi_k\rangle.
$$

They are symmetric in $i,j$ because mixed partials commute. The second kind
is obtained by raising the first index with the inverse metric:

$$
\Gamma^k_{\;ij}=g^{k\ell}\Gamma_{\ell ij}.
$$

### Deriving the metric formula

Start from the metric and differentiate:

$$
\partial_i g_{jk}
=\partial_i\langle \Phi_j,\Phi_k\rangle
=\langle \Phi_{ij},\Phi_k\rangle+\langle \Phi_j,\Phi_{ik}\rangle
=\Gamma_{kij}+\Gamma_{jik}.
$$

Now write the same identity with the indices cyclically permuted:

$$
\partial_i g_{jk}=\Gamma_{kij}+\Gamma_{jik},
$$

$$
\partial_j g_{ik}=\Gamma_{kji}+\Gamma_{ijk},
$$

$$
\partial_k g_{ij}=\Gamma_{jki}+\Gamma_{ijk}.
$$

Because $\Gamma_{kij}=\Gamma_{kji}$ (symmetry in the last two indices),
add the first two and subtract the third:

$$
\partial_i g_{jk}+\partial_j g_{ik}-\partial_k g_{ij}
=2\Gamma_{kij}.
$$

Therefore

$$
\boxed{\;
\Gamma_{kij}
=\frac12\bigl(\partial_i g_{jk}+\partial_j g_{ik}-\partial_k g_{ij}\bigr)
\;}
$$

and

$$
\boxed{\;
\Gamma^k_{\;ij}
=\frac12\,g^{k\ell}
\bigl(\partial_i g_{j\ell}+\partial_j g_{i\ell}-\partial_\ell g_{ij}\bigr).
\;}
$$

This is the formula in `Curvature.lean`. It is the only formula you need to
compute geodesics from a metric alone. The code computes
$\partial_i g_{jk}$ **exactly** from the stored second partials:

$$
\partial_i g_{jk}
=\langle \Phi_{ij},\Phi_k\rangle+\langle \Phi_j,\Phi_{ik}\rangle,
$$

via `dMetricEntry`.

### Why the tangent projection is the right correction

If $V=V^k e_k$ is a vector field, then

$$
D_i V
=\partial_i(V^k e_k)
=(\partial_i V^k)e_k+V^k\Phi_{ik}
=\bigl(\partial_i V^k+\Gamma^k_{\;ij}V^j\bigr)e_k
+II_{ik}V^k n.
$$

The **covariant derivative** is the tangent part:

$$
\boxed{\;\nabla_i V^k=\partial_i V^k+\Gamma^k_{\;ij}V^j\;}
$$

The extra $\Gamma$ term subtracts the rotation of the basis. This is exactly
the “connection” that lets you compare vectors in neighbouring tangent
planes.

The metric is compatible with this derivative:

$$
\nabla_k g_{ij}=0,
$$

which is another way to derive the formula above: impose
$\nabla_k g_{ij}=0$ and solve for $\Gamma$.

---

## 2.4 The Jacobian/directional-derivative view

Let $J=[\Phi_u\ \Phi_v]$. Moving in the parameter direction $j$ changes the
Jacobian:

$$
D_j J
=\bigl[\Phi_{uj}\ \ \Phi_{vj}\bigr].
$$

The Gauss formula says each column of this “derivative of the Jacobian” is

$$
D_j\Phi_i
=\Gamma^k_{\;ij}\Phi_k+II_{ij}n.
$$

So:

- **$\Gamma$** is the tangential part of $D J$;
- **$II$** is the normal part of $D J$.

If you only remember one sentence: *the Christoffel symbols are the
directional derivative of the tangent basis, projected back onto the tangent
plane.*

---

## 2.5 Example: the flat plane in polar coordinates

This is the most important example, because it separates **coordinate
twist** from **curvature**. The plane is flat ($K=0$), yet the polar
Christoffel symbols are not zero.

$$
\Phi(r,\theta)=(r\cos\theta,r\sin\theta,0).
$$

Basis:

$$
e_r=(\cos\theta,\sin\theta,0),\qquad
e_\theta=(-r\sin\theta,r\cos\theta,0).
$$

Derivatives:

$$
\partial_r e_r=0,\qquad
\partial_\theta e_r=(-r\sin\theta,r\cos\theta,0)=\frac1r e_\theta,
$$

$$
\partial_r e_\theta=(-\sin\theta,\cos\theta,0)=\frac1r e_\theta,
\qquad
\partial_\theta e_\theta=(-r\cos\theta,-r\sin\theta,0)=-r\,e_r.
$$

Therefore the nonzero symbols are

$$
\Gamma^\theta_{\;\theta r}
=\Gamma^\theta_{\;r\theta}=\frac1r,
\qquad
\Gamma^r_{\;\theta\theta}=-r.
$$

Metric and inverse:

$$
g=\begin{bmatrix}1&0\\0&r^2\end{bmatrix},
\qquad
g^{-1}=\begin{bmatrix}1&0\\0&1/r^2\end{bmatrix}.
$$

Check with the metric formula: the only nonzero metric derivative is
$\partial_r g_{\theta\theta}=2r$, hence

$$
\Gamma^r_{\;\theta\theta}
=\tfrac12 g^{rr}(-\partial_r g_{\theta\theta})
=-r,
$$

and

$$
\Gamma^\theta_{\;r\theta}
=\tfrac12 g^{\theta\theta}(\partial_r g_{\theta\theta})
=\frac1r.
$$

### Geodesic equations and the centrifugal force

The geodesic equation is

$$
\ddot q^k+\Gamma^k_{\;ij}\dot q^i\dot q^j=0.
$$

In polar coordinates:

$$
\ddot r-r\dot\theta^2=0,
\qquad
\ddot\theta+\frac{2}{r}\dot r\dot\theta=0.
$$

The first equation is the familiar **centrifugal term**: to stay on a
straight line, a particle in polar coordinates needs an inward coordinate
acceleration $-r\dot\theta^2$. The second equation says

$$
\frac{d}{ds}(r^2\dot\theta)
=r^2\Bigl(\ddot\theta+\frac{2}{r}\dot r\dot\theta\Bigr)=0,
$$

i.e. **conservation of angular momentum**. Straight lines in the plane are
the geodesics; the Christoffel terms are coordinate effects, not curvature.
Indeed the Riemann curvature of the plane is zero.

> **Lesson.** Nonzero $\Gamma$ does **not** mean the surface is curved. It
> means the coordinates are curved. Curvature is a combination of derivatives
> of $\Gamma$ (the Riemann tensor).

---

## 2.6 Example: the unit sphere

Use polar coordinates on the sphere:

$$
\Phi(\theta,\varphi)
=(\sin\theta\cos\varphi,\ \sin\theta\sin\varphi,\ \cos\theta).
$$

Metric:

$$
g=
\begin{bmatrix}
1&0\\
0&\sin^2\theta
\end{bmatrix}.
$$

The nonzero metric derivative is $\partial_\theta g_{\varphi\varphi}
=2\sin\theta\cos\theta$. The nonzero symbols are

$$
\Gamma^\theta_{\;\varphi\varphi}
=-\sin\theta\cos\theta,
\qquad
\Gamma^\varphi_{\;\theta\varphi}
=\Gamma^\varphi_{\;\varphi\theta}
=\cot\theta.
$$

The geodesic equations are

$$
\ddot\theta-\sin\theta\cos\theta\,\dot\varphi^2=0,
$$

$$
\ddot\varphi+2\cot\theta\,\dot\theta\dot\varphi=0.
$$

The second equation gives

$$
\frac{d}{ds}\bigl(\sin^2\theta\,\dot\varphi\bigr)=0,
$$

conservation of angular momentum about the polar axis. Solutions are great
circles: set $\varphi=\text{const}$ to get a meridian
($\ddot\theta=0$), or set $\theta=\pi/2$ to get the equator
($\ddot\varphi=0$). The sphere is genuinely curved: $K=1$.

---

## 2.7 Example: the paraboloid

For $z=0.3(x^2+y^2)$, the metric at $(x,y)$ is

$$
I=
\begin{bmatrix}
1+0.36x^2 & 0.36xy\\
0.36xy & 1+0.36y^2
\end{bmatrix}.
$$

At $(1,0)$ the code prints

$$
I=\begin{bmatrix}1.36&0\\0&1\end{bmatrix},
\qquad
\Gamma^u_{\;uu}=0.2647,\quad
\Gamma^u_{\;uv}=0,\quad
\Gamma^u_{\;vv}=0.2647,
$$

and all $\Gamma^v_{\;ij}=0$. The nonzero $\Gamma^u_{\;vv}$ is a good check of
the geometric picture: the $v$-coordinate line $x=\text{const}$ does not look
“straight” in the tilted tangent basis at $(1,0)$, so its second derivative
has a component along $\Phi_u$. The normal part of that same second
derivative is the bending $II_{vv}$.

---

## 2.8 How the code computes $\Gamma$

```lean
Surface.dPhi  s i x y   -- Φ_i
Surface.d2Phi s i j x y -- Φ_ij
Surface.dMetricEntry s u i j x y
                        -- ∂_u g_ij, exact
Surface.christoffel s x y
                        -- Γ^k_ij as a Fin 2 → Fin 2 → Fin 2 → Float
Surface.metricI s x y   -- I
```

`christoffel` is a literal transcription of

$$
\Gamma^k_{\;ij}
=\tfrac12 g^{k0}
\bigl(\partial_i g_{j0}+\partial_j g_{i0}-\partial_0 g_{ij}\bigr)
+\tfrac12 g^{k1}
\bigl(\partial_i g_{j1}+\partial_j g_{i1}-\partial_1 g_{ij}\bigr).
$$

Because `Smooth2` stores exact partials, no finite differences are involved.
`Geodesic.accel` then forms
$-\Gamma^k_{\;ij}\dot q^i\dot q^j$.

---

## 2.9 Exercises

1. **Polar from the metric.** Derive
   $\Gamma^r_{\;\theta\theta}=-r$ and
   $\Gamma^\theta_{\;r\theta}=1/r$ from
   $g=\mathrm{diag}(1,r^2)$ alone. Then check the basis-derivative
   computation in §2.5.
2. **Cartesian flatness.** Show that all $\Gamma^k_{\;ij}=0$ for the plane in
   Cartesian coordinates. Why does this not contradict the nonzero polar
   symbols?
3. **First-kind identity.** Prove
   $\partial_k g_{ij}=\Gamma_{ijk}+\Gamma_{jik}$ and use it to re-derive the
   first-kind formula.
4. **Cylinder.** Compute the Christoffel symbols for the cylinder
   $\Phi(\theta,z)=(R\cos\theta,R\sin\theta,z)$. Which are nonzero? Is the
   cylinder curved in the Gaussian sense?
5. **Sphere check.** Starting from
   $g=\mathrm{diag}(1,\sin^2\theta)$, derive
   $\Gamma^\theta_{\;\varphi\varphi}=-\sin\theta\cos\theta$ and
   $\Gamma^\varphi_{\;\theta\varphi}=\cot\theta$.
6. **Angular momentum.** Show that the polar geodesic equation implies
   $r^2\dot\theta=\text{const}$, and interpret this physically.
7. **Coordinate artifact.** Explain why $\Gamma$ is not a tensor. (Hint: in
   Cartesian coordinates it vanishes, but a tensor that vanishes in one
   coordinate system must vanish in all.)
8. **Code experiment.** Use `Surface.christoffel` on the paraboloid at
   $(0,0)$ and at $(2,0)$. Which symbols vanish/change? Relate the change to
   the metric entries.
