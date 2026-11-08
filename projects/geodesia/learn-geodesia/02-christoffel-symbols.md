# 2️⃣ Christoffel symbols — the tangent part of the second derivatives

> **Goal.** Understand $\Gamma^k_{\;ij}$ as pure vector algebra: the
> coefficients you get when you express each second-derivative arrow
> $\Phi_{uu}, \Phi_{uv}, \Phi_{vv}$ in the tangent basis
> $\{\Phi_u,\Phi_v\}$, plus the dot-product formula that computes them.

---

## 2.1 The problem: the basis arrows change from point to point

In the $xy$-plane the basis arrows $e_1=(1,0)$, $e_2=(0,1)$ are the same
everywhere, so differentiating a vector field is a component-wise operation.
On a surface — or even in the plane using *other* coordinates — the basis
arrows themselves depend on where you are.

Recall from Chapter 1 that the tangent basis is made of the two columns of
the Jacobian, $\Phi_u$ and $\Phi_v$. Differentiate these arrows: each
derivative is again a 3-D arrow, and we can split it into a part inside the
tangent plane and a part along the unit normal $n$:

```math
\boxed{\;
\Phi_{uu} = \Gamma^u_{\;\color{red}{uu}}\,\Phi_u
          + \Gamma^v_{\;uu}\,\Phi_v + II_{11}\,n,
\qquad
\Phi_{uv} = \Gamma^u_{\;uv}\,\Phi_u + \Gamma^v_{\;uv}\,\Phi_v + II_{12}\,n,
\qquad
\Phi_{vv} = \Gamma^u_{\;vv}\,\Phi_u + \Gamma^v_{\;vv}\,\Phi_v + II_{22}\,n.
\;}
```

The normal pieces are exactly the entries of $II$ from Chapter 1 (just dot
with $n$ to see it). The tangent pieces, expressed in the basis
$\{\Phi_u,\Phi_v\}$, are the **Christoffel symbols** $\Gamma^k_{\;ij}$:

> $\Gamma^k_{\;ij}$ = the coefficient of the basis arrow $\Phi_k$ in the
> tangent-plane part of the second-derivative arrow $\Phi_{ij}$.

In words: "how much of the $\Phi_u$ arrow reproduces itself when I move
$\Phi_u$ in the $u$-direction", and so on. It measures the turning of the
basis arrows (a coordinate effect), *not* the bending of the surface (that
is $II$).

---

## 2.2 Computing the coefficients by dot products

How do you find the coefficients of an arrow expressed in a basis? You take
dot products with the basis arrows and invert the Gram matrix — the same
least-squares move as in Chapter 1.

Write the unknown expansion $\Phi_{ij} = a\,\Phi_u + b\,\Phi_v + II_{ij}\,n$.
Dot with $\Phi_u$ and $\Phi_v$ (the $n$ term drops out because
$n\perp\Phi_u,\Phi_v$):

```math
\Phi_u\cdot\Phi_{ij} = a\,(\Phi_u\cdot\Phi_u) + b\,(\Phi_v\cdot\Phi_u),\qquad
\Phi_v\cdot\Phi_{ij} = a\,(\Phi_u\cdot\Phi_v) + b\,(\Phi_v\cdot\Phi_v).
```

That is a $2\times 2$ linear system whose matrix is $I = J^{\mathsf T}J$:

```math
\boxed{\;
\Gamma^k_{\;ij} = \bigl(I^{-1}\bigr)^{k\ell}\,
\bigl(\Phi_\ell\cdot\Phi_{ij}\bigr).\;}
```

This is the cleanest definition: *invert the tangent-basis Gram matrix and
multiply by the dot products of the second-derivative arrows with the basis
arrows.*

### The Gram-derivative formula

There is a second formula, computable from $I$ alone. Differentiate the dot
products that make up $I$:

```math
\frac{\partial}{\partial u}\,(\Phi_j\cdot\Phi_k)
= \Phi_{uj}\cdot\Phi_k + \Phi_j\cdot\Phi_{uk}.
```

Every $\Phi_{\bullet\bullet}\cdot\Phi_\bullet$ on the right is a
"first-kind" coefficient. Combine the three such differentiated identities
(cyclically in $i,j,k$) and eliminate the shared terms:

```math
\boxed{\;
\Phi_m\cdot\Phi_{ij}
= \tfrac12\Bigl(
\frac{\partial}{\partial u_i}\,(\Phi_j\cdot\Phi_m)
+ \frac{\partial}{\partial u_j}\,(\Phi_i\cdot\Phi_m)
- \frac{\partial}{\partial u_m}\,(\Phi_i\cdot\Phi_j)\Bigr)\;}
```

where $u_1=u$, $u_2=v$. Combining with the boxed formula above gives the
all-metric version

```math
\Gamma^k_{\;ij}
= \tfrac12 \bigl(I^{-1}\bigr)^{k\ell}
\Bigl(
\partial_i\, I_{j\ell} + \partial_j\, I_{i\ell}
- \partial_\ell\, I_{ij}\Bigr),
```

with $I_{ij} = \Phi_i\cdot\Phi_j$. This is the formula implemented in
`Curvature.lean`, and it needs no normals at all — only the Gram matrix and
its derivatives. The code computes the $\partial I$ terms **exactly** from
stored second partials:

```math
\partial_u\,(\Phi_j\cdot\Phi_k)
= \Phi_{uj}\cdot\Phi_k + \Phi_j\cdot\Phi_{uk},
```

see `Surface.dMetricEntry`.

### Both formulas agree

Because $\Phi_\ell \cdot (\Phi_{ij} - (\text{tangent part})) = 0$ by
construction, the two routes — "project the arrow" and "differentiate the
Gram matrix" — produce the same numbers. Use whichever is convenient;
Chapter 1's `Surface.christoffel` can be read through either lens.

---

## 2.3 How symbols turn into a corrected derivative

Take a surface vector field $V = V^u\,\Phi_u + V^v\,\Phi_v$ (two scalar
coefficient functions plus basis arrows). Differentiate in the $i$-th
parameter direction:

```math
\frac{\partial V}{\partial u_i}
= \frac{\partial V^k}{\partial u_i}\,\Phi_k
+ V^k\,\Phi_{ki}
= \Bigl(\frac{\partial V^k}{\partial u_i}
+ \Gamma^k_{\;ij}\,V^j\Bigr)\Phi_k
+ (\text{normal terms}).
```

The normal terms are discarded (they are $II$'s business). What remains,
with the correction term added back in, is the **covariant tangent part**:

```math
\boxed{\;\nabla_i V^k = \frac{\partial V^k}{\partial u_i}
+ \Gamma^k_{\;ij}\,V^j.\;}
```

Read it as: the honest rate of change of the arrows, minus the part caused
merely by the basis arrows themselves rotating. On the flat plane in
standard coordinates the $\Gamma$ all vanish — and indeed no correction is
needed there. In curved coordinates the same plane needs corrections; a
sphere needs them everywhere. This is the correction that lets
neighboring tangent planes be compared.

---

## 2.4 The Jacobian / directional derivative view

Collect everything into matrices. The derivative of the Jacobian in the
$i$-th parameter direction is a new $3\times 2$ matrix
$[\,\Phi_{ui}\;\ \Phi_{vi}\,]$, and each column of it splits:

```math
(\text{column } j)\ = \
\Gamma^u_{\;ij}\,\Phi_u + \Gamma^v_{\;ij}\,\Phi_v + II_{j i}\,n.
```

So

- $I$, $II$ and $\Gamma$ are all read off the **same object**, the derivative
  of the Jacobian;
- $\Gamma$ is its tangent-plane projection written back in the tangent basis,
  $II$ is its normal part.

If you remember one sentence: *the Christoffel symbols are the tangent-plane
shadow of the second derivatives of the surface, written in the tangent
basis.*

---

## 2.5 Example: the flat plane in polar coordinates

This example separates **coordinate twist** from genuine curvature, because
the surface is just a flat plane.

```math
\Phi(r,\theta) = (r\cos\theta,\ r\sin\theta,\ 0).
```

The basis arrows are

```math
\Phi_r = (\cos\theta,\ \sin\theta,\ 0),
\qquad
\Phi_\theta = (-r\sin\theta,\ r\cos\theta,\ 0),
```

and their Gram matrix is diagonal:
$I = \begin{bmatrix}1&0\\0&r^2\end{bmatrix}$. Differentiate the arrows:

```math
\frac{\partial\Phi_r}{\partial\theta} = \frac1r\,\Phi_\theta,
\qquad
\frac{\partial\Phi_\theta}{\partial r} = \frac1r\,\Phi_\theta,
\qquad
\frac{\partial\Phi_\theta}{\partial\theta} = -r\,\Phi_r.
```

Reading off the coefficients in the basis $\{\Phi_r,\Phi_\theta\}$:

```math
\Gamma^\theta_{\;\theta r} = \Gamma^\theta_{\;r\theta} = \frac1r,
\qquad
\Gamma^r_{\;\theta\theta} = -r.
```

Cross-check with the Gram-derivative formula: the only nonzero Gram
derivative is $\partial_r\, (\Phi_\theta\cdot\Phi_\theta) = 2r$, and
$\tfrac12 (I^{-1})^{rr} (\partial_r I_{\theta\theta}) \ldots$ reproduces
exactly the same two symbols. Good — the formulas agree, as they must.

As a first sanity check that this is coordinate machinery, not bending: the
plane is flat ($K=0$ from Chapter 1), yet its polar symbols are nonzero.
The symbols record that the *basis arrows* rotate as you walk; they say
nothing about the *plane itself* being bent — it isn't.

---

## 2.6 Example: the sphere

For the unit sphere
$\Phi(\theta,\varphi) = (\sin\theta\cos\varphi,\ \sin\theta\sin\varphi,\ \cos\theta)$,
the tangent
arrows have Gram matrix
$I = \begin{bmatrix}1&0\\0&\sin^2\theta\end{bmatrix}$, and computing the
second derivatives gives

```math
\Gamma^\theta_{\;\varphi\varphi} = -\sin\theta\cos\theta,
\qquad
\Gamma^\varphi_{\;\theta\varphi}
= \Gamma^\varphi_{\;\varphi\theta} = \cot\theta.
```

Their role: a curve that keeps $\varphi$ fixed is automatically a geodesic
(its second-derivative arrow is purely along $n$, no tangent part); a curve
with $\theta = \pi/2$ is one too; general great circles follow from combining
the two.

---

## 2.7 Example: the paraboloid

For $z = 0.3(x^2+y^2)$ the Gram matrix of the tangent arrows is

```math
I = \begin{bmatrix}1+0.36x^2 & 0.36xy\\ 0.36xy & 1+0.36y^2\end{bmatrix}.
```

At $(1,0)$ the code prints

```math
I=\begin{bmatrix}1.36&0\\0&1\end{bmatrix},
\qquad
\Gamma^u_{\;uu}=0.2647,\quad
\Gamma^u_{\;uv}=0,\quad
\Gamma^u_{\;vv}=0.2647,
```

and all $\Gamma^v_{\;ij}=0$. Why is $\Gamma^u_{\;vv}$ nonzero? The
$v$-coordinate curve $x = \text{const}$ has second-derivative arrow
$(0,0,0.6)$ at that point, which is not parallel to $n=(\,-0.6,0,1)/\sqrt{1.36}$;
its tangent-plane component therefore has a piece along the
$\Phi_u=(1,0,0.6)$ arrow — the coefficient of that piece is
$\Gamma^u_{\;vv} = 0.2647$. The normal part of the very same arrow is the
$II_{22}$ entry from Chapter 1. One arrow, two projections — that repeating
pattern is the whole subject.

---

## 2.8 How the code computes $\Gamma$

```lean
Surface.dPhi  s i x y   -- the basis arrows Φ_i
Surface.d2Phi s i j x y -- the second-derivative arrows Φ_ij
Surface.dMetricEntry s u i j x y
                        -- exact Gram-matrix derivative ∂_u (Φ_i·Φ_j)
Surface.christoffel s x y
                        -- Γ^k_ij : Fin 2 → Fin 2 → Fin 2 → Float
Surface.metricI s x y   -- the Gram matrix I
```

`christoffel` is a literal transcription of the Gram-derivative formula:

```math
\Gamma^k_{\;ij}
= \tfrac12 \bigl(I^{-1}\bigr)^{k0}
\bigl(\partial_i I_{j0} + \partial_j I_{i0} - \partial_0 I_{ij}\bigr)
+ \tfrac12 \bigl(I^{-1}\bigr)^{k1}
\bigl(\partial_i I_{j1} + \partial_j I_{i1} - \partial_1 I_{ij}\bigr).
```

Because `Smooth2` stores exact partials, no finite differences are
involved. `Geodesic.accel` then forms
$-\Gamma^k_{\;ij}\,d^i d^j$ for a direction vector $d=(\dot u,\dot v)$ —
see the next chapter.

---

## 2.9 Exercises

1. **Polar, both ways.** Compute $\Gamma^r_{\;\theta\theta} = -r$ and
   $\Gamma^\theta_{\;r\theta}=1/r$ from the basis-derivative arrows in
   §2.5, and check them using the Gram-derivative formula with
   $I=\mathrm{diag}(1,r^2)$.
2. **Standard coordinates.** Show that every
   $\Gamma^k_{\;ij}=0$ for the plane in standard $(x,y)$ coordinates. Why
   does this not contradict the nonzero polar symbols?
3. **Symmetry.** Use $\Phi_{uv}=\Phi_{vu}$ to explain why
   $\Gamma^k_{\;ij}=\Gamma^k_{\;ji}$ always holds.
4. **First-kind identity.** Prove that
   $\partial_k\,(\Phi_i\cdot\Phi_j) = \Phi_{ik}\cdot\Phi_j + \Phi_i\cdot\Phi_{jk}$,
   and use it to derive the boxed "Gram-derivative formula" of §2.2.
5. **Cylinder.** Compute the basis-arrow derivatives for
   $\Phi(\theta,z)=(R\cos\theta,R\sin\theta,z)$ and show that only
   $\Gamma^z_{\;\theta\theta}$-type symbols could exist — and actually all
   of them vanish. Is that consistent with the cylinder being flat?
6. **Sphere check.** With $I=\mathrm{diag}(1,\sin^2\theta)$, derive
   $\Gamma^\theta_{\;\varphi\varphi}=-\sin\theta\cos\theta$ and
   $\Gamma^\varphi_{\;\theta\varphi}=\cot\theta$ from the
   Gram-derivative formula.
7. **Conserved quantity.** Show that on the sphere the combination
   $\sin^2\theta\,\dot\varphi$ stays constant along any curve obeying the
   geodesic equations of Chapter 3. (Hint: differentiate
   $\sin^2\theta\,\dot\varphi$ with respect to $s$ and substitute the two
   geodesic equations.)
8. **Code experiment.** Evaluate `Surface.christoffel` on the paraboloid at
   $(0,0)$ and at $(2,0)$. Which symbols change, and can you predict the
   change from the Gram-matrix entries?
