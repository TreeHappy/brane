# 1️⃣ The first and second fundamental forms — the ruler and the bending

> **Goal.** Understand exactly what `metricI`, `metricII`, `gaussK`, `meanH`
> and `principalK` in `Curvature.lean` are computing, and why the Jacobian
> $J$ is the object behind all of them.

![The tangent plane, the first fundamental form as ruler, and the split of a curve's acceleration into bending and geodesic-curvature parts.](images/fundamental-forms.png)

---

## 1.1 The parameterization and its Jacobian

A surface is a map

$$
\Phi : \mathbb{R}^2 \to \mathbb{R}^3,
\qquad \Phi(u,v) = \big(x(u,v),\,y(u,v),\,z(u,v)\big).
$$

Think of $(u,v)$ as **coordinates painted on the surface**. A curve in the
parameter plane, $s \mapsto (u(s),v(s))$, becomes a curve on the surface:

$$
\gamma(s) = \Phi\big(u(s),v(s)\big).
$$

Differentiate with the chain rule:

$$
\gamma'(s)
= \frac{\partial\Phi}{\partial u}\dot u
 + \frac{\partial\Phi}{\partial v}\dot v
= \underbrace{\begin{bmatrix}\Phi_u & \Phi_v\end{bmatrix}}_{J}
  \begin{bmatrix}\dot u\\ \dot v\end{bmatrix}.
$$

The **Jacobian** is the $3\times 2$ matrix

$$
J = \begin{bmatrix}\Phi_u & \Phi_v\end{bmatrix},
\qquad
\Phi_u=\frac{\partial\Phi}{\partial u},\quad
\Phi_v=\frac{\partial\Phi}{\partial v}.
$$

Two interpretations, both important:

- **Push-forward of a velocity.** $J$ maps parameter velocity
  $\xi=(\dot u,\dot v)$ to the ambient 3-D velocity $\gamma'=J\xi$.
- **A basis of the tangent plane.** The columns $\Phi_u,\Phi_v$ span the
  tangent plane $T_{\Phi(u,v)}S$. The tangent plane is literally the image of
  $J$.

> **Directional derivative.** If you walk in the parameter direction
> $\xi=(a,b)$, your 3-D direction on the surface is
> $$
> D_\xi\Phi = J\xi = a\,\Phi_u + b\,\Phi_v.
> $$
> This is the directional derivative of the embedding, and it is the only
> “directional derivative” you need for the first fundamental form.

### A graph is a special case

For a height field $z=f(x,y)$,

$$
\Phi(x,y)=(x,y,f(x,y)),
\qquad
J=
\begin{bmatrix}
1 & 0\\
0 & 1\\
f_x & f_y
\end{bmatrix}.
$$

This is the most common surface in the code: `Surface.ofHeight`.

---

## 1.2 The first fundamental form $I$ (the metric)

The **first fundamental form** is

$$
\boxed{\;I = J^{\mathsf T}J\;}
\qquad
I_{ij}=\langle \Phi_i,\Phi_j\rangle.
$$

Written out:

$$
I=
\begin{bmatrix}
E & F\\
F & G
\end{bmatrix}
=
\begin{bmatrix}
\langle\Phi_u,\Phi_u\rangle & \langle\Phi_u,\Phi_v\rangle\\
\langle\Phi_v,\Phi_u\rangle & \langle\Phi_v,\Phi_v\rangle
\end{bmatrix}.
$$

For a parameter velocity $\xi=(\dot u,\dot v)$:

$$
|\gamma'|^2
= (J\xi)^{\mathsf T}(J\xi)
= \xi^{\mathsf T}J^{\mathsf T}J\,\xi
= \xi^{\mathsf T}I\,\xi
= E\,\dot u^2 + 2F\,\dot u\dot v + G\,\dot v^2.
$$

So $I$ is the **ruler**: it converts coordinate velocity into actual length.
It also gives angles between tangent vectors $v,w$ via
$\langle v,w\rangle = v^{\mathsf T}I w$, and area via

$$
dA = \sqrt{\det I}\;du\,dv = |\Phi_u\times\Phi_v|\,du\,dv.
$$

### Recovering a 3-D tangent vector’s coordinates

Suppose someone gives you an ambient tangent vector $w\in T_pS$ (a 3-vector)
and asks for its parameter coordinates $\xi=(\dot u,\dot v)$. Since
$w=J\xi$ and $J^{\mathsf T}J=I$,

$$
\boxed{\;\xi = I^{-1}J^{\mathsf T}w\;}
$$

This is exactly the formula Eberly uses to convert an initial 3-D velocity
$y'(0)$ into initial parameter velocity $(du/ds,dv/ds)$:

$$
\begin{bmatrix}du/ds\\ dv/ds\end{bmatrix}
= G^{-1}J^{\mathsf T}y'(s),
\qquad G=I.
$$

It is implemented as `Geodesic.tangentVecToParams`.

### Graph metric

For $\Phi=(x,y,f)$,

$$
I=
\begin{bmatrix}
1+f_x^2 & f_xf_y\\
f_xf_y & 1+f_y^2
\end{bmatrix}.
$$

Check: a step $(dx,dy)$ in the plane lifts to $(dx,dy,f_xdx+f_ydy)$, whose
squared length is $dx^2+dy^2+(f_xdx+f_ydy)^2$, precisely the quadratic form
above.

---

## 1.3 The second fundamental form $II$ (the bending)

The first fundamental form sees only *intrinsic* geometry: distances and
angles measured by a bug living on the surface. To see how the surface bends
in the surrounding 3-D space, we need the **unit normal**

$$
n = \frac{\Phi_u\times\Phi_v}{|\Phi_u\times\Phi_v|}.
$$

The **second fundamental form** is

$$
\boxed{\;II_{ij}=\langle n,\Phi_{ij}\rangle\;}
\qquad
II=
\begin{bmatrix}
e & f\\
f & g_2
\end{bmatrix}
=
\begin{bmatrix}
\langle n,\Phi_{uu}\rangle & \langle n,\Phi_{uv}\rangle\\
\langle n,\Phi_{vu}\rangle & \langle n,\Phi_{vv}\rangle
\end{bmatrix}.
$$

### Why this measures bending

Take a curve $\gamma(s)=\Phi(q(s))$ on the surface, parameterized by arc
length so that $|\gamma'|=1$. Differentiate twice:

$$
\gamma' = q^i{}\Phi_i,
$$

$$
\gamma'' = \ddot q^i\,\Phi_i + \dot q^i\dot q^j\,\Phi_{ij}.
$$

Dot with the normal $n$. Since $n\cdot\Phi_i=0$,

$$
\boxed{\;n\cdot\gamma'' = II_{ij}\,\dot q^i\dot q^j\;}
$$

The normal component of the curve’s acceleration is the **normal curvature**

$$
\kappa_n
= \frac{n\cdot\gamma''}{|\gamma'|^2}
= \frac{II(\xi,\xi)}{I(\xi,\xi)},
\qquad \xi=(\dot u,\dot v).
$$

So $II$ answers: *how much does the surface leave its tangent plane as I walk
in direction $\xi$?* Positive/negative signs depend on the choice of $n$, but
the product of the two principal curvatures (Gaussian curvature) does not.

### Second equivalent formula

Differentiate $n\cdot\Phi_i=0$:

$$
n_j\cdot\Phi_i + n\cdot\Phi_{ij}=0
\quad\Longrightarrow\quad
II_{ij}=-n_i\cdot\Phi_j.
$$

This is the **Weingarten** view: $II$ is the derivative of the normal field,
with a minus sign. The shape operator below makes this precise.

### Local height above the tangent plane

Choose coordinates at a point so that $\Phi(0,0)=0$, the tangent plane is
$z=0$, and $n=(0,0,1)$. Taylor’s theorem gives

$$
z(u,v)=\tfrac12\,II_{ij}(0,0)\,u^i u^j + O(3).
$$

So the second fundamental form is literally the **quadratic approximation of
the surface above its tangent plane**. A bowl has definite $II$; a saddle has
mixed signs; a flat plane or a cylinder (in one direction) has a zero
eigenvalue.

---

## 1.4 The shape operator, principal curvatures, $K$ and $H$

The two forms are not independent: both are bilinear forms on the same
tangent plane. The linear map that compares them is the **shape operator**

$$
S = I^{-1}II,
\qquad
II(\xi,\eta)=I(S\xi,\eta).
$$

Its eigenvalues are the **principal curvatures** $\kappa_1,\kappa_2$ and its
eigenvectors are the **principal directions** — the directions of maximal and
minimal bending. In an orthonormal basis, $S$ is just the matrix of $II$.

The two classical scalar curvatures are

$$
K=\det S=\frac{\det II}{\det I}
\qquad\text{(Gaussian curvature)},
$$

$$
H=\tfrac12\operatorname{tr}S
=\tfrac12\operatorname{tr}(I^{-1}II)
\qquad\text{(mean curvature)}.
$$

They determine the principal curvatures:

$$
\kappa_{1,2}=H\pm\sqrt{H^2-K}.
$$

- $K>0$: bowl/ellipsoid (both directions bend the same way).
- $K<0$: saddle (directions bend oppositely).
- $K=0$: developable (plane, cylinder, cone); one principal curvature is zero.
- $H=0$: minimal surface (soap film), e.g. a catenoid.

### Weingarten map

The derivative of the normal in direction $\xi$ is

$$
D_\xi n = -S\xi
$$

(with our outward-normal sign convention). So the shape operator *is* the
derivative of the Gauss map $n:S\to S^2$. That is the cleanest geometric
meaning of $II$: “how fast does the normal turn when I move?”

---

## 1.5 The Gauss formula — the bridge to Christoffel symbols

Decompose the second derivative $\Phi_{ij}$ into a tangent part and a normal
part:

$$
\boxed{\;\Phi_{ij}
= \Gamma^k_{\;ij}\,\Phi_k + II_{ij}\,n\;}
$$

This is the **Gauss formula**. Dotting with $n$ recovers
$II_{ij}=\langle n,\Phi_{ij}\rangle$. Dotting with $\Phi_k$ gives the
Christoffel symbols of the first kind:

$$
\Gamma_{kij}
= \langle \Phi_{ij},\Phi_k\rangle.
$$

The next chapter is entirely about $\Gamma$: it is the tangent part of the
“second Jacobian” $\Phi_{ij}$.

> **One sentence summary.** $J$ gives the tangent plane; $I=J^{\mathsf T}J$
> measures inside it; $II$ measures how the tangent plane tilts; and $\Gamma$
> measures how the coordinate basis itself twists.

---

## 1.6 Worked examples

### Plane

$$
\Phi(u,v)=(u,v,0),\quad
I=\begin{bmatrix}1&0\\0&1\end{bmatrix},\quad
II=0.
$$

$K=H=0$. The plane neither stretches nor bends.

### Cylinder of radius $R$

$$
\Phi(\theta,z)=(R\cos\theta,R\sin\theta,z).
$$

$$
\Phi_\theta=(-R\sin\theta,R\cos\theta,0),\qquad
\Phi_z=(0,0,1),
$$

$$
I=\begin{bmatrix}R^2&0\\0&1\end{bmatrix},\qquad
n=(\cos\theta,\sin\theta,0).
$$

$$
\Phi_{\theta\theta}=(-R\cos\theta,-R\sin\theta,0)=-R\,n,
\qquad
\Phi_{zz}=0,
$$

so

$$
II=\begin{bmatrix}-R&0\\0&0\end{bmatrix}.
$$

One principal curvature is $-1/R$, the other is $0$; $K=0$ (a cylinder is
developable — you can unroll it without stretching). The nonzero curvature is
the bending of the circle, the zero one is the straight generator.

### Sphere of radius $R$

$$
\Phi(\theta,\varphi)
=R(\sin\theta\cos\varphi,\ \sin\theta\sin\varphi,\ \cos\theta).
$$

With the outward normal $n=\Phi/R$,

$$
I=
\begin{bmatrix}
R^2&0\\
0&R^2\sin^2\theta
\end{bmatrix},
\qquad
II=
\begin{bmatrix}
-R&0\\
0&-R\sin^2\theta
\end{bmatrix},
$$

so $S=-(1/R)\,\mathrm{Id}$ and

$$
\kappa_1=\kappa_2=-\frac1R,\qquad
K=\frac1{R^2},\qquad
H=-\frac1R.
$$

The sign is negative because the outward normal points *away* from the centre
while the surface bends toward it. Using the inward normal flips all signs
but leaves $K$ unchanged. The code prints $K=1$ for the unit sphere.

### Paraboloid graph $z=0.3(x^2+y^2)$

At $(x,y)=(1,0)$:

$$
\nabla f=(0.6,0),\qquad
\operatorname{Hess}f=
\begin{bmatrix}0.6&0\\0&0.6\end{bmatrix},
$$

$$
I=
\begin{bmatrix}1.36&0\\0&1\end{bmatrix},
\qquad
II=
\frac{1}{\sqrt{1.36}}
\begin{bmatrix}0.6&0\\0&0.6\end{bmatrix}
=
\begin{bmatrix}0.5145&0\\0&0.5145\end{bmatrix}.
$$

Notice $II$ is a multiple of the identity here, but $I$ is not. Therefore

$$
S=I^{-1}II
=
\begin{bmatrix}
0.5145/1.36&0\\
0&0.5145
\end{bmatrix}
=
\begin{bmatrix}
0.3783&0\\
0&0.5145
\end{bmatrix}.
$$

So the principal curvatures are $0.3783$ and $0.5145$, even though the
“bending matrix” $II$ looked isotropic. The **metric** decides which direction
is more curved. These are exactly the numbers printed by `lake exe geodesia`.
Also

$$
K=\frac{\det II}{\det I}
=\frac{0.5145^2}{1.36}
=0.1946,
\qquad
H=\frac{0.3783+0.5145}{2}=0.4464.
$$

---

## 1.7 How the code represents this (exact derivatives)

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
the main path. A paraboloid is built with exact combinators:

```lean
def paraboloid : Smooth2 :=
  Smooth2.smul 0.3
    (Smooth2.add (Smooth2.mul Smooth2.coordX Smooth2.coordX)
                 (Smooth2.mul Smooth2.coordY Smooth2.coordY))
```

`Smooth2.mul` implements the product rule in both variables, `Smooth2.sin`
implements the chain rule twice, etc. If you only have a bare
`Float → Float → Float`, use `Smooth2.ofFiniteDiff` as a fallback; that is the
only place an `ε` appears.

Then:

```lean
Surface.metricI s x y       -- I
Surface.metricII s x y      -- II
Surface.gaussK s x y        -- K = det II / det I
Surface.meanH s x y         -- H = ½ tr(I⁻¹ II)
Surface.principalK s x y    -- (κ_min, κ_max)
Surface.shapeOperator s x y w
Surface.normalCurvature s x y w
Surface.firstForm s x y w
Surface.secondForm s x y w
```

`Surface.dPhi` and `Surface.d2Phi` read the exact partials; `Curvature.lean`
never differentiates numerically.

---

## 1.8 Exercises

1. **Cone.** Take $\Phi(r,\theta)=(r\cos\theta,r\sin\theta,r)$ (a cone).
   Compute $I$ and $II$. Which principal curvature is zero? Why is $K=0$
   inevitable for a cone?
2. **Determinant identity.** Show
   $\det I=\det(J^{\mathsf T}J)=|\Phi_u\times\Phi_v|^2$.
3. **Graph curvature.** Starting from
   $I=\begin{bmatrix}1+f_x^2&f_xf_y\\f_xf_y&1+f_y^2\end{bmatrix}$ and
   $II=\operatorname{Hess}f/\sqrt{1+|\nabla f|^2}$, derive
   $$
   K=\frac{f_{xx}f_{yy}-f_{xy}^2}{(1+f_x^2+f_y^2)^2}.
   $$
4. **Weingarten sign.** Differentiate $n\cdot\Phi_i=0$ to prove
   $II_{ij}=-n_i\cdot\Phi_j$. What changes if you choose $-n$ instead of $n$?
5. **Cylinder vs sphere.** Explain in words why the cylinder has $K=0$ but the
   sphere has $K=1/R^2$, using the Gauss formula and the fact that the
   cylinder can be unrolled onto a plane.
6. **Code experiment.** In `Main.lean`, move the paraboloid point from
   $(1,0)$ to $(2,0)$ and predict how $\Gamma^u_{uu}$, $K$ and the principal
   curvatures change before running. (Hint: $K$ decreases as you move out.)
7. **Shape operator.** Verify numerically that
   $II(\xi,\eta)=I(S\xi,\eta)$ for a few vectors on the paraboloid.
