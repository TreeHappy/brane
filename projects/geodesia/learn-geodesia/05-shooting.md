# 5️⃣ Shooting methods — hitting a target with geodesics

> **Goal.** Understand the shooting method as *“integrate an initial-value
> problem, measure the miss, adjust the initial direction, repeat”*, see where
> the Jacobian comes from, and learn why the unknown geodesic length makes the
> naive version tricky (Eberly §2).

The shooting method is the numerical version of archery:

1. Stand at the source $a$.
2. Choose an initial direction and speed (the *shot*).
3. Integrate the geodesic ODE to see where you land.
4. Measure the miss relative to the target $b$.
5. Adjust the shot and repeat.

Everything else is bookkeeping: what “direction”, “miss” and “adjust” mean
precisely.

---

## 5.1 The boundary-value problem

The geodesic equation is a second-order ODE. An **initial-value problem**
(IVP) specifies the starting point and starting direction:

```math
q(0)=a,\qquad q'(0)=v.
```

A **boundary-value problem** (BVP) specifies two positions:

```math
q(0)=a,\qquad q(L)=b,
```

with the total length $L$ generally **unknown**. The BVP is harder: you are
not allowed to choose the starting direction; it is determined by the
requirement of
hitting $b$. Shooting turns the BVP into a root-finding problem for the
initial direction.

Write the solution of the IVP as

```math
q(s;a,v),
```

emphasising its dependence on the initial data. For a fixed $L$, the miss is

```math
\boxed{\;R(v)=q(L;a,v)-b.\;}
```

We want $R(v)=0$. This is a nonlinear equation in $v\in\mathbb{R}^2$ (for a
surface). The function $R$ is only available through numerical integration,
but that is enough for Newton’s method.

---

## 5.2 Standard shooting with a known length

Newton’s method is

```math
\boxed{\;v_{k+1}=v_k-\bigl[DR(v_k)\bigr]^{-1}R(v_k).\;}
```

The only hard part is the **Jacobian** $DR=\partial R/\partial v$. Formally,

```math
DR(v)=\frac{\partial q(L;a,v)}{\partial v}.
```

Let $Y(s)=\partial q(s;a,v)/\partial v$ be the $2\times2$ matrix of
derivatives with respect to the 2-component initial direction $v$. Differentiate the ODE
$q''=F(q,q')$ with respect to $v$:

```math
\boxed{\;
Y''=F_q(q,q')\,Y+F_{q'}(q,q')\,Y'.
\;}
\tag{first variational equation}
```

The initial conditions follow from $q(0)=a$ (independent of $v$) and
$q'(0)=v$:

```math
Y(0)=0,\qquad Y'(0)=I.
```

Then $DR(v)=Y(L)$. This is a **linear matrix ODE**; integrate it together
with the geodesic (or after it) using RK4, and you have an exact (up to RK4
error) Jacobian. This is the standard algorithm in Eberly §2.1.

### The geodesic variational equation explicitly

Recall what $F$ is doing, from Chapter 3's vector statement
$I\ddot q + g = 0$ with $g_k = \Gamma_{kij}\,d^id^j$ the coordinate recipe of
the tangent parts: $F$ is that equation solved for $\ddot q$ (then
multiplied through by $I^{-1}$). In coordinates:

```math
F^i(q,q')=-\Gamma^i_{\;jk}(q)\,q'^j q'^k,
```

we have

```math
\frac{\partial F^i}{\partial q'^j}
=-2\Gamma^i_{\;jk}q'^k,
\qquad
\frac{\partial F^i}{\partial q^m}
=-\partial_m\Gamma^i_{\;jk}\,q'^j q'^k.
```

Therefore each column $Y_{\cdot j}$ of $Y$ satisfies

```math
\boxed{\;
(Y^i)''+2\Gamma^i_{\;jk}q'^k\,(Y^j)'
+\partial_m\Gamma^i_{\;jk}q'^j q'^k\,Y^m=0.
\;}
```

Geometrically this is the *linearization of the tangential-projection
equation*: perturb the curve slightly and ask how the residual $T\gamma''$
then reads in coordinates — the Jacobi equation is what the arrow-and-
projection condition says about a pencil of nearby geodesics. Its solutions
are **Jacobi fields**: they measure how a nearby geodesic separates from the
original one. The columns of $DR$ are precisely Jacobi fields with initial
data $Y(0)=0$, $Y'(0)=I$.

---

## 5.3 The unknown length problem

The standard algorithm assumes $L$ is known. For a geodesic between two
points, $L$ is exactly what we are trying to find. Eberly §2.2 explains why
you cannot simply guess it:

- If $L$ is **too small**, no shot can reach $b$: every curve of length $L$
  falls short.
- If $L$ is **too large**, the geodesic may pass through $b$ before the end;
  forcing the endpoint at $s=L$ to equal $b$ is then impossible.

### Fix 1: overestimate and minimise the miss

Choose an overestimate $L'>L$ (for example, the length of the straight chord
from $a$ to $b$). For each initial direction $v$, integrate to
$c(v)=q(L';a,v)$. The true geodesic passes through $b$ at some $s=L<L'$, so
$b$ lies on the curve even though the endpoint is not $b$. Minimise the
distance from $b$ to the endpoint:

```math
\ell(v)=\operatorname{dist}\bigl(b,c(v)\bigr).
```

Eberly uses a line-segment distance weighted by the metric along the segment
(his eq. 18):

```math
\ell(v)=\int_0^1
\sqrt{(c-b)^{\mathsf T}G\bigl((1-t)b+tc\bigr)(c-b)}\,dt.
```

The gradient $\ell_v$ involves $c_v=\partial c/\partial v$, i.e. the same
variational matrix $Y$ as before. Use steepest descent or conjugate gradients
to minimise $\ell$. This is robust but requires a good overestimate and a
minimiser.

### Fix 2: solve for angle and length together

On a 2-D surface there is a cleaner formulation. Parameterise the initial
initial data by a **unit** surface direction and an unknown length:

```math
v=\cos\alpha\,e_1+\sin\alpha\,e_2,
\qquad (e_1,e_2)\text{ orthonormal in }T_aS.
```

The unknowns are $(\alpha,L)$, and the equations are

```math
q(L;a,\alpha)-b=0
\qquad (2\text{ equations}).
```

So $2$ unknowns for $2$ equations. This is the version used in the worked
example below and in the figure script. It is the $n=2$ case of the general
statement “unit direction on $S^{n-1}$ ($n-1$ angles) plus length $L$ gives
$n$ unknowns, matching $n$ boundary equations.”

Initial guesses:

- $\alpha_0$: direction of the chord $b-a$, projected into the tangent plane;
- $L_0$: the chord length (or $1.1$–$1.3$ times it).

Then run damped Newton on $(\alpha,L)$.

### Fix 3: reparameterise to $t\in[0,1]$

Eberly §2.3 writes $q(s)=y(t)$ with $t\in[0,1]$ and derives an ODE for
$y(t)$ involving $s'=ds/dt$ and $s''=d^2s/dt^2$. The boundary conditions
become $y(0)=a$, $y(1)=b$, so the standard shooting method can be applied
with a known parameter interval. The price is a more complicated ODE (eq. 24
in the PDF) because the unknown length is hidden in $s'$.

---

## 5.4 Which unknowns, which equations?

| Situation | Unknowns | Equations |
|---|---|---|
| target is a point, $L$ known | starting direction $v\in\mathbb{R}^2$ | $q(L)-b=0$ (2) |
| target is a point, unit speed | angle $\alpha$, length $L$ | $q(L)-b=0$ (2) |
| general $n$-D manifold, unit speed | $n-1$ angles + $L$ | $n$ endpoint eqs |
| target is a curve (e.g. boundary) | angle $\alpha$ (or $v$ modulo tangent to curve) | one scalar condition |

For a point target on a surface, the ambient residual $q(L)-b\in\mathbb{R}^3$
has only two independent components because the endpoint is constrained to
the surface. Working in parameter coordinates avoids the redundancy.

---

## 5.5 The full Newton algorithm (angle + length)

```
shoot(alpha, L):
    e1, e2 = orthonormal tangent frame at a
    w = cos(alpha) e1 + sin(alpha) e2
    xi = I(a)^{-1} J(a)^T w        # tangent-basis coords of w
    q = integrate_geodesic(a, xi, L)
    return q

residual(x = (alpha, L)):
    return shoot(alpha, L) - b

Newton:
    x = (alpha0, L0)
    repeat:
        r = residual(x)
        if ||r|| < tol: break
        J = finite-difference Jacobian of residual
        solve J * dx = -r
        choose step t (backtracking line search)
        x = x + t * dx
```

The finite-difference Jacobian costs two extra integrations per iteration
(one per unknown). For larger systems, use the variational equation instead:

```
augmented state = (q, p, Y, Y')      # Y is n x n
integrate the geodesic and Y'' = F_q Y + F_p Y' together
DR = Y(L)
```

This costs one integration of a larger system, and is more accurate.

### Damping and safeguards

- **Backtracking line search:** accept a step only if it reduces $\|R\|$.
- **Bounds:** keep $\alpha\in[0,2\pi)$; require $L>0$.
- **Multiple seeds:** geodesics are not unique. Try several angles, especially
  for targets near the cut locus.
- **Continuation:** solve a sequence of targets moving gradually from an easy
  one to the desired one.
- **Rescale:** work in units where the surface’s typical length is $O(1)$.

---

## 5.6 Why shooting fails: conjugate points and caustics

The Jacobian $DR=Y(L)$ can be singular. A nonzero Jacobi field $Y$ with
$Y(0)=0$ and $Y(L)=0$ is a **conjugate point** of $a$ at parameter $L$.
Geometrically, nearby geodesics launched from $a$ refocus at $b$; the map
“initial direction $\mapsto$ endpoint” is not invertible.

Examples:

- **Sphere, antipodal points.** Every great circle from the north pole reaches
  the south pole. The target is reachable by a continuum of initial
  directions, and $DR$ is singular. Newton has no unique answer.
- **Ellipsoid.** There can be several geodesics between two points; shooting
  converges to the one selected by the initial guess.
- **Near the cut locus.** Even when $DR$ is invertible, it is ill-conditioned;
  small changes in the shot produce large changes in the endpoint. Use
  damped Newton, better initial guesses, or switch to relaxation.

This is not a bug in the implementation; it is geometry. Relaxation methods
(Chapter 6) are often more robust precisely because they do not rely on
inverting the endpoint map.

---

## 5.7 Worked example: paraboloid

Take

```math
\Phi(u,v)=\bigl(u,v,0.3(u^2+v^2)\bigr),
\qquad
a=(0.4,0.2),\quad b=(1.2,0.9).
```

The straight chord has length $1.2281$. Start with the chord direction and
$L_0=1.15\times\text{chord}=1.4123$. Damped Newton on $(\alpha,L)$ gives:

| iteration | $\alpha$ | $L$ | $\|R\|$ |
|---:|---:|---:|---:|
| 0 | 0.67857208 | 1.41231284 | $1.26\times10^{-1}$ |
| 1 | 0.67976725 | 1.23649445 | $2.44\times10^{-3}$ |
| 2 | 0.67974028 | 1.23976961 | $8.87\times10^{-7}$ |
| 3 | 0.67974028 | 1.23977080 | $2.02\times10^{-14}$ |

The geodesic length is $1.23977$, slightly longer than the chord, as it must
be on a curved surface. The final unit tangent direction is about
$\alpha=0.6797$ rad ($38.9^\circ$).

![Shooting on the paraboloid: grey dashed shots miss the target at the initial length; the red path is the converged geodesic. The right panel shows the residual collapsing under Newton iteration.](images/shooting.png)

The Lean function `Geodesic.shoot` computes one shot:

```lean
Geodesic.shoot s u v angle length steps   -- ambient endpoint Vec3
```

It is the IVP half of the method. A root finder around it (in Lean or in a
host script) is the other half.

---

## 5.8 A fan of shots: the Jacobi metric

The metric need not be inherited from an embedding. A large family comes from
mechanics. For a unit-mass particle at energy $E$ moving in a potential
$V(q)$, the spatial trajectories of

```math
q''=-\nabla V(q)
```

are the unparameterised geodesics of the **Jacobi metric**

```math
\boxed{\,g_{ij}=2\bigl(E-V(q)\bigr)\,\delta_{ij}\,}.
```

The parameter is different, but the curves are the same; energy conservation
fixes the speed along them:

```math
\tfrac12|q'|^2+V(q)=E
\quad\Longrightarrow\quad
|q'(0)|=\sqrt{2\bigl(E-V(q_0)\bigr)}.
```

Take a Gaussian hill

```math
V(x,y)=A\,e^{-(x^2+y^2)/2\sigma^2},
```

with height $A$, width $\sigma$, and energy $E>V_{\max}=A$. Then $E-V>0$
everywhere and every shot is defined.

Shooting is now a fan. Fix a source $q_0$ and $N$ launch angles

```math
\alpha_k=\frac{2\pi k}{N},\qquad k=0,1,\dots,N-1.
```

For each $k$, integrate the IVP

```math
q_k(0)=q_0,\qquad
q_k'(0)=\sqrt{2\bigl(E-V(q_0)\bigr)}\,
(\cos\alpha_k,\sin\alpha_k),
```

and draw the spatial curve $q_k(t)$ until it leaves the plotting square. The
initial speed is not a free parameter: it is determined by the energy.

![Shooting a fan of geodesics of the Jacobi metric $g=2(E-V)\delta$ near a Gaussian hill $V=A e^{-(x^2+y^2)/2\sigma^2}$ with $A=1$, $\sigma=0.70$, $E=1.10>V_{\max}$. Forty shots leave $q_0=(1.50,1.50)$ at equal angles, with speed fixed by the energy. The hill deflects most shots; the single shot aimed at the summit (mauve) is the only one that reaches the lower-left shadow behind the hill.](images/jacobi-fan.png)

The hill bends every shot that passes near it, and because $q_0$ is off-centre
the fan is not symmetric. Most shots are deflected around the hill. The region
behind the hill (lower-left in the figure) receives only the shot aimed at the
summit (mauve). That empty region is a **shadow**, not a caustic: the forward
map $\alpha\mapsto q_k(\text{exit})$ is still monotone here, but the hill has
swept a whole interval of launch angles away from it. A caustic appears only
when that map folds, which is the failure mode of §5.6.

---

## 5.9 🥾 Walking a hill: the tangent-plane step, projection, and why a straight path bends

> **Picture.** §5.8 shot geodesics of the Jacobi metric in the *parameter
> plane*. This section stays in the room instead: the hill becomes an
> honest embedded surface, and we find a geodesic the way a hiker would —
> take a straight step, fall back onto the ground, re-aim along the new
> ground, repeat. No geodesic ODE is written down anywhere in the walk —
> and the walk lands on the same curve the ODE produces.

The hill is now carried by the embedding, not by a metric:

```math
\Phi(x,y) = \bigl(x,\,y,\,h(x,y)\bigr),
\qquad h(x,y) = A\,e^{-(x^2+y^2)/2\sigma^2},
\qquad A = 1,\ \sigma = 0.7,
```

the same Gaussian as §5.8, lifted into $\mathbb{R}^3$.

(Watch the notation: the hill's height is $h(x,y)$ *and* the step length
below is $h$ too. Subscripted symbols like $h_{xx}$, $h_{xy}$, $h_{yy}$
always belong to the height; the bare $h$ is a step length.)

> **What does "ambient" mean?** The *ambient space* is just the surrounding
> $\mathbb{R}^3$ the surface lives in. Positions (`Surface.point : Vec3`),
> directions, and dot products all live there. The surface's ruler is
> ambient too: $I = J^{\mathsf T}J$ is the plain ambient dot product pulled
> back along $J$, so a surface length is an ambient chord length. And
> `Geodesic.shoot` returns an *ambient* endpoint `Vec3` — which is why the
> shooting residual $q(L) - b$ of §5.1 is an ambient vector. But $b$ and
> $q(L)$ both sit on the surface, so only two of the residual's three
> components are independent; that redundancy is exactly why §5.4 prefers
> parameter coordinates.

### The walker

Stand at a point $p$ on the surface holding a *unit tangent* direction $w$
(a unit ambient tangent: $I(w,w) = 1$, so straight ambient steps of size
$h$ are honest length-$h$ surface steps, up to curvature). The loop has
three moves:

```
walk(p, w, h, L):                    # p on the surface, unit tangent w at p
    repeat round(L / h) times:
        tip = p + h*w                # 1. move straight: a chord in ambient R^3
        q   = project(tip)           # 2. pull the chord tip back onto the surface
        n   = unit normal at q
        w   = (w - (w·n) n) / |w - (w·n) n|   # 3. re-tangent + renormalize
        p   = q
    return p
```

Step 1 is the hiker's stride: a perfectly straight segment in $\mathbb{R}^3$
that — this is the whole point — *leaves the surface*, because the surface
curves away from its own tangent plane. Step 2 is the fall back to the
ground; for a graph the cheapest exact projection is the vertical one: keep
the tip's $(x, y)$, replace its $z$ by the surface height there. (A
closest-point projection works too; the numbers below use the vertical
one.) Step 3 is the re-aim, and it is a Gauss split you already own: the
new tangent plane is tilted relative to the old one, so the old $w$ picks
up a normal component; kill it — that is exactly
`Surface.tangentialPart` — and renormalize to a unit vector again.

> **Step 3 is parallel transport.** In the $h \to 0$ limit the
> re-tangent step is exactly the parallel transport of a tangent vector
> along the walked curve. Carrying any tangent vector $W$ along a curve
> with unit tangent $t$, parallel transport solves
> $\tfrac{dW}{ds} = -\langle W,\, Dn\,t\rangle\, n$ — the normal tips as
> you move, and the vector must shed exactly that much normal residue to
> stay tangent. The walker's update
> $\Delta W = -h\,\langle W,\, Dn\,t\rangle\, n + O(h^2)$ is that ODE's
> Euler step, and renormalization only matters at $O(h^2)$. And "$w$
> parallel along itself" — transport a tangent vector along the very
> curve it spans — *is* the geodesic condition. That is why the walk
> converges to a geodesic: each step parallel-transports, and
> self-parallel curves are geodesics.

### The walk, in numbers

Launch data, all fixed:

- start $a = \Phi(1.40,\,0.90) = (1.40,\,0.90,\,0.059218)$;
- launch angle $\alpha = 215^\circ$ in the orthonormal tangent frame
  $(e_1, e_2)$ at $a$ (the `Geodesic.tangentFrame` construction:
  $e_1$ along $\Phi_x$, $e_2$ equal to $\Phi_y$ with its $e_1$-part
  removed), so $w = \cos\alpha\,e_1 + \sin\alpha\,e_2$;
- the parameter-space velocity the ODE route needs comes from the chapter's
  conversion $\xi = I^{-1}J^{\mathsf T}w = (-0.797470,\,-0.570306)$ — the
  same line as §5.5's pseudocode;
- total walk length $L = 1.6$.

The reference curve is the geodesic from the same launch data, integrated
with RK4 on the identity $I\,(u'', v'') = -g$ — the one `Geodesic.accelVector`
implements and [`Proofs.lean`](Proofs.lean) verifies. Its endpoint is

```math
\Phi(0.315158,\,0.107340) = (0.315158,\,0.107340,\,0.893054),
```

so the walk has to climb from $z = 0.059$ up the shoulder of the hill to
$z \approx 0.893$, just below the summit plateau. Here is the walker's
convergence — the ambient distance between its endpoint and the ODE
endpoint:

| $h$ | steps | ambient endpoint error |
|---:|---:|---:|
| 0.2 | 8 | 0.014587 |
| 0.1 | 16 | 0.006906 |
| 0.05 | 32 | 0.003350 |
| 0.025 | 64 | 0.001649 |
| 0.0125 | 128 | 0.000818 |

Every halving of $h$ halves the error: **first-order convergence**,
error $\propto h$. The first three steps at $h = 0.1$ show the two
per-step measurements the next section explains:

| step | chord-tip miss (gap to the surface) | turn of $w$ |
|---:|---:|---:|
| 1 | 0.0029 | $3.272^\circ$ |
| 2 | 0.0033 | $3.650^\circ$ |
| 3 | 0.0036 | $3.856^\circ$ |

At the finer $h = 0.025$ the largest turn of the whole walk is
$1.262^\circ$ per step, and the total bend accumulated over the 1.6-length
walk is $39.78^\circ$.

![The hill walk on the embedded Gaussian surface $\Phi(x,y)=(x,y,e^{-(x^2+y^2)/2\sigma^2})$: at the start point the tangent plane $T_pS$, the unit direction $w$ and the normal $n$; the first chord (dashed) is perfectly straight in ℝ³ but floats a gap of 0.0029 above the curving surface — the side-view inset shows chord vs surface to scale; the piecewise-linear walker (peach, 16 chords of $h=0.1$) and the RK4 geodesic (red) of the same surface are visually one curve, climbing together from $z=0.059$ to the shared endpoint at $z=0.893$.](images/hill-walk.png)

### Why the straight path bends

Each chord is perfectly straight in $\mathbb{R}^3$, yet it floats *off* the
surface. The floating is not an accident of the scheme: because the chord
direction $w$ is tangent, the surface's height above the tangent plane has
**no linear term** along the chord — first order cancels — and what
survives is the surface's own second-order departure,

```math
\operatorname{gap}(h) \;=\; \frac{h^2}{2}\,
\bigl(h_{xx}w_x^2 + 2\,h_{xy}\,w_x w_y + h_{yy}w_y^2\bigr) \;+\; O(h^3).
```

At the launch point the quadratic factor is $0.269$, so the formula
predicts $0.0027$ where the walk measured $0.0029$ — the misses in the
table above are $h^2$-sized, exactly as the formula says.

The projection snaps the point back down, but it cannot undo the *tilt*:
the normal $n$ at the new point has turned, so the old $w$ now has a
normal component, and the re-tangenting step rotates $w$ by a purely
normal correction. Per unit length, that rotation rate is exactly

```math
-\langle w,\,Dn\,w\rangle \;=\; \kappa_n(w),
```

the normal curvature of the hill in the walking direction — the shape
operator of Chapter 2 at work. The numbers agree: at the launch point
$\kappa_n(w) = 0.5283$ rad per unit length ($30.25^\circ$/unit). The
$h = 0.1$ walk turns its first step $3.272^\circ$, a rate of
$32.7^\circ$/unit; the $h = 0.025$ walk turns its first step $0.773^\circ$,
a rate of $30.9^\circ$/unit — closing in on $\kappa_n$ as $h$ shrinks. The
turn also *grows* along the walk ($3.272^\circ \to 3.650^\circ \to
3.856^\circ$ in three steps) as the shoulder steepens under the walker's
feet.

Compound hundreds of these tiny straight moves on successively tilted
tangent planes and the zig-zag melts into a smooth curve whose ambient
acceleration is **purely normal**:

```math
\boxed{\;\frac{dw}{ds} \;=\; -\langle w,\,Dn\,w\rangle\,n
\;=\; II(w,w)\,n,
\qquad T\,\frac{dw}{ds} \;=\; 0.\;}
```

The bend you *see* in 3D is the surface bending under your feet (normal
curvature), not the path choosing to steer sideways: the tangential part
of the bending vanishes. That is precisely the geodesic condition —
"straightest inside the surface" — and its formal version is the
machine-checked statement `GeodesicProof.tangential_gamma2_eq_zero` in
[`Proofs.lean`](Proofs.lean): with the geodesic acceleration,
$\langle \Phi_k, \Phi''\rangle = 0$ for $k = u, v$, i.e. the ambient
acceleration is purely normal.

In the limit $h \to 0$ the walker converges to the true geodesic, and the
convergence table says the accuracy is first-order. One sentence of why:
each step lands within $O(h^3)$ of the true geodesic arc *along the
normal*, but the vertical projection discards the $O(h^2)$ *tangential*
drift that the geodesic's own normal bending produces in the $(x,y)$
plane; summed over $L/h$ steps, that per-step $O(h^2)$ becomes a global
$O(h)$.

> **Why the Lean code doesn't walk.** The walk is the *understanding*
> device: three plain moves, no Christoffel symbols, no ODE — and
> first-order accurate. The library integrates the ODE instead:
> `Geodesic.integrate` runs RK4, which is fourth-order accurate, so
> refining the grid buys entire digits at once rather than one at a time.
> The two routes meet at the same curve — the convergence table above is
> that meeting, watched from the walker's side, and the figure shows it:
> peach and red are the same path.

### The same walk, read in parameter space

Everything the walker did above has a second life in $(u, v)$
coordinates — and read there, the walk stops being a trick and becomes
the geodesic ODE in disguise.

**The parameter velocity is the tangent, seen from above.** For the
graph $\Phi(x,y) = (x, y, h(x,y))$, the parameter velocity
$\xi = I^{-1}J^{\mathsf T}w$ of the launch data (the
`Geodesic.tangentVecToParams` conversion) is literally
$\xi = (w_x, w_y)$: tangency forces $w_z = h_x w_x + h_y w_y$, which is
exactly the statement $J\xi = w$ with $\xi = (w_1, w_2)$. And the
vertical projection of step 2 does not touch $(x, y)$ at all. So read in
parameter coordinates, the walker takes perfectly straight steps

```math
q_{k+1} \;=\; q_k + h\,\xi_k,
```

and the *only* thing the re-tangent step does is update $\xi$. This
graph identity is the kind of statement [`Proofs.lean`](Proofs.lean)
machine-checks <!-- LEMMA-REF: graph identity -->.

**The re-tangent step IS the Christoffel term.** Compare the re-tangent
update on $\xi$ against one Euler step of Chapter 3's geodesic ODE
$I\,\ddot q = -g$, i.e. $\xi \mapsto \xi + h\cdot\mathrm{accel}$ with
$\mathrm{accel} = -I^{-1}g$ — or, in §5.2's Christoffel recipe
$g_k = \Gamma_{kij}\,\xi^i\xi^j$, the update
$\xi \mapsto \xi - h\,\Gamma(\xi,\xi)$. From the chapter's launch data
$\xi_0 = (-0.797470, -0.570306)$:

| $h$ | $\lVert\,\text{re-tangent} - \text{Euler}\,\rVert$ |
|---:|---:|
| 0.1 | $2.45\times 10^{-3}$ |
| 0.05 | $5.88\times 10^{-4}$ |
| 0.025 | $1.44\times 10^{-4}$ |

Each halving of $h$ quarters the difference: the two updates agree up to
$O(h^2)$. The ambient normal $n$ and the Christoffel symbols $\Gamma$
are two coordinate systems for the same bending.

**The pure parameter-space walker.** That reading collapses the loop to
two lines:

```
walk_params(q, xi, h, L):
    repeat round(L / h) times:
        xi = xi + h * accel(q, xi)     # the "re-tangent" step: ξ -= h Γ(ξ,ξ)
        q  = q  + h * xi               # the "straight chord" step
    return q
```

No projection, no tip, no normal — the chart enforces the constraint for
free: the coordinates ARE the surface. Endpoint errors against the RK4
reference endpoint $(0.315158,\ 0.107340)$, same launch data,
$L = 1.6$, measured in parameter space (the ambient-distance version of
this comparison is the table in "The walk, in numbers"):

| $h$ | ambient walker (§5.9's loop) | parameter Euler walker |
|---:|---:|---:|
| 0.2 | 0.013029 | 0.027849 |
| 0.1 | 0.006169 | 0.014040 |
| 0.05 | 0.002994 | 0.007044 |

Both columns halve with $h$: first order, and the same curve — but the
ambient walker has the better constant, because the straight chord in
$\mathbb{R}^3$ does part of the curvature correction for free.

> **Straight in parameter space is not a thing.** Walk straight in
> parameter space with *no* correction — constant $\xi$,
> $q_{k+1} = q_k + h\,\xi_k$ forever — and you do *not* walk a geodesic:
> from the same launch data you end at $(0.124048,\ -0.012490)$, a miss
> of $0.2256$ — two orders of magnitude worse than either walker above.
> In ambient $\mathbb{R}^3$ "straight" is unambiguous (a chord); in
> parameter space it is meaningless until you say which ruler — and the
> Christoffel term is exactly what "straight in $\mathbb{R}^3$" looks
> like written in $(u, v)$ coordinates. This is Chapter 3's
> $I\ddot q = -g$ seen from the walker's side: the velocity bends so
> that the ambient path stays straight.

**Three viewpoints, one walk:**

| Viewpoint | How to step | What measures the bending |
|---|---|---|
| ambient | kill the normal part of $w$ (`Surface.tangentialPart`) | turn rate $\kappa_n(w)$ |
| parameter | $\xi \mapsto \xi - h\,\Gamma(\xi,\xi)$ | Christoffel symbols of the graph metric |
| invariant | tangential part of the acceleration $= 0$ | the geodesic condition — `GeodesicProof.tangential_gamma2_eq_zero` in [`Proofs.lean`](Proofs.lean) |

### Where it lives in the code

Every ingredient of the walk is already in the library — only the loop
itself lives in the figure script, because the Lean side deliberately
takes the accurate ODE route (that is the point of the callout above):

| What | Code | File |
|---|---|---|
| the Gauss split the re-tangent step runs on: $T w = w - (w\cdot n)\,n$ | `Surface.tangentialPart`, `Surface.normalPart`, `Surface.normalAt` | `Geodesia/Surface.lean` |
| ambient direction → parameter velocity $\xi = I^{-1}J^{\mathsf T}w$ | `Geodesic.tangentVecToParams` | `Geodesia/Geodesic.lean` |
| the ODE shot — the accurate route to the same curve | `Geodesic.shoot`, `Geodesic.integrate` | `Geodesia/Geodesic.lean` |
| the walker of this section (figure code, not Lean) | `hill_walk` / `fig_hill_walk` | `images/make_figs.py` |
| the parameter-space walker of the subsection above, plus the printed comparison (figure code, not Lean) | `hill_walk_param_space` / `hill_walk_param_report` | `images/make_figs.py` |
| the walk as a game: arrow keys steer the §5.9 walker (and its straight-in-$\mathbb{R}^2$ twin) on the hill, no dependencies — open the file in a browser | the playable walker | [`WALKER.html`](WALKER.html) |
| the formal backbone: no tangential residue | `GeodesicProof.tangential_gamma2_eq_zero` | `Proofs.lean` |

> **Try it.** Rerun the walk with $\alpha$ aimed straight at the summit,
> then tangentially across the shoulder, and watch each walk's total bend
> track $\kappa_n(w)$: the summit-directed launch spends its 1.6 length
> where the normal curvature is largest, the tangential launch mostly does
> not. Second experiment: change $h$ and watch the total bend settle
> ($39.98^\circ$ at $h = 0.2$, $39.78^\circ$ at $h = 0.025$) even as the
> endpoint error collapses — the bend is geometry; the error is numerics.
> Or skip the scripts and *play* it: open [`WALKER.html`](WALKER.html) in a
> browser, walk with the arrow keys, and press G to switch between the
> geodesic walker and its straight-in-$\mathbb{R}^2$ twin — the HUD reads
> out $\kappa_n(w)$ and the chord gap live.

---

## 5.10 Exercises

1. **Variational equation.** Derive
   $Y''=F_qY+F_{q'}Y'$ from $q''=F(q,q')$ by differentiating with respect to
   the initial direction. What are $Y(0)$ and $Y'(0)$?
2. **Jacobi equation.** Substitute the coordinate recipe
   $F^i=-\Gamma^i_{\;jk}q'^jq'^k$ (Chapter 2 §2.8 / Chapter 3: the
   coordinates of $T\gamma''$'s acceleration) and derive the explicit
   second-order equation for $Y$ in §5.2. Then sanity-check it geometrically:
   the linearization must say "vanishing tangent projection stays vanishing".
3. **1-D shooting.** Fix $L$ and target a curve (e.g. the circle
   $u^2+v^2=1$). Reduce the residual to one scalar equation and implement
   bisection or secant. Why is no 2-D Jacobian needed?
4. **Underestimate.** Explain geometrically why an underestimated $L$ can
   never hit a target farther than $L$ away. What happens numerically?
5. **Overestimate objective.** Write the gradient of Eberly’s line-segment
   distance $\ell(v)$ in terms of $c_v=Y(L')$. (You may differentiate under
   the integral.)
6. **Antipodal sphere.** Show that for $a=$ north pole and
   $b=$ south pole, every unit initial direction reaches $b$ at $s=\pi R$.
   What does $DR$ look like?
7. **Nearby geodesics.** Integrate two geodesics from $a$ with initial angles
   differing by $10^{-3}$ and plot their separation vs $s$. Where does it
   grow fastest?
8. **Code.** Use `Geodesic.shoot` to tabulate the residual for
   $\alpha\in[0.6,0.8]$ at the converged $L$, and implement a secant method
   to find the zero.
