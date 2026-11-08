# 4️⃣ Ordinary differential equations and Runge–Kutta

> **Goal.** Read the geodesic equation as a first-order ODE, understand what
> an integrator actually computes, and know why RK4 is the default in
> `Geodesic.lean`.

---

## 4.1 From one ODE to a state vector

A first-order ODE is

```math
y'(s)=f\bigl(s,y(s)\bigr),
\qquad y(s_0)=y_0.
```

Here $y$ can be a vector. The **state** is everything you need to know to
continue the solution; $f$ is the **vector field** that says how the state
changes.

The geodesic equation is second order:

```math
q''=F(q,q'),
\qquad F^k(q,q')=-\Gamma^k_{\;ij}(q)\,q'^i q'^j.
```

Introduce velocities $p=q'$:

```math
\boxed{
\begin{aligned}
q' &= p,\\
p' &= F(q,p).
\end{aligned}}
```

In components, with $q=(u,v)$ and $p=(\dot u,\dot v)$:

```math
u'=\dot u,\quad
v'=\dot v,\quad
\dot u'=F^u(u,v,\dot u,\dot v),\quad
\dot v'=F^v(u,v,\dot u,\dot v).
```

This is the `GeoState` in the code:

```lean
structure GeoState where
  u  : Float
  v  : Float
  du : Float
  dv : Float
```

and `Geodesic.rhs` returns `(du, dv, u'', v'')`.

> **Order reduction** is the standard trick: an $n$-th order ODE becomes a
> first-order system of dimension $n\times(\text{number of components})$.
> A geodesic in two parameters has a 4-dimensional state.

---

## 4.2 Euler’s method

The simplest integrator is

```math
y_{n+1}=y_n+h\,f(s_n,y_n),
\qquad s_{n+1}=s_n+h.
```

It walks along the tangent of the vector field for a time $h$.

- Local truncation error: $O(h^2)$.
- Global error after a fixed interval: $O(h)$ (first order).
- It is explicit and cheap, but it drifts badly and is not accurate enough
  for geodesics except as a teaching tool.

---

## 4.3 Midpoint and fourth-order Runge–Kutta

The midpoint method samples the vector field halfway along the Euler step:

```math
k_1=f(s_n,y_n),
```

```math
k_2=f\bigl(s_n+\tfrac h2,\ y_n+\tfrac h2 k_1\bigr),
```

```math
y_{n+1}=y_n+h\,k_2.
```

It is second order: global error $O(h^2)$.

**Classical RK4** uses four samples:

```math
\begin{aligned}
k_1&=f(s_n,y_n),\\
k_2&=f\bigl(s_n+\tfrac h2,\ y_n+\tfrac h2 k_1\bigr),\\
k_3&=f\bigl(s_n+\tfrac h2,\ y_n+\tfrac h2 k_2\bigr),\\
k_4&=f\bigl(s_n+h,\ y_n+h k_3\bigr),
\end{aligned}
```

```math
\boxed{\;
y_{n+1}=y_n+\frac h6\bigl(k_1+2k_2+2k_3+k_4\bigr).
\;}
```

It is fourth order: global error $O(h^4)$. The weights
$(1,2,2,1)/6$ are Simpson’s rule in disguise: RK4 matches the Taylor
expansion of the exact solution through $h^4$.

### A quick accuracy check

For $y'=y$, $y(0)=1$, integrating to $s=2$ (exact $e^2=7.389056$):

| step $h$ | Euler | midpoint | RK4 |
|---:|---:|---:|---:|
| $0.5$ | $5.0625$ | $6.9729$ | $7.38397$ |
| $0.25$ | $5.9605$ | $7.2622$ | $7.38867$ |

Halving $h$ reduces the Euler error by about a factor $2$, the midpoint
error by about $4$, and the RK4 error by about $16$ — exactly the predicted
orders $1,2,4$.

---

## 4.4 Stability and step size

Apply a method to the test equation

```math
y'=\lambda y,
\qquad \lambda\in\mathbb{C}.
```

The numerical solution is $y_n=R(z)^n y_0$ with $z=\lambda h$. For RK4,

```math
R(z)=1+z+\frac{z^2}{2}+\frac{z^3}{6}+\frac{z^4}{24}.
```

The method is stable when $|R(z)|\le 1$. For real negative $\lambda$, RK4 is
stable for roughly $-2.8\le \lambda h\le 0$. This is a **conditional**
stability limit: too large a step makes the numerical solution explode even
though the exact solution decays.

For geodesics:

- The ODE is usually not stiff, but near high curvature, a nearly singular
  metric, or a coordinate singularity (e.g. the poles of spherical
  coordinates), the effective $\lambda$ can be large.
- The step must resolve the *smallest* length scale of the surface, not just
  the total path length.
- Check convergence by halving the step; a well-resolved RK4 solution changes
  by a factor of about $16$ in the error.
- For very long integrations, monitor conserved quantities (energy, angular
  momentum). A drift is a warning sign.

### Adaptive steps and symplectic methods

Production codes often use an embedded RK45 pair: two estimates of different
order give a local error estimate, and the step is accepted/rejected
accordingly. For Hamiltonian systems, **symplectic** integrators (Störmer–
Verlet, implicit midpoint, Gauss–Legendre) preserve the phase-space
structure and show much better long-term energy behaviour.

---

## 4.5 The geodesic ODE in code

`Geodesia/Geodesic.lean` contains the whole pipeline:

```lean
def accel (s : Surface) (st : GeoState) : Vec2 :=
  -- u'' = -Γ^u_ij du^i du^j,  v'' = -Γ^v_ij du^i du^j

def rhs (s : Surface) (st : GeoState) : GeoState :=
  -- (du, dv, u'', v'')

def rk4Step (s : Surface) (h : Float) (st : GeoState) : GeoState :=
  -- the four k's above

def integrate (s : Surface) (h : Float) (n : Nat) (st : GeoState) : GeoState :=
  -- n steps of size h
```

Because `Surface.christoffel` uses exact derivatives, the only numerical
error is the RK4 discretisation. The demo reports

```
E(0), E(1) = 0.680000, 0.680000
```

for the paraboloid, i.e. energy is conserved to the printed precision.

### Choosing initial conditions

To integrate a geodesic you need a point and a direction. A convenient way to
specify the direction is a unit vector in the tangent plane:

```math
w=\cos\alpha\,e_1+\sin\alpha\,e_2,
```

where $(e_1,e_2)$ is an orthonormal tangent frame (`Geodesic.tangentFrame`).
Then convert to parameter velocity:

```math
\xi=I^{-1}J^{\mathsf T}w
\qquad\text{(`Geodesic.tangentVecToParams`)}.
```

This is exactly the conversion in Eberly’s *Geodesics on a Parameterized
Surface*: $J^{\mathsf T}J=I$, so $I^{-1}J^{\mathsf T}J\xi=\xi$.

`Geodesic.shoot` does this for a given angle and length. Chapter 5 turns it
into a boundary-value solver.

---

## 4.6 Preview: the variational equation

For shooting we need to know how the final point changes when the initial
velocity changes. If $y(s;v)$ solves

```math
y'=f(s,y),
\qquad y(0)=y_0,
```

then the matrix

```math
Y(s)=\frac{\partial y(s;v)}{\partial v}
```

solves the **variational equation**

```math
Y'=f_y(s,y(s;v))\,Y,
\qquad Y(0)=I
```

(where $I$ is the identity on the velocity components, and zeros on the
position components). This is just the chain rule applied to the ODE. It is a
linear matrix ODE that can be integrated alongside the original one. Chapter 5
shows how it produces the Jacobian for Newton’s method.

---

## 4.7 Exercises

1. **Order reduction.** Write the third-order ODE $y'''=f(y,y',y'')$ as a
   first-order system. How many state components?
2. **Euler error.** For $y'=y$, $y(0)=1$, show by induction that Euler gives
   $y_n=(1+h)^n$. Expand to prove first-order accuracy.
3. **RK4 weights.** Verify that the RK4 weights $(1,2,2,1)/6$ reproduce the
   Taylor expansion of $e^h$ through $h^4$.
4. **Stability boundary.** For $y'=\lambda y$ with real $\lambda<0$, find the
   approximate boundary of the RK4 stability region. What happens for
   $\lambda h=-3$?
5. **Energy check.** Integrate the paraboloid geodesic with $h=0.1,0.05,0.025$
   and compare the energy drift. Estimate the observed order.
6. **Adaptive idea.** Suppose you have two RK4 solutions with steps $h$ and
   $h/2$. How could you estimate the error of the finer one? (This is the
   idea behind Richardson extrapolation.)
7. **Tangent frame.** Show that `tangentFrame` produces an orthonormal basis
   when $\Phi_u,\Phi_v$ are linearly independent.
8. **Initial velocity.** Derive $\xi=I^{-1}J^{\mathsf T}w$ from $w=J\xi$ and
   $I=J^{\mathsf T}J$.
