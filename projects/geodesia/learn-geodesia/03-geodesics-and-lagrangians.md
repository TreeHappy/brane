# 3️⃣ Geodesics — curves whose bending is purely normal

> *Chapter 2 ended at a door: the condition $T\gamma'' = 0$ that the hill
> walk's winner satisfied — and the conserved ruler-speed it walked away
> with, $E = 0.680000$ from start to finish. This chapter opens the door:
> "straightest" gets three equivalent definitions, the condition becomes an
> ordinary ODE the library integrates, and the conservation becomes a
> theorem with a symmetry behind it.*

> **What this chapter is.** One equation, three derivations, many surfaces:
> a **geodesic** is a curve whose acceleration $\gamma''$ has *no tangent
> part* — the surface's own bending supplies all of it. Dotting that
> condition with the basis vectors turns it into the ODE
> $I\ddot q + g = 0$ (§3.2), which the variational route re-derives from
> the length functional (§3.3) and which carries two conserved quantities
> for free: the speed $E$ always, and a pairing $p_\ell$ whenever the
> surface has a symmetry (§3.4). The chapter is written like the Lean
> documentation: every boxed formula is followed by the **Lean
> declarations that compute it** — a **Signature** block — and every
> numeric claim is a **Checked example**: the Lean snippet together with
> the output it really prints. Nothing in an `Output:` block is
> hand-computed; each was produced by running the snippet against the
> library. The math notation is kept in full ($I$, $g$, $p_\ell$,
> $\Gamma$, …); the Lean identifiers are the computational shadows of the
> same objects.
>
> Sections are numbered `§3.1–3.7`, matching the cross-references from
> Chapters 4–8, `index.md` and `GA_VIEW.md`. The deep geometric-algebra
> and category-theory treatments live in `GA_VIEW.md`; this chapter ends
> with a coda connecting to them.
>
> **A note on words.** This chapter never says "arrow" for a mathematical
> object: each vector is called by what it *is* — acceleration vector,
> velocity, basis vector, curvature vector, normal vector, ambient vector
> — reserving "arrow" for the drawn glyphs in the figures.

**Reproducing the checked examples.** All snippets run against this project:

```bash
cd projects/geodesia
lake build
lake exe geodesia                   # the demo printed in Ch. 1 §1.6
lake env lean <your-scratch>.lean   # any snippet below (it begins with `import Geodesia`)
```

The dictionary this chapter uses (Ch. 1–2 tabulated the rest):

| Math object | Lean declaration | File |
|---|---|---|
| $\Phi(u,v)$, $\Phi_i$, $\Phi_{ij}$, $n$ | `Surface.point`, `dPhi`, `d2Phi`, `normalAt` | `Geodesia/Surface.lean` |
| $J\,d$ — the tangent vector of a parameter direction | `Surface.tangentVec` | `Geodesia/Surface.lean` |
| $I = J^{\mathsf T}J$ — the ruler | `Surface.metricI`, `Surface.firstForm` | `Geodesia/Curvature.lean` |
| $T w = w - (w\cdot n)\,n$ — the tangent part of the split | `Surface.tangentialPart` | `Geodesia/Surface.lean` |
| $\Gamma^k_{\;ij}$ (coordinates of $T\Phi_{ij}$) | `Surface.christoffel` | `Geodesia/Curvature.lean` |
| $(u'', v'') = -\,\Gamma^k_{\;ij}\dot q^i\dot q^j$ — the geodesic ODE | `Geodesic.accel` | `Geodesia/Geodesic.lean` |
| first-order right-hand side $(\dot u, \dot v, \ddot u, \ddot v)$ | `Geodesic.rhs` | `Geodesia/Geodesic.lean` |
| one RK4 step / $n$ steps | `Geodesic.rk4Step` / `Geodesic.integrate` | `Geodesia/Geodesic.lean` |
| $E = \tfrac12\,\dot q^{\mathsf T}I\,\dot q$ — the conserved speed squared | `Geodesic.energy` | `Geodesia/Geodesic.lean` |
| $\lvert\gamma'\rvert = \sqrt{2E}$ | `Geodesic.speed` | `Geodesia/Geodesic.lean` |
| orthonormal frame $(e_1, e_2)$ of the tangent plane | `Geodesic.tangentFrame` | `Geodesia/Geodesic.lean` |
| shoot and land (Ch. 5's subject) | `Geodesic.shoot` | `Geodesia/Geodesic.lean` |
| exact derivatives (stored as functions) | `Smooth2` + its combinators | `Geodesia/Diff.lean` |

**The notation bridge.** The Lean code uses short identifiers where the text
uses math symbols. Every occurrence is translatable by this table:

| Lean notation | Vector-calculus meaning |
|---|---|
| `GeoState` fields `u v du dv` | $q(s) = (u,v)$ and $q'(s) = (u',v')$ along the curve |
| `Vec2` fields `.x`, `.y` | the two components of a 2-vector |
| `Vec3` fields `.x`, `.y`, `.z` | the three components of an ambient vector |
| `Mat2` with fields `a`, `b`, `c` | the symmetric matrix with rows $(a,\,b)$ and $(b,\,c)$ |
| `Surface.metricI s u v` | the Gram matrix $I(u,v)$ at the point $(u,v)$ |
| `Surface.firstForm s u v d` | the quadratic form $I(d, d) = d^{\mathsf T}I\,d$ |
| `Surface.dPhi s 0 u v` / `s 1 u v` | the basis vectors $\Phi_u$ / $\Phi_v$ |
| `Surface.tangentVec s u v ⟨du, dv⟩` | the surface velocity $J\,\dot q$ |
| `Geodesic.accel s st` | the parameter acceleration $(u'', v'') = -\,I^{-1}g$ |
| `Geodesic.energy s st` | $E = \tfrac12\,\dot q^{\mathsf T}I\,\dot q$ |
| `Geodesic.speed s st` | the ruler-speed $\lvert J\,\dot q\rvert = \sqrt{2E}$ |
| `Geodesic.integrate s h n st` | $n$ RK4 steps of size $h$ from the state `st` |
| `Geodesic.tangentFrame s u v` | the orthonormal frame $(e_1, e_2)$ at $(u,v)$ |
| `Geodesic.shoot s u v angle length steps` | where a shot lands (initial-value aim; Ch. 5) |
| `pConserved s st` (helper defined in §3.4) | the pairing $p_\ell = \Phi_\ell\cdot\gamma'$ at $\ell = v$ |

(Implementations are shown verbatim from the sources, with the fully
qualified names they get from their `namespace`.)

---

## 3.1 Three equivalent definitions

Chapter 2 §2.10 left one condition on the table: $T\gamma'' = 0$, the split
of Ch. 1 §1.5 applied to an acceleration. Here is the same condition stated
three ways — all the same equation:

1. **Geometric (vector).** $\;T\gamma'' = 0$: the acceleration vector of the
   surface curve points *only* along the normal — zero component inside the
   tangent plane.
2. **Variational.** A critical point of the length functional
   $\ell[\gamma] = \int \lvert\gamma'\rvert\, ds$ (or of the quadratic
   functional $\int \tfrac12\,\lvert J\,\dot q\rvert^2 ds$).
3. **Intrinsic straightness.** Cut the surface locally and lay it flat: the
   geodesic becomes a straight line — the straightest possible curve on the
   surface.

> **Caveat.** A geodesic is *locally* length-minimising, not globally
> shortest: on a sphere the long arc of a great circle is a geodesic too,
> and between antipodal points every great circle ties for the same length.

The chapter proves the vector form and the variational form produce the same
ODE (§3.2 = §3.3), then cashes the equation in twice: once as conserved
quantities (§3.4), once as concrete surfaces to walk (§3.5).

---

## 3.2 The vector condition, and the equation it produces

Let $\gamma(s) = \Phi(q(s))$ be a surface curve with parameter velocity
$\dot q = (\dot u, \dot v)$. Differentiating twice with the chain rule:

```math
\gamma' = \dot u\,\Phi_u + \dot v\,\Phi_v,
```

```math
\gamma'' = \ddot u\,\Phi_u + \ddot v\,\Phi_v
+ \Phi_{uu}\,\dot u^2 + 2\,\Phi_{uv}\,\dot u\dot v + \Phi_{vv}\,\dot v^2.
```

Everything on the right is an ordinary 3-D vector combination — no
coordinate machinery yet. The geometric definition says its **tangential
projection must vanish**, $T\gamma'' = 0$, i.e.
$\gamma'' = \bigl(n\cdot\gamma''\bigr)n$. To get a *differential equation
for the two scalar functions* $(u,v)$, take it the Chapter-2 way: dot each
basis vector into $\gamma''$ (the normal part drops out of these dot
products, because $\Phi_u$ and $\Phi_v$ are both orthogonal to $n$):

```math
\Phi_k \cdot \gamma'' \;=\; 0, \qquad k \in \{u, v\}.
```

Expand — orthogonality to $n$ kills the second-half terms of the split:

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
2-component vector assembled purely by *dotting basis vectors with
curvature vectors* — vectors, twice, no symbols. Solved for the
acceleration:

```math
\ddot q = -\,I^{-1}\,g.
```

That is the **geodesic equation** in this course's native currency: compare
by eye against the boxed form the $\Gamma$'s usually occupy,

```math
\underbrace{\Bigl(\ddot q^k + \Gamma^k_{\;ij}\dot q^i\dot q^j\Bigr)}_{\text{tangent-plane part}} = 0
\quad\Longleftrightarrow\quad
I\,\ddot q + g = 0,
```

— the two right-hand sides are identical because
$g_k = \Gamma_{kij}\,\dot q^i\dot q^j$, where
$\Gamma_{kij} = \Phi_k\cdot\Phi_{ij}$ are the "first-kind" coordinates of
the same tangent parts (Chapter 2 §2.8). The vector form is nicer for two
reasons: it does not need the $\Gamma$ array at all — just vectors and dot
products — and it *exposes* what the equation is: *make the parameter
acceleration exactly cancel the turning of the basis vectors, so the
ambient acceleration vector $\gamma''$ is left purely normal.*

Note what does *not* appear: the normal part of $\gamma''$ never vanishes on
a curved surface. A geodesic is the curve that lets the surface's own
bending supply *all* of its acceleration.

**Signature.** The library's ODE right-hand side for $(u'', v'')$ — the
coefficient form of the boxed equation, with the tangent-part coordinates
$\Gamma^k_{\;ij}$ contracted against the velocity:

```lean
/-- Parameter-space acceleration `(u'', v'')` from the geodesic equation. -/
def Geodesic.accel (s : Surface) (st : GeoState) : Vec2
```

**Implementation.** Verbatim from `Geodesia/Geodesic.lean`:

```lean
def Geodesic.accel (s : Surface) (st : GeoState) : Vec2 :=
  let Γ := Surface.christoffel s st.u st.v
  let g1uu := Γ 0 0 0
  let g1uv := Γ 0 0 1
  let g1vv := Γ 0 1 1
  let g2uu := Γ 1 0 0
  let g2uv := Γ 1 0 1
  let g2vv := Γ 1 1 1
  let du := st.du
  let dv := st.dv
  ⟨-(g1uu * du * du + 2 * g1uv * du * dv + g1vv * dv * dv),
   -(g2uu * du * du + 2 * g2uv * du * dv + g2vv * dv * dv)⟩
```

**Vector-calculus reading.** Each `Γ k i j` is the coefficient of $\Phi_k$
inside $T\Phi_{ij}$ (Chapter 2 §2.8), so each bracket is one component of
$-\,g$ with the first-kind coordinates swapped for the second-kind ones via
$I^{-1}$ — the boxed equation solved for $\ddot q$, term by term. Nothing
here is new machinery: it is the projection $T$ of Chapter 1, wearing
coefficients.

**The surfaces of this chapter.** The three walked below, in one fence so
the checked examples can share a scratch file — the paraboloid hill of
§2.9, the polar plane of §2.5, the unit sphere of §2.6, plus the pairing
helper of §3.4:

```lean
import Geodesia
open Geodesia

def paraboloid : Smooth2 :=
  Smooth2.smul 0.3
    (Smooth2.add (Smooth2.mul Smooth2.coordX Smooth2.coordX)
      (Smooth2.mul Smooth2.coordY Smooth2.coordY))
def surf : Surface := Surface.ofHeight paraboloid

def polarPlane : Surface where
  sx := Smooth2.mul Smooth2.coordX (Smooth2.cos Smooth2.coordY)
  sy := Smooth2.mul Smooth2.coordX (Smooth2.sin Smooth2.coordY)
  sz := Smooth2.const 0

def sphere : Surface where
  sx := Smooth2.mul (Smooth2.sin Smooth2.coordX) (Smooth2.cos Smooth2.coordY)
  sy := Smooth2.mul (Smooth2.sin Smooth2.coordX) (Smooth2.sin Smooth2.coordY)
  sz := Smooth2.cos Smooth2.coordX

/-- The pairing `p_ℓ = Φ_ℓ · γ′` for `ℓ = 1` (the second coordinate):
conserved whenever the metric is independent of that coordinate. -/
def pConserved (s : Surface) (st : GeoState) : Float :=
  (Surface.dPhi s 1 st.u st.v).dot
    (Surface.tangentVec s st.u st.v ⟨st.du, st.dv⟩)
```

**Checked example.** The boxed equation, verified at a single point with no
integrator involved — on the paraboloid at $(1,0)$ with parameter velocity
$(1,0)$. There the velocity-vector sum in the $g$ formula is just
$\Phi_{uu}$, so $g_u = \Phi_u\cdot\Phi_{uu}$ and
$g_v = \Phi_v\cdot\Phi_{uu}$:

```lean
#eval show IO Unit from do
  let s := surf
  -- q̇ = (1,0): the vector sum in the boxed equation is just Φ_uu
  let gq := (Surface.dPhi s 0 1 0).dot (Surface.d2Phi s 0 0 1 0)
  let gv := (Surface.dPhi s 1 1 0).dot (Surface.d2Phi s 0 0 1 0)
  IO.println s!"g at (1,0), q̇ = (1,0):   g_u = {gq}   g_v = {gv}"
  let st : GeoState := ⟨1.0, 0.0, 1.0, 0.0⟩
  let a := Geodesic.accel s st
  IO.println s!"accel = (u″, v″)          = ({a.x}, {a.y})"
  -- the boxed equation, checked: I · accel + g  must be (0, 0)
  let I := Surface.metricI s 1 0
  let r1 := I.a * a.x + I.b * a.y + gq
  let r2 := I.b * a.x + I.c * a.y + gv
  IO.println s!"I·accel + g               = ({r1}, {r2})"
```

**Output:**
```text
g at (1,0), q̇ = (1,0):   g_u = 0.360000   g_v = 0.000000
accel = (u″, v″)          = (-0.264706, -0.000000)
I·accel + g               = (0.000000, 0.000000)
```

**How to read the numbers.** Three facts line up. First,
$g_u = 0.360000 = \Phi_u\cdot\Phi_{uu}$ — the $\Phi_u$-component of the
velocity-vector sum, exactly the box's definition. Second,
$u'' = -0.264706$ is minus the very number §2.7 printed as
$\Gamma^u_{\;uu}$ — the coefficient form and the vector form agree on the
same point. Third, the residual is $I\cdot\text{accel} + g = (0.000000,\,0.000000)$:
the two forms of the equation are the same equation, to the last float
digit, with no integrator in sight. (The printed $-0.000000$ in
the $v''$ slot is a negative zero from the float arithmetic — magnitude
below $10^{-15}$, kept verbatim as the program emitted it.)

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
"a tangent part and a normal part" of Chapters 1–2 is all it takes to go
between them.

Like any Lagrangian system, this one automatically transports its own
"energy":

```math
\frac{d}{ds}\Bigl(\tfrac12\,\dot q^{\mathsf T} I\, \dot q\Bigr) = 0,
```

i.e. the surface velocity $\gamma'$ keeps constant length along the curve —
the conserved quadratic form used as a numerical health check in the code.

---

## 3.4 The conserved quadratic form, and what symmetry adds

Direct check straight from the boxed equation — notice it needs *no*
$\Gamma$'s, only $I\ddot q + g = 0$ and the dot-product expansion of $g$:

```math
\frac{d}{ds}\Bigl(\tfrac12\,\dot q^{\mathsf T} I\,\dot q\Bigr)
= \dot q^{\mathsf T}\bigl(I\,\ddot q + g\bigr) = 0.
```

Consequences:

- $\lvert J\dot q\rvert = \sqrt{2 \times (\text{that value})}$ is
  **constant**: a geodesic parameterized this way is automatically
  constant-speed.
- If $s \mapsto as+b$ is an affine reparameterization, the equation is
  unchanged; arc length is a conventional choice.
- For the raw length functional
  $\ell = \int \sqrt{\dot q^{\mathsf T} I\,\dot q}\, ds$, the equations only
  determine the *image* of the curve; picking constant speed reduces them
  to the same geodesic equation. The quadratic functional is the
  computationally convenient version of the same geometry.

**Checked example — the speed never changes.** Integrate the paraboloid
geodesic of §2.9 (start $\langle 1, 0, 1, 0\rangle$, step $h = 0.01$,
$100$ RK4 steps) and print the state, the ruler-speed and $E$ at five
checkpoints:

```lean
#eval show IO Unit from do
  let s := surf
  let mut st : GeoState := ⟨1.0, 0.0, 1.0, 0.0⟩
  let ss := ["0.000", "0.250", "0.500", "0.750", "1.000"]
  for i in [0:5] do
    IO.println s!"paraboloid  s={ss[i]!}  (u,v) = ({st.u}, {st.v})  speed = \
      {Geodesic.speed s st}  E = {Geodesic.energy s st}"
    if i < 4 then
      st := Geodesic.integrate s 0.01 25 st
```

**Output:**
```text
paraboloid  s=0.000  (u,v) = (1.000000, 0.000000)  speed = 1.166190  E = 0.680000
paraboloid  s=0.250  (u,v) = (1.241823, 0.000000)  speed = 1.166190  E = 0.680000
paraboloid  s=0.500  (u,v) = (1.467985, 0.000000)  speed = 1.166190  E = 0.680000
paraboloid  s=0.750  (u,v) = (1.679912, 0.000000)  speed = 1.166190  E = 0.680000
paraboloid  s=1.000  (u,v) = (1.879182, 0.000000)  speed = 1.166190  E = 0.680000
```

$E = 0.680000$ at every checkpoint — the same value §2.9 measured at the
endpoints, now held across the whole walk. And the speed
$1.166190 = \sqrt{2 \times 0.68} = \sqrt{1.36}$ is exactly the length of
$J\,(1,0)$ at the start: $\Phi_u\cdot\Phi_u = 1.36$ there, so the
ruler-speed starts and stays at $\sqrt{1.36}$. The geodesic's parameter
speed $\lvert\dot q\rvert$ does sink along the way (the $u$-column:
$1 \to 1.879$ in unit parameter time) — but only because the ruler grows
outward; the *ruler-speed*, what the surface's own measuring stick reads,
never budges.

![The height map of the paraboloid hill. Left: contours of equal height h = 0.3(u²+v²) with three paths leaving (1,0) with the same initial parameter velocity (1,0) — the constant-radius ring rides one contour (equal height = equal I-cost ring), the parameter-straight line runs out to (2,0), and the geodesic (green) ends at u = 1.879182, hoarding its speed where the contours are expensive. Right: along the geodesic the ruler-speed |J q̇| stays flat (E = 0.680000 conserved) while its parameter speed sinks; the straight line's ruler-cost per parameter second rises as the contours spread apart.](images/height-map.png)

**Conserved quantities from symmetry.** Whenever the Gram matrix $I$ does
not depend on a coordinate $q^\ell$, the pairing

```math
p_\ell = \bigl(I\,\dot q\bigr)_\ell = \Phi_\ell\cdot\gamma'
```

is conserved — i.e. the projection of the surface velocity along a symmetry
direction of the surface cannot change:

- **Polar plane:** $\Phi_\theta\cdot\gamma' = r^2\dot\theta$ is constant —
  the straightness of the geodesic, seen in rotating coordinates.
- **Sphere:** $\sin^2\theta\,\dot\varphi$ is constant.
- **Surfaces of revolution** $\Phi(r,\theta) = (r\cos\theta, r\sin\theta, f(r))$:
  the same $\Phi_\theta$-pairing is constant (historically "Clairaut's
  relation").

This is a practical numerical check: a correct integrator must conserve
these pairings, and any measurable drift is a red flag. The helper
`pConserved` in the setup fence computes $p_\ell$ at $\ell = v$ — the
second coordinate of whichever surface it is handed.

**Checked example — the polar plane's pairing.** Start at
$(r,\theta) = (1.5,\,0.7)$ with velocity $(0.8,\,0.6)$; at every
checkpoint the pairing equals its starting value
$1.5^2 \times 0.6 = 1.350000$, and the direct product $r^2\dot\theta$
agrees digit for digit:

```lean
#eval show IO Unit from do
  let s := polarPlane
  let mut st : GeoState := ⟨1.5, 0.7, 0.8, 0.6⟩
  let ss := ["0.000", "0.250", "0.500", "0.750", "1.000"]
  for i in [0:5] do
    let direct : Float := st.u * st.u * st.dv
    IO.println s!"polarPlane  s={ss[i]!}  (r,θ) = ({st.u}, {st.v})  \
      p_θ = {pConserved s st}  r²θ̇ = {direct}  speed = {Geodesic.speed s st}"
    if i < 4 then
      st := Geodesic.integrate s 0.01 25 st
```

**Output:**
```text
polarPlane  s=0.000  (r,θ) = (1.500000, 0.700000)  p_θ = 1.350000  r²θ̇ = 1.350000  speed = 1.204159
polarPlane  s=0.250  (r,θ) = (1.714825, 0.831588)  p_θ = 1.350000  r²θ̇ = 1.350000  speed = 1.204159
polarPlane  s=0.500  (r,θ) = (1.952562, 0.932557)  p_θ = 1.350000  r²θ̇ = 1.350000  speed = 1.204159
polarPlane  s=0.750  (r,θ) = (2.205816, 1.010998)  p_θ = 1.350000  r²θ̇ = 1.350000  speed = 1.204159
polarPlane  s=1.000  (r,θ) = (2.469818, 1.072988)  p_θ = 1.350000  r²θ̇ = 1.350000  speed = 1.204159
```

Read the drift: $r$ walks from $1.5$ out to $2.47$ — the geodesic is a
straight line in the plane, and straight lines leave the origin. Since
$r^2\dot\theta$ is pinned at $1.350000$ while $r$ grows, the angular speed
$\dot\theta$ must fall like $1/r^2$: the rotating coordinates see the same
straight line slowing down azimuthally. The speed column is constant too —
the plane's ruler is the plain $\delta$, so this is the ordinary
constant velocity of a free straight line, $1.204159$ throughout.

**Checked example — the sphere's pairing.** Start at
$(\theta,\varphi) = (1.0,\,0.2)$ with velocity $(0.8,\,0.6)$ — deliberately
*not* a meridian, so the great circle has to drift. The pairing
$p_\varphi = \sin^2\theta\,\dot\varphi$ holds at $0.424844$ from start to
finish:

```lean
#eval show IO Unit from do
  let s := sphere
  let mut st : GeoState := ⟨1.0, 0.2, 0.8, 0.6⟩
  let ss := ["0.000", "0.250", "0.500", "0.750", "1.000"]
  for i in [0:5] do
    let direct : Float := (Float.sin st.u) * (Float.sin st.u) * st.dv
    IO.println s!"sphere      s={ss[i]!}  (θ,φ) = ({st.u}, {st.v})  \
      p_φ = {pConserved s st}  sin²θ·φ̇ = {direct}  \
      speed = {Geodesic.speed s st}"
    if i < 4 then
      st := Geodesic.integrate s 0.01 25 st
```

**Output:**
```text
sphere      s=0.000  (θ,φ) = (1.000000, 0.200000)  p_φ = 0.424844  sin²θ·φ̇ = 0.424844  speed = 0.945995
sphere      s=0.250  (θ,φ) = (1.204101, 0.334357)  p_φ = 0.424844  sin²θ·φ̇ = 0.424844  speed = 0.945995
sphere      s=0.500  (θ,φ) = (1.413343, 0.448735)  p_φ = 0.424844  sin²θ·φ̇ = 0.424844  speed = 0.945995
sphere      s=0.750  (θ,φ) = (1.624478, 0.555634)  p_φ = 0.424844  sin²θ·φ̇ = 0.424844  speed = 0.945995
sphere      s=1.000  (θ,φ) = (1.834986, 0.665015)  p_φ = 0.424844  sin²θ·φ̇ = 0.424844  speed = 0.945995
```

The great-circle drift, live: $\theta$ climbs from $1.0$ toward
$1.834986$ — the curve is walking toward the equator, where $\sin\theta$
is largest. With $\sin^2\theta\,\dot\varphi$ pinned at $0.424844$, the
azimuthal speed $\dot\varphi$ must shrink as $\sin^2\theta$ grows: near the
equator a fixed amount of eastward momentum buys fewer radians per second,
because the latitude circles are wider there. That is the whole content of
"the projection of the velocity along a symmetry direction cannot change"
— the symmetry here is the sphere's rotation about its polar axis, and
$p_\varphi$ is the momentum of that rotation.

Vector framing once more: for the coordinate function $q^\ell$, the pairing
$\Phi_\ell\cdot\gamma'$ is "the surface velocity's $\ell$-component against
the symmetric-direction basis vector" — a coordinate-free statement of *the
projection of the velocity along a symmetry cannot change*.

---

## 3.5 Worked examples

### Plane in Cartesian coordinates

$I$ is the identity; $g = 0$ because every curvature vector $\Phi_{ij}$
vanishes; so $\ddot q = 0$: geodesics are straight lines at constant
parameter speed.

### Plane in polar coordinates

With the polar basis vectors of Chapter 2 (§2.5), the two boxed component
equations $I\ddot q = -g$ become

```math
\ddot r - r\,\dot\theta^2 = 0,
\qquad
\ddot\theta + \frac{2}{r}\,\dot r\,\dot\theta = 0.
```

These are straight lines written in rotating coordinates. The conserved
pairing $r^2\dot\theta$ — verified to five checkpoints in §3.4 — is the
projection of the straight line's constant velocity onto the rotating basis
vector $\Phi_\theta$; the $\sim 1/r$ in the second equation comes from that
vector's changing length $r$.

### Sphere

With the sphere basis vectors of Chapter 2 (§2.6):

```math
\ddot\theta - \sin\theta\cos\theta\,\dot\varphi^2 = 0,
\qquad
\ddot\varphi + 2\cot\theta\,\dot\theta\dot\varphi = 0.
```

Solutions are great circles. Meridians ($\varphi = \text{const}$) give
$\ddot\theta = 0$; the equator ($\theta = \pi/2$) gives
$\ddot\varphi = 0$. The conserved pairing is
$\sin^2\theta\,\dot\varphi$ — verified at $0.424844$ in §3.4 — the
projection of the velocity against the rotation-symmetry vector
$\Phi_\varphi$. The §3.4 table shows the general (non-meridian) case
drifting toward the equator as the pairing holds fixed.

### Paraboloid

For $z = f(x,y)$ the Gram matrix is

```math
I = \begin{bmatrix}1+f_x^2 & f_xf_y \\ f_xf_y & 1+f_y^2\end{bmatrix};
```

the $g$-vector is assembled from the curvature vectors
$\Phi_{xx}, \Phi_{xy}, \Phi_{yy}$ exactly as in §3.2's box; and no
closed-form solution exists in general — which is precisely why the
numerical methods of Chapters 4–6 exist. The §3.4 table (and the height map
figure there) shows the integrated geodesic holding
$E = 0.680000$ across the whole walk — a first check that the ODE and the
integrator agree.

### Cone

A cone is developable: cut it along a generator and unroll it to a flat
sector, where geodesics become straight lines. Beautiful check of the
intrinsic nature of geodesics: $K = 0$, and the "shortest path" is a
straight line after unrolling.

---

## 3.6 How the code implements it

```lean
Geodesic.accel s st      -- (u'', v''): the boxed equation of §3.2
Geodesic.rhs s st        -- first-order system (u, v, u', v')
Geodesic.rk4Step s h st  -- one RK4 step
Geodesic.integrate s h n st
Geodesic.energy s st     -- ½ · dot qᵀ I dot q  (the conserved quadratic form)
Geodesic.speed s st      -- |J·dotq|  =  √(2 · that value)
```

`accel` (implementation in §3.2) computes exactly the boxed equation of
§3.2, up to coordinates. In the code the tangent-part coefficients come
pre-computed as the `Surface.christoffel` array (those are, again, purely
the coordinates of $T\Phi_{ij}$ in the basis — the routine *is* the vector
formula of Chapter 2 §2.8, in coefficient clothing), and `accel` contracts
them with the direction numbers — the coefficient-image of the vector
recipe "form the vector sum
$\dot u^2\Phi_{uu} + 2\dot u\dot v\Phi_{uv} + \dot v^2\Phi_{vv}$, dot it
with each basis vector to get $g$, solve $I\ddot q = -g$" (the boxed
equation of §3.2).

The demo (`lake exe geodesia`) integrates a geodesic on the paraboloid and
prints the conserved quadratic form holding to the printed precision — the
discrete analogue of §3.4's identity, and the same health check the three
tables there put on display.

---

## 3.7 Exercises

1. **Normal acceleration.** Show that a geodesic satisfies
   $n\cdot\gamma'' = II(\dot q,\dot q)$ and that its *tangential projection*
   vanishes — then verify both numerically on the sphere (the §3.4 setup
   fence has `sphere` ready; `Surface.tangentialPart` and
   `Surface.normalPart` do the split).
2. **Variational algebra.** Fill in the symmetrization step in §3.3 and
   explain why contracting with $\dot q^i \dot q^j$ lets you replace
   $\partial_i g_{kj}$ by
   $\tfrac12(\partial_i g_{kj} + \partial_j g_{ki})$ inside $I$-derivative
   terms.
3. **Length vs quadratic functional.** Derive Euler–Lagrange for
   $\ell = \int\sqrt{\dot q^{\mathsf T}I\,\dot q}\,ds$ and show it reduces
   to the geodesic equation when $\lvert\gamma'\rvert$ is constant.
4. **Polar straight lines.** Show that
   $r\cos(\theta - \theta_0) = d$ satisfies both polar component equations.
5. **Sphere great circles.** Verify that
   $\gamma(s) = \cos s\,U + \sin s\,V$ with two orthonormal vectors $U,V$
   satisfies the sphere component equations — *first* by checking that its
   surface second derivative points purely radially (the vector test!),
   *then* from the components.
6. **Symmetry pairing.** Prove that $\partial_\ell\,I = 0$ implies
   $p_\ell = \Phi_\ell\cdot\gamma'$ is constant along a geodesic. *(The
   §3.4 tables are the numerical shadow of exactly this statement.)*
7. **Surface of revolution.** Derive the Gram matrix for
   $\Phi(r,\theta) = (r\cos\theta, r\sin\theta, f(r))$ and the corresponding
   conserved pairing.
8. **Affine parameter.** Show that if $q(s)$ solves the geodesic equation,
   then $q(as+b)$ also solves it. What happens to the conserved quadratic
   form?
9. **Code check.** Run `lake exe geodesia` and verify that the conserved
   quadratic form is unchanged after 100 RK4 steps. Then change the step
   size and see how the error behaves. *(Compare with §3.4's tables: at
   $h = 0.01$, $E = 0.680000$ at every checkpoint.)*

---

## 🧭 Optional coda — the same chapter through geometric algebra and category theory

> **Optional.** Everything above is complete without this coda. Read on if
> you want to see the same objects — the geodesic condition, the conserved
> speed, the symmetry pairing — wearing two other coats: multivectors
> ($\mathrm{Cl}_3$) and universal properties (category theory). The deep GA
> treatment lives in `GA_VIEW.md`; each chapter of this tutorial ends with
> a coda like this one.

### The GA view: the acceleration is parallel to the normal

Chapter 1–2's codas used one unit 2-blade
$\hat B = (\Phi_u \wedge \Phi_v)/\lvert\Phi_u \wedge \Phi_v\rvert$,
$\hat B^2 = -1$, with the tangency test $w \wedge \hat B = 0$. The geodesic
condition is the *mirror image* of that test: the acceleration passes the
test for being **normal**,

```math
\gamma'' \wedge n = 0
\qquad\Longleftrightarrow\qquad
\gamma'' = (\gamma''\cdot n)\,n
\qquad\Longleftrightarrow\qquad
\underbrace{(\gamma''\cdot\hat B)\,\hat B^{-1}}_{\text{the tangent part } T\gamma''} = 0.
```

The sandwich projection onto the blade vanishes, and the whole acceleration
lives on the normal 1-blade. And the size of that normal part is not
incidental — exercise 1 shows $n\cdot\gamma'' = II(\dot q,\dot q)$: the
acceleration's magnitude is the second fundamental form evaluated on the
velocity, so "purely normal acceleration" and "the bending table drives the
curve" are one statement.

| This chapter | GA form |
|---|---|
| geodesic condition $T\gamma'' = 0$ | $\gamma'' \wedge n = 0$ — acceleration parallel to the normal |
| the tangent part $T\gamma''$ | the blade projection $(\gamma''\cdot\hat B)\,\hat B^{-1}$, here $0$ |
| the normal part $(\gamma''\cdot n)\,n$ | the *whole* acceleration: $\gamma'' = (\gamma''\cdot n)\,n$ |
| $n\cdot\gamma'' = II(\dot q,\dot q)$ (exercise 1) | the acceleration's magnitude is the second form on the velocity |
| conserved pairing $p_\ell = \Phi_\ell\cdot\gamma'$ | an inner product of two vectors in the blade's plane |

### The CT view: the geodesic is what the symmetries conserve

The variational form of §3.3 turns the chapter into a showcase for
Noether's principle: every continuous symmetry of the action yields a
conserved quantity. Both conserved quantities of §3.4 are of this kind:

- $E = \tfrac12\,\dot q^{\mathsf T}I\,\dot q$ is conserved because the
  Lagrangian has no explicit $s$-dependence — "shifting the whole walk in
  time" is a symmetry of the action.
- $p_\ell = \Phi_\ell\cdot\gamma'$ is conserved when $\partial_\ell I = 0$:
  the surface is invariant under the one-parameter family of "slides" in
  the coordinate $\ell$ — the rotations of the polar plane about the
  origin, the longitudes of the sphere — a group action on the object
  $\mathbb{R}^2$ that the ruler $I = \Phi^{*}\delta$ does not notice. The
  pairing is that action's momentum.

| This chapter | CT form |
|---|---|
| $L = \tfrac12\,I(\dot q,\dot q)$ | a Lagrangian on the tangent bundle of the object $\mathbb{R}^2$ |
| geodesic equation $I\ddot q + g = 0$ | Euler–Lagrange: critical paths of the action |
| $E$ conserved | Noether: the action is invariant under $s$-translations |
| $p_\ell$ conserved when $\partial_\ell I = 0$ | Noether: the metric $\Phi^{*}\delta$ is invariant under the $\ell$-slides; $p_\ell$ is the momentum of that symmetry |
| $r^2\dot\theta$, $\sin^2\theta\,\dot\varphi$ (§3.4 tables) | the momenta of the respective rotation ($S^1$) actions |

Why the pairing *had* to be conserved: the symmetry is a relabeling that
leaves $I$ unchanged; the geodesic equation is written entirely in $I$, so
relabeling a geodesic gives a geodesic; and $p_\ell$ measures how the walk
drifts along the relabeling — a quantity nothing in the equation can
change.

### One object, three coats

| Object | Vector chapter | GA | CT |
|---|---|---|---|
| the geodesic condition | $T\gamma'' = 0$, i.e. $I\ddot q + g = 0$ | $\gamma'' \wedge n = 0$ — acceleration on the normal blade | Euler–Lagrange of the action |
| the conserved speed | $E = \tfrac12\,\dot q^{\mathsf T}I\,\dot q$ | $\lvert J\dot q\rvert$, read from the blade's inner product | Noether charge of the $s$-translations |
| the conserved pairing | $p_\ell = \Phi_\ell\cdot\gamma'$ | a dot product inside the blade's plane | momentum of the $\ell$-slide symmetry |

> **By the way.** Why end chapters with two extra lenses at all? Because
> the same fact can sit at three depths: a computation (this chapter), a
> product structure (GA), and a universal property (CT). Seeing all three
> is the difference between knowing the formula and knowing the formula
> *could not have been otherwise*. When a coda feels redundant, skip it —
> it is written to be skippable.
