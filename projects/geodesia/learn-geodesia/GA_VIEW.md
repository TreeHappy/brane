# 🌀 Geometric Algebra view — surfaces in Cl₃

*A supplement to Chapters 1–3. Nothing in the code changes; this re-reads the
same objects (tangent arrows, second-derivative arrows, geodesics, curvatures)
so that every "projection" becomes a *product* and every coordinate story has
a geometric double.*

> Why this fits geodesia so well: our whole track runs on one surface
> $\Phi : \mathbb{R}^2 \to \mathbb{R}^3$, i.e. on vectors in $\mathbb{R}^3$.
> Geometric algebra (GA, the algebra $\mathrm{Cl}_3$) is *vectors closed under
> a product that also makes directed planes and volumes*. Exactly the
> objects the surface theory keeps wanting: a tangent plane as one element, a
> normal as its dual, and "geodesic" as one commutation law.

---

## §1 The cast, in grades

GA upgrades $\mathbb{R}^3$ by closing it under the **geometric product**.
Basis arrows $e_1,e_2,e_3$ with $e_ie_i=1$, $e_ie_j=e_{ij}=-e_{ji}$
($i\neq j$), $e_{123}=I_3$ and $I_3^2=-1$ give subspaces "by grade":

| Object of the chapter | GA native |
|---|---|
| ambient arrow $w$ | grade-1 vector |
| tangent plane at a point | grade-2 **bivector** $B \propto \Phi_u\wedge\Phi_v$ |
| unit normal $n$ | grade-1 vector, obtained from $B$ by duality (below) |
| second-derivative arrows $\Phi_{ij}$ | grade-1 vectors — *nothing to upgrade* |
| "tangent part vs normal part" split | **geometric product parts** of $nw$ |
| Gaussian curvature $K$ | ratio of a *pushed blade* to its source blade |
| moving frame along a curve | **rotor** $R(s)$, one grade-0+grade-2 element |
| Darboux scalars $(\kappa_g,\kappa_n,\tau_g)$ | the **three bivector coefficients of $R^{-1}\dot R$** |

```mermaid
flowchart LR
    B["Φᵤ ∧ Φᵥ<br/>tangent blade"] -->|dualize| N["n (unit normal)"]
    W["ambient arrow w"] -->|"nwn"| REF["reflection in n"]
    REF -->|½(w − nwn)| T["T(w) = tangential part"]
    G["γ″ of a curve"] -->|"γ″∧n = 0 ?"| GEO["geodesic test"]
    R["rotor R(s)<br/>moving frame"] -->|"R⁻¹Ṙ"| OMEGA["bivector Ω:<br/>κ_g, κ_n, τ_g coefficients"]
```

---

## §2 The split, as one product (vs Chapter 2 §2.4)

Chapter 2's workhorse was the tangential projection
$T w = w - (w\cdot n)\,n$. Take a **unit** normal $n$; reflecting $w$ in the
normal *line* gives (GA reflection identity) $n w n = 2\,(w\cdot n)n - w$.
Half the difference from $w$ is exactly the tangential part, half the sum is
exactly the normal part:

```math
\boxed{\;T\,w \;=\; \tfrac12\bigl(w - n\,w\,n\bigr),\qquad
(w\cdot n)\,n \;=\; \tfrac12\bigl(w + n\,w\,n\bigr).\;}
```

One geometric product, both projections; no linear-algebra machinery. The
Gauss split of $\Phi_{ij}$ — "tangent part $+$ bending part" — becomes:

```math
\Phi_{ij} \;=\; \tfrac12(\Phi_{ij} - n\Phi_{ij}n) \;+\; \tfrac12(\Phi_{ij}+ n\Phi_{ij}n).
```

Everything Chapter 2 calls "covariant derivative" — $T(\partial V/\partial
u_i)\,$ — is then literally: differentiate the arrows (ordinary calculus),
*and halve a sandwich product*:

```math
\nabla_i\,V \;=\; \tfrac12\Bigl(\partial_i V - n\,\partial_i V\,n\Bigr).
```

The sandwich survives differentiation unchanged (it's linear), so the whole
connection becomes "project with the sandwich at each point".

---

## §3 The normal is the tangent blade, dualized

In Chapter 1 the normal comes from the cross product, which is available
only in dimension 3. GA generalizes it: the *tangent plane* is the blade

```math
B \;=\; \Phi_u\wedge\Phi_v \qquad (\text{a bivector, magnitude } \sqrt{\det I}),
```

and its **Hodge dual** — multiplication by the pseudoscalar $I_3$ (division
is well-defined since $I_3^{-1} = -I_3$, and here by the magnitude too) —
lands in grade 1:

```math
n \;=\; \frac{B\,I_3^{\,\leftarrow}}{\lvert B\rvert},
\qquad I_3 = e_{123},\; I_3^2 = -1.
```

Concretely: dualizing turns the *oriented plane* into the *oriented normal
arrow*, automatically with the right orientation, and the magnitude division
is the normalizing step that makes $n$ unit. So the pair
$(t_1,t_2 \to B \to n)$ is the same content as Chapter 1's Gram-vs-cross
jump — with the plane, not a pair of arrows, as the carrier object. (Blades
also make the "two projections of one arrow" motto of Ch. 2 a generative
principle: project *any* arrow into/out of the plane by sandwiching with
$n$, or grade-wedge to test membership.)

---

## §4 Geodesics: one wedge-test

Chapter 3 defined geodesics via $T\gamma''=0$, $\gamma''$ purely normal. In
GA, "parallel to $n$" has a certificate: two vectors are parallel exactly
when their wedge vanishes. So:

```math
\boxed{\;\gamma'' \wedge n \;=\; 0\;}
```

—is the entire geodesic condition. It is the *same content* as
$\gamma'' = (n\cdot\gamma'')n$, but as an oriented element: the acceleration
has no component in the tangent *blade*. The wedge makes "no tangential
part" a property of one element rather than a subtraction performed twice.

The Eulerian variational route (Chapter 3 §3.3) still ends at
$I\ddot q + g = 0$: GA doesn't simplify the *numerics* (the ODE still lives
in coordinate functions and is integrated with RK4), it simplifies the
*statement* and the checks. See §7 for a confession about cost.

---

## §5 Curvature as blade maps

The shape operator remains a linear map $S(v) = -T\,\partial_v n$ (its
tangential image is automatically tangent, since $\partial_v|n|^2 = 0$). Two
GA readings:

**(a) Normal via the blade, again.** With $n = B/|B|$ the normal is the
*dualized blade*; differentiating the *unit* blade keeps its derivative
*within* the tangent plane — the very content of Rodrigues'
$\partial_v n = -S(v)$: $\partial_v n$ is tangent, so it lies inside $B$
itself.

**(b) Gaussian curvature as a pushed blade.** For a linear map $S$,
"blade stays blade" gives

```math
S(e_u \wedge e_v) = S(e_u)\wedge S(e_v) = K\,\bigl(e_u\wedge e_v\bigr),
```

— the area-map coefficient. **Gaussian curvature is exactly how the tangent
plane's area element rescales under $S$** — the same "push the top blade,
read the scalar" pattern as any Jacobian determinant story: one blade, one
number.

Mean curvature $H = \tfrac12\,\text{trace}(S)$ up to sign; principal
directions diagonalize $S$ *within the tangent blade* — a symmetric map
restricted to a 2-D boost: honest eigen-stuff, no coordinates implied.

---

## §6 Moving frames are rotors (the payoff)

Take an orthonormal frame $(T, h, n)$ along a curve ($h$ the tangent-plane
normal, $h = n\times T$; Chapter 5 uses it for shooting). One rotor
$R(s)$ — a grade-$0/2$ element with $R\tilde R = 1$ — carries a *reference*
frame $(E_1, E_2, E_3)$ onto it:

```math
T = R\,E_1\tilde R,\qquad h = R\,E_2\tilde R,\qquad n = R\,E_3\tilde R.
```

The *rate of the whole frame* is a single bivector:

```math
\dot R\,\tilde R \;=\; -\tfrac12\bigl(
\kappa_g\, T\wedge h
\;+\; \kappa_n\, T\wedge n
\;+\; \tau_g\, h\wedge n \bigr),
\quad\text{(up to the sign/}\tfrac12\text{ convention of §7)}
```

whose three bivector coefficients are **exactly the Darboux scalars**
(geodesic curvature, normal curvature, geodesic torsion) — Chapter 3's
conserved pairings, Chapter 5's shooting normal, all facets of one element:

| Darboux scalar | meaning | where the chapters use it |
|---|---|---|
| $\kappa_g$ | in-plane turning | **geodesic test**: $\kappa_g = 0$ ⟺ geodesic |
| $\kappa_n = II(T,T)$ | bending out of the plane | normal acceleration (Ch. 3 exercise 1) |
| $\tau_g$ | frame twist about $T$ | torsion bookkeeping in shooting (Ch. 5) |

And the punchline lands cleanly:

```math
\gamma \text{ geodesic} \;\Longleftrightarrow\; \kappa_g = 0
\;\Longleftrightarrow\; \tilde R\dot R \text{ has no } T\wedge h \text{ component}
\;\Longleftrightarrow\; \gamma''\wedge n = 0.
```

Same geometry, three certificates: *blade test*, *bivector coefficient*,
*scalar condition*. The rotor form is also how you'd implement a frame ODE
(without trig-free $3\times3$ bookkeeping): evolve $R$ with the bivector
$\dot R = R\,(R^{-1}\dot R)$, normalize, and never re-orthogonalize.

---

## §7 Honest audit

- **Nothing gets cheaper for the ODE.** The shooting/relaxation numerics of
  Chapters 4–6 live in ordinary `Float` coordinates; GA is a reading layer.
  The one place it might earn runtime later: rotor-based frame propagation
  (evolve 4 components instead of re-Gram-Schmidt-ing 3 arrows).
- **Lean support.** Core Lean has nothing GA; Mathlib has exterior algebra
  (the grade-2 skeleton). A Cl₃ structure is 8 `Float`s, a product table,
  and the sandwich/reversion conventions — very buildable (see exercise).
- **Sign discipline.** $I_3^2 = -1$ and rotor conventions hide minus signs
  that a coordinate treatment exposes. State conventions once and derive
  everything from the reflection identity — otherwise it's the same
  "where did that $\Gamma$ minus come from" game with less vocabulary.

---

## §8 Exercise: build $\mathrm{Cl}_3$ in core Lean

A multivector is $(a_0, a_1e_1, a_2e_2, a_3e_3, a_{12}e_{12},
a_{13}e_{13}, a_{23}e_{23}, a_{123}I_3)$ — 8 `Float`s:

```lean
structure Cl3 where
  s : Float; e1 : Float; e2 : Float; e3 : Float
  e12 : Float; e13 : Float; e23 : Float
  i : Float   -- e₁e₂e₃ coefficient (I₃)
```

Geometric product = reindexing with the $\mathrm{Cl}_3$ law
(sign $-$ per swapped basis pair, $I_3^2 = -1$). Then check the chapter's
headline claims numerically:

1. **Sandwich split.** For $w$ with $|n|=1$: check in a `#eval` that
   $\tfrac12(w - nwn)$ is tangent (zero dot with $n$) and equals
   $w - (w\cdot n)\,n$, while $\tfrac12(w + nwn)$ equals $(w\cdot n)\,n$
   — both projections of Chapter 2 come from one product.
2. **Wedge test.** Great circles on the sphere, `γ(s) = cos s·U + sin s·V`
   with orthonormal `U, V`: here $\gamma'' = -\gamma$, so `γ″ ∧ n = 0` holds;
   a latitude circle at $\theta \ne \pi/2$ fails it: `γ″` has a tangent part
   (Chapter 2 §2.6's rows, in one `#eval`).
3. **Rotor frame.** Evolve $R(s) = \exp(\Omega s/2)$ for the sphere
   parallel curve with $\Omega$'s known Darboux values; reconstruct
   $(T, h, n)$ by sandwiching the reference frame and check orthonormality
   numerically.

That's the whole algebra Chapters 1–3 ask for, with zero dependencies.
