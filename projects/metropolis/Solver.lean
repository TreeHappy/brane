import Std

-- A model's log-density (higher = more likely). Unnormalized is fine.
abbrev Model := Array Float → Float

structure StepResult where
  params   : Array Float
  gen      : StdGen
  accepted : Bool

abbrev StepFn := Model → Array Float → StdGen → StepResult

structure Solver where
  name : String
  step : StepFn

-- Generic runner: burn-in + thinning-free collection. Solver-agnostic.
def runChain (solver : Solver) (model : Model) (init : Array Float)
    (nBurn nSamples : Nat) (gen : StdGen) : Array (Array Float) × StdGen :=
  let rec loop (k : Nat) (θ : Array Float) (g : StdGen) (acc : Array (Array Float)) :
      Array (Array Float) × StdGen :=
    if k = 0 then (acc, g)
    else
      let r := solver.step model θ g
      let acc := if k ≤ nSamples then acc.push r.params else acc
      loop (k - 1) r.params r.gen acc
    termination_by k
  loop (nBurn + nSamples) init gen #[]

def runChainDetailed (solver : Solver) (model : Model) (init : Array Float)
    (nBurn nSamples : Nat) (gen : StdGen) :
    Array (Array Float) × Nat × Nat × StdGen :=
  let rec loop (k : Nat) (θ : Array Float) (g : StdGen) (acc : Array (Array Float)) (accA : Nat) :
      Array (Array Float) × Nat × StdGen :=
    if k = 0 then (acc, accA, g)
    else
      let r := solver.step model θ g
      let acc := if k ≤ nSamples then acc.push r.params else acc
      let accA := if r.accepted then accA + 1 else accA
      loop (k - 1) r.params r.gen acc accA
    termination_by k
  let (s, a, g) := loop (nBurn + nSamples) init gen #[] 0
  (s, a, nBurn + nSamples, g)

def randFloat (g : StdGen) : Float × StdGen :=
  let (n, g) := randNat g 0 1000000
  (Float.ofNat n / 1000000.0, g)

def randNormal (g : StdGen) : Float × StdGen :=
  let (u1, g) := randFloat g
  let (u2, g) := randFloat g
  let u1 := if u1 == 0.0 then 0.000001 else u1
  let r := Float.sqrt (-2.0 * Float.log u1)
  let θ := 2.0 * 3.141592653589793 * u2
  (r * Float.cos θ, g)

def metropolis (σ : Float) : Solver := {
  name := "metropolis"
  step := fun model θ gen =>
    let n := θ.size
    if n == 0 then { params := θ, gen, accepted := false }
    else
      let (iNat, gen) := randNat gen 0 (n - 1)
      let (eps, gen) := randNormal gen
      let θ' := θ.set! iNat (θ[iNat]! + eps * σ)
      let logp  := model θ
      let logp' := model θ'
      let (u, gen) := randFloat gen
      let u := if u == 0.0 then 0.000001 else u
      let accept := Float.log u < (logp' - logp)
      if accept then { params := θ', gen, accepted := true }
      else            { params := θ,  gen, accepted := false }
}
