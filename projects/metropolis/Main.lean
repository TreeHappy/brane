import Solver

def data : Array Float := #[
  2.8, 3.4, 2.9, 3.1, 3.6, 2.7, 3.2, 3.0, 3.5, 2.6,
  3.3, 3.1, 2.8, 3.4, 3.0, 3.2, 2.9, 3.6, 3.1, 2.7,
  3.4, 3.0, 3.3, 2.8, 3.2, 3.5, 2.9, 3.1, 3.0, 3.4,
  2.7, 3.3, 3.2, 2.9, 3.1, 3.0, 3.4, 2.8, 3.2, 3.1,
  3.0, 3.3, 2.9, 3.2, 3.1, 3.0, 2.8, 3.4, 3.2, 2.9
]

def model : Model := fun θ =>
  let μ := θ[0]!
  let prior := -(μ * μ) / (2.0 * 25.0)
  let ll := data.foldl (fun acc x => acc - 0.5 * (x - μ) * (x - μ)) 0.0
  prior + ll

def main : IO Unit := do
  let solver := metropolis 0.3
  let init := #[0.0]
  let gen := StdGen.mk 12345 67890
  let (samples, nAcc, nTot, _) := runChainDetailed solver model init 2000 5000 gen
  let (sum, n) := samples.foldl (fun (s, c) θ => (s + θ[0]!, c + 1)) (0.0, 0)
  let mean := if n == 0 then 0.0 else sum / Float.ofNat n
  -- posterior std
  let sq := samples.foldl (fun acc θ => acc + (θ[0]! - mean) * (θ[0]! - mean)) 0.0
  let std := Float.sqrt (sq / Float.ofNat (n - 1))
  IO.println s!"true μ = 3.0"
  IO.println s!"posterior mean = {mean}"
  IO.println s!"posterior std  = {std}"
  IO.println s!"acceptance rate = {Float.ofNat nAcc / Float.ofNat nTot}"
  IO.println s!"num samples = {n}"
