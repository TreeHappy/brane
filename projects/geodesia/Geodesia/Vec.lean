/-!
# Geodesia — vectors

Small `Float` vector layer so the rest of the library reads like differential
geometry, not like index manipulations.
-/

namespace Geodesia

/-- A 2-vector: a point of the parameter domain, or a tangent there. -/
structure Vec2 where
  x : Float
  y : Float

/-- A 3-vector: an ambient point, tangent, or normal. -/
structure Vec3 where
  x : Float
  y : Float
  z : Float

namespace Vec2

instance : ToString Vec2 := ⟨fun v => s!"({v.x}, {v.y})"⟩

def add (u v : Vec2) : Vec2 := ⟨u.x + v.x, u.y + v.y⟩
def sub (u v : Vec2) : Vec2 := ⟨u.x - v.x, u.y - v.y⟩
def smul (a : Float) (u : Vec2) : Vec2 := ⟨a * u.x, a * u.y⟩
def dot (u v : Vec2) : Float := u.x * v.x + u.y * v.y
def norm (u : Vec2) : Float := Float.sqrt (u.dot u)
def normalize (u : Vec2) : Vec2 := smul (1 / u.norm) u

/-- Explicit coordinate access as an indexed family (nice for Christoffel code). -/
def comp (v : Vec2) : Fin 2 → Float :=
  fun i => match i with
    | 0 => v.x
    | 1 => v.y

end Vec2

namespace Vec3

instance : ToString Vec3 := ⟨fun v => s!"({v.x}, {v.y}, {v.z})"⟩

def add (u v : Vec3) : Vec3 := ⟨u.x + v.x, u.y + v.y, u.z + v.z⟩
def sub (u v : Vec3) : Vec3 := ⟨u.x - v.x, u.y - v.y, u.z - v.z⟩
def smul (a : Float) (u : Vec3) : Vec3 := ⟨a * u.x, a * u.y, a * u.z⟩
def dot (u v : Vec3) : Float := u.x * v.x + u.y * v.y + u.z * v.z
def norm (u : Vec3) : Float := Float.sqrt (u.dot u)
def normalize (u : Vec3) : Vec3 := smul (1 / u.norm) u

/-- The cross product — the one nonplanar operation we need, for normals. -/
def cross (u v : Vec3) : Vec3 :=
  ⟨u.y * v.z - u.z * v.y, u.z * v.x - u.x * v.z, u.x * v.y - u.y * v.x⟩

def comp (v : Vec3) : Fin 3 → Float :=
  fun i => match i with
    | 0 => v.x
    | 1 => v.y
    | 2 => v.z

end Vec3

end Geodesia
