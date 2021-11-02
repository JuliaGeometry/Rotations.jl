using StaticArrays
using Rotations
using LinearAlgebra
using Test

import Rotations: RotationError, rotation_error, add_error, params

@testset "Rotation Error" begin
q1 = rand(UnitQuaternion)
dq = expm(0.1 * normalize(@SVector rand(3)))
q2 = q1*dq
@test q1\q2 ≈ dq
dg = RodriguesParam(dq)

emap = CayleyMap()
e1 = rotation_error(q2, q1, emap)
@test e1 ≈ q2 ⊖ q1
@test e1 ≈ params(dg)
@test add_error(q1, e1) ≈ q2
@test q1 ⊕ e1 isa UnitQuaternion
@test rotation_angle(UnitQuaternion(e1)) ≈ 0.1

aa1 = AngleAxis(q1)
@test add_error(aa1, e1) ≈ q2

e1 = rotation_error(q2, aa1, emap)
@test add_error(aa1, e1) ≈ AngleAxis(q2)
@test aa1 ⊕ e1 ≈ q2

# Test IdentityMap
g1 = rand(RodriguesParam)
g2 = RodriguesParam(g1*dg)
@test rotation_angle(g1\g2) ≈ 0.1

e1 = rotation_error(g2, g1, IdentityMap())
@test e1 ≈ g2 ⊖ g1
@test e1 ≈ params(dg)

p1 = rand(MRP)
dp = MRP(dg)
p2 = p1*dp
@test rotation_angle(p1\p2) ≈ 0.1

e1 = rotation_error(p2, p1, IdentityMap())
@test e1 ≈ p2 ⊖ p1
@test e1 ≈ params(dp)

@test_throws ArgumentError rotation_error(q2, q1, IdentityMap())
@test_throws ArgumentError q1 ⊕ e1

# Test inv
emap = CayleyMap()
imap = Rotations.InvCayleyMap()
g = inv(emap)(q1)
@test g ≈ imap(q1)
@test emap(g) ≈ q1
@test inv(imap)(g) ≈ q1
@test inv(inv(emap))(g) ≈ q1
@test inv(inv(imap))(q1) ≈ g

g = RodriguesParam(g)
@test UnitQuaternion(g) ≈ q1
@test imap(g) ≈ params(g)
end
