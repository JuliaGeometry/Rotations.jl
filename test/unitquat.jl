using ForwardDiff

import Rotations: jacobian, ∇rotate, ∇composition1, ∇composition2
import Rotations: kinematics, pure_quaternion, params
import Rotations: vmat, rmult, lmult, hmat, tmat

@testset "Unit Quaternions" begin
    q1 = rand(QuatRotation)
    q2 = rand(QuatRotation)
    r = @SVector rand(3)
    ω = @SVector rand(3)

    # Constructors
    @test QuatRotation(1.0, 0.0, 0.0, 0.0) isa QuatRotation{Float64}
    @test QuatRotation(1.0, 0, 0, 0) isa QuatRotation{Float64}
    @test QuatRotation(1.0, 0, 0, 0) isa QuatRotation{Float64}
    @test QuatRotation(1, 0, 0, 0) isa QuatRotation{Float64}
    @test QuatRotation(1.0f0, 0, 0, 0) isa QuatRotation{Float32}

    # indexing
    qi = QuatRotation(1, 0, 0, 0)
    @test qi[1] == 1.0
    @test qi[6] == 0.0
    @test qi[5] == 1.0
    @test_throws BoundsError qi[10]

    q = normalize(@SVector rand(4))
    q32 = SVector{4,Float32}(q)
    @test QuatRotation(q) isa QuatRotation{Float64}
    @test QuatRotation(q32) isa QuatRotation{Float32}
    r = QuatRotation(q)
    @test r.w == r.q.s
    @test r.x == r.q.v1
    @test r.y == r.q.v2
    @test r.z == r.q.v3

    r = @SVector rand(3)
    r32 = SVector{3,Float32}(r)
    @test pure_quaternion(r) isa Quaternion{Float64}
    @test pure_quaternion(r32) isa Quaternion{Float32}
    @test real(pure_quaternion(r)) == 0

    @test QuatRotation{Float64}(1, 0, 0, 0) isa QuatRotation{Float64}
    @test QuatRotation{Float32}(1, 0, 0, 0) isa QuatRotation{Float32}
    @test QuatRotation{Float32}(1.0, 0, 0, 0) isa QuatRotation{Float32}
    @test QuatRotation{Float64}(1f0, 0, 0, 0) isa QuatRotation{Float64}

    # normalization
    @test QuatRotation(2.0, 0, 0, 0, true) == one(QuatRotation)
    @test QuatRotation(2q, true) ≈ QuatRotation(q)

    # Copy constructors
    q = rand(QuatRotation)
    @test QuatRotation(q) === q
    @test QuatRotation{Float32}(q) isa QuatRotation{Float32}
    QuatRotation{Float32}(q)

    # rand
    @test rand(QuatRotation) isa QuatRotation{Float64}
    @test rand(QuatRotation{Float32}) isa QuatRotation{Float32}

    # Test math
    @test QuatRotation(I) isa QuatRotation{Float64}

    ϕ = inv(ExponentialMap())(q1)
    @test expm(ϕ * 2) ≈ q1
    q = Rotations.pure_quaternion(ϕ)
    @test QuatRotation(exp(q)) ≈ q1
    @test exp(q) isa Quaternion

    q = QuatRotation((@SVector [1, 2, 3, 4.0]), false)
    @test 2 * q == 2 * Matrix(q)
    @test q * 2 == 2 * Matrix(q)

    @test_throws ErrorException QuatRotation(SVector(1, 2, 3, 4, 5))

    # Axis-angle
    ϕ = 0.1 * @SVector [1, 0, 0]
    q = expm(ϕ)
    @test q isa QuatRotation
    @test logm(expm(ϕ)) ≈ ϕ
    @test expm(logm(q1)) ≈ q1
    @test rotation_angle(q) ≈ 0.1
    @test rotation_axis(q) == [1, 0, 0]

    @test norm(q1 * ExponentialMap()(ϕ)) ≈ √3
    @test q1 ⊖ q2 isa StaticVector{3}
    @test (q1 * CayleyMap()(ϕ)) ⊖ q1 ≈ ϕ


    # Test inverses
    q3 = q2 * q1
    @test q2 \ q3 ≈ q1
    @test q3 / q1 ≈ q2
    @test inv(q1) * r ≈ q1 \ r
    @test r ≈ q3 \ (q2 * q1 * r)
    @test q3 ⊖ q2 ≈ inv(CayleyMap())(q1)

    q = q1
    rhat = Rotations.pure_quaternion(r)
    @test q * r ≈ vmat() * lmult(q) * rmult(q)' * vmat()'r
    @test q * r ≈ vmat() * lmult(q) * rmult(q)' * hmat(r)
    @test q * r ≈ vmat() * lmult(q) * lmult(rhat) * tmat() * params(q)
    @test q * r ≈ vmat() * rmult(q)' * rmult(rhat) * params(q)
    @test q * r ≈ hmat()' * rmult(q)' * rmult(rhat) * params(q)

    @test rmult(params(q)) == rmult(q)
    @test lmult(params(q)) == lmult(q)
    r_pure = pure_quaternion(r)
    @test hmat(r) == SVector(r_pure.s, r_pure.v1, r_pure.v2, r_pure.v3)

    # Test Jacobians
    @test ForwardDiff.jacobian(q -> QuatRotation(q, false) * r, params(q)) ≈
          ∇rotate(q, r)

    @test ForwardDiff.jacobian(q -> params(q2 * QuatRotation(q, false)), params(q1)) ≈
          ∇composition1(q2, q1)
    @test ForwardDiff.jacobian(q -> params(QuatRotation(q, false) * q1), params(q2)) ≈
          ∇composition2(q2, q1)

    b = @SVector rand(4)
    qval = params(q1)
    ForwardDiff.jacobian(q -> ∇composition1(q2, QuatRotation(q))'b,
        @SVector [1, 0, 0, 0.0,])
    diffcomp = ϕ -> params(q2 * CayleyMap()(ϕ))
    ∇diffcomp(ϕ) = ForwardDiff.jacobian(diffcomp, ϕ)
    @test ∇diffcomp(@SVector zeros(3)) ≈ Rotations.∇differential(q2)
    @test ForwardDiff.jacobian(ϕ -> ∇diffcomp(ϕ)'b, @SVector zeros(3)) ≈
          Rotations.∇²differential(q2, b)

    @test lmult(q) ≈ ∇composition1(q, q2)

    ϕ = @SVector zeros(3)
    @test Rotations.∇differential(q) ≈ lmult(q) * jacobian(QuatVecMap(), ϕ)
    @test Rotations.∇differential(q) ≈ lmult(q) * jacobian(ExponentialMap(), ϕ)
    @test Rotations.∇differential(q) ≈ lmult(q) * jacobian(CayleyMap(), ϕ)
    @test Rotations.∇differential(q) ≈ lmult(q) * jacobian(MRPMap(), ϕ)


    # Check ops with Float32
    ϕ = SA_F32[1, 2, 3]
    @test expm(SA_F32[1, 2, 3]) isa QuatRotation{Float32}

    q32 = rand(QuatRotation{Float32})
    @test Rotations._log_as_quat(q32) isa Quaternion{Float32}
    @test log(q32) isa RotationGenerator{3}
    @test eltype(logm(q32)) == Float32
    @test expm(logm(q32)) ≈ q32

    @test normalize(q32) isa SMatrix{3,3,Float32}

    ω = @SVector rand(3)
    ω32 = Float32.(ω)
    @test Rotations.kinematics(q, ω) isa SVector{4,Float64}
    @test Rotations.kinematics(q32, ω32) isa SVector{4,Float32}
    @test Rotations.kinematics(q32, ω) isa SVector{4,Float32}
    @test Rotations.kinematics(q32, [1, 2, 3]) isa SVector{4,Float32}

    @test eltype(lmult(q32)) == Float32
    @test eltype(lmult(q)) == Float64

    @test eltype(tmat()) == Float64
    @test eltype(tmat(Int)) == Int
    @test eltype(vmat(Float32)) == Float32
    @test eltype(hmat(Float32)) == Float32
end
