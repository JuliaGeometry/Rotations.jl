using ForwardDiff

import Rotations: jacobian, ∇rotate, ∇composition1, ∇composition2
import Rotations: kinematics, pure_quaternion, params
import Rotations: vmat, rmult, lmult, hmat, tmat

@testset "Unit Quaternions" begin
    q1 = rand(UnitQuaternion)
    q2 = rand(UnitQuaternion)
    r = @SVector rand(3)
    ω = @SVector rand(3)

    # Constructors
    @test UnitQuaternion(1.0, 0.0, 0.0, 0.0) isa UnitQuaternion{Float64}
    @test UnitQuaternion(1.0, 0, 0, 0) isa UnitQuaternion{Float64}
    @test UnitQuaternion(1.0, 0, 0, 0) isa UnitQuaternion{Float64}
    @test UnitQuaternion(1, 0, 0, 0) isa UnitQuaternion{Int}
    @test UnitQuaternion(1.0f0, 0, 0, 0) isa UnitQuaternion{Float32}

    q = normalize(@SVector rand(4))
    q32 = SVector{4,Float32}(q)
    @test UnitQuaternion(q) isa UnitQuaternion{Float64}
    @test UnitQuaternion(q32) isa UnitQuaternion{Float32}

    r = @SVector rand(3)
    r32 = SVector{3,Float32}(r)
    @test pure_quaternion(r) isa UnitQuaternion{Float64}
    @test pure_quaternion(r32) isa UnitQuaternion{Float32}
    @test pure_quaternion(r).w == 0

    @test UnitQuaternion{Float64}(1, 0, 0, 0) isa UnitQuaternion{Float64}
    @test UnitQuaternion{Float32}(1, 0, 0, 0) isa UnitQuaternion{Float32}
    @test UnitQuaternion{Float32}(1.0, 0, 0, 0) isa UnitQuaternion{Float32}
    @test UnitQuaternion{Float64}(1f0, 0, 0, 0) isa UnitQuaternion{Float64}

    # normalization
    @test UnitQuaternion(2.0, 0, 0, 0, true) == one(UnitQuaternion)
    @test UnitQuaternion(2q, true) ≈ UnitQuaternion(q)

    # Copy constructors
    q = rand(UnitQuaternion)
    @test UnitQuaternion(q) === q
    @test UnitQuaternion{Float32}(q) isa UnitQuaternion{Float32}
    UnitQuaternion{Float32}(q)

    # rand
    @test rand(UnitQuaternion) isa UnitQuaternion{Float64}
    @test rand(UnitQuaternion{Float32}) isa UnitQuaternion{Float32}

    # Test math
    @test UnitQuaternion(I) isa UnitQuaternion{Float64}

    ϕ = inv(ExponentialMap())(q1)
    @test expm(ϕ * 2) ≈ q1
    q = Rotations.pure_quaternion(ϕ)
    @test exp(q) ≈ q1

    q = UnitQuaternion((@SVector [1, 2, 3, 4.0]), false)
    @test 2 * q == UnitQuaternion((@SVector [2, 4, 6, 8.0]), false)
    @test q * 2 == UnitQuaternion((@SVector [2, 4, 6, 8.0]), false)

    # Axis-angle
    ϕ = 0.1 * @SVector [1, 0, 0]
    q = expm(ϕ)
    @test logm(expm(ϕ)) ≈ ϕ
    @test expm(logm(q1)) ≈ q1
    @test rotation_angle(q) ≈ 0.1
    @test rotation_axis(q) == [1, 0, 0]

    @test norm(q1 * ExponentialMap()(ϕ)) ≈ 1
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
    @test hmat(r) == params(pure_quaternion(r))

    # Test Jacobians
    @test ForwardDiff.jacobian(q -> UnitQuaternion(q, false) * r, params(q)) ≈
          ∇rotate(q, r)

    @test ForwardDiff.jacobian(q -> params(q2 * UnitQuaternion(q, false)), params(q1)) ≈
          ∇composition1(q2, q1)
    @test ForwardDiff.jacobian(q -> params(UnitQuaternion(q, false) * q1), params(q2)) ≈
          ∇composition2(q2, q1)

    b = @SVector rand(4)
    qval = params(q1)
    ForwardDiff.jacobian(q -> ∇composition1(q2, UnitQuaternion(q))'b,
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
    @test expm(SA_F32[1, 2, 3]) isa UnitQuaternion{Float32}

    q32 = rand(UnitQuaternion{Float32})
    @test Rotations._log_as_quat(q32) isa UnitQuaternion{Float32}
    @test log(q32) isa SMatrix
    @test eltype(logm(q32)) == Float32
    @test expm(logm(q32)) ≈ q32

    @test normalize(q32) isa UnitQuaternion{Float32}

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
