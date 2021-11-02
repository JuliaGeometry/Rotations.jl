using ForwardDiff
import Rotations: ∇rotate, ∇composition1, ∇composition2, skew, params


@testset "$R basic tests" for R in (RodriguesParam, MRP)

    # Constructors
    @test R(1.0, 0.0, 0.0) isa R{Float64}
    @test R(1.0, 0, 0) isa R{Float64}
    @test R(1.0f0, 0f0, 0f0) isa R{Float32}
    @test R(1.0f0, 0, 0) isa R{Float32}
    @test R(1.0, 0f0, 0) isa R{Float64}
    @test R(1, 0, 0) isa R{Int}
    @test R{Float64}(1, 0, 0) isa R{Float64}
    @test R{Float64}(1f0, 0f0, 0f0) isa R{Float64}
    @test R{Float32}(1.0, 0, 0) isa R{Float32}

    # Copy constructors
    g = rand(R)
    @test R{Float32}(g) isa R{Float32}
    @test R{Float64}(rand(R{Float32})) isa R{Float64}

    # initializers
    @test rand(R) isa R{Float64}
    @test rand(R{Float32}) isa R{Float32}
    @test one(R) isa R{Float64}
    @test one(R{Float32}) isa R{Float32}
    @test params(one(R)) === @SVector [0,0,0.]


    # Math operations
    g = rand(R)
    @test norm(g) == sqrt(g.x^2 + g.y^2 + g.z^2)

    # Test Jacobians
    R = RodriguesParam
    g1 = rand(R)
    g2 = rand(R)
    r = @SVector rand(3)
    @test ForwardDiff.jacobian(g->R(g)*r, params(g1)) ≈ ∇rotate(g1, r)

    function compose(g2,g1)
        params(R(g2)*R(g1))
    end
    @test ForwardDiff.jacobian(g->compose(params(g2),g), params(g1)) ≈ ∇composition1(g2,g1)
    @test ForwardDiff.jacobian(g->compose(g,params(g1)), params(g2)) ≈ ∇composition2(g2,g1)

    g0 = R{Float64}(0,0,0)
    @test ∇composition1(g2, g0) ≈ Rotations.∇differential(g2)

    gval = params(g1)
    b = @SVector rand(3)
    @test ForwardDiff.jacobian(g->∇composition1(g2,R(g))'b, gval) ≈
        Rotations.∇²composition1(g2,g1,b)
    @test Rotations.∇²differential(g2, b) ≈
        Rotations.∇²composition1(g2, g0, b)

end

@testset "kinematics" begin
    # UnitQuaternion
    q = rand(UnitQuaternion)
    ω = @SVector rand(3)
    q_ = Rotations.params(q)
    qdot = Rotations.kinematics(q,ω)
    @test qdot ≈ 0.5*lmult(q)*hmat()*ω
    @test qdot ≈ 0.5*lmult(q)*hmat(ω)
    @test ω ≈ 2*vmat()*lmult(q)'qdot
    @test ω ≈ 2*vmat()*lmult(inv(q))*qdot
    @test qdot ≈ 0.5 * params(q*pure_quaternion(ω))

    # MRPs
    ω = @SVector rand(3)
    g = rand(MRP)
    p = Rotations.params(g)
    A = Diagonal(I,3) + 2*(skew(p)^2 + skew(p))/(1+p'p)
    @test Rotations.kinematics(g, ω) ≈ 0.25*(1 + p'p) * A*ω
    @test ω ≈ 4/(1+p'p) * A'Rotations.kinematics(g,ω)

    # RPs
    g = rand(RodriguesParam)
    p = Rotations.params(g)
    gdot = Rotations.kinematics(g, ω)
    @test gdot ≈ 0.5*(Diagonal(I,3) + skew(p) + p*p')*ω
    @test ω ≈ 2/(1+p'p)*(gdot - p × gdot)
end
