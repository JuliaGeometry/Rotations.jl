using ForwardDiff
import Rotations: jacobian, ErrorMap
import Rotations: CayleyMap, ExponentialMap, MRPMap, IdentityMap, QuatVecMap

@testset "Quaternion Maps" begin
    ϕ = @SVector rand(3)
    v = 0.1 * @SVector rand(3)
    g = @SVector rand(3)
    p = params(rand(MRP))

    @testset "Forward Maps" begin
        check_forward_jacobian(emap::ErrorMap, ϕ) =
            ForwardDiff.jacobian(x -> params(emap(x)), ϕ) ≈ jacobian(emap, ϕ)

        # Exponential
        @test check_forward_jacobian(ExponentialMap(), ϕ)

        ϕ = 1e-6 * @SVector rand(3)
        @test check_forward_jacobian(ExponentialMap(), ϕ)

        # Vector Part
        @test check_forward_jacobian(QuatVecMap(), v)

        # Gibbs Vectors
        @test check_forward_jacobian(CayleyMap(), g)

        # MRPs
        @test check_forward_jacobian(MRPMap(), p)

        μ0 = 1 / Rotations.scaling(QuatVecMap)
        jac_eye = [@SMatrix zeros(1, 3); μ0 * Diagonal(@SVector ones(3))]
        @test jacobian(ExponentialMap(), p * 1e-10) ≈ jac_eye
        @test jacobian(MRPMap(), p * 1e-10) ≈ jac_eye
        @test jacobian(CayleyMap(), p * 1e-10) ≈ jac_eye
        @test jacobian(QuatVecMap(), p * 1e-10) ≈ jac_eye
    end



    ############################################################################################
    #                                 INVERSE RETRACTION MAPS
    ############################################################################################
    @testset "Inverse Maps" begin

        # Exponential Map
        Random.seed!(1)
        function test_inverse_map(emap, invmap, q, ϕ)
            μ = Rotations.scaling(emap)
            @test inv(emap)(q) ≈ μ * invmap(params(q))
            @test ForwardDiff.jacobian(invmap, params(q)) * μ ≈ jacobian(inv(emap), q)
            @test emap(inv(emap)(q)) ≈ q
            @test inv(emap)(emap(ϕ)) ≈ ϕ
            b = @SVector rand(3)
            @test ForwardDiff.jacobian(
                q -> jacobian(inv(emap), UnitQuaternion(q, false))'b,
                params(q),
            ) ≈ Rotations.∇jacobian(inv(emap), q, b)
        end

        q = rand(UnitQuaternion)

        function invmap(q)
            v = @SVector [q[2], q[3], q[4]]
            s = q[1]
            θ = norm(v)
            M = 2 * atan(θ, s) / θ
            return M * v
        end
        test_inverse_map(ExponentialMap(), invmap, q, ϕ)

        qI = ExponentialMap()(v * 1e-5)
        @test ForwardDiff.jacobian(invmap, params(q)) * Rotations.scaling(ExponentialMap) ≈
              jacobian(inv(ExponentialMap()), q)
        @test ForwardDiff.jacobian(invmap, params(qI)) *
              Rotations.scaling(ExponentialMap) ≈ jacobian(inv(ExponentialMap()), qI)

        # Vector Part
        invmap(q) = sign(q[1]) * @SVector [q[2], q[3], q[4]]
        test_inverse_map(QuatVecMap(), invmap, q, v)

        # Cayley
        invmap(q) = 1 / q[1] * @SVector [q[2], q[3], q[4]]
        test_inverse_map(CayleyMap(), invmap, q, g)

        # MRP
        invmap(q) = 1 / (1 + q[1]) * @SVector [q[2], q[3], q[4]]
        test_inverse_map(MRPMap(), invmap, q, p)

        # Test near origin
        μ0 = Rotations.scaling(CayleyMap)
        jacT_eye = [@SMatrix zeros(1, 3); μ0 * Diagonal(@SVector ones(3))]'
        @test isapprox(jacobian(inv(ExponentialMap()), qI), jacT_eye, atol = 1e-5)
        @test isapprox(jacobian(inv(QuatVecMap()), qI), jacT_eye, atol = 1e-5)
        @test isapprox(jacobian(inv(CayleyMap()), qI), jacT_eye, atol = 1e-5)
        @test isapprox(jacobian(inv(MRPMap()), qI), jacT_eye, atol = 1e-5)

    end
end
