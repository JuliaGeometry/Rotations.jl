@testset "Eigen_3D" begin
    all_types = (RotMatrix{3}, AngleAxis, RotationVec,
                 UnitQuaternion, RodriguesParam, MRP,
                 RotXYZ, RotYZX, RotZXY, RotXZY, RotYXZ, RotZYX,
                 RotXYX, RotYZY, RotZXZ, RotXZX, RotYXY, RotZYZ,
                 RotX, RotY, RotZ,
                 RotXY, RotYZ, RotZX, RotXZ, RotYX, RotZY)
    oneaxis_types = (RotX, RotY, RotZ)

    @testset "$(T)" for T in all_types, F in (one, rand)
        R = F(T)
        λs = eigvals(R)
        vs = eigvecs(R)
        E = eigen(R)
        v1 = vs[:,1]
        v2 = vs[:,2]
        v3 = vs[:,3]
        @test R * vs ≈ transpose(λs) .* vs
        @test norm(v1) ≈ 1
        @test norm(v2) ≈ 1
        @test norm(v3) ≈ 1
        @test E.values == λs
        @test E.vectors == vs
        if !(T in oneaxis_types) && VERSION ≥ v"1.2"
            # If the rotation angle is in [0°, 180°], then the eigvals will be equal.
            # Note that the randomized RotX (and etc.) have rotation angle in [0°, 360°].
            # This needs Julia(≥1.2) to get sorted eigenvalues in a canonical order
            # See https://github.com/JuliaLang/julia/pull/21598
            @test eigvals(R) ≈ eigvals(collect(R))
        end
    end
end

@testset "Eigen_2D" begin
    all_types = (RotMatrix{2}, Angle2d)

    @testset "$(T)" for T in all_types, θ in 0.0:0.1:π
        R = T(Angle2d(θ))
        λs = eigvals(R)
        vs = eigvecs(R)
        E = eigen(R)
        v1 = vs[:,1]
        v2 = vs[:,2]
        @test R * vs ≈ transpose(λs) .* vs
        @test norm(v1) ≈ 1
        @test norm(v2) ≈ 1
        @test E.values == λs
        @test E.vectors == vs
        if VERSION ≥ v"1.2"
            # If the rotation angle θ is in [0°, 180°], then the eigvals will be equal.
            # Note that the randomized RotX (and etc.) have rotation angle in [0°, 360°].
            # This needs Julia(≥1.2) to get sorted eigenvalues in a canonical order
            # See https://github.com/JuliaLang/julia/pull/21598
            @test eigvals(R) ≈ eigvals(collect(R))
        end
    end
end
