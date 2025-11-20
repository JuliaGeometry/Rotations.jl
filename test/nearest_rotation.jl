@testset "nearest_rotation" begin
    @testset "$(N)-dim" for N in [2,3,4]
        for _ in 1:100
            M = randn(N,N)
            R = nearest_rotation(M)
            V = M/R
            U = R\M

            # R is rotation matrix
            @test isrotation(R)
            # U, V are (approximately) symmetric matrices
            # See [polar decomposition](https://en.wikipedia.org/wiki/Polar_decomposition)
            @test U ≈ U'
            @test V ≈ V'
            # Eigen values are non-negative, except for the first value.
            # We just need sort for Julia 1.0
            @test all(sort(eigvals(U))[2:end] .≥ 0)
            @test all(sort(eigvals(V))[2:end] .≥ 0)
        end

        for _ in 1:100
            M = @SMatrix randn(N,N)
            R = nearest_rotation(M)
            V = M/R
            U = R\M

            @test isrotation(R)
            @test U ≈ U'
            @test V ≈ V'
            # Can't calculate eigen values of general SMatrix.
            @test all(sort(eigvals(Symmetric(U)))[2:end] .≥ 0)
            @test all(sort(eigvals(Symmetric(V)))[2:end] .≥ 0)
        end
    end
    @testset "return eltype" begin
        for T in (Float16, Float32, Float64)
            M = @SMatrix randn(T, 3, 3)
            R = nearest_rotation(M)
            @test eltype(M) == eltype(R)
        end
    end
end
