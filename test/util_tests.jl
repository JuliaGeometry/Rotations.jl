@testset "Util" begin
    @testset "Perpendicular vector" begin
        for i = 1 : 100
            vec = rand(SVector{3, Float64})
            perp = Rotations.perpendicular_vector(vec)
            @test norm(perp) >= maximum(abs.(vec))
            @test isapprox(dot(vec, perp), 0.; atol = 1e-10)
        end
    end

    @testset "Angle difference" begin
        for i = 1 : 100
            a = rand(-5π : 5π)
            b = rand(-5π : 5π)
            c = angle_difference(a, b)
            @test isapprox(cos(b + c), cos(a); atol = 1e-12)
            @test isapprox(sin(b + c), sin(a); atol = 1e-12)
            @test c >= -π
            @test c < π
        end
    end
end
