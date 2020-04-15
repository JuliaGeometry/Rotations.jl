using Test

@testset "Deprecations" begin
    # RodriguesVec
    @test_deprecated RodriguesVec(1,2,3)
    @test_deprecated RodriguesVec(rand(UnitQuaternion))
    @test_deprecated RodriguesVec(Tuple(rand(UnitQuaternion)))
    @test_deprecated RodriguesVec{Float32}(1,2,3)
    @test_deprecated one(RodriguesVec)
    @test_deprecated one(RodriguesVec{Float64})
    @test_deprecated rand(RodriguesVec)
    @test_deprecated rand(RodriguesVec{Float64})

    # Quat
    @test_deprecated Quat(1,2,3,4.)
    @test_deprecated Quat{Float64}(1,2,3,4)
    @test_deprecated Quat(rand(MRP))
    @test_deprecated one(Quat{Float32})
    @test_deprecated one(Quat)
end
