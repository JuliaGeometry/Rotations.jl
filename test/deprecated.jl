
@testset "Deprecations" begin
    # RodriguesVec
    @test_deprecated RodriguesVec(1,2,3)
    @test_deprecated RodriguesVec(rand(UnitQuaternion))
    @test_deprecated RodriguesVec(Tuple(rand(UnitQuaternion)))
    @test_deprecated RodriguesVec{Float32}(1,2,3)
    @test_deprecated one(RodriguesVec)
    @test_deprecated one(RodriguesVec{Float64})
end
