
@testset "Deprecations" begin
    # RodriguesVec
    @test Base.isdeprecated(Rotations, :RodriguesVec)
end
