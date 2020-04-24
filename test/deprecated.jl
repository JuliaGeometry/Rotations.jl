
@testset "Deprecations" begin
    # RodriguesVec
    @test Base.isdeprecated(Rotations, :RodriguesVec)
    @test Base.isdeprecated(Rotations, :Quat)
    @test Base.isdeprecated(Rotations, :SPQuat)
end
