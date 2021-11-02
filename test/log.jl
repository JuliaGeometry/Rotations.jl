@testset "log" begin
    all_types = (RotMatrix{3}, AngleAxis, RotationVec,
                 UnitQuaternion, RodriguesParam, MRP,
                 RotXYZ, RotYZX, RotZXY, RotXZY, RotYXZ, RotZYX,
                 RotXYX, RotYZY, RotZXZ, RotXZX, RotYXY, RotZYZ,
                 RotX, RotY, RotZ,
                 RotXY, RotYZ, RotZX, RotXZ, RotYX, RotZY,
                 RotMatrix{2}, Angle2d)

    @testset "$(T)" for T in all_types, F in (one, rand)
        R = F(T)
        @test R ≈ exp(log(R))
        @test log(R) isa SMatrix
    end
end
