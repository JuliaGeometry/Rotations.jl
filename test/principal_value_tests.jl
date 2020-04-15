@testset "Principal Value (SPQuat)" begin
    for i = 1:1000
        spq = SPQuat(5.0 * randn(), 5.0 * randn(), 5.0 * randn())
        spq_prin = principal_value(spq)
        @test spq_prin ≈ spq
        @test (spq_prin.x^2 + spq_prin.y^2 + spq_prin.z^2) < (1.0 + 1.0e-13)
        # end
    end
end

@testset "Principal Value (Quat)" begin
    for i = 1:1000
        qq = rand(Quat)
        qq_prin = principal_value(qq)
        @test 0.0 <= qq_prin.w
        @test qq_prin ≈ qq
    end
end

@testset "Principal Value (Angle Axis)" begin
    for i = 1:1000
        aa = AngleAxis(100.0 * randn(), randn(), randn(), randn())
        aa_prin = principal_value(aa)
        @test 0.0 <= aa_prin.theta
        @test aa_prin ≈ aa
    end
end

@testset "Principal Value (Rotation Vector)" begin
    for i = 1:1000
        rv = RotationVec(100.0 * randn(), 100.0 * randn(), 100.0 * randn())
        rv_prin = principal_value(rv)
        @test rotation_angle(rv_prin) <= pi
        @test rv_prin ≈ rv
    end
    rv = RotationVec(0.0, 0.0, 0.0)
    rv_prin = principal_value(rv)
    @test rotation_angle(rv_prin) <= pi
    @test rv_prin ≈ rv
end

@testset "Principal Value ($(rot_type))" for rot_type in [:RotX, :RotY, :RotZ] begin
        @eval begin
            for i = 1:1000
                r = $(rot_type)(100.0 * randn())
                r_prin = principal_value(r)
                @test -pi <= r_prin.theta <= pi
                @test r_prin ≈ r
            end
        end
    end
end

@testset "Principal Value ($(rot_type))" for rot_type in [:RotXY, :RotYX, :RotZX, :RotXZ, :RotYZ, :RotZY] begin
        @eval begin
            for i = 1:1000
                r = $rot_type(100.0 * randn(), 100.0 * randn())
                r_prin = principal_value(r)
                @test -pi <= r_prin.theta1 <= pi
                @test -pi <= r_prin.theta2 <= pi
                @test r_prin ≈ r
            end
        end
    end
end

@testset "Principal Value ($(rot_type))" for rot_type in [:RotXYX, :RotYXY, :RotZXZ, :RotXZX, :RotYZY, :RotZYZ, :RotXYZ, :RotYXZ, :RotZXY, :RotXZY, :RotYZX, :RotZYX] begin
        @eval begin
            for i = 1:1000
                r = $(rot_type)(100.0 * randn(), 100.0 * randn(), 100.0 * randn())
                r_prin = principal_value(r)
                @test -pi <= r_prin.theta1 <= pi
                @test -pi <= r_prin.theta2 <= pi
                @test -pi <= r_prin.theta3 <= pi
                @test r_prin ≈ r
            end
        end
    end
end
