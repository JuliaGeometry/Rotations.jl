# function to perform tests of the rotation functions in the Rotations module

#=
########################
# Define helper methods
########################

import Rotations: numel
numel{T <: Rotations.RotationTypes}(X::T) = Rotations.numel(Rotations.strip_eltype(T))
numel{T}(X::T) = (T <: FixedSizeArrays.FixedArray) || (T <: FixedSizeArrays.AbstractArray)  ? length(X) : error("numel undefined for input of type $(T)")

# a macro to test if the contents of two containers are approximatley equal
macro contents_approx_eq(a, b)
    quote
        n = numel($(esc(a)))
        @test n == numel($(esc(b)))
        for i in 1:n
            ai, bi = getindex($(esc(a)), i), getindex($(esc(b)), i)
            @test typeof(ai) == typeof(bi)
            @test_approx_eq ai bi
        end
    end
end
macro contents_approx_eq_eps(a, b, eps)
    quote
        n = numel($(esc(a)))
        @test n == numel($(esc(b)))
        for i in 1:n
            ai, bi = getindex($(esc(a)), i), getindex($(esc(b)), i)
            @test typeof(ai) == typeof(bi)
            @test_approx_eq_eps ai bi $(esc(eps))
        end
    end
end


# a macro to test if two types are approximatey equal
macro types_approx_eq(a, b)
    quote
        @test typeof($(esc(a))) == typeof($(esc(b)))
        @contents_approx_eq($(esc(a)), $(esc(b)))
    end
end
macro types_approx_eq_eps(a, b, eps)
    quote
        @test typeof($(esc(a))) == typeof($(esc(b)))
        @contents_approx_eq_eps($(esc(a)), $(esc(b)), $(esc(eps)))
    end
end

# a macro to test if the contents of two containers are approximatley equal, without examining the element types
macro contents_approx_eq_notype(a, b)
    quote
        n = numel($(esc(a)))
        @test n == numel($(esc(b)))
        for i in 1:n
            ai, bi = getindex($(esc(a)), i), getindex($(esc(b)), i)
            @test_approx_eq ai bi
        end
    end
end
=#

#####################################################################################
# build a full list of rotation types including the different angle ordering schemas
#####################################################################################

rot_types = (RotMatrix3, RotMatrix{3}, AngleAxis, RotationVec,
             QuatRotation, RodriguesParam, MRP,
             RotXYZ, RotYZX, RotZXY, RotXZY, RotYXZ, RotZYX,
             RotXYX, RotYZY, RotZXZ, RotXZX, RotYXY, RotZYZ)

one_types = (RotX, RotY, RotZ)
two_types = (RotXY, RotYZ, RotZX, RotXZ, RotYX, RotZY)
taitbyran_types = (RotXYZ, RotYZX, RotZXY, RotXZY, RotYXZ, RotZYX)
all_types = (RotMatrix3, RotMatrix{3}, AngleAxis, RotationVec,
             QuatRotation, RodriguesParam, MRP,
             RotXYZ, RotYZX, RotZXY, RotXZY, RotYXZ, RotZYX,
             RotXYX, RotYZY, RotZXZ, RotXZX, RotYXY, RotZYZ,
             RotX, RotY, RotZ,
             RotXY, RotYZ, RotZX, RotXZ, RotYX, RotZY)

###############################
# Start testing
###############################

@testset "Rotations Tests" begin
    # Ensure we're testing all 3D rotation types
    @test setdiff(subtypes(Rotation), all_types) == [Angle2d, RotMatrix]
    @test setdiff(all_types, subtypes(Rotation)) == [RotMatrix3, RotMatrix{3}]

    ###############################
    # Check fixed relationships
    ###############################

    @testset "Identity rotation checks" begin
        I = one(SMatrix{3,3,Float64})
        I32 = one(SMatrix{3,3,Float32})
        @testset "$(R)" for R in all_types
            # one(R) should always return something of type R (#114)
            @test one(R)::R == I
            @test one(one(R))::R == I
            @test one(R{Float32})::R{Float32} == I32
            @test one(one(R{Float32}))::R{Float32} == I32
        end
    end

    ###############################
    # Check zero function
    ###############################

    @testset "zero checks" begin
        @testset "$(R)" for R in all_types
            # zero
            @test zero(R) == zero(R{Float64}) == zero(one(R))
            @test zero(R) isa SMatrix
            @test zero(R{Float64}) isa SMatrix
            @test zero(one(R)) isa SMatrix
            # zeros
            @test zeros(R)[1] == zeros(R,3)[1] == zeros(R,3,3)[1] == zeros(R,(3,3,3))[1] == zero(R)
            @test zeros(R) isa Array{<:SMatrix,0}
            @test zeros(R,3) isa Array{<:SMatrix,1}
            @test zeros(R,3,3) isa Array{<:SMatrix,2}
            @test zeros(R,(3,3,3)) isa Array{<:SMatrix,3}
        end

        @test_throws ErrorException zero(Rotation)
        @test_throws ErrorException zero(RotMatrix)
    end

    ###############################
    # Check real function
    ###############################

    @testset "real checks" begin
        @testset "$(R)" for R in all_types
            @test real(R) === R
            @test real(one(R)) === one(R)
        end
        @test real(Rotation) === Rotation
        @test real(RotMatrix) === RotMatrix
    end

    ################################
    # check on the inverse function
    ################################

    @testset "Testing inv()" begin
        repeats = 100
        I = one(RotMatrix{3,Float64})
        @testset "$(R)" for R in all_types
            Random.seed!(0)
            for i in 1:repeats
                r1 = rand(R)
                r2 = rand(R)
                v = @SVector rand(3)
                @test inv(r1) === adjoint(r1) === transpose(r1)
                @test inv(r1)*r1 ≈ I
                @test r1*inv(r1) ≈ I
                @test r1/r1 ≈ I
                @test r1\r1 ≈ I
                @test r1/r2 ≈ r1*inv(r2)
                @test r1\r2 ≈ inv(r1)*r2
                @test r1\v ≈ inv(r1)*v ≈ Matrix(r1)\Vector(v)
            end
        end
    end

    ################################
    # check on the norm functions
    ################################

    @testset "Testing norm() and normalize()" begin
        repeats = 100
        for R in all_types
            I = one(R)
            Random.seed!(0)
            for i in 1:repeats
                r = rand(R)
                @test norm(r) ≈ norm(Matrix(r))
                @test normalize(r) ≈ normalize(Matrix(r))
                @test normalize(r) isa SMatrix
            end
        end
    end

    #########################################################################
    # Rotate some stuff
    #########################################################################

    # a random rotation of a random point
    @testset "Rotate Points" begin
        repeats = 100
        @testset "$(R)" for R in all_types
            I = one(R)::R

            Random.seed!(0)
            for i in 1:repeats
                r = rand(R)
                m = SMatrix(r)
                v1 = randn(SVector{3})
                v2 = randn(3)

                @test I*v1 ≈ v1
                @test I*v2 ≈ v2
                @test r*v1 ≈ m*v1
                @test r*v2 ≈ m*v2
            end
        end
    end

    @testset "Quaternion double cover" begin
        repeats = 100
        Q = QuatRotation
        for i in 1 : repeats
            q = rand(QuatRotation)

            q2 = QuatRotation(-q.q) # normalize: need a tolerance
            @test Rotations.params(q2) ≈ -Rotations.params(q) atol = 100 * eps()
            @test q ≈ q2 atol = 100 * eps()

            q3 = QuatRotation(-real(q.q), (-).(imag_part(q.q))..., false) # don't normalize: everything is exact
            @test Rotations.params(q3) == -Rotations.params(q)
            @test q == q3

            Δq = q \ q3
            @test Δq ≈ one(QuatRotation) atol = 100 * eps()
        end
    end

    # compose two random rotations
    @testset "Compose rotations" begin
        repeats = 100
        @testset "$(R1) * $(R2)" for R1 in all_types, R2 in all_types
            Random.seed!(0)
            for i in 1:repeats
                r1 = rand(R1)
                m1 = SMatrix(r1)

                r2 = rand(R2)
                m2 = SMatrix(r2)

                @test r1*r2 ≈ m1*m2
            end
        end
    end

    @testset "Slerp rotations" begin
        repeats = 100
        @testset "slerp($(R1), $(R2), rand())" for R1 in all_types, R2 in all_types
            Random.seed!(0)
            for i in 1:repeats
                r1 = rand(R1)
                q1 = QuatRotation(r1)
                r2 = rand(R2)
                q2 = QuatRotation(r2)
                t = rand()
                @test slerp(r1, r2, t) ≈ slerp(q1, q2, t)
            end
        end
    end

    #########################################################################
    # Test conversions between rotation types
    #########################################################################
    @testset "Convert rotations" begin
        repeats = 100
        @testset "convert $(R1) -> $(R2)" for R1 in all_types, R2 in rot_types
            Random.seed!(0)
            for i in 1:repeats
                r1 = rand(R1)
                m1 = SMatrix(r1)

                r2 = R2(r1)

                @test r2 ≈ m1
            end
        end
    end

    @testset "Convert rotations 1-axis -> 2-axis" begin
        repeats = 100
        @testset "convert $(R1) -> $(R2)" for R1 in one_types, R2 in two_types
            # Check if the two-axis include one-axis
            if string(R1)[end] in string(R2)[end-1:end]
                Random.seed!(0)
                for i in 1:repeats
                    r1 = rand(R1)
                    m1 = SMatrix(r1)

                    r2 = R2(r1)

                    @test r2 ≈ m1
                end
            end
        end
    end

    #########################################################################
    # Test robustness of DCM to QuatRotation function
    #########################################################################
    @testset "DCM to QuatRotation" begin
        pert = 1e-3
        for type_rot in all_types
            for _ = 1:100
                not_orthonormal = rand(type_rot) + rand(3, 3) * pert
                quat_ill_cond = QuatRotation(not_orthonormal)
                @test 0 <= real(quat_ill_cond.q)
                @test norm(quat_ill_cond - nearest_rotation(not_orthonormal)) < 10 * pert
            end
        end
    end

    #########################################################################
    # Check angle and axis and inv work as expected
    #########################################################################

    @testset "Testing angle / axis extraction" begin
        repeats = 100
        @testset "$(R)" for R in rot_types
            Random.seed!(0)
            for i in 1:repeats
                r1 = rand(AngleAxis)

                angle = rotation_angle(r1)
                axis = rotation_axis(r1)

                r2 = R(r1)

                @test rotation_angle(r2) ≈ angle
                @test rotation_axis(r2) ≈ axis
            end
        end
        @test norm(rotation_axis(QuatRotation(1.0, 0.0, 0.0, 0.0))) ≈ 1.0


        # TODO RotX, RotXY?
    end

    @testset "RotationVec->QuatRotation, especially near zero" begin
        for i = 1:10
            p = rand(RotationVec)
            all(iszero, Rotations.params(p)) && continue
            @test QuatRotation(p) ≈ QuatRotation(AngleAxis(p))
        end

        @test rotation_angle(RotationVec(0.0, 0.0, 0.0)) ≈ 0.0
        @test rotation_axis(RotationVec(0.0, 0.0, 0.0)) ≈ [1.0, 0.0, 0.0]
        @test QuatRotation(RotationVec(0.0, 0.0, 0.0)) == QuatRotation(1.0, 0.0, 0.0, 0.0)
        qr = QuatRotation(RotationVec(2e-30, 0.0, 0.0))
        @test qr.w ≈ 1
        @test qr.x ≈ 1e-30
    end

    #########################################################################
    # Check construction of QuatRotation given two vectors
    #########################################################################

    @testset "Testing construction of QuatRotation given two vectors" begin
        angle_axis_test(from, to, rot, atol) = isapprox(rot * from * norm(to) / norm(from), to; atol = atol)

        for i in 1 : 10000
            from = randn(SVector{3, Float64})
            to = rand(SVector{3, Float64})
            rot = rotation_between(from, to)
            @test angle_axis_test(from, to, rot, 1e-10)
        end

        # degenerate cases
        for i in 1 : 10000
            from = randn(SVector{3, Float64})
            to = randn() * from # either from and to are aligned, or in opposite directions
            rot = rotation_between(from, to)
            # @show rot
            @test angle_axis_test(from, to, rot, 1e-7)
        end
        for direction in 1 : 3
            for i in 1 : 10000
                from = @SVector [ifelse(i == direction, 1., 0.) for i in 1 : 3] # unit vector in direction 'direction'
                to = randn() * from
                rot = rotation_between(from, to)
                @test angle_axis_test(from, to, rot, 1e-7)
            end
        end
        @test_throws ArgumentError rotation_between(zero(SVector{3}), rand(SVector{3}))
        @test_throws ArgumentError rotation_between(rand(SVector{3}), zero(SVector{3}))
    end

    #########################################################################
    # Check roll, pitch, yaw constructors
    #########################################################################

    @testset "Testing roll / pitch / yaw constructors" begin
        repeats = 100
        s = 1e-4
        @testset "$(R)" for R in taitbyran_types
            Random.seed!(0)
            for i in 1:repeats
                roll = s*(rand()-0.5)
                pitch = s*(rand()-0.5)
                yaw = s*(rand()-0.5)

                # This tests whether the rotations are the same to first order
                # in roll, pitch and yaw. Second-order terms are small enough
                # to pass the isapprox() test, but not first-order terms.
                r1 = RotXYZ(roll, pitch, yaw)
                r2 = R(roll=roll, pitch=pitch, yaw=yaw)

                @test r1 ≈ r2
            end
        end
    end

    @testset "$(N)-dimensional rotation_between" for N in 2:7
        @testset "random check" begin
            for _ in 1:100
                u = randn(SVector{N})
                v = randn(SVector{N})
                R = rotation_between(u,v)
                @test isrotation(R)
                @test R isa Rotation
                @test normalize(v) ≈ R * normalize(u)
            end
        end

        @testset "mixed types" begin
            u = randn(MVector{N})
            v = randn(SizedVector{N})
            R = rotation_between(u,v)
            @test isrotation(R)
            @test R isa Rotation
            @test normalize(v) ≈ R * normalize(u)
        end

        @testset "zero-vector" begin
            @test_throws ArgumentError rotation_between(zero(SVector{N}), rand(SVector{N}))
            @test_throws ArgumentError rotation_between(rand(SVector{N}), zero(SVector{N}))
        end
    end

#########################################################################
    # Check that the eltype is inferred in Rot constructors
    @testset "Rot constructor eltype promotion" begin
        @test eltype(RotX(10)) == Float64
        @test eltype(RotX(10.0f0)) == Float32
        @test eltype(RotX(BigInt(10))) == BigFloat

        @test eltype(RotXY(10, 20)) == Float64
        @test eltype(RotXY(10.0, 20.0)) == Float64
        @test eltype(RotXY(10.0f0, 20.0f0)) == Float32
        # Mixing ints + floats promotes to narrowest float type
        @test eltype(RotXY(10.0f0, 20)) == Float32
        @test eltype(RotXY(20.0, BigInt(10))) == BigFloat

        @test eltype(RotXYZ(10, 20, 30)) == Float64
        @test eltype(RotXYZ(10.0, 20.0, 30.0)) == Float64
        @test eltype(RotXYZ(10.0f0, 20.0f0, 30.0f0)) == Float32
        @test eltype(RotXYZ(10.0f0, 20, 30)) == Float32

        # Promotion is correct with dimensionless Unitful types
        ° = Unitful.°
        rad = Unitful.rad
        @test eltype(RotX(10°)) == Float64
        @test eltype(RotX(10.0f0°)) == Float32
        @test eltype(RotX(10rad)) == Float64
        @test eltype(RotX(BigInt(10)*rad)) == BigFloat

        @test RotX(10°) ≈ RotX(deg2rad(10.0))
        @test RotXY(10°,20°) ≈ RotXY(deg2rad(10.0), deg2rad(20.0))
        @test RotXYZ(10°,20°,30°) ≈ RotXYZ(deg2rad(10.0), deg2rad(20.0), deg2rad(30.0))
    end

    #########################################################################
    # Check that isrotation works
    #########################################################################
    @testset "Testing isrotation" begin
        # Rotate 90° around x-axis
        a=[1.0 0.0 0.0
           0.0 0.0 -1.0
           0.0 1.0 0.0]
        @test isrotation(a)

        # Random generated Rotation must be rotation
        foreach(rot_types) do rot_type
            foreach(1:20) do idx
            @test isrotation(rand(rot_type))
            end
        end

        # Scaling in x-axis is not rotation
        a=[4.0 0.0 0.0
           0.0 1.0 0.0
           0.0 0.0 1.0]
        @test !isrotation(a)

        # Flipping y-axis and z-axis is not rotation
        a=[1.0 0.0 0.0
           0.0 0.0 1.0
           0.0 1.0 0.0]
        @test !isrotation(a)

        # Non-square matrix is not rotation
        @test !isrotation(zeros(2,3))
        @test !isrotation(@SMatrix zeros(2,3))

        # isrotation should work for integer (or boolean) matrices (issue #94)
        @test isrotation([0 1 0; -1 0 0; 0 0 1])

        # identity matrix is rotation
        @test isrotation(I(1))
        @test isrotation(I(3))
        @test isrotation(I(4))
        @test isrotation(Matrix(I,1,1))
        @test isrotation(Matrix(I,3,3))
        @test isrotation(Matrix(I,4,4))
        @test isrotation(one(SMatrix{1,1}))
        @test isrotation(one(SMatrix{3,3}))
        @test isrotation(one(SMatrix{4,4}))
        @test isrotation(one(RotMatrix{1}))
        @test isrotation(one(RotMatrix{3}))
        @test isrotation(one(RotMatrix{4}))

        # Rotation matrix can be AbstractMatrix{<:Complex}
        @test isrotation(I(3) .+ 0im)
        @test isrotation(one(SMatrix{4,4}) .+ 0im)

        # Including NaNs are not rotaion
        @test !isrotation(MRP(0, 0, NaN))
        @test !isrotation(RotXYZ(0, NaN, NaN))
        @test !isrotation(RotX(NaN))
        @test !isrotation(AngleAxis(1.,0.,0.,0.))
        @test !isrotation(QuatRotation(0.,0.,0.,0.))

        # Complex unitary matrix is not rotation
        M = Hermitian(randn(Complex{Float64},3,3))
        A = eigvecs(M)
        @test A*A' ≈ one(A)
        @test !isrotation(A)
    end

    @testset "Testing type aliases" begin
        @test one(RotMatrix{2, Float64}) isa RotMatrix2{Float64}
        @test one(RotMatrix{3, Float64}) isa RotMatrix3{Float64}
    end

    @testset "Testing normalization" begin
        θ, x, y, z = 1., 2., 3., 4.
        aa = AngleAxis(θ, x, y, z)
        @test norm([aa.axis_x, aa.axis_y, aa.axis_z]) ≈ 1.
        aa = AngleAxis(θ, x, y, z, false)
        @test norm([aa.axis_x, aa.axis_y, aa.axis_z]) ≈ norm([x, y, z])

        aa = AngleAxis(0.0, 0.0, 0.0, 0.0, false)
        @test convert(typeof(aa), aa) === aa

        w, x, y, z = 1., 2., 3., 4.
        quat = QuatRotation(w, x, y, z)
        @test norm(Rotations.params(quat)) ≈ 1.
        quat = QuatRotation(w, x, y, z, false)
        @test norm(Rotations.params(quat)) ≈ norm([w, x, y, z])

        w, x, y, z = 1., 2., 3., 4.
        quat = QuatRotation(w, x, y, z)
        @test norm(Rotations.params(quat)) ≈ 1.
        quat = QuatRotation(w, x, y, z, false)
        @test norm(Rotations.params(quat)) ≈ norm([w, x, y, z])
    end

    @testset "Testing RotMatrix conversion to Tuple" begin
        rot = one(RotMatrix{3, Float64})
        @inferred Tuple(rot)
    end

    @testset "Testing show" begin
        io = IOBuffer()
        r = rand(RotMatrix{2})
        show(io, MIME("text/plain"), r)
        str = String(take!(io))
        @test startswith(str, "2×2 RotMatrix2{Float64}")

        rxyz = RotXYZ(1.0, 2.0, 3.0)
        show(io, MIME("text/plain"), rxyz)
        str = String(take!(io))
        @test startswith(str, "3×3 RotXYZ{Float64}") && occursin("(1.0, 2.0, 3.0):", str)
    end

    #########################################################################
    # Check that 137 is solved
    #########################################################################
    @testset "Regression test issue 137" begin
        @test det(RotationVec(1e19, 0.0, 0.0)) ≈ 1.
    end

    @testset "params" begin
        p1, p2, p3, p4 = randn(4)
        @test Rotations.params(RotX(p1)) == [p1]
        @test Rotations.params(RotXY(p1,p2)) == [p1,p2]
        @test Rotations.params(RotXYZ(p1,p2,p3)) == [p1,p2,p3]
        @test Rotations.params(AngleAxis(p1,p2,p3,p4)) ≈ pushfirst!(normalize([p2,p3,p4]),p1)
        @test Rotations.params(RotationVec(p1,p2,p3)) == [p1,p2,p3]
        @test Rotations.params(QuatRotation(p1,p2,p3,p4)) ≈ normalize([p1,p2,p3,p4])
        @test Rotations.params(MRP(p1,p2,p3)) == [p1,p2,p3]
        @test Rotations.params(RodriguesParam(p1,p2,p3)) == [p1,p2,p3]
    end
end
