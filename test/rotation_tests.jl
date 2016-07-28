# function to perform tests of the rotation functions in the Rotations module
using Compat
using FixedSizeArrays
using Rotations



###############################
# Check fixed relationships
###############################

# Do no rotation
@testset "Identity rotation checks" begin
    I = eye(RotMatrix{Float64})
    @testset "$(rT)" for rT in rot_types
        rot_var = rT(I)
        null_var = eye(rT)
        @types_approx_eq(rot_var, null_var)
    end
end

# Constructor tests
@testset "Testing constructors" begin
    @testset "Testing integer construction (should convert to Float64)" begin
        @testset "$(rT)" for rT in setdiff(rot_types, [RotMatrix])
            rot_var1 = rT(zeros(Int, Rotations.numel(rT))...)
            rot_var2 = rT(zeros(Rotations.numel(rT))...)
            @types_approx_eq(rot_var1, rot_var2)
        end
    end

    @testset "Testing mixed element construction" begin
        @testset "$(rT)" for rT in rot_types
            rot_var1 = rT(vcat(zero(Float32), zeros(Rotations.numel(rT)-1))...)
            rot_var2 = rT(zeros(Rotations.numel(rT))...)
            @types_approx_eq(rot_var1, rot_var2)
        end
    end

    # test the axis vector constructor for angle axis
    Iaa = eye(AngleAxis)
    @testset "Testing Angle Axis vector axis construction" begin
        @test Iaa == AngleAxis(rotation_angle(Iaa), rotation_axis(Iaa))
        @test Iaa == AngleAxis(rotation_angle(Iaa), Vector(rotation_axis(Iaa)))
        @test_throws DimensionMismatch AngleAxis(0.0, [1.0, 0, 0, 0])
    end
end


################################
# check on the inverse function
################################

@testset "Testing inverse()" begin
    repeats, thresh = 100, 1e-6
    I = eye(RotMatrix{Float64})
    @testset "$(rT)" for rT in [AngleAxis, Quaternion, SpQuat, RodriguesVec]
        srand(0)
        for i = 1:repeats
            X = rT(nquatrand())
            Rout = RotMatrix(X) * RotMatrix(inv(X))
            rd = vecnorm(I - Rout)
            @test rd <= thresh
        end
    end
end


#########################################################################
# Rotate some stuff
#########################################################################

# a random rotation
@testset "Rotate Points" begin
    thresh, repeats = 1e-10, 1000
    @testset "$(rT)" for rT in rot_types
        srand(0)
        for i = 1:repeats
            q = i == 1 ? Quaternion(1.0, 0.0, 0.0, 0.0) : nquatrand()  # the identity is a bit special for most parametrizations
            Rm = RotMatrix(q)
            X = Vec{3, Float64}(randn(), randn(), randn())  # a point to rotate
            Xo = Rm * X  # Fixed size arrays better get this right

            R = rT(Rm)  # convert R to this formulation
            Xo_t = rotate(R, X)
            @types_approx_eq_eps(Xo_t, Xo, thresh)
        end

        # rotate a Vector as well
        Rm = RotMatrix(nquatrand())
        X = randn(Vec{3, Float64})  # a point to rotate
        Xo = Rm * X
        Xo_t = rotate(rT(Rm), [Float32(X[i]) for i in 1:3])
        @test reduce(&, [abs(Xo[i] - Xo_t[i]) <= 1e-6 for i in 1:3])
    end


    #
    # both Quaternions and rotation matrices should be able to do rotations with scaling
    #
    @testset "Rotations with scaling: $(rT)" for rT in [Quaternion, RotMatrix]
        qhat = nquatrand() # a unit quaternion
        s = .5 + rand()    # add a scaling
        qs = qhat * s      # non-unit quaternion (i.e. rotate and scale)
        X = randn(Vec{3, Float64})  # a point to rotate

        # rotate then scale as a base line
        Xo = s * rotate(qhat, X)

        # do the rotation and scaling together
        Xo_t = rotate(rT(qs), X)
        @test reduce(&, [abs(Xo[i] - Xo_t[i]) <= 1e-6 for i in 1:3])
    end
end


#########################################################################
# Test conversions between rotation types
#########################################################################

# test random round trip conversions
@testset "Rotation parameterization conversion checks" begin
    repeats, thresh = 100, 1e-12
    I = eye(RotMatrix{Float64})
    for rT_in in rot_types
        srand(0)

        # straight conversions
        @testset "$(rT_in) -> $(rT_out)" for rT_out in rot_types
            for i = 1:repeats                               # and each test

                # start with a random quaternion
                q = i == 1 ? Quaternion(1.0, 0.0, 0.0, 0.0) : nquatrand()  # the identity is a bit special for most parametrizations
                X = convert(rT_in, q)

                # round trip conversion
                Xd = convert(rT_in, convert(rT_out, X))

                # compare rotations before and after the round trip
                Rout = RotMatrix(X) * RotMatrix(Xd)'  # should be the identity
                rd = vecnorm(I - Rout)
                @test rd <= thresh
            end
        end

        # element type conversion
        if (length(rT_in.parameters) == 1) || (TypeVar != typeof(rT_in.parameters[end-1]))          # dont try when there'a missing order parameter
            @testset "$(rT_in) -> $(rT_out) with element type conversion" for rT_out in rot_types[!ordered_type]

                # start with a random quaternion
                q = Quaternion(1.0, 0.0, 0.0, 0.0) # the identity is a bit special for most parametrizations
                X = convert(rT_in{Float64}, q)

                # check type inference without any element type conversions
                @inferred convert(rT_out, X)
                @inferred convert(rT_in, convert(rT_out, X))

                # check type inference with element type conversions
                @inferred convert(rT_out{Float32}, X)
                @inferred convert(rT_in{Float64}, convert(rT_out{Float32}, X))
                Xd = convert(rT_in{Float64}, convert(rT_out{Float32}, X))

                # compare rotations before and after the round trip
                Rout = RotMatrix(X) * RotMatrix(Xd)'  # should be the identity
                rd = vecnorm(I - Rout)
                @test rd <= thresh
            end
        end
    end
end



#####################################################################################
# Check that Quaternions can import a 3 vector correctly
# (this is to make sure the "point" import functions from Quaternions.jl still work)
#####################################################################################

@testset "Quaternion 3 vector import checks" begin
    @types_approx_eq(Quaternion([1.0,2.0,3.0]), Quaternion(0.0, 1.0, 2.0, 3.0))
    @types_approx_eq(Quaternion(Vec(1.0,2.0,3.0)), Quaternion(0.0, 1.0, 2.0, 3.0))
end


#########################################################################
# Check angle and axis and inv work as expected
#########################################################################

@testset "Testing angle / axis extraction" begin
    theta_v, avec = linspace(-2*pi, 2*pi, 17), Vec(1.0, 0.0, 0.0)
    thresh = 1e-9
    @testset "$(rT)" for rT in [AngleAxis, RodriguesVec, Quaternion, SpQuat]
        for theta in theta_v
            aa = AngleAxis(theta, avec...)
            rot_var = rT(aa)
            theta_ex, avec_ex = rotation_angle(rot_var), rotation_axis(rot_var)

            # N.B. its OK of the axis is reversed and theta = 2*pi - theta
            dp = dot(avec_ex, avec)
            # @test_approx_eq_eps abs(dp)  1 thresh
            if (dp > 0)
               @test abs(Rotations.wrap_angle(theta_ex - theta)) <= thresh
            else
               @test abs(Rotations.wrap_angle((2*pi - theta_ex) - theta)) <= thresh
            end
        end
    end

    # make sure a slightly unnormalized Quaternion doesn't throw when extracting the angle
    q = Quaternion(1.0000000000000002, -5.040577330528428e-13, 5.663457441733299e-13, 2.452760217153127e-13)
    @test abs(rotation_angle(q)) < thresh

end

#########################################################################
# Check conversion to and from mutable and immutable vectors,
# as well as element conversions
#########################################################################

@testset "Vector import / export tests" begin
    eltypes = subtypes(AbstractFloat)  # test different AbstractFloats
    @testset "$(rT)" for rT in rot_types

        # create one
        rot_var = eye(rT)

        # export to mutable
        @inferred vec(rot_var)
        mvu = vec(rot_var)
        @contents_approx_eq(rot_var, mvu)

        # import from mutable
        if (rT != RotMatrix)
            @inferred rT(mvu)   # Julia v0.4 (v0.5 can) can't infer the type here for RotMatrix and I don't know why.  whatevs
        end
        rot_mvu = rT(mvu)
        @types_approx_eq(rot_var, rot_mvu)

        # export to immutable
        @inferred Vec(rot_var) # Julia v0.4 again
        ivu = Vec(rot_var)
        @contents_approx_eq(rot_var, ivu)

        # import from immutable
        if (rT != RotMatrix)
            @inferred rT(ivu)  # Julia v0.4 (v0.5 can) can't infer the type here for RotMatrix and I don't know why.  whatevs
        end
        rot_ivu = rT(ivu)
        @types_approx_eq(rot_var, rot_ivu)


        # test typed stuff
        order_param = (length(rT.parameters) > 1) && (TypeVar == typeof(rT.parameters[end-1]))  # does it have a missing order template parameter?
        @testset "$(order_param ? rT{rT.parameters[1], eT} : rT{eT})" for eT in [Float32]#eltypes

            # export to mutable
            mvt = Vector{eT}(rot_var)
            mvc = convert(Vector{eT}, mvu)
            @types_approx_eq(mvt, mvc)

            # import from mutable
            rot_mvt = rT(mvt)
            @contents_approx_eq_notype(rot_var, rot_mvt)

            # export to immutable
            ivt = Vec{numel(rT), eT}(rot_var)
            ivc = convert(Vec{numel(rT), eT}, ivt)
            @types_approx_eq(ivt, ivc)

            # import from immutable
            rot_ivt = rT(ivt)
            @contents_approx_eq_notype(rot_var, rot_ivt)

            # and test the element conversion on the rotation parameterization directly
            if (!order_param)
                rot_c = convert(rT{eT}, rot_var)
                @contents_approx_eq(rot_c, rot_mvt)
            end
        end
    end
end

