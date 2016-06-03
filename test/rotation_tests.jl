# function to perform tests of the rotation functions in the Rotations module
using FixedSizeArrays
using Quaternions
using Rotations


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
        for i = 1:n
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
        for i = 1:n
            ai, bi = getindex($(esc(a)), i), getindex($(esc(b)), i)
            @test typeof(ai) == typeof(bi)
            @test_approx_eq_eps ai bi $(esc(eps))
        end
    end
end


# a macro to test if tow types are approximatey equal
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
        for i = 1:n
            ai, bi = getindex($(esc(a)), i), getindex($(esc(b)), i)
            @test_approx_eq ai bi
        end
    end
end


#####################################################################################
# build a full list of rotation types including the different angle ordering schemas
#####################################################################################

rot_types = Vector{Any}(0)
for rT in Rotations.RotTypeList
    if (length(rT.parameters) == 2)

        # get the super type for the order parameter
        order_type = super(Rotations.euler_order(rT))
        for order in subtypes(order_type)
            push!(rot_types, rT{order})
        end
    end

    # ordered ones should have defaults so leave them in
    push!(rot_types, rT)
end





###############################
# Start testing
###############################

@testset "Rotations Tests" begin

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

    # convert Int -> Float on construction
    @testset "Testing integer construction (should convert to Float64)" begin
        @testset "$(rT)" for rT in setdiff(rot_types, [RotMatrix])
            rot_var1 = rT(zeros(Int, Rotations.numel(rT))...)
            rot_var2 = rT(zeros(Rotations.numel(rT))...)
            @types_approx_eq(rot_var1, rot_var2)
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
    @testset "Rotation Points" begin
        repeats = 1000
        @testset "$(rT)" for rT in rot_types
            srand(0)
            for i = 1:repeats
                Rm = RotMatrix(nquatrand())
                X = Vec{3, Float64}(randn(), randn(), randn())  # a point to rotate
                Xo = Rm * X  # Fixed size arrays better get this right

                R = rT(Rm)  # convert R to this formulation
                Xo_t = rotate(R, X)
                @types_approx_eq_eps(Xo_t, Xo, 1e-10)
            end
        end
    end


    #########################################################################
    # Test conversions between rotation types
    #########################################################################

    # test random round trip conversions
    @testset "Rotation parameterization conversion checks" begin
        repeats, thresh = 100, 1e-6
        I = eye(RotMatrix{Float64})
        for rT_in in rot_types
            srand(0)
            @testset "$(rT_in) -> $(rT_out)" for rT_out in rot_types
                for i = 1:repeats                               # and each test

                    # start with a random quaternion
                    q = nquatrand()
                    X = convert(rT_in, q)

                    # round trip conversion
                    Xd = convert(rT_in, convert(rT_out, X))

                    # compare rotations before and after the round trip
                    Rout = RotMatrix(X) * RotMatrix(Xd)'  # should be the identity
                    rd = vecnorm(I - Rout)
                    @test rd <= thresh
                end
            end
        end
    end


    #########################################################################
    # Check that Quaternions can import a 3 vector correctly
    #########################################################################

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
            mvu = vec(rot_var)
            @contents_approx_eq(rot_var, mvu)

            # import from mutable
            rot_mvu = rT(mvu)
            @types_approx_eq(rot_var, rot_mvu)

            # export to immutable
            ivu = Vec(rot_var)
            @contents_approx_eq(rot_var, ivu)

            # import from immutable
            rot_ivu = rT(ivu)
            @types_approx_eq(rot_var, rot_ivu)

            # test typed stuff
            order_param = (length(rT.parameters) > 1) && (TypeVar == typeof(rT.parameters[end-1]))  # doe it have a missing order template parameter?
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
end

