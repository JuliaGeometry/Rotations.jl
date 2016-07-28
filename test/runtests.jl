if VERSION >= v"0.5.0-dev+7720"
    using Base.Test
else
    using BaseTestNext
    const Test = BaseTestNext
end

using Compat
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
ordered_type = Vector{Bool}(0)
for rT in Rotations.RotTypeList
    if (length(rT.parameters) == 2)

        # get the super type for the order parameter
        @compat order_type = supertype(Rotations.euler_order(rT))
        for order in subtypes(order_type)
            push!(rot_types, rT{order})
            push!(ordered_type, false)
        end
    end

    # ordered ones should have defaults so leave them in
    push!(rot_types, rT)
    push!(ordered_type, (length(rT.parameters) == 2))
end



###############################
# Start testing
###############################

@testset "Rotations Tests" begin

    # include your own tests here
    include("rotation_tests.jl")
    include("derivative_tests.jl")

end

