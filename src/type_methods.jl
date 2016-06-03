# create a union for all of the types
RotationTypes  = Union{RotTypeList...}

# build expression for tuples of the parameterization's data
function expr_gen(typed, fields)
    expr = :(())                            # build an expression for a tuple
    acc = typed ? typed_acc : untyped_acc   # accessor method
    append!(expr.args, [acc(field) for field in fields])
    return expr
end
untyped_acc(i) = isa(i, Int) ? :(X[$(i)]) : :(X.$(i))
typed_acc(i) = isa(i, Int) ? :(T(X[$(i)])) : :(T(X.$(i)))


#
# Define special import rules for Quaternions because we want to enable 3 vectors and 4 vectors to be promoted to Quaternions
#
@inline convert{T <: Quaternion, U <: Real}(::Type{T}, X::Vector{U}) =
        (length(X) == 3) ? Quaternion(zero(U), X[1], X[2], X[3], false) : Quaternion(X[1], X[2], X[3], X[4], false)  # allow 3 Vec -> Quat
@inline call{T <: Quaternion, U <: Real}(::Type{T}, X::Vector{U}) = convert(T, X)
@inline convert{T <: Quaternion, U <: Real}(::Type{T}, X::FixedVector{3, U}) = T(0.0, X[1], X[2], X[3], false)



#
# function to create a code block to convert to and from Vec and Vector rtypes
#
function add_methods(rot_type)

    #
    # Build expressions for building typed and untyped tuples using an Integer Index
    #
    untyped_indexed_expr = expr_gen(false, 1:numel(rot_type))
    typed_indexed_expr = expr_gen(true, 1:numel(rot_type))

    #
    # Build expressions for building typed and untyped tuples using a field
    #

    fields = fieldnames(rot_type)
    untyped_field_expr = (rot_type <: Mat) ? untyped_indexed_expr : expr_gen(false, fields[1:numel(rot_type)])
    typed_field_expr = (rot_type <: Mat) ? typed_indexed_expr : expr_gen(true, fields[1:numel(rot_type)])


    # the code block to return
    qb = quote end

    #
    # convert to / from a mutable vector
    #
    qn = quote

        # convert to a mutable vector
        @inline vec(X::$(rot_type)) = vcat($(untyped_field_expr.args...))
        @inline convert(::Type{Vector}, X::$(rot_type)) = vec(X)
        @inline convert{T <: Real}(::Type{Vector{T}}, X::$(rot_type)) = vcat($(typed_field_expr.args...))
        @inline call{T <: Vector}(::Type{T}, X::$(rot_type)) = convert(T, X)

        # convert from a mutable vector
        if ($(rot_type) != Quaternion)  # Quaternions has this defined already
            @inline convert{T <: $(rot_type), U <: Real}(::Type{T}, X::Vector{U}) = T($(untyped_indexed_expr.args...))
            @inline call{T <: $(rot_type), U <: Real}(::Type{T}, X::Vector{U}) = convert(T, X)
        end

        # convert to an immutable vector
        @inline convert(::Type{Vec}, X::$(rot_type)) = Vec($(untyped_field_expr.args...))
        @inline convert{T <: Real}(::Type{Vec{$(numel(rot_type)), T}}, X::$(rot_type)) = Vec($(typed_field_expr.args...))
        @inline call(::Type{Vec}, X::$(rot_type)) = convert(Vec, X)
        @inline call{T <: Real}(::Type{Vec{$(numel(rot_type)), T}}, X::$(rot_type)) = convert(Vec{$(numel(rot_type)), T}, X)

        # convert from an immutable vector
        @inline convert{T <: $(rot_type), U <: Real}(::Type{T}, X::FixedVector{$(numel(rot_type)), U}) = T($(untyped_indexed_expr.args...))
        @inline call{T <: $(rot_type), U <: Real}(::Type{T}, X::FixedVector{$(numel(rot_type)), U}) = convert(T, X)

        # convert to a tuple
        # @inline tuple(X::$(rot_type)) = tuple($(untyped_field_expr.args...))

        # convert from a tuple
        @inline convert{T <: $(rot_type), U <: Real}(::Type{T}, X::NTuple{$(numel(rot_type)), U}) = T($(untyped_indexed_expr.args...))
        @inline call{T <: $(rot_type), U <: Real}(::Type{T}, X::NTuple{$(numel(rot_type)), U}) = convert(T, X)


    end
    append!(qb.args, qn.args)

    #
    # NaN checking
    #
    nan_expr = Expr(:||, :(isnan($(untyped_field_expr.args[1]))),  :(isnan($(untyped_field_expr.args[2]))))
    for i = 3:numel(rot_type)
        nan_expr = Expr(:||, nan_expr, :(isnan($(untyped_field_expr.args[i]))))
    end
    push!(qb.args, :( isnan(X::$(rot_type)) = $(nan_expr)))

    #
    # Sneak in an extra constructor that casts Ints to Floats (who wants Ints for these thing anyways)
    #
    int_cons1 = [:( $(symbol("x$(i)"))::Integer ) for i in 1:numel(rot_type)]
    int_cons2 = [:( Float64($(symbol("x$(i)"))) ) for i in 1:numel(rot_type)]
    int_expr = :(
                    call(::Type{$(rot_type)}, $(int_cons1...)) = call($(rot_type), $(int_cons2...))
                )
    push!(qb.args, int_expr)
    if (rot_type in [EulerAngles, ProperEulerAngles])
        int_expr = :(
                    call{ORDER}(::Type{$(rot_type){ORDER}}, $(int_cons1...)) = call($(rot_type){ORDER}, $(int_cons2...))
                    )
        push!(qb.args, int_expr)
    end

    return qb
end

#
# Now go through and add everything
#
for rt in RotTypeList

    # println("type: $(rt), vars: $(numel(rt))")

    # add import / export
    eval(add_methods(rt))
end

