# define maths function for rotations types
import Base: (*), .*, (-), (.-), (+), (.+), (/), (./)

const enable_mult = false   # enable RotType * Vec to rotate a point.  This may not be real maths

function add_maths(rot_type)

    #
    # Define the operations to build
    #

    # unary expression
    unary_ops = [:(-)]

    # binary expressions below
    scalar_ops = [:(*), :(+), :(-), :(/)]

    vector_ops = [:(.*), :(.+), :(.-), :(./)]

    type_ops = [:(.+), :(.-), :(.*), :(./)]  # return a Vec

    # create an expression for a tuple X[1], X[2], etc...
    fields = fieldnames(rot_type)

    # the code block we are creating
    qb = quote; end

    #
    # allow rotating a Vec{3, ...} using *
    # N.B. this is now legacy
    if (enable_mult)
        qn = quote
            (*){T <: Real}(R::$(rot_type), X::Vec{3,T}) = rotate(R, X)
        end
        append!(qb.args, qn.args)
    end


    #
    # Build unary ops
    #
    for op in unary_ops

        # build the rhs expression
        exprs = [:($(op)(X.$(fields[i]))) for i in 1:numel(rot_type)]

        # form the whole expression
        qn = quote
            ($(op)){T <: $(rot_type)}(X::T) = T($(exprs...))
        end
        append!(qb.args, qn.args)

    end

    #
    # Build operations with scalars
    #
    for op in scalar_ops


        # build the rhs epxression
        exprs1 = [Expr(:call, :($(op)), :(s), :(X.$(field))) for field in fields[1:numel(rot_type)]]
        exprs2 = [Expr(:call, :($(op)), :(X.$(field)), :(s)) for field in fields[1:numel(rot_type)]]

        # form the whole expression
        qn = quote
            ($(op)){T <: $(rot_type), U <: Real}(s::U, X::T) = strip_eltype(T)($(exprs1...))
            ($(op)){T <: $(rot_type), U <: Real}(X::T, s::U) = strip_eltype(T)($(exprs2...))
        end
        append!(qb.args, qn.args)

    end



    #
    # Build operations with vectors
    #
    for op in vector_ops

        # build the rhs epxression
        exprs1 = [Expr(:call, :($(op)), :(s[$(i)]), :(X.$(fields[i]))) for i in 1:numel(rot_type)]
        exprs2 = [Expr(:call, :($(op)), :(X.$(fields[i])), :(s[$(i)])) for i in 1:numel(rot_type)]

        # form the whole expression
        qn = quote
            ($(op)){T <: $(rot_type), U <: Real}(s::Vec{$(numel(rot_type)), U}, X::T) = strip_eltype(T)($(exprs1...))
            ($(op)){T <: $(rot_type), U <: Real}(X::T, s::Vec{$(numel(rot_type)), U}) = strip_eltype(T)($(exprs2...))
        end
        append!(qb.args, qn.args)

    end

    #
    # Build operations on themselves
    #
    for op in type_ops

        # build the rhs epxression
        exprs = [Expr(:call, :($(op)), :(X1.$(field)), :(X2.$(field))) for field in fields[1:numel(rot_type)]]

        # form the whole expression
        qn = quote
            function ($(op)){T1 <: $(rot_type), T2 <: $(rot_type)}(X1::T1, X2::T2)
                @assert strip_eltype(T1)==strip_eltype(T2) "Different euler angle orders"
                Vec($(exprs...))
            end
        end
        append!(qb.args, qn.args)

    end


    return qb
end


for rt in RotTypeList

    # println("type: $(rt), vars: $(numel(rt))")

    # add some basic maths
    if !(rt <: RotMatrix) && (!(rt <: Quaternion))
        eval(add_maths(rt))
    elseif (rt <: Quaternion) && (enable_mult)
        # allow rotation via * with quaternions
        eval(:((*){T <: Real}(R::$(rt), X::Vec{3,T}) = rotate(R, X)))
    end
end

