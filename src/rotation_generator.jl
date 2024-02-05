"""
    abstract type RotationGenerator{N,T} <: StaticMatrix{N,N,T}

An abstract type representing `N`-dimensional rotation generator. More abstractly, they represent
skew-symmetric real `N`×`N` matrices.
"""
abstract type RotationGenerator{N,T} <: StaticMatrix{N,N,T} end

Base.adjoint(r::RotationGenerator) = -r
Base.transpose(r::RotationGenerator{N,T}) where {N,T<:Real} = -r
Base.real(R::Type{<:RotationGenerator}) = R

# Generate identity-matrix with SMatrix
# Note that ones(RotationGenerator3,dims...) is not Array{<:RotationGenerator} but Array{<:StaticMatrix{3,3}}
Base.one(::RotationGenerator{N,T}) where {N,T} = one(SMatrix{N,N,T})
Base.one(::Type{RotationGenerator}) = error("The dimension of rotation is not specified.")
Base.one(::Type{<:RotationGenerator{N}}) where N = one(SMatrix{N,N})
Base.one(::Type{<:RotationGenerator{N,T}}) where {N,T} = one(SMatrix{N,N,T})
Base.ones(::Type{R}) where {R<:RotationGenerator} = ones(R, ()) # avoid StaticArray constructor
Base.ones(::Type{R}, dims::Base.DimOrInd...) where {R<:RotationGenerator} = ones(typeof(one(R)),dims...)
Base.ones(::Type{R}, dims::NTuple{N, Integer}) where {R<:RotationGenerator, N} = ones(typeof(one(R)),dims)
Base.ones(::Type{R}, dims::Tuple{}) where {R<:RotationGenerator} = ones(typeof(one(R)),dims)

# Generate zero rotation matrix
Base.zero(r::RotationGenerator) = zero(typeof(r))

# difference
Base.:-(r1::RotationGenerator, r2::RotationGenerator) = r1 + (-r2)

# `convert` goes through the constructors, similar to e.g. `Number`
Base.convert(::Type{R}, rot::RotationGenerator{N}) where {N,R<:RotationGenerator{N}} = R(rot)

# RotationGenerator matrices should be skew-symmetric. Only the operations we define,
# like addition, will stay as RotationGenerators, otherwise users will get an
# SMatrix{3,3} (e.g. rot1 * rot2 -> SMatrix)
Base.@pure StaticArrays.similar_type(::Union{R,Type{R}}) where {R <: RotationGenerator} = SMatrix{size(R)..., eltype(R), prod(size(R))}
Base.@pure StaticArrays.similar_type(::Union{R,Type{R}}, ::Type{T}) where {R <: RotationGenerator, T} = SMatrix{size(R)..., T, prod(size(R))}

################################################################################
################################################################################
"""
    struct RotMatrixGenerator{N,T} <: RotationGenerator{N,T}

A statically-sized, N×N skew-symmetric matrix.

Note: the skew-symmetricity of the input matrix is *not* checked by the constructor.
"""
struct RotMatrixGenerator{N,T,L} <: RotationGenerator{N,T} # which is <: AbstractMatrix{T}
    mat::SMatrix{N, N, T, L} # The final parameter to SMatrix is the "length" of the matrix, 3 × 3 = 9
    RotMatrixGenerator{N,T,L}(x::AbstractArray) where {N,T,L} = new{N,T,L}(convert(SMatrix{N,N,T,L}, x))
    # fixes #49 ambiguity introduced in StaticArrays 0.6.5
    RotMatrixGenerator{N,T,L}(x::StaticArray) where {N,T,L} = new{N,T,L}(convert(SMatrix{N,N,T,L}, x))
end
RotMatrixGenerator(x::SMatrix{N,N,T,L}) where {N,T,L} = RotMatrixGenerator{N,T,L}(x)

Base.zero(::Type{RotMatrixGenerator}) = error("The dimension of rotation is not specified.")

# These functions (plus size) are enough to satisfy the entire StaticArrays interface:
for N in 2:3
    L = N*N
    RotMatrixGeneratorN = Symbol(:RotMatrixGenerator, N)
    @eval begin
        @inline RotMatrixGenerator(t::NTuple{$L})  = RotMatrixGenerator(SMatrix{$N,$N}(t))
        @inline (::Type{RotMatrixGenerator{$N}})(t::NTuple{$L}) = RotMatrixGenerator(SMatrix{$N,$N}(t))
        @inline RotMatrixGenerator{$N,T}(t::NTuple{$L}) where {T} = RotMatrixGenerator(SMatrix{$N,$N,T}(t))
        @inline RotMatrixGenerator{$N,T,$L}(t::NTuple{$L}) where {T} = RotMatrixGenerator(SMatrix{$N,$N,T}(t))
        const $RotMatrixGeneratorN{T} = RotMatrixGenerator{$N, T, $L}
    end
end
Base.@propagate_inbounds Base.getindex(r::RotMatrixGenerator, i::Int) = r.mat[i]
@inline Base.Tuple(r::RotMatrixGenerator) = Tuple(r.mat)

@inline RotMatrixGenerator(v::Number) = RotMatrixGenerator{2}(v)
@inline function (::Type{RotMatrixGenerator{2}})(v::Number)
    RotMatrixGenerator(@SMatrix [zero(v) -v; v zero(v)])
end
@inline function RotMatrixGenerator{2,T}(v::Number) where T
    RotMatrixGenerator(@SMatrix T[zero(v) -v; v zero(v)])
end

# A rotation is more-or-less defined as being an orthogonal (or unitary) matrix
Base.:-(r::RotMatrixGenerator) = RotMatrixGenerator(r.mat')

# By default, composition of rotations will go through RotMatrixGenerator, unless overridden
@inline Base.:+(r1::RotationGenerator, r2::RotationGenerator) = RotMatrixGenerator(r1) + RotMatrixGenerator(r2)
@inline Base.:+(r1::RotMatrixGenerator, r2::RotationGenerator) = r1 + RotMatrixGenerator(r2)
@inline Base.:+(r1::RotationGenerator, r2::RotMatrixGenerator) = RotMatrixGenerator(r1) + r2
@inline Base.:+(r1::RotMatrixGenerator, r2::RotMatrixGenerator) = RotMatrixGenerator(r1.mat + r2.mat)

@inline Base.:*(t::Number, r::RotMatrixGenerator) = RotMatrixGenerator(t*r.mat)
@inline Base.:*(r::RotMatrixGenerator, t::Number) = t*r
@inline Base.:/(r::RotMatrixGenerator, t::Number) = RotMatrixGenerator(r.mat/t)

# Special case multiplication of 2×2 rotation matrices: speedup using skew-symmetricity
@inline function Base.:+(r1::RotMatrixGenerator{2}, r2::RotMatrixGenerator{2})
    v = r1[2,1]+r2[2,1]
    s = @SMatrix [0  -v
                  v   0]
    return RotMatrixGenerator(s)
end

# Special case multiplication of 3×3 rotation matrices: speedup using skew-symmetricity
@inline function Base.:+(r1::RotMatrixGenerator{3}, r2::RotMatrixGenerator{3})
    x = r1[6]+r2[6]
    y = r1[7]+r2[7]
    z = r1[2]+r2[2]
    s = @SMatrix [ 0 -z  y
                   z  0 -x
                  -y  x  0]

    return RotMatrixGenerator(s)
end

"""
    struct Angle2dGenerator{T} <: RotationGenerator{2,T}
        v::T
    end

A 2×2 rotation generator matrix (i.e. skew-symmetric matrix).
[ 0 -v
  v  0 ]
"""
struct Angle2dGenerator{T} <: RotationGenerator{2,T}
    v::T
    Angle2dGenerator{T}(r::Number) where T = new{T}(r)
end

@inline function Angle2dGenerator(r::T) where T <: Number
    Angle2dGenerator{T}(r)
end

@inline Angle2dGenerator(r::RotationGenerator{2}) = Angle2dGenerator(r[2,1])
@inline Angle2dGenerator{T}(r::RotationGenerator{2}) where {T} = Angle2dGenerator{T}(r[2,1])

@inline Base.zero(::Type{A}) where {A<: Angle2dGenerator} = A(0.0)

@inline Base.:+(r1::Angle2dGenerator, r2::Angle2dGenerator) = Angle2dGenerator(r1.v + r2.v)
@inline Base.:-(r1::Angle2dGenerator, r2::Angle2dGenerator) = Angle2dGenerator(r1.v - r2.v)
@inline Base.:*(t::Number, r::Angle2dGenerator) = Angle2dGenerator(t*r.v)
@inline Base.:*(r::Angle2dGenerator, t::Number) = t*r
@inline Base.:/(r::Angle2dGenerator, t::Number) = Angle2dGenerator(r.v/t)
@inline Base.:-(r::Angle2dGenerator) = Angle2dGenerator(-r.v)

@inline function Base.getindex(r::Angle2dGenerator{T}, i::Int) where T
    if i == 1
        zero(T)
    elseif i == 2
        r.v
    elseif i == 3
        -r.v
    elseif i == 4
        zero(T)
    else
        throw(BoundsError(r,i))
    end
end

@inline params(r::Angle2dGenerator) = SVector{1}(r.v)

"""
    struct RotationVecGenerator{T} <: RotationGenerator{2,T}
        x::T
        y::T
        z::T
    end

A 3×3 rotation generator matrix (i.e. skew-symmetric matrix).
[ 0 -z  y
  z  0 -x
 -y  x  0]
"""
struct RotationVecGenerator{T} <: RotationGenerator{3,T}
    x::T
    y::T
    z::T
    RotationVecGenerator{T}(x,y,z) where T = new{T}(x,y,z)
end

@inline function RotationVecGenerator(x::X,y::Y,z::Z) where {X<:Number, Y<:Number, Z<:Number}
    RotationVecGenerator{promote_type(X,Y,Z)}(x, y, z)
end

@inline RotationVecGenerator(r::RotationGenerator{3}) = RotationVecGenerator(r[6], r[7], r[2])
@inline RotationVecGenerator{T}(r::RotationGenerator{3}) where {T} = RotationVecGenerator{T}(r[6], r[7], r[2])

@inline Base.zero(::Type{R}) where {R<: RotationVecGenerator} = R(0.0,0.0,0.0)

@inline Base.:+(r1::RotationVecGenerator, r2::RotationVecGenerator) = RotationVecGenerator(r1.x + r2.x, r1.y + r2.y, r1.z + r2.z)
@inline Base.:-(r1::RotationVecGenerator, r2::RotationVecGenerator) = RotationVecGenerator(r1.x - r2.x, r1.y - r2.y, r1.z - r2.z)
@inline Base.:*(t::Number, r::RotationVecGenerator) = RotationVecGenerator(t*r.x, t*r.y, t*r.z)
@inline Base.:*(r::RotationVecGenerator, t::Number) = t*r
@inline Base.:/(r::RotationVecGenerator, t::Number) = RotationVecGenerator(r.x/t, r.y/t, r.z/t)
@inline Base.:-(r::RotationVecGenerator) = RotationVecGenerator(-r.x, -r.y, -r.z)

@inline function Base.getindex(r::RotationVecGenerator{T}, i::Int) where T
    if i == 1
        zero(T)
    elseif i == 2
        r.z
    elseif i == 3
        -r.y
    elseif i == 4
        -r.z
    elseif i == 5
        zero(T)
    elseif i == 6
        r.x
    elseif i == 7
        r.y
    elseif i == 8
        -r.x
    elseif i == 9
        zero(T)
    else
        throw(BoundsError(r,i))
    end
end

@inline params(r::RotationVecGenerator) = SVector{3}(r.x,r.y,r.z)

################################################################################
################################################################################


# A simplification and specialization of the Base.show function for AbstractArray makes
# everything sensible at the REPL.
function Base.show(io::IO, ::MIME"text/plain", X::RotationGenerator)
    if !haskey(io, :compact)
        io = IOContext(io, :compact => true)
    end
    summary(io, X)
    if !isa(X, RotMatrixGenerator)
        n_fields = length(fieldnames(typeof(X)))
        print(io, "(")
        for i in 1:n_fields
            print(io, getfield(X, i))
            if i < n_fields
                print(io, ", ")
            end
        end
        print(io, ")")
    end
    print(io, ":")
    println(io)
    io = IOContext(io, :typeinfo => eltype(X))
    Base.print_array(io, X)
end

"""
    isrotationgenerator(r)

Check whether `r` is a rotation generator matrix (skew-symmetric matrix).
"""
isrotationgenerator

function isrotationgenerator(r::StaticMatrix{N,N,T}) where {N,T<:Real}
    # This method can be removed, but is a little faster than isrotationgenerator(r::AbstractMatrix{T}).
    m = SMatrix(r)
    return m == -m'
end

function isrotationgenerator(r::AbstractMatrix{T}) where T<:Real
    indsm, indsn = axes(r)
    if indsm != indsn
        return false
    end
    for i in indsn, j in i:last(indsn)
        if r[i,j] != -r[j,i]
            return false
        end
    end
    return true
end

function isrotationgenerator(r::AbstractMatrix{T}) where T<:Number
    if !isreal(r)
        return false
    end
    return isrotationgenerator(real(r))
end
