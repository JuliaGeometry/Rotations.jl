"""
    abstract type Rotation{N,T} <: StaticMatrix{N,N,T}

An abstract type representing `N`-dimensional rotations. More abstractly, they represent
unitary (orthogonal) `N`×`N` matrices.
"""
abstract type Rotation{N,T} <: StaticMatrix{N,N,T} end

Base.@pure StaticArrays.Size(::Type{Rotation{N}}) where {N} = Size(N,N)
Base.@pure StaticArrays.Size(::Type{Rotation{N,T}}) where {N,T} = Size(N,N)
Base.@pure StaticArrays.Size(::Type{R}) where {R<:Rotation} = Size(supertype(R))
Base.adjoint(r::Rotation) = inv(r)
Base.transpose(r::Rotation{N,T}) where {N,T<:Real} = inv(r)
Base.real(R::Type{<:Rotation}) = R

# Generate zero-matrix with SMatrix
# Note that zeros(Rotation3,dims...) is not Array{<:Rotation} but Array{<:StaticMatrix{3,3}}
Base.zero(::Rotation{N,T}) where {N,T} = @SMatrix zeros(T,N,N)
Base.zero(::Type{Rotation}) = error("The dimension of rotation is not specified.")
Base.zero(::Type{<:Rotation{N}}) where N = @SMatrix zeros(N,N)
Base.zero(::Type{<:Rotation{N,T}}) where {N,T} = @SMatrix zeros(T,N,N)
Base.zeros(::Type{R}) where {R<:Rotation} = zeros(R, ()) # avoid StaticArray constructor
Base.zeros(::Type{R}, dims::Base.DimOrInd...) where {R<:Rotation} = zeros(typeof(zero(R)),dims...)
Base.zeros(::Type{R}, dims::NTuple{N, Integer}) where {R<:Rotation, N} = zeros(typeof(zero(R)),dims)
Base.zeros(::Type{R}, dims::Tuple{}) where {R<:Rotation} = zeros(typeof(zero(R)),dims) # avoid ambiguity

# Generate identity rotation matrix
Base.one(r::Rotation) = one(typeof(r))

# Rotation angles and axes can be obtained by converting to the AngleAxis type
rotation_angle(r::Rotation{3}) = rotation_angle(AngleAxis(r))
rotation_axis(r::Rotation{3}) = rotation_axis(AngleAxis(r))

# `convert` goes through the constructors, similar to e.g. `Number`
Base.convert(::Type{R}, rot::R) where {N,R<:Rotation{N}} = rot
Base.convert(::Type{R}, rot::Rotation{N}) where {N,R<:Rotation{N}} = R(rot)

# Rotation matrices should be orthoginal/unitary. Only the operations we define,
# like multiplication, will stay as Rotations, otherwise users will get an
# SMatrix{3,3} (e.g. rot1 + rot2 -> SMatrix)
Base.@pure StaticArrays.similar_type(::Union{R,Type{R}}) where {R <: Rotation} = SMatrix{size(R)..., eltype(R), prod(size(R))}
Base.@pure StaticArrays.similar_type(::Union{R,Type{R}}, ::Type{T}) where {R <: Rotation, T} = SMatrix{size(R)..., T, prod(size(R))}

@inline function Base.:/(r1::Rotation, r2::Rotation)
    r1 * inv(r2)
end

@inline function Base.:\(r1::Rotation, r2::Rotation)
    inv(r1) * r2
end

@inline function Base.:\(r::Rotation, v::AbstractVector)
    inv(r) * v
end

# This definition is for avoiding anbiguity
@inline function Base.:\(r::Rotation, v::StaticVector)
    inv(r) * v
end

################################################################################
################################################################################
"""
    struct RotMatrix{N,T} <: Rotation{N,T}

A statically-sized, N×N unitary (orthogonal) matrix.

Note: the orthonormality of the input matrix is *not* checked by the constructor.
"""
struct RotMatrix{N,T,L} <: Rotation{N,T} # which is <: AbstractMatrix{T}
    mat::SMatrix{N, N, T, L} # The final parameter to SMatrix is the "length" of the matrix, 3 × 3 = 9
    RotMatrix{N,T,L}(x::AbstractArray) where {N,T,L} = new{N,T,L}(convert(SMatrix{N,N,T,L}, x))
    # fixes #49 ambiguity introduced in StaticArrays 0.6.5
    RotMatrix{N,T,L}(x::StaticArray) where {N,T,L} = new{N,T,L}(convert(SMatrix{N,N,T,L}, x))
end
RotMatrix(x::SMatrix{N,N,T,L}) where {N,T,L} = RotMatrix{N,T,L}(x)

Base.zero(::Type{RotMatrix}) = error("The dimension of rotation is not specified.")

# These functions (plus size) are enough to satisfy the entire StaticArrays interface:
@inline function RotMatrix(t::Tuple)
    n = sqrt(length(t))
    if !isinteger(n)
        throw(DimensionMismatch("The length of input tuple $(t) must be square number."))
    end
    N = Int(n)
    RotMatrix(SMatrix{N,N}(t))
end
@inline (::Type{RotMatrix{N}})(t::Tuple) where N = RotMatrix(SArray{Tuple{N,N}}(t))
@inline RotMatrix{N,T}(t::Tuple) where {N,T} = RotMatrix(SArray{Tuple{N,N},T}(t))
@inline RotMatrix{N,T,L}(t::Tuple) where {N,T,L} = RotMatrix(SArray{Tuple{N,N},T}(t))

# Create aliases RotMatrix2{T} = RotMatrix{2,T,4} and RotMatrix3{T} = RotMatrix{3,T,9}
for N in 2:3
    L = N*N
    RotMatrixN = Symbol(:RotMatrix, N)
    @eval begin
        const $RotMatrixN{T} = RotMatrix{$N, T, $L}
        @inline $RotMatrixN(t::Tuple) = RotMatrix(SArray{Tuple{$N,$N}}(t))
    end
end

Base.@propagate_inbounds Base.getindex(r::RotMatrix, i::Int) = r.mat[i]
@inline Base.Tuple(r::RotMatrix) = Tuple(r.mat)

@inline RotMatrix(θ::Number) = RotMatrix{2}(θ)
@inline function (::Type{<:RotMatrix{2}})(θ::Number)
    s, c = sincos(θ)
    RotMatrix(@SMatrix [c -s; s c])
end
@inline function RotMatrix{2,T}(θ::Number) where T
    s, c = sincos(θ)
    RotMatrix(@SMatrix T[c -s; s c])
end

Base.one(::Type{R}) where {N,R<:RotMatrix{N}} = R(I)

# A rotation is more-or-less defined as being an orthogonal (or unitary) matrix
Base.inv(r::RotMatrix) = RotMatrix(r.mat')

# By default, composition of rotations will go through RotMatrix, unless overridden
@inline Base.:*(r1::Rotation, r2::Rotation) = RotMatrix(r1) * RotMatrix(r2)
@inline Base.:*(r1::RotMatrix, r2::Rotation) = r1 * RotMatrix(r2)
@inline Base.:*(r1::Rotation, r2::RotMatrix) = RotMatrix(r1) * r2
@inline Base.:*(r1::RotMatrix, r2::RotMatrix) = RotMatrix(r1.mat * r2.mat) # TODO check that this doesn't involve extra copying.

# Special case multiplication of 3×3 rotation matrices: speedup using cross product
@inline function Base.:*(r1::RotMatrix{3}, r2::RotMatrix{3})
    ret12 = r1 * r2[:, SVector(1, 2)]
    ret3 = ret12[:, 1] × ret12[:, 2]
    RotMatrix([ret12 ret3])
end

"""
    struct Angle2d{T} <: Rotation{2,T}
        theta::T
    end

A 2×2 rotation matrix parameterized by a 2D rotation by angle.
Only the angle is stored inside the `Angle2d` type, values
of `getindex` etc. are computed on the fly.
"""
struct Angle2d{T} <: Rotation{2,T}
    theta::T
    Angle2d{T}(theta::Number) where T = new{T}(theta)
end

@inline function Angle2d(theta::Number)
    Angle2d{rot_eltype(typeof(theta))}(theta)
end

params(r::Angle2d) = SVector{1}(r.theta)

Angle2d(r::Rotation{2}) = Angle2d(rotation_angle(r))
Angle2d{T}(r::Rotation{2}) where {T} = Angle2d{T}(rotation_angle(r))
@inline (::Type{R})(t::Tuple) where R<:Angle2d = convert(R, RotMatrix2(t))

Base.one(::Type{A}) where {A<: Angle2d} = A(0)

rotation_angle(rot::Angle2d) = rot.theta
function rotation_angle(rot::Rotation{2})
    c = @inbounds rot[1,1]
    s = @inbounds rot[2,1]
    atan(s, c)
end

@inline function Base.:*(r::Angle2d, v::StaticVector)
    if length(v) != 2
        throw(DimensionMismatch("Cannot rotate a vector of length $(length(v))"))
    end
    x,y = v
    s,c = sincos(r.theta)
    T = Base.promote_op(*, typeof(s), eltype(v))
    similar_type(v,T)(c*x - s*y, s*x + c*y)
end

Base.:*(r1::Angle2d, r2::Angle2d) = Angle2d(r1.theta + r2.theta)
Base.:^(r::Angle2d, t::Real) = Angle2d(r.theta*t)
Base.:^(r::Angle2d, t::Integer) = Angle2d(r.theta*t)
Base.inv(r::Angle2d) = Angle2d(-r.theta)

@inline function Base.getindex(r::Angle2d, i::Int)
    if i == 1
        cos(r.theta)
    elseif i == 2
        sin(r.theta)
    elseif i == 3
        -sin(r.theta)
    elseif i == 4
        cos(r.theta)
    else
        throw(BoundsError(r,i))
    end
end

################################################################################
################################################################################

_isrotation_eps(T::Type{<:Integer}) = zero(T)
_isrotation_eps(T) = eps(T)

"""
    isrotation(r)
    isrotation(r, tol)

Check whether `r` is a rotation matrix, where `r * r'` is within
`tol` of the identity matrix (using the Frobenius norm). `tol` defaults to
`1000 * eps(float(eltype(r)))` or `zero(T)` for integer `T`.
"""
isrotation

function isrotation(r::StaticMatrix{N,N,T}, tol::Real = 1000 * _isrotation_eps(eltype(T))) where {N,T<:Real}
    m = SMatrix(r)
    d = norm(m*m'-one(SMatrix{N,N}))
    return d ≤ tol && det(m) > 0
end

function isrotation(r::AbstractMatrix{T}, tol::Real = 1000 * _isrotation_eps(eltype(T))) where T<:Real
    s = size(r)
    if s[1] != s[2]
        return false
    end
    d = norm(r*r'-one(r))
    return d ≤ tol && det(r) > 0
end

function isrotation(r::AbstractMatrix{T}, tol::Real = 1000 * _isrotation_eps(eltype(real(T)))) where T<:Number
    if !isreal(r)
        return false
    end
    return isrotation(real(r), tol)
end

"""
    nearest_rotation(M) -> RotMatrix

Get the nearest special orthonormal matrix from given matrix `M`.
See [Wahba's problem](https://en.wikipedia.org/wiki/Wahba%27s_problem) for more information.
"""
function nearest_rotation(M::StaticMatrix{N,N}) where N
    u, _, v = svd(M)
    s = sign(det(u * v'))
    d = @SVector ones(N-1)
    R = u * Diagonal(push(d,s)) * v'
    return RotMatrix{N}(R)
end

function nearest_rotation(M::AbstractMatrix{T}) where T
    N = size(M,1)
    L = N^2
    M_ = convert(SMatrix{N,N,T,L}, M)
    return nearest_rotation(M_)
end

# A simplification and specialization of the Base.show function for AbstractArray makes
# everything sensible at the REPL.
function Base.show(io::IO, ::MIME"text/plain", X::Rotation)
    if !haskey(io, :compact)
        io = IOContext(io, :compact => true)
    end
    summary(io, X)
    if !isa(X, RotMatrix)
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
