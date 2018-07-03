"""
    abstract type Rotation{N,T} <: StaticMatrix{N,N,T}

An abstract type representing `N`-dimensional rotations. More abstractly, they represent
unitary (orthogonal) `N`×`N` matrices.
"""
abstract type Rotation{N,T} <: StaticMatrix{N,N,T} end

Base.@pure StaticArrays.Size(::Type{Rotation{N}}) where {N} = Size(N,N)
Base.@pure StaticArrays.Size(::Type{Rotation{N,T}}) where {N,T} = Size(N,N)
Base.@pure StaticArrays.Size(::Type{R}) where {R<:Rotation} = Size(supertype(R))
Compat.adjoint(r::Rotation) = inv(r)
Compat.transpose(r::Rotation{N,T}) where {N,T<:Real} = inv(r)

# Rotation angles and axes can be obtained by converting to the AngleAxis type
rotation_angle(r::Rotation) = rotation_angle(AngleAxis(r))
rotation_axis(r::Rotation) = rotation_axis(AngleAxis(r))

# Rotation matrices should be orthoginal/unitary. Only the operations we define,
# like multiplication, will stay as Rotations, otherwise users will get an
# SMatrix{3,3} (e.g. rot1 + rot2 -> SMatrix)
Base.@pure StaticArrays.similar_type(::Union{R,Type{R}}) where {R <: Rotation} = SMatrix{size(R)..., eltype(R), prod(size(R))}
Base.@pure StaticArrays.similar_type(::Union{R,Type{R}}, ::Type{T}) where {R <: Rotation, T} = SMatrix{size(R)..., T, prod(size(R))}

function Base.rand(::Type{R}) where R <: Rotation{2}
    T = eltype(R)
    if T == Any
        T = Float64
    end

    R(2π * rand(T))
end

# A random rotation can be obtained easily with unit quaternions
# The unit sphere in R⁴ parameterizes quaternion rotations according to the
# Haar measure of SO(3) - see e.g. http://math.stackexchange.com/questions/184086/uniform-distributions-on-the-space-of-rotations-in-3d
function Base.rand(::Type{R}) where R <: Rotation{3}
    T = eltype(R)
    if T == Any
        T = Float64
    end

    q = Quat(randn(T), randn(T), randn(T), randn(T))
    return R(q)
end

# Useful for converting arrays of rotations to another rotation eltype, for instance.
# Only works because parameters of all the rotations are of a similar form
# Would need to be more sophisticated if we have arbitrary dimensions, etc
@inline function Base.promote_op(::Type{R1}, ::Type{R2}) where {R1 <: Rotation, R2 <: Rotation}
    size(R1) == size(R2) || throw(DimensionMismatch("cannot promote rotations of $(size(R1)[1]) and $(size(R2)[1]) dimensions"))
    if isleaftype(R1)
        return R1
    else
        return R1{eltype(R2)}
    end
end

@inline function Base.:/(r1::Rotation, r2::Rotation)
    r1 * inv(r2)
end

@inline function Base.:\(r1::Rotation, r2::Rotation)
    inv(r1) * r2
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

# These functions (plus size) are enough to satisfy the entire StaticArrays interface:
# @inline (::Type{R}){R<:RotMatrix}(t::Tuple)  = error("No precise constructor found. Length of input was $(length(t)).")
for N = 2:3
    L = N*N
    RotMatrixN = Symbol(:RotMatrix, N)
    @eval begin
        @inline RotMatrix(t::NTuple{$L})  = RotMatrix(SMatrix{$N,$N}(t))
        @inline (::Type{RotMatrix{$N}})(t::NTuple{$L}) = RotMatrix(SMatrix{$N,$N}(t))
        @inline RotMatrix{$N,T}(t::NTuple{$L}) where {T} = RotMatrix(SMatrix{$N,$N,T}(t))
        @inline RotMatrix{$N,T,$L}(t::NTuple{$L}) where {T} = RotMatrix(SMatrix{$N,$N,T}(t))
        const $RotMatrixN{T} = RotMatrix{$N, T, $L}
    end
end
Base.@propagate_inbounds Base.getindex(r::RotMatrix, i::Int) = r.mat[i]
@inline Tuple(r::RotMatrix) = Tuple(r.mat)

@inline RotMatrix(θ::Real) = RotMatrix(@SMatrix [cos(θ) -sin(θ); sin(θ) cos(θ)])
@inline (::Type{RotMatrix{2}})(θ::Real)      = RotMatrix(@SMatrix [cos(θ) -sin(θ); sin(θ) cos(θ)])
@inline RotMatrix{2,T}(θ::Real) where {T} = RotMatrix(@SMatrix T[cos(θ) -sin(θ); sin(θ) cos(θ)])

# A rotation is more-or-less defined as being an orthogonal (or unitary) matrix
inv(r::RotMatrix) = RotMatrix(r.mat')

# A useful constructor for identity rotation (eye is already provided by StaticArrays, but needs an eltype)
@inline eye(::Type{RotMatrix{N}}) where {N} = RotMatrix((eye(SMatrix{N,N,Float64})))
@inline eye(::Type{RotMatrix{N,T}}) where {N,T} = RotMatrix((eye(SMatrix{N,N,T})))

# By default, composition of rotations will go through RotMatrix, unless overridden
@inline *(r1::Rotation, r2::Rotation) = RotMatrix(r1) * RotMatrix(r2)
@inline *(r1::RotMatrix, r2::Rotation) = r1 * RotMatrix(r2)
@inline *(r1::Rotation, r2::RotMatrix) = RotMatrix(r1) * r2
@inline *(r1::RotMatrix, r2::RotMatrix) = RotMatrix(r1.mat * r2.mat) # TODO check that this doesn't involve extra copying.

# Special case multiplication of 3×3 rotation matrices: speedup using cross product
@inline function *(r1::RotMatrix{3}, r2::RotMatrix{3})
    ret12 = r1 * r2[:, SVector(1, 2)]
    ret3 = ret12[:, 1] × ret12[:, 2]
    RotMatrix([ret12 ret3])
end

################################################################################
################################################################################

"""
    isrotation(r)
    isrotation(r, tol)

Check whether `r` is a 3×3 rotation matrix, where `r * r'` is within `tol` of
the identity matrix (using the Frobenius norm). (`tol` defaults to
`1000 * eps(eltype(r))`.)
"""
function isrotation(r::AbstractMatrix{T}, tol::Real = 1000 * eps(eltype(T))) where T
    if size(r) == (2,2)
        # Transpose is overloaded for many of our types, so we do it explicitly:
        r_trans = @SMatrix [conj(r[1,1])  conj(r[2,1]);
                            conj(r[1,2])  conj(r[2,2])]
        d = Compat.norm((r * r_trans) - eye(SMatrix{2,2}))
    elseif size(r) == (3,3)
        r_trans = @SMatrix [conj(r[1,1])  conj(r[2,1])  conj(r[3,1]);
                            conj(r[1,2])  conj(r[2,2])  conj(r[2,3]);
                            conj(r[1,3])  conj(r[2,3])  conj(r[3,3])]
        d = Compat.norm((r * r_trans) - eye(SMatrix{3,3}))
    else
        return false
    end

    return d < tol
end

@static if VERSION < v"0.7-"
    # A simplification and specialization of the Base.showarray() function makes
    # everything sensible at the REPL.
    function Base.showarray(io::IO, X::Rotation, repr::Bool = true; header = true)
        if !haskey(io, :compact)
            io = IOContext(io, compact=true)
        end
        if repr
            if isa(X, RotMatrix)
                Base.print_matrix_repr(io, X)
            else
                print(io, typeof(X).name.name)
                n_fields = length(fieldnames(typeof(X)))
                print(io, "(")
                for i = 1:n_fields
                    print(io, getfield(X, i))
                    if i < n_fields
                        print(io, ", ")
                    end
                end
                print(io, ")")
            end
        else
            if header
                print(io, summary(X))
                if !isa(X, RotMatrix)
                    n_fields = length(fieldnames(typeof(X)))
                    print(io, "(")
                    for i = 1:n_fields
                        print(io, getfield(X, i))
                        if i < n_fields
                            print(io, ", ")
                        end
                    end
                    print(io, ")")
                end
                println(io, ":")
            end
            punct = (" ", "  ", "")
            Base.print_matrix(io, X, punct...)
        end
    end
else
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
            for i = 1:n_fields
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
end

# Removes module name from output, to match other types
function Base.summary(r::Rotation{N,T}) where {T,N}
    inds = indices(r)
    typestring = last(split(string(typeof(r)), '.'; limit = 2))
    string(Base.dims2string(length.(inds)), " ", typestring)
end
