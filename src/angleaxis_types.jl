################################################################################
################################################################################
"""
    struct AngleAxis{T} <: Rotation{3,T}
    AngleAxis(Θ, x, y, z)

A 3×3 rotation matrix parameterized by a 3D rotation by angle θ about an
arbitrary axis `[x, y, z]`.

Note that the axis is not unique for θ = 0, and that this parameterization does
not continuously map the neighbourhood of the identity rotation (and therefore
might not be suitable for autodifferentation and optimization purposes).

Note: by default, the constructor will renormalize the input so that the axis
has length 1 (x² + y² + z² = 1).

Renormalization can be skipped by passing `false` as an additional constructor
argument, in which case the user provides the guarantee that the input arguments
represent a normalized rotation axis. Operations on an `AngleAxis` with a rotation
axis that does not have unit norm, created by skipping renormalization in this fashion,
are not guaranteed to do anything sensible.
"""
struct AngleAxis{T} <: Rotation{3,T}
    theta::T
    axis_x::T
    axis_y::T
    axis_z::T

    @inline function AngleAxis{T}(θ, x, y, z, normalize::Bool = true) where {T}
        if normalize
            # Not sure what to do with theta?? Should it become theta * norm ?
            norm = sqrt(x*x + y*y + z*z)
            new(θ, x/norm, y/norm, z/norm)
        else
            new(θ, x, y, z)
        end
    end
end

params(aa::AngleAxis) = SVector{4}(aa.theta, aa.axis_x, aa.axis_y, aa.axis_z)

# StaticArrays will take over *all* the constructors and put everything in a tuple...
# but this isn't quite what we mean when we have 4 inputs (not 9).
@inline function AngleAxis(θ::Θ, x::X, y::Y, z::Z, normalize::Bool = true) where {Θ,X,Y,Z}
    AngleAxis{promote_type(promote_type(promote_type(Θ, X), Y), Z)}(θ, x, y, z, normalize)
end

# These functions are enough to satisfy the entire StaticArrays interface:
@inline (::Type{AA})(t::NTuple{9}) where {AA <: AngleAxis} = AA(QuatRotation(t)) # TODO: consider going directly from tuple (RotMatrix) to AngleAxis
@inline Base.getindex(aa::AngleAxis, i::Int) = QuatRotation(aa)[i]

@inline function Base.Tuple(aa::AngleAxis{T}) where T
    # Rodrigues' rotation formula.
    s, c = sincos(aa.theta)
    c1 = one(T) - c

    c1x2 = c1 * aa.axis_x^2
    c1y2 = c1 * aa.axis_y^2
    c1z2 = c1 * aa.axis_z^2

    c1xy = c1 * aa.axis_x * aa.axis_y
    c1xz = c1 * aa.axis_x * aa.axis_z
    c1yz = c1 * aa.axis_y * aa.axis_z

    sx = s * aa.axis_x
    sy = s * aa.axis_y
    sz = s * aa.axis_z

    # Note that the RotMatrix constructor argument order makes this look transposed:
    (one(T) - c1y2 - c1z2, c1xy + sz, c1xz - sy,
        c1xy - sz, one(T) - c1x2 - c1z2, c1yz + sx,
        c1xz + sy, c1yz - sx, one(T) - c1x2 - c1y2)
end

@inline function (::Type{Q})(aa::AngleAxis) where Q <: QuatRotation
    s, c = sincos(aa.theta / 2)
    return Q(c, s * aa.axis_x, s * aa.axis_y, s * aa.axis_z, false)
end

@inline function (::Type{AA})(q::QuatRotation) where AA <: AngleAxis
    w = real(q.q)
    x, y, z = imag_part(q.q)
    s2 = x * x + y * y + z * z
    sin_t2 = sqrt(s2)
    theta = 2 * atan(sin_t2, w)
    num_pert = eps(typeof(theta))^4
    inv_sin_t2 = 1 / (sin_t2 + num_pert)
    return principal_value(AA(theta, inv_sin_t2 * (x + num_pert), inv_sin_t2 * y, inv_sin_t2 * z, false))
end

# Trivial type conversions for RotX, RotY and RotZ
@inline function (::Type{AA})(r::RotX) where AA <: AngleAxis
    return AA(r.theta, 1, 0, 0)
end
@inline function (::Type{AA})(r::RotY) where AA <: AngleAxis
    return AA(r.theta, 0, 1, 0)
end
@inline function (::Type{AA})(r::RotZ) where AA <: AngleAxis
    return AA(r.theta, 0, 0, 1)
end

# Using Rodrigues formula on an AngleAxis parameterization (assume unit axis length) to do the rotation
# (implementation from: https://ceres-solver.googlesource.com/ceres-solver/+/1.10.0/include/ceres/rotation.h)
function Base.:*(aa::AngleAxis, v::StaticVector)
    if length(v) != 3
        throw("Dimension mismatch: cannot rotate a vector of length $(length(v))")
    end

    w = rotation_axis(aa)
    st, ct = sincos(aa.theta)
    w_cross_pt = cross(w, v)
    m = dot(v, w) * (one(w_cross_pt[1]) - ct)
    T = Base.promote_op(*, eltype(aa), eltype(v))
    return similar_type(v,T)(v[1] * ct + w_cross_pt[1] * st + w[1] * m,
                             v[2] * ct + w_cross_pt[2] * st + w[2] * m,
                             v[3] * ct + w_cross_pt[3] * st + w[3] * m)
end

@inline Base.:*(aa::AngleAxis, r::Rotation) = QuatRotation(aa) * r
@inline Base.:*(aa::AngleAxis, r::RotMatrix) = QuatRotation(aa) * r
@inline Base.:*(aa::AngleAxis, r::MRP) = QuatRotation(aa) * r
@inline Base.:*(r::Rotation, aa::AngleAxis) = r * QuatRotation(aa)
@inline Base.:*(r::RotMatrix, aa::AngleAxis) = r * QuatRotation(aa)
@inline Base.:*(r::MRP, aa::AngleAxis) = r * QuatRotation(aa)
@inline Base.:*(aa1::AngleAxis, aa2::AngleAxis) = QuatRotation(aa1) * QuatRotation(aa2)

@inline Base.inv(aa::AngleAxis) = AngleAxis(-aa.theta, aa.axis_x, aa.axis_y, aa.axis_z)

# define identity rotations for convenience
@inline Base.one(::Type{AngleAxis}) = AngleAxis(0.0, 1.0, 0.0, 0.0)
@inline Base.one(::Type{AngleAxis{T}}) where {T} = AngleAxis{T}(zero(T), one(T), zero(T), zero(T))

# accessors
@inline rotation_angle(aa::AngleAxis) = aa.theta #  - floor((aa.theta+pi) / (2*pi)) * 2*pi
@inline rotation_axis(aa::AngleAxis) = SVector(aa.axis_x, aa.axis_y, aa.axis_z)


################################################################################
################################################################################
"""
    struct RotationVec{T} <: Rotation{3,T}
    RotationVec(sx, sy, sz)

Rodrigues vector parameterization of a 3×3 rotation matrix. The direction of the
vector [sx, sy, sz] defines the axis of rotation, and the rotation angle is
given by its norm.
"""
struct RotationVec{T} <: Rotation{3,T}
    sx::T
    sy::T
    sz::T
end

params(aa::RotationVec) = SVector{3}(aa.sx, aa.sy, aa.sz)

# StaticArrays will take over *all* the constructors and put everything in a tuple...
# but this isn't quite what we mean when we have 3 inputs (not 9).
@inline RotationVec(x::X, y::Y, z::Z) where {X,Y,Z} = RotationVec{promote_type(promote_type(X, Y), Z)}(x, y, z)

# These functions are enough to satisfy the entire StaticArrays interface:
@inline (::Type{RV})(t::NTuple{9}) where {RV <: RotationVec} = RV(QuatRotation(t)) # TODO: go through AngleAxis once it's faster
@inline Base.getindex(aa::RotationVec, i::Int) = QuatRotation(aa)[i]
@inline Base.Tuple(rv::RotationVec) = Tuple(QuatRotation(rv))

function (::Type{AA})(rv::RotationVec) where AA <: AngleAxis
    # TODO: consider how to deal with derivative near theta = 0. There should be a first-order expansion here.
    theta = rotation_angle(rv)
    return theta > 0 ? AA(theta, rv.sx / theta, rv.sy / theta, rv.sz / theta, false) : AA(zero(theta), one(theta), zero(theta), zero(theta), false)
end

function (::Type{RV})(aa::AngleAxis) where RV <: RotationVec
    return RV(aa.theta * aa.axis_x, aa.theta * aa.axis_y, aa.theta * aa.axis_z)
end

@inline function (::Type{Q})(rv::RotationVec) where Q <: QuatRotation
    theta = rotation_angle(rv)
    if theta < sqrt(eps(typeof(theta)))
        ϕ = theta / π / 2
        sc = sinc(ϕ) / 2    # this form gracefully handles theta = 0
        qx, qy, qz = sc * rv.sx, sc * rv.sy, sc * rv.sz
        return Q(cos(theta / 2), qx, qy, qz, false)
    else
        s, c = sincos(theta / 2)
        sθ = s / theta
        return Q(c, sθ * rv.sx, sθ * rv.sy, sθ * rv.sz, false)
    end
end

(::Type{RV})(q::QuatRotation) where {RV <: RotationVec} = RV(AngleAxis(q))

Base.:*(rv::RotationVec, v::StaticVector{3}) = QuatRotation(rv) * v

@inline Base.:*(rv::RotationVec, r::Rotation) = QuatRotation(rv) * r
@inline Base.:*(rv::RotationVec, r::RotMatrix) = QuatRotation(rv) * r
@inline Base.:*(rv::RotationVec, r::MRP) = QuatRotation(rv) * r
@inline Base.:*(rv::RotationVec, r::AngleAxis) = QuatRotation(rv) * r
@inline Base.:*(r::Rotation, rv::RotationVec) = r * QuatRotation(rv)
@inline Base.:*(r::RotMatrix, rv::RotationVec) = r * QuatRotation(rv)
@inline Base.:*(r::MRP, rv::RotationVec) = r * QuatRotation(rv)
@inline Base.:*(r::AngleAxis, rv::RotationVec) = r * QuatRotation(rv)
@inline Base.:*(rv1::RotationVec, rv2::RotationVec) = QuatRotation(rv1) * QuatRotation(rv2)

@inline Base.inv(rv::RotationVec) = RotationVec(-rv.sx, -rv.sy, -rv.sz)

# rotation properties
@inline rotation_angle(rv::RotationVec) = hypot(rv.sx, rv.sy, rv.sz)
function rotation_axis(rv::RotationVec)     # what should this return for theta = 0?
    theta = rotation_angle(rv)
    return (theta > 0 ? SVector(rv.sx / theta, rv.sy / theta, rv.sz / theta) : SVector(one(theta), zero(theta), zero(theta)))
end

# define identity rotations for convenience
@inline Base.one(::Type{RotationVec}) = RotationVec(0.0, 0.0, 0.0)
@inline Base.one(::Type{RotationVec{T}}) where {T} = RotationVec{T}(zero(T), zero(T), zero(T))
