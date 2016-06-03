import Base: convert, call, getindex, *, eltype, eye, norm, inv
import Quaternions: axis, angle

@deprecate(rot_angle, rotation_angle)

# defines common methods for the various rotation types
RotTypeList = Vector{Type}(0)

# build a dictionary to build transformation paths from one representation to the next
# I hate this but trying to make a generic overload for convert causes no end of problems...
conversion_path = Dict()
defined_conversions = Vector{NTuple{2,Type}}(0)  # use this to keep track of all conversions that have been defined (fill in the blanks later)


numel{T}(::Type{T}) = length(fieldnames(T))  # for getting the number of elements in the representation
                                             # e.g. length(Quaternion) == 1                # its a number
                                             # e.g. length(fieldnames(Quaternion)) == 5    # it has a normalized field
                                             # e.g. numel(Quaternion) == 4                 # what we want

###################################################
# Rotation matrix (typeofalias of Mat{3,3,Float64})
# N.B. all other parameterizations T should implement
# a convert_rotation(T, RotMatrix) method
###################################################

"""
A type alias for rotation matrices (alias of Mat{3,3,T})
"""
typealias RotMatrix{eT}   Mat{3,3,eT}

# add to the type list
push!(RotTypeList, RotMatrix)

# Fixes for FixedSizeArray stuff (based on version 0.1.0)
@inline call{T}(::Type{RotMatrix}, mat::RotMatrix{T}) = mat          # weird result otherwise
@inline convert{T}(::Type{RotMatrix}, mat::Matrix{T}) = mat          # this hangs othewise
@inline call{T}(::Type{RotMatrix}, mat::Matrix{T}) = convert(RotMatrix{T}, mat)    # this hangs othewise

@inline call{T}(::Type{RotMatrix}, xt::NTuple{9,T}) = RotMatrix{T}(xt[1:3], xt[4:6], xt[7:9])

@inline call{T}(::Type{RotMatrix}, c1::NTuple{3,T}, c2::NTuple{3,T}, c3::NTuple{3,T}) = RotMatrix{T}(c1, c2, c3)
@inline call{T1,T2,T3}(::Type{RotMatrix}, c1::NTuple{3,T1}, c2::NTuple{3,T2}, c3::NTuple{3,T3}) = RotMatrix{promote_type(T1,T2,T3)}(c1, c2, c3)

@inline call{T}(::Type{RotMatrix}, x1::T, x2::T, x3::T, x4::T, x5::T, x6::T, x7::T, x8::T, x9::T) = RotMatrix{T}((x1,x2,x3), (x4,x5,x6), (x7,x8,x9))
@inline call{T}(::Type{RotMatrix{T}}, x1::Real, x2::Real, x3::Real, x4::Real, x5::Real, x6::Real, x7::Real, x8::Real, x9::Real) =
                    RotMatrix{T}((T(x1),T(x2),T(x3)), (T(x4),T(x5),T(x6)), (T(x7),T(x8),T(x9)))

conversion_path[RotMatrix] = nothing  # the blessed type, everything goes through RotMatrixs by default

# number of numeric elements
numel(::Type{RotMatrix}) = 9

strip_eltype{T <: RotMatrix}(::Type{T}) = RotMatrix

# define null rotations for convenience
@inline eye(::Type{RotMatrix}) = eye(RotMatrix{Float64})



###################################################
# Quaternions
# (Quaternions are defined in Quaternions.jl)
###################################################

# add to the type list
push!(RotTypeList, Quaternion)

# an extra constructor so we don't need to specify whether its normalized
@inline call{T}(::Type{Quaternion{T}}, a::Real, b::Real, c::Real, d::Real) = Quaternion(T(a), T(b), T(c), T(d))

# define its interaction with RotMatrixs
@inline convert(::Type{RotMatrix}, q::Quaternion) = quat_to_rot(q)
@inline convert{T}(::Type{RotMatrix{T}}, q::Quaternion) = convert(RotMatrix{T}, quat_to_rot(q))

@inline convert(::Type{Quaternion}, R::RotMatrix) = rot_to_quat(R)
@inline convert{T}(::Type{Quaternion{T}}, R::RotMatrix) = convert(Quaternion{T}, rot_to_quat(R))

# go from a Quaternion to any other representation via the Rotation Matrix
# @inline convert{T}(::Type{T}, q::Quaternion) = convert(T, quat_to_rot(q))  # go via the rotation representation
conversion_path[Quaternion] = RotMatrix
append!(defined_conversions, [(Quaternion, RotMatrix), (RotMatrix, Quaternion)])

# define its interaction with other angle representations
@inline eltype{T}(::Type{Quaternion{T}}) = T
@inline eltype{T}(::Quaternion{T}) = T

# Quaternions have an extra bool field that type_methods.jl doesn't want to know about
numel(::Type{Quaternion}) = 4
@inline getindex(X::Quaternion, i::Integer) = X.(i)

# angle and axis functions
@inline rotation_angle(q::Quaternion) = 2 * acos(q.s)
@inline function rotation_axis(q::Quaternion)
    aa = AngleAxis(q)
    Vec(aa.axis_x, aa.axis_y, aa.axis_z)
end

strip_eltype{T <: Quaternion}(::Type{T}) = Quaternion

# define null rotations for convenience
@inline eye(::Type{Quaternion}) = Quaternion(1.0, 0.0, 0.0, 0.0)
@inline eye{T}(::Type{Quaternion{T}}) = Quaternion{T}(T(1), T(0), T(0), T(0))



###################################################
# Stereographic projection of Quaternions
###################################################

"""

A 3 element parametrization of a unit quaternion Q formed by the intersection of a line from [-1,0,0,0] to Q, with a plane containing the origin and with normal direction [1,0,0,0]. This
is a compact representation of rotations where the derivitives of the rotation matrix's elements wrt the SpQuat parameters are rational functions (good for optimization)

See:

    [1] "A Recipe on the Parameterization of Rotation Matrices for Non-Linear Optimization using Quaternions",
    Terzakis, Culverhouse, Bugmann, Sharma, Sutton,
    MIDAS technical report, 2011
    http://www.tech.plymouth.ac.uk/sme/springerusv/2011/publications_files/Terzakis%20et%20al%202012,%20A%20Recipe%20on%20the%20Parameterization%20of%20Rotation%20Matrices...MIDAS.SME.2012.TR.004.pdf

    Note 1: the singularity (origin) has been moved from [0,0,0,-1] in Ref[1] to [-1,0,0,0], so the 0 rotation quaternion [1,0,0,0] maps to [0,0,0] as opposed of to [1,0,0].
    Note 2: the quaternion rotation ambiguity q = -q means there will be a rotation with ||SpQuat|| <= 1 and an equivilent rotation with ||SpQuat|| > 1.  This package
            uses the solution with ||SpQuat|| <= 1

"""
immutable SpQuat{T}
    x::T
    y::T
    z::T
end

# add to the type list
push!(RotTypeList, SpQuat)

# enable mixed element construction...
@inline call{T <: Real}(::Type{SpQuat}, x1::T, x2::T, x3::T) = SpQuat{T}(promote(x1,x2,x3)...)
@inline call(::Type{SpQuat}, x1::Real, x2::Real, x3::Real) = SpQuat(promote(x1,x2,x3)...)

@inline getindex(X::SpQuat, i::Integer) = X.(i)

# element conversions
@inline convert{T <: Real}(::Type{SpQuat{T}}, spq::SpQuat{T}) = spq
@inline convert{T <: Real}(::Type{SpQuat{T}}, spq::SpQuat) = SpQuat{T}(T(spq.x), T(spq.y), T(spq.z))

# convert to and from quaternions
@inline convert(::Type{Quaternion}, spq::SpQuat) = spquat_to_quat(spq)
@inline convert{T}(::Type{Quaternion{T}}, spq::SpQuat) = convert(Quaternion{T}, spquat_to_quat(spq))

@inline convert(::Type{SpQuat}, q::Quaternion) = quat_to_spquat(q)
@inline convert{T}(::Type{SpQuat{T}}, q::Quaternion) = convert(SpQuat{T}, quat_to_spquat(q))

# go from a SpQuat to any other representation via Quaternions
# @inline convert{T}(::Type{T}, spq::SpQuat) = convert(T, spquat_to_quat(spq))
conversion_path[SpQuat] = Quaternion
append!(defined_conversions, [(Quaternion, SpQuat), (SpQuat, Quaternion)])

# define an inverse for the SPQuat since its trivial.  Inverse in the sense of the corresponding inverse rotation
@inline inv(X::SpQuat) = SpQuat(-X.x, -X.y, -X.z)

# define multiplication for SpQuat corresponding to combining rotations
@inline *(lhs::SpQuat, rhs::SpQuat) = SpQuat(Quaternion(lhs) * Quaternion(rhs))

# rotation properties
@inline rotation_angle(X::SpQuat) = rotation_angle(Quaternion(X))
@inline function rotation_axis(spq::SpQuat)
    s = sqrt(sum(spq.x .* spq.x + spq.y .* spq.y + spq.z .* spq.z))
    Vec(spq.x / s, spq.y / s, spq.z / s)
end

# element type is handy
@inline eltype{T}(::Type{SpQuat{T}}) = T
@inline eltype{T}(::SpQuat{T}) = T
strip_eltype{T <: SpQuat}(::Type{T}) = SpQuat

# define null rotations for convenience
@inline eye(::Type{SpQuat}) = SpQuat(0.0, 0.0, 0.0)
@inline eye{T}(::Type{SpQuat{T}}) = SpQuat{T}(T(0), T(0), T(0))



###################################################
# Arbitrary axis rotation representations
# (limited support for now)
###################################################

"""
Rotation by theta about an arbitrary axis [axis_x, axis_y, axis_z]
"""
immutable AngleAxis{T}
    theta::T
    axis_x::T
    axis_y::T
    axis_z::T
end

# add to the type list
push!(RotTypeList, AngleAxis)

@inline getindex(X::AngleAxis, i::Integer) = X.(i)

# enable mixed element construction...
@inline call{T <: Real}(::Type{AngleAxis}, x1::T, x2::T, x3::T, x4::T) = AngleAxis{T}(promote(x1,x2,x3,x4)...)
@inline call(::Type{AngleAxis}, x1::Real, x2::Real, x3::Real, x4::Real) = AngleAxis(promote(x1,x2,x3,x4)...)

@inline convert{T <: Real}(::Type{AngleAxis{T}}, aa::AngleAxis{T}) = aa
@inline convert{T <: Real}(::Type{AngleAxis{T}}, aa::AngleAxis) = AngleAxis{T}(T(aa.theta), T(aa.axis_x), T(aa.axis_y), T(aa.axis_z))

# define its interaction with other angle representations
@inline convert(::Type{Quaternion}, aa::AngleAxis) = angleaxis_to_quat(aa)
@inline convert{T}(::Type{Quaternion{T}}, aa::AngleAxis) = convert(Quaternion{T}, angleaxis_to_quat(aa))

@inline convert(::Type{AngleAxis}, q::Quaternion) = quat_to_angleaxis(q)
@inline convert{T}(::Type{AngleAxis{T}}, q::Quaternion) = convert(AngleAxis{T}, quat_to_angleaxis(q))

@inline inv(aa::AngleAxis) = AngleAxis(-aa.theta, aa.axis_x, aa.axis_y, aa.axis_z)

# go from an AngleAxis to any other representation via Quaternions
# @inline convert{T}(::Type{T}, aa::AngleAxis) = convert(T, angleaxis_to_quat(spq))
conversion_path[AngleAxis] = Quaternion
append!(defined_conversions, [(Quaternion, AngleAxis), (AngleAxis, Quaternion)])

# accessors
@inline rotation_angle(aa::AngleAxis) = aa.theta #  - floor((aa.theta+pi) / (2*pi)) * 2*pi
@inline rotation_axis(aa::AngleAxis) = Vec(aa.axis_x, aa.axis_y, aa.axis_z)

# element type is handy
@inline eltype{T}(::Type{AngleAxis{T}}) = T
@inline eltype{T}(::AngleAxis{T}) = T

strip_eltype{T <: AngleAxis}(::Type{T}) = AngleAxis


# define null rotations for convenience
@inline eye(::Type{AngleAxis}) = AngleAxis(0.0, 1.0, 0.0, 0.0)
@inline eye{T}(::Type{AngleAxis{T}}) = AngleAxis{T}(T(0), T(1), T(0), T(0))


###################################################
# Rodrigues Vector theta *axis
# Don't trust autodiff on these!
###################################################

"""
Rodrigues vector parameterization (stored as theta * axis)
"""
immutable RodriguesVec{T}
    sx::T
    sy::T
    sz::T
end

# add to the type list
push!(RotTypeList, RodriguesVec)

@inline getindex(X::RodriguesVec, i::Integer) = X.(i)

# enable mixed element construction...
@inline call{T <: Real}(::Type{RodriguesVec}, x1::T, x2::T, x3::T) = RodriguesVec{T}(promote(x1,x2,x3)...)
@inline call(::Type{RodriguesVec}, x1::Real, x2::Real, x3::Real) = RodriguesVec(promote(x1,x2,x3)...)

@inline convert{T <: Real}(::Type{RodriguesVec{T}}, rv::RodriguesVec{T}) = rv
@inline convert{T <: Real}(::Type{RodriguesVec{T}}, rv::RodriguesVec) = RodriguesVec{T}(T(rv.sx), T(rv.sy), T(rv.sz))

# define its interaction with other angle representations
@inline convert(::Type{Quaternion}, rv::RodriguesVec) = rodrigues_to_quat(rv)
@inline convert{T}(::Type{Quaternion{T}}, rv::RodriguesVec) = convert(Quaternion{T}, rodrigues_to_quat(rv))
@inline convert(::Type{AngleAxis}, rv::RodriguesVec) = rodrigues_to_angleaxis(rv)
@inline convert{T}(::Type{AngleAxis{T}}, rv::RodriguesVec) = convert(AngleAxis{T}, rodrigues_to_angleaxis(rv))

@inline convert(::Type{RodriguesVec}, q::Quaternion) = quat_to_rodrigues(q)
@inline convert{T}(::Type{RodriguesVec{T}}, q::Quaternion) = convert(RodriguesVec{T}, quat_to_rodrigues(q))
@inline convert(::Type{RodriguesVec}, aa::AngleAxis) = angleaxis_to_rodrigues(aa)
@inline convert{T}(::Type{RodriguesVec{T}}, aa::AngleAxis) = convert(RodriguesVec{T}, angleaxis_to_rodrigues(aa))

# get here from a rotation matrix
@inline inv(rv::RodriguesVec) = RodriguesVec(-rv.sx, -rv.sy, -rv.sz)

# go from an RodriguesVec to any other representation via Quaternions
conversion_path[RodriguesVec] = Quaternion
append!(defined_conversions, [(Quaternion, RodriguesVec), (RodriguesVec, Quaternion)])
append!(defined_conversions, [(AngleAxis, RodriguesVec), (RodriguesVec, AngleAxis)])

# rotation properties
@inline norm(rv::RodriguesVec) = sqrt(rv.sx * rv.sx + rv.sy * rv.sy + rv.sz * rv.sz) # norm is meaningful for these things
@inline rotation_angle(rv::RodriguesVec) = norm(rv)
function rotation_axis{T}(rv::RodriguesVec{T})     # what should this return for theta = 0?
    theta = rotation_angle(rv)
    return (theta > eps(T) ? Vec(rv.sx / theta, rv.sy / theta, rv.sz / theta) : Vec(one(theta), zero(theta), zero(theta)))
end

# element type is handy
@inline eltype{T}(::Type{RodriguesVec{T}}) = T
@inline eltype{T}(::RodriguesVec{T}) = T

strip_eltype{T <: RodriguesVec}(::Type{T}) = RodriguesVec

# define null rotations for convenience
@inline eye(::Type{RodriguesVec}) = RodriguesVec(0.0, 0.0, 0.0)
@inline eye{T}(::Type{RodriguesVec{T}}) = RodriguesVec{T}(T(0), T(0), T(0))


###################################################
# Euler Angles
###################################################

# this file defines the types:
# EulerAngles       - for tait bryan ordering
# ProperEulerAngles - for proper Euler ordering
include("euler_types.jl")

# add to the type list
append!(RotTypeList, [EulerAngles, ProperEulerAngles])

#
# Euler Angles
#
@inline convert(::Type{EulerAngles}, x::Real, y::Real, z::Real) = EulerAngles{DefaultEulerOrder()}(x,y,z)

@inline getindex(X::EulerAngles, i::Integer) = X.(i)

@inline convert(::Type{EulerAngles}, ea::EulerAngles) = ea
@inline convert{ORD <: TaitByranOrder,T <: Real}(::Type{EulerAngles{ORD}}, ea::EulerAngles{ORD,T}) = ea
@inline convert{ORD <: TaitByranOrder,T <: Real}(::Type{EulerAngles{ORD,T}}, ea::EulerAngles{ORD,T}) = ea
@inline convert{ORD <: TaitByranOrder,T <: Real, U <: Real}(::Type{EulerAngles{ORD,T}}, ea::EulerAngles{ORD,U}) = EulerAngles{ORD,T}(T(ea.theta_x), T(ea.theta_y), T(ea.theta_z))


# define its interaction with RotMatrixs
@inline convert(::Type{RotMatrix}, ea::EulerAngles) = euler_to_rot(ea)
@inline convert{T}(::Type{RotMatrix{T}}, ea::EulerAngles) = convert(RotMatrix{T}, euler_to_rot(ea))

@inline convert{T <: EulerAngles}(::Type{T}, R::RotMatrix) = rot_to_euler(EulerAngles{euler_order(T)}, R)
@inline convert{ORD,T}(::Type{EulerAngles{ORD,T}}, R::RotMatrix) = rot_to_euler(EulerAngles{ORD,T}, R)

# go from a EulerAngles to any other representation via a rotation matrix
# @inline convert{T}(::Type{T}, ea::EulerAngles) = convert(T, euler_to_rot(ea))
conversion_path[EulerAngles] = RotMatrix
append!(defined_conversions, [(RotMatrix, EulerAngles), (EulerAngles, RotMatrix)])

# element type is handy
@inline eltype{ORD,T}(::Type{EulerAngles{ORD,T}}) = T
@inline eltype{ORD,T}(::EulerAngles{ORD,T}) = T

# order as well
@inline euler_order(::Type{EulerAngles}) = DefaultEulerOrder()
@inline euler_order{ORD}(::Type{EulerAngles{ORD}}) = ORD
@inline euler_order{ORD,T}(::Type{EulerAngles{ORD,T}}) = ORD
@inline euler_order{ORD,T}(::EulerAngles{ORD,T}) = ORD

# need special constructors to fill in the template parametera
@inline call(::Type{EulerAngles}, x::Real, y::Real, z::Real) = EulerAngles{DefaultEulerOrder()}(promote(x, y, z)...)
@inline call{ORD, T <: Real}(::Type{EulerAngles{ORD}}, x::T, y::T, z::T) = EulerAngles{ORD,T}(x, y, z)
@inline call{ORD}(::Type{EulerAngles{ORD}}, x::Real, y::Real, z::Real) = EulerAngles{ORD}(promote(x, y, z)...)

# convert from one order to another
@inline convert{ORD1, ORD2, eT}(::Type{EulerAngles{ORD1}}, ea::EulerAngles{ORD2, eT}) = convert(EulerAngles{ORD1}, euler_to_rot(ea))
@inline convert{ORD1, ORD2, eT1, eT2}(::Type{EulerAngles{ORD1, eT1}}, ea::EulerAngles{ORD2, eT2}) = convert(EulerAngles{ORD1, eT1}, euler_to_rot(ea))

strip_eltype{T <: EulerAngles}(::Type{T}) = EulerAngles{euler_order(T)}


#
# Proper Euler Angles
#
@inline convert(::Type{ProperEulerAngles}, x::Real, y::Real, z::Real) = EulerAngles{DefaultProperEulerOrder()}(x,y,z)

@inline getindex(X::ProperEulerAngles, i::Integer) = X.(i)

@inline convert(::Type{ProperEulerAngles}, ea::ProperEulerAngles) = ea
@inline convert{ORD <: ProperEulerOrder,T <: Real}(::Type{ProperEulerAngles{ORD}}, ea::ProperEulerAngles{ORD,T}) = ea
@inline convert{ORD <: ProperEulerOrder,T <: Real}(::Type{ProperEulerAngles{ORD,T}}, ea::ProperEulerAngles{ORD,T}) = ea
@inline convert{ORD <: ProperEulerOrder,T <: Real, U <: Real}(::Type{ProperEulerAngles{ORD,T}}, ea::ProperEulerAngles{ORD,U}) = ProperEulerAngles{ORD,T}(T(ea.theta_1), T(ea.theta_2), T(ea.theta_3))

# define its interaction with RotMatrixs
@inline convert(::Type{RotMatrix}, ea::ProperEulerAngles) = euler_to_rot(ea)
@inline convert{T}(::Type{RotMatrix{T}}, ea::ProperEulerAngles) = convert(RotMatrix{T}, euler_to_rot(ea))

@inline convert{T <: ProperEulerAngles}(::Type{T}, R::RotMatrix) = rot_to_euler(ProperEulerAngles{euler_order(T)}, R)
@inline convert{ORD,T}(::Type{ProperEulerAngles{ORD,T}}, R::RotMatrix) = rot_to_euler(ProperEulerAngles{ORD,T}, R)

# go from a ProperEulerAngles to any other representation via a rotation matrix
# @inline convert{T}(::Type{T}, ea::ProperEulerAngles) = convert(T, euler_to_rot(ea))
conversion_path[ProperEulerAngles] = RotMatrix
append!(defined_conversions, [(RotMatrix, ProperEulerAngles), (ProperEulerAngles, RotMatrix)])

# element type is handy
@inline eltype{ORD,T}(::Type{ProperEulerAngles{ORD,T}}) = T
@inline eltype{ORD,T}(::ProperEulerAngles{ORD,T}) = T

# order as well
@inline euler_order(::Type{ProperEulerAngles}) = DefaultProperEulerOrder()
@inline euler_order{ORD}(::Type{ProperEulerAngles{ORD}}) = ORD
@inline euler_order{ORD,T}(::Type{ProperEulerAngles{ORD,T}}) = ORD
@inline euler_order{ORD,T}(::ProperEulerAngles{ORD,T}) = ORD

# need special constructors to fill in the template parametera
@inline call(::Type{ProperEulerAngles}, x::Real, y::Real, z::Real) = ProperEulerAngles{DefaultProperEulerOrder()}(promote(x, y, z)...)
@inline call{ORD, T <: Real}(::Type{ProperEulerAngles{ORD}}, x::T, y::T, z::T) = ProperEulerAngles{ORD,T}(x, y, z)
@inline call{ORD}(::Type{ProperEulerAngles{ORD}}, x::Real, y::Real, z::Real) = ProperEulerAngles{ORD}(promote(x, y, z)...)

# convert from one order to another
@inline convert{ORD1, ORD2, eT}(::Type{ProperEulerAngles{ORD1}}, ea::ProperEulerAngles{ORD2, eT}) = convert(ProperEulerAngles{ORD1}, euler_to_rot(ea))
@inline convert{ORD1, ORD2, eT1, eT2}(::Type{ProperEulerAngles{ORD1, eT1}}, ea::ProperEulerAngles{ORD2, eT2}) = convert(ProperEulerAngles{ORD1, eT1}, euler_to_rot(ea))

strip_eltype{T <: ProperEulerAngles}(::Type{T}) = ProperEulerAngles{euler_order(T)}

# define null rotations for convenience
@inline eye(::Type{EulerAngles}) = EulerAngles(0.0, 0.0, 0.0)
@inline eye{ORD}(::Type{EulerAngles{ORD}}) = EulerAngles{ORD}(0.0, 0.0, 0.0)
@inline eye{ORD, T}(::Type{EulerAngles{ORD, T}}) = EulerAngles{ORD, T}(T(0), T(0), T(0))

# define null rotations for convenience
@inline eye(::Type{ProperEulerAngles}) = ProperEulerAngles(0.0, 0.0, 0.0)
@inline eye{ORD}(::Type{ProperEulerAngles{ORD}}) = ProperEulerAngles{ORD}(0.0, 0.0, 0.0)
@inline eye{ORD, T}(::Type{ProperEulerAngles{ORD, T}}) = ProperEulerAngles{ORD, T}(T(0), T(0), T(0))

