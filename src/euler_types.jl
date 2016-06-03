# File created: 08-Mar-2016 15:00:09
# using Matlab 2013b with Symbolic Math Toolbox version 5.11


##############################
# Euler Representations
##############################

"""
EulerOrder specifies the extrinsic (axes of rotation are fixed) order of rotations
"""
abstract EulerOrder     # Defines the order of rotations


"""
An abstract type for EulerAngles and ProperEulerAngles
"""
abstract AbstractEulerAngles{ORDER <: EulerOrder, T}     # An abstract type for Euler parameterizations




#
# tait-byran ordering
#
abstract TaitByranOrder   <: EulerOrder


"""
EulerZYX corresponds to:

    Extrinsic: R = Rz(theta_z) * Ry(theta_y) * Rx(theta_x)

or equivalently,

    Intrinsic: R = Rx''(theta_x) * Ry'(theta_y) * Rz(theta_z) (where the axis of rotation inherits from previous rotations)
"""
type EulerZYX <: TaitByranOrder end


"""
EulerZXY corresponds to:

    Extrinsic: R = Rz(theta_z) * Rx(theta_x) * Ry(theta_y)

or equivalently,

    Intrinsic: R = Ry''(theta_y) * Rx'(theta_x) * Rz(theta_z) (where the axis of rotation inherits from previous rotations)
"""
type EulerZXY <: TaitByranOrder end


"""
EulerYZX corresponds to:

    Extrinsic: R = Ry(theta_y) * Rz(theta_z) * Rx(theta_x)

or equivalently,

    Intrinsic: R = Rx''(theta_x) * Rz'(theta_z) * Ry(theta_y) (where the axis of rotation inherits from previous rotations)
"""
type EulerYZX <: TaitByranOrder end


"""
EulerYXZ corresponds to:

    Extrinsic: R = Ry(theta_y) * Rx(theta_x) * Rz(theta_z)

or equivalently,

    Intrinsic: R = Rz''(theta_z) * Rx'(theta_x) * Ry(theta_y) (where the axis of rotation inherits from previous rotations)
"""
type EulerYXZ <: TaitByranOrder end


"""
EulerXYZ corresponds to:

    Extrinsic: R = Rx(theta_x) * Ry(theta_y) * Rz(theta_z)

or equivalently,

    Intrinsic: R = Rz''(theta_z) * Ry'(theta_y) * Rx(theta_x) (where the axis of rotation inherits from previous rotations)
"""
type EulerXYZ <: TaitByranOrder end


"""
EulerXZY corresponds to:

    Extrinsic: R = Rx(theta_x) * Rz(theta_z) * Ry(theta_y)

or equivalently,

    Intrinsic: R = Ry''(theta_y) * Rz'(theta_z) * Rx(theta_x) (where the axis of rotation inherits from previous rotations)
"""
type EulerXZY <: TaitByranOrder end


#
# proper-euler ordering
#
abstract ProperEulerOrder   <: EulerOrder


"""
EulerZYZ corresponds to:

    Extrinsic: R = Rz(theta_1) * Ry(theta_2) * Rz(theta_3)

or equivalently,

    Intrinsic: R = Rz''(theta_3) * Ry'(theta_2) * Rz(theta_1) (where the axis of rotation inherits from previous rotations)
"""
type EulerZYZ <: ProperEulerOrder end


"""
EulerZXZ corresponds to:

    Extrinsic: R = Rz(theta_1) * Rx(theta_2) * Rz(theta_3)

or equivalently,

    Intrinsic: R = Rz''(theta_3) * Rx'(theta_2) * Rz(theta_1) (where the axis of rotation inherits from previous rotations)
"""
type EulerZXZ <: ProperEulerOrder end


"""
EulerYZY corresponds to:

    Extrinsic: R = Ry(theta_1) * Rz(theta_2) * Ry(theta_3)

or equivalently,

    Intrinsic: R = Ry''(theta_3) * Rz'(theta_2) * Ry(theta_1) (where the axis of rotation inherits from previous rotations)
"""
type EulerYZY <: ProperEulerOrder end


"""
EulerYXY corresponds to:

    Extrinsic: R = Ry(theta_1) * Rx(theta_2) * Ry(theta_3)

or equivalently,

    Intrinsic: R = Ry''(theta_3) * Rx'(theta_2) * Ry(theta_1) (where the axis of rotation inherits from previous rotations)
"""
type EulerYXY <: ProperEulerOrder end


"""
EulerXYX corresponds to:

    Extrinsic: R = Rx(theta_1) * Ry(theta_2) * Rx(theta_3)

or equivalently,

    Intrinsic: R = Rx''(theta_3) * Ry'(theta_2) * Rx(theta_1) (where the axis of rotation inherits from previous rotations)
"""
type EulerXYX <: ProperEulerOrder end


"""
EulerXZX corresponds to:

    Extrinsic: R = Rx(theta_1) * Rz(theta_2) * Rx(theta_3)

or equivalently,

    Intrinsic: R = Rx''(theta_3) * Rz'(theta_2) * Rx(theta_1) (where the axis of rotation inherits from previous rotations)
"""
type EulerXZX <: ProperEulerOrder end


##############################
# Euler (Tait Byran) Angles
##############################

@inline DefaultEulerOrder() = EulerZXY

"""
EulerAngles{Order <: TaitByranOrder, T <: AbstractFloat} where order specifies the extrinsic (axes of rotation are fixed) order of rotations.
See subtypes(Rotations.TaitByranOrder) for a list of possible orders

Fields:
    theta_x   - Angle to rotate about the x axis (radians)
    theta_y   - Angle to rotate about the y axis (radians)
    theta_z   - Angle to rotate about the z axis (radians)
"""
immutable EulerAngles{Order <: TaitByranOrder, T} <: AbstractEulerAngles{Order, T}
    theta_x::T
    theta_y::T
    theta_z::T
end


##############################
# Proper Euler Angles
##############################

@inline DefaultProperEulerOrder() = EulerXZX

"""
ProperEulerAngles{Order <: ProperEulerOrder, T <: AbstractFloat} where order specifies the extrinsic (axes of rotation are fixed) order of rotations.
See subtypes(Rotations.ProperEulerOrder) for a list of possible orders

Fields:
    theta_1   - angle to rotate about the first axis in the order (radians)
    theta_2   - angle to rotate about the second axis in the order (radians)
    theta_3   - angle to rotate about the final axis in the order (radians)
"""
immutable ProperEulerAngles{Order <: ProperEulerOrder, T} <: AbstractEulerAngles{Order, T}
    theta_1::T
    theta_2::T
    theta_3::T
end
