# we generate code in this module, so precompile where possible
__precompile__(true)

module Rotations

using FixedSizeArrays
using Quaternions  # TODO: Switch to the UnitQuaternions package? https://github.com/kam3k/UnitQuaternions.jl

import Base: convert, mean, eltype, isnan, length, call, getindex, inv, .*, tuple, vec
import Quaternions.Quaternion

include("rotation_types.jl")
include("rotation_conversions.jl")
include("euler_conversions.jl")
include("type_methods.jl")
include("rotate.jl")
include("quaternion_ops.jl")
include("math_funcs.jl")
include("rot_covariances.jl")


export

        # Euler angle ordering
        EulerOrder,
        EulerZXY,

        # representations
        RotMatrix,
        SpQuat,
        AngleAxis,
        RodriguesVec,
        EulerAngles,         # order is in the type specification
        ProperEulerAngles,   # order is in the type specification

        # Quaternions as well
        Quaternion,

        # We use Vec for points so export it from fixed size arrays
        Vec,

        # rotate stuff
        rotate,

        # retrieve the euler order
        euler_order,

        # check validity of the rotation
        valid_rotation,

        # extra quaternion operations
        rotation_angle,  # N.B. now depricated in favour of angle
        rotation_axis

end # end rotations
