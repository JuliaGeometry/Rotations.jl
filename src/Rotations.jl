# we generate code in this module, so precompile where possible
__precompile__(true)

module Rotations

using StaticArrays

import Base: convert, eltype, size, length, getindex, inv, *, Tuple, eye

include("util.jl")
include("core_types.jl")
include("quaternion_types.jl")
include("angleaxis_types.jl")
include("euler_types.jl")
include("mean.jl")
include("derivatives.jl")


export
    Rotation, RotMatrix, RotMatrix2, RotMatrix3,
    Quat, SPQuat,
    AngleAxis, RodriguesVec,
    RotX, RotY, RotZ,
    RotXY, RotYX, RotZX, RotXZ, RotYZ, RotZY,
    RotXYX, RotYXY, RotZXZ, RotXZX, RotYZY, RotZYZ,
    RotXYZ, RotYXZ, RotZXY, RotXZY, RotYZX, RotZYX,

    # check validity of the rotation (is it close to unitary?)
    isrotation,

    # angle and axis introspection
    rotation_angle,
    rotation_axis,

    # quaternion from two vectors
    rotation_between,

    # difference between two angles
    angle_difference

    # derivatives (names clash with ForwarDiff?)
    #jacobian, hessian

end # module
