# we generate code in this module, so precompile where possible
__precompile__(true)

module Rotations

using LinearAlgebra
using StaticArrays
using Random
using Quaternions

import Statistics

include("util.jl")
include("core_types.jl")
include("unitquaternion.jl")
include("mrps.jl")
include("euler_types.jl")
include("angleaxis_types.jl")
# TODO: Add method `mean_rotation` instead of `mean`
# include("mean.jl")
include("derivatives.jl")
include("principal_value.jl")

include("rodrigues_params.jl")
include("error_maps.jl")
include("rotation_error.jl")
include("rotation_generator.jl")
include("logexp.jl")
include("eigen.jl")
include("rand.jl")
include("deprecated.jl")

export
    # Rotation types
    Rotation, RotMatrix, RotMatrix2, RotMatrix3,
    Angle2d,
    QuatRotation,
    AngleAxis, RotationVec, RodriguesParam, MRP,
    RotX, RotY, RotZ,
    RotXY, RotYX, RotZX, RotXZ, RotYZ, RotZY,
    RotXYX, RotYXY, RotZXZ, RotXZX, RotYZY, RotZYZ,
    RotXYZ, RotYXZ, RotZXY, RotXZY, RotYZX, RotZYX,

    # Deprecated, but export for compatibility
    UnitQuaternion,

    # rotation generators
    RotationGenerator,
    RotMatrixGenerator,
    RotMatrixGenerator2,
    RotMatrixGenerator3,
    Angle2dGenerator,
    RotationVecGenerator,

    # Quaternion math ops
    logm, expm, ⊖, ⊕, slerp,

    # Quaternion maps
    ExponentialMap, QuatVecMap, CayleyMap, MRPMap, IdentityMap,

    # check validity of the rotation (is it close to orthonormal?)
    isrotation,
    # check validity of the rotation (is it skew-symmetric?)
    isrotationgenerator,

    # Get nearest rotation matrix
    nearest_rotation,

    # angle and axis introspection
    rotation_angle,
    rotation_axis,

    # quaternion from two vectors
    rotation_between,

    # principal value of a rotation
    principal_value

    # derivatives (names clash with ForwarDiff?)
    #jacobian, hessian

end # module
