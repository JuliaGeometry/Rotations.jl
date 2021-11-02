
"""
    RotationError{T<:Real, D<:ErrorMap} <: StaticVector{3,T}

A three-parameter rotation error, converted to/from a `UnitQuaternion` using the
`ErrorMap` `D`.

# Usage
A `RotationError` is typically created using one of the following methods

    rotation_error(R1::Rotation, R2::Rotation, error_map::ErrorMap)
    R1 ⊖ R2

which compute the difference between the rotations `R1` and `R2` and convert the
result to a 3D rotation error using `error_map`.

The error can be "added" back to a rotation using the inverse operation:

    add_error(R1::Rotation, e::RotationError)
    R1::Rotation ⊕ e::RotationError
"""
struct RotationError{T,D} <: StaticVector{3,T}
    err::SVector{3,T}
    map::D
    @inline function RotationError(err::SVector{3,T}, map::D) where {T,D <: ErrorMap}
        new{T,D}(err, map)
    end
end

@inline UnitQuaternion(e::RotationError)::Rotation = e.map(e.err)

"""
    rotation_error(R1::Rotation, R2::Rotation, error_map::ErrorMap)

Compute the `RotationError` by calculating the "difference" between `R1` and `R2`, i.e.
`R2\\R1`, then mapped to a three-parameter error using `error_map`.

Can be equivalently called using the default map with `R1 ⊖ R2`

If `error_map::IdentityMap`, then `SVector(R1\\R2)::SVector{3}` is used as the error. Note
this only works for three-parameter rotation representations such as `RodriguesParam` or `MRP`.
"""
function rotation_error(R1::Rotation, R2::Rotation, error_map::ErrorMap)
    return RotationError(inv(error_map)(R2\R1), error_map)
end

function rotation_error(R1::Rotation, R2::Rotation, error_map::IdentityMap)
    err = params(R2\R1)
    if length(err) != 3
        throw(ArgumentError("R2\\R1 must be a three-dimensional parameterization, got $(length(err))"))
    end
    return RotationError(err, error_map)
end

#   set the default error map to the CayleyMap
@inline ⊖(R1::Rotation, R2::Rotation) = rotation_error(R1, R2, CayleyMap())

#   default to the identity map for Rodrigues Params and MRPs
@inline ⊖(R1::RodriguesParam, R2::RodriguesParam) = rotation_error(R1, R2, IdentityMap())
@inline ⊖(R1::MRP,            R2::MRP           ) = rotation_error(R1, R2, IdentityMap())


# Static Arrays interface
function (::Type{E})(t::NTuple{3}) where E <: RotationError
    E(t[1], t[2], t[3])
end
Base.@propagate_inbounds Base.getindex(e::RotationError, i::Int) = e.err[i]
@inline Base.Tuple(e::RotationError) = Tuple(e.err)

# Compose a rotation with an error
"""
    add_error(R1::Rotation, e::RotationError)

"Adds" the rotation error `e` to the rotation `R1` by converting `e` to a quaternion via
its `ErrorMap` and then composing with `R1`.

Equivalent to

    R1 ⊕ e
"""
function add_error(R1::Rotation, e::RotationError)
    R1 * UnitQuaternion(e)
end

function add_error(R1::R, e::RotationError{<:Any, IdentityMap}) where R <: Rotation
    # must be able to construct R from a SVector{3}
    if length(params(R1)) != 3
        throw(ArgumentError("Can't construct a rotation of type $R from three parameters"))
    end
    R1 * R(e.err)
end

@inline ⊕(R1::Rotation, e::RotationError) = add_error(R1, e)
