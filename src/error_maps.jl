
"""
    ErrorMap

A nonlinear mapping between the space of unit quaternions and three-dimensional
rotation errors.

These mappings are extremely useful for converting from globally nonsingular 3D rotation
representations such as `UnitQuaternion` or `RotMatrix3` to a three-parameter error that can
be efficiently used in gradient-based optimization methods that optimize deviations about
a current iterate using first or second-order information.

# Usage
    errmap(v::AbstractVector)  # "forward" map from a 3D error to a `UnitQuaternion`
    inv(errmap)(R::Rotation)   # "inverse" map from a rotation (via `UnitQuaternion`) to a 3D error

where `errmap <: ErrorMap`

# Implemented Maps
- `CayleyMap`:  Uses `RodriguesParam` as the error representation (default).
    Goes singular at 180° and does not have a sign ambiguity.
- `ExponentialMap`: Uses the canonical exponential map from Lie Group theory. Computationally
    expensive to compute. Exhibits kinematic singularities.
- `MRPMap`: Uses a scaled `MRP` as the error representation. Singular at 360° but has
    a sign ambiguity (with the "shadow" MRP set).
- `QuatVecMap`: Uses the vector part of the quaternions. Cheapest map to compute, but
    goes singular at 180° and suffers from sign ambiguity.
- `IdentityMap`: Maps values through directly. Only works with three-parameter rotation
    representations with the following methods: `R(::SVector{3})` and `SVector(::R)::SVector{3}`
"""
abstract type ErrorMap end
struct CayleyMap <: ErrorMap end
struct ExponentialMap <: ErrorMap end
struct MRPMap <: ErrorMap end
struct QuatVecMap <: ErrorMap end
struct IdentityMap <: ErrorMap end

"""
    InvErrorMap

The nonlinear mapping from unit quaternions to a three-dimensional error state. Obtained by
inverting an `ErrorMap`, i.e.

    InvCayleyMap() = inv(CayleyMap())

# Usage
    imap(R::Rotation)             # "inverse" map from a rotation to a 3D error
    inv(imap)(v::AbstractVector)  # "forward" map from a 3D error to a `UnitQuaternion`

where `imap <: InvErrorMap`.

See `ErrorMap` for documentation on the implemented maps.
"""
abstract type InvErrorMap end
struct InvCayleyMap <: InvErrorMap end
struct InvExponentialMap <: InvErrorMap end
struct InvMRPMap <: InvErrorMap end
struct InvQuatVecMap <: InvErrorMap end

inv(::CayleyMap) = InvCayleyMap()
inv(::ExponentialMap) = InvExponentialMap()
inv(::QuatVecMap) = InvQuatVecMap()
inv(::MRPMap) = InvMRPMap()
inv(::IdentityMap) = IdentityMap()

inv(::InvCayleyMap) = CayleyMap()
inv(::InvExponentialMap) = ExponentialMap()
inv(::InvQuatVecMap) = QuatVecMap()
inv(::InvMRPMap) = MRPMap()


# Scalings
@inline scaling(::Type{ExponentialMap}) = 0.5
@inline scaling(::Type{QuatVecMap}) = 1.0
@inline scaling(::Type{CayleyMap}) = 1.0
@inline scaling(::Type{MRPMap}) = 2.0
scaling(m::M) where M <: ErrorMap = scaling(M)

# Quaternion Maps
(::ExponentialMap)(ϕ::AbstractVector) = expm(ϕ/scaling(ExponentialMap))

function (::QuatVecMap)(v::AbstractVector)
    check_length(v, 3)
    μ = 1/scaling(QuatVecMap)
    UnitQuaternion(sqrt(1-μ^2*v'v), μ*v[1], μ*v[2], μ*v[3])
end

function (::CayleyMap)(g::AbstractVector)
    check_length(g, 3)
    g /= scaling(CayleyMap)
    M = 1/sqrt(1+g'g)
    UnitQuaternion(M, M*g[1], M*g[2], M*g[3])
end

function (::MRPMap)(p::AbstractVector)
    check_length(p, 3)
    p /= scaling(MRPMap)
    n2 = p'p
    M = 2/(1+n2)
    UnitQuaternion((1-n2)/(1+n2), M*p[1], M*p[2], M*p[3])
end

function (::IdentityMap)(q::AbstractVector)
    check_length(q, 3)
    UnitQuaternion(q[1], q[2], q[3], q[4])
end


# Quaternion Map Jacobians
"""
    jacobian(::ErrorMap, ϕ)

Jacobian of the quaternion map that takes a three-dimensional vector `ϕ` and returns a
    unit quaternion.
Returns a 4x3 Static Matrix

For all the maps (except the `IdentityMap`)
`jacobian(::ErrorMap, zeros(3)) = [0; I] = Hmat()'`
"""
function jacobian(::ExponentialMap, ϕ, eps=1e-5)
    μ = 1/scaling(ExponentialMap)
    θ = norm(ϕ)
    cθ = cos(μ*θ/2)
    sincθ = sinc(μ*θ/2π)
    if θ < eps
        0.5*μ*[-0.5*μ*sincθ*ϕ'; sincθ*I + (cθ - sincθ)*ϕ*ϕ']
    else
        0.5*μ*[-0.5*μ*sincθ*ϕ'; sincθ*I + (cθ - sincθ)*ϕ*ϕ'/(ϕ'ϕ)]
    end
end

function jacobian(::QuatVecMap, v)
    μ = 1/scaling(QuatVecMap)
    μ2 = μ*μ
    M = -μ2/sqrt(1-μ2*v'v)
    @SMatrix [v[1]*M v[2]*M v[3]*M;
              μ 0 0;
              0 μ 0;
              0 0 μ]
end

function jacobian(::CayleyMap, g)
    μ = 1/scaling(CayleyMap)
    μ2 = μ*μ
    n = 1+μ2*g'g
    ni = 1/n
    μ*[-μ*g'; -μ2*g*g' + I*n]*ni*sqrt(ni)
end

function jacobian(::MRPMap, p)
    μ = 1/scaling(MRPMap)
    μ2 = μ*μ
    n = 1+μ2*p'p
    2*[-2*μ2*p'; I*μ*n - 2*μ*μ2*p*p']/n^2
end

jacobian(::IdentityMap, q) = I



############################################################################################
#                             INVERSE RETRACTION MAPS
############################################################################################
# (emap::ErrorMap)(R::Rotation) =
#     throw(ArgumentError("Must use the inverse map to convert a rotation to an error state, e.g. inv(emap)(q)"))
# (emap::InvErrorMap)(R::Rotation) = emap(UnitQuaternion(R))  # automatically call the inverse map for rotations
# NOTE: Julia v1.0 doesn't allow these methods on abstract types. Leave out for now.

for imap in (InvExponentialMap,InvCayleyMap,InvMRPMap,InvQuatVecMap)
    @eval begin
        @inline (imap::$imap)(R::Rotation) = imap(UnitQuaternion(R))
    end
end

(::InvExponentialMap)(q::UnitQuaternion) = scaling(ExponentialMap)*logm(q)
(::InvCayleyMap)(q::UnitQuaternion) = scaling(CayleyMap) * vector(q)/q.w
(::InvMRPMap)(q::UnitQuaternion) = scaling(MRPMap)*vector(q)/(1+q.w)
(::InvQuatVecMap)(q::UnitQuaternion) = scaling(QuatVecMap)*vector(q) * sign(q.w)


# ~~~~~~~~~~~~~~~ Inverse map Jacobians ~~~~~~~~~~~~~~~ #
"""
    jacobian(::InvErrorMap, q::UnitQuaternion)

Jacobian of the inverse quaternion map, returning a 3x4 matrix.
For all maps: `jacobian(::InvErrorMap, UnitQuaternion(I)) = [0 I] = Hmat()'`
"""
function jacobian(::InvExponentialMap, q::UnitQuaternion, eps=1e-5)
    μ = scaling(ExponentialMap)
    s = scalar(q)
    v = vector(q)
    θ2 = v'v
    θ = sqrt(θ2)
    datan = 1/(θ2 + s^2)
    ds = -datan*v

    if θ < eps
        return 2*μ*[ds (v*v' + I)/s]
    else
        atanθ = atan(θ,s)
        dv = ((s*datan - atanθ/θ)v*v'/θ + atanθ*I )/θ
        d0 = ((s*datan - atanθ/θ)v*v'/θ^2 + atanθ/θ*I )
        d1 = (s*datan - atanθ/θ)
        d2 = v*v'/θ2
        d3 = atanθ/θ * I
        return 2*μ*[ds dv]
    end
end


function jacobian(::InvQuatVecMap, q::UnitQuaternion)
    μ = scaling(QuatVecMap)
    return sign(q.w) * SA[
                0. μ 0 0;
                0. 0 μ 0;
                0. 0 0 μ]
end


function jacobian(::InvCayleyMap, q::UnitQuaternion)
    μ = scaling(CayleyMap)
    si = 1/q.w
    return μ*@SMatrix [-si^2*q.x si 0 0;
                       -si^2*q.y 0 si 0;
                       -si^2*q.z 0 0 si]
end


function jacobian(::InvMRPMap, q::UnitQuaternion)
    μ = scaling(MRPMap)
    si = 1/(1+q.w)
    return μ*@SMatrix [-si^2*q.x si 0 0;
                       -si^2*q.y 0 si 0;
                       -si^2*q.z 0 0 si]
end


jacobian(::IdentityMap, q::UnitQuaternion) = I


# ~~~~~~~~~~~~~~~ Inverse map Jacobian derivative ~~~~~~~~~~~~~~~ #
"""
    ∇jacobian(::InvErrorMap, q::UnitQuaternion, b::SVector{3})

Jacobian of G(q)'b, where G(q) = jacobian(::InvErrorMap, q),
    b is a 3-element vector
"""
function ∇jacobian(::InvExponentialMap, q::UnitQuaternion, b::SVector{3}, eps=1e-5)
    μ = scaling(ExponentialMap)
    s = scalar(q)
    v = vector(q)
    θ2 = v'v
    θ = sqrt(θ2)
    datan = 1/(θ2 + s^2)
    ds = -datan*v

    if θ < eps
        # return 2*μ*[b'ds; (v*v'b + b)/s]
        return 2*μ*[b'*(datan^2*2s*v) -b'datan*I;
                    -(v*v'b +b)/s^2 (I*(v'b) + v*b')/s]
    else
        dsds = 2b's*datan^2*v
        dsdv = b'*(-datan*I + 2datan^2*v*v')

        atanθ = atan(θ,s)
        d1 = (s*datan - atanθ/θ)
        d2 = v*v'b/θ2
        d3 = atanθ/θ*b
        d1ds = (datan - 2s^2*datan^2 + datan)
        dvds = d1ds*d2 - datan*b

        d1dv =  (-2s*datan^2*v' - s*datan*v'/θ^2 + atanθ/θ^3*v')
        d2dv = (I*(v'b) + v*b')/θ2 - 2(v*v'b)/θ^4 * v'
        d3dv = b*(s*datan*v'/θ^2 - atanθ/θ^3*v')
        dvdv = d2*d1dv + d1*d2dv + d3dv

        # return 2*μ*[ds'b; dv'b]
        return 2*μ*@SMatrix [
            dsds    dsdv[1] dsdv[2] dsdv[3];
            dvds[1] dvdv[1] dvdv[4] dvdv[7];
            dvds[2] dvdv[2] dvdv[5] dvdv[8];
            dvds[3] dvdv[3] dvdv[6] dvdv[9];
        ]
    end
end

function ∇jacobian(::InvCayleyMap, q::UnitQuaternion, b::SVector{3})
    μ = scaling(CayleyMap)
    si = 1/q.w
    v = vector(q)
    μ*@SMatrix [
        2*si^3*(v'b) -si^2*b[1] -si^2*b[2] -si^2*b[3];
       -si^2*b[1] 0 0 0;
       -si^2*b[2] 0 0 0;
       -si^2*b[3] 0 0 0;
    ]
end

function ∇jacobian(::InvMRPMap, q::UnitQuaternion, b::SVector{3})
    μ = scaling(MRPMap)
    si = 1/(1+q.w)
    v = vector(q)
    μ * @SMatrix [
        2*si^3*(v'b) -si^2*b[1] -si^2*b[2] -si^2*b[3];
       -si^2*b[1] 0 0 0;
       -si^2*b[2] 0 0 0;
       -si^2*b[3] 0 0 0;
    ]
end

function ∇jacobian(::InvQuatVecMap, q::UnitQuaternion, b::SVector{3})
    μ = scaling(QuatVecMap)
    @SMatrix zeros(4,4)
end
