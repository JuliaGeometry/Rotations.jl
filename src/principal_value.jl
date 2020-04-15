"""
    principal_value(R::Rotation{3})

**Background:** All non `RotMatrix` rotation types can represent the same `RotMatrix` in two or more ways. Sometimes a
particular set of numbers is better conditioned (e.g. `SPQuat`) or obeys a particular convention (e.g. `AngleAxis` has
non-negative rotation). In order to preserve differentiability it is necessary to allow rotation representations to
travel slightly away from the nominal domain; this is critical for applications such as optimization or dynamics.

This function takes a rotation type (e.g. `Quat`, `RotXY`) and outputs a new rotation of the same type that corresponds
to the same `RotMatrix`, but that obeys certain conventions or is better conditioned. The outputs of the function have
the following properties:

- all angles are between between `-pi` to `pi` (except for `AngleAxis` which is between `0` and `pi`).
- all `Quat` have non-negative real part
- the components of all `SPQuat` have a norm that is at most 1.
- the `RotationVec` rotation is at most `pi`

"""
principal_value(r::RotMatrix) = r
principal_value(q::Quat{T}) where {T} = q.w < zero(T) ? Quat{T}(-q.w, -q.x, -q.y, -q.z, false) : q
principal_value(q::Q) where Q <: UnitQuaternion = q.w < zero(eltype(q)) ? Q(-q.w, -q.x, -q.y, -q.z, false) : q
function principal_value(spq::SPQuat{T}) where {T}
    # A quat with positive real part: Quat( qw,  qx,  qy,  qz)
    #
    # A spq corresponding to the Quat with a positive real part:
    # SPQuat( qx / (1 + qw),  qy / (1 + qw),  qz / (1 + qw)) ≡ SPQuat(spx, spy, spz)
    #
    # A spq corresponding to the Quat with a negative real part:
    # SPQuat(-qx / (1 - qw), -qy / (1 - qw), -qz / (1 - qw)) ≡ SPQuat(snx, sny, snz)
    #
    # Claim:         spx / snx = -1 / (spx^2 + spy^2 + spz^2)
    #     -(1 + qw) / (1 - qw) = -1 / (spx^2 + spy^2 + spz^2)
    #      (1 - qw) / (1 + qw) = (spx^2 + spy^2 + spz^2)
    #      (1 - qw) / (1 + qw) = (qx^2 + qy^2 + qz^2) / (1 + qw)^2
    #      (1 - qw) * (1 + qw) =  qx^2 + qy^2 + qz^2
    #                 1 - qw^2 =  qx^2 + qy^2 + qz^2 (Q.E.D.)
    alpha_2 = spq.x^2 + spq.y^2 + spq.z^2
    if one(T) < alpha_2
        scale = -one(T) / alpha_2
        return SPQuat(scale * spq.x, scale * spq.y, scale * spq.z)
    else
        return spq
    end
end

function principal_value(aa::AngleAxis{T}) where {T}
    theta = rem2pi(aa.theta, RoundNearest)
    if theta < zero(T)
        return AngleAxis(-theta, -aa.axis_x, -aa.axis_y, -aa.axis_z, false)
    else
        return AngleAxis( theta,  aa.axis_x,  aa.axis_y,  aa.axis_z, false)
    end
end

function principal_value(rv::RotationVec{T}) where {T}
    theta = rotation_angle(rv)
    if pi < theta
        re_s = rem2pi(theta, RoundNearest) / theta
        return RotationVec(re_s * rv.sx, re_s * rv.sy, re_s * rv.sz)
    else
        return rv
    end
end

for rot_type in [:RotX, :RotY, :RotZ]
    @eval begin
        function principal_value(r::$rot_type{T}) where {T}
            return $(rot_type){T}(rem2pi(r.theta, RoundNearest))
        end
    end
end

for rot_type in [:RotXY, :RotYX, :RotZX, :RotXZ, :RotYZ, :RotZY]
    @eval begin
        function principal_value(r::$rot_type{T}) where {T}
            theta1 = rem2pi(r.theta1, RoundNearest)
            theta2 = rem2pi(r.theta2, RoundNearest)
            return $(rot_type){T}(theta1, theta2)
        end
    end
end

for rot_type in [:RotXYX, :RotYXY, :RotZXZ, :RotXZX, :RotYZY, :RotZYZ, :RotXYZ, :RotYXZ, :RotZXY, :RotXZY, :RotYZX, :RotZYX]
    @eval begin
        function principal_value(r::$rot_type{T}) where {T}
            theta1 = rem2pi(r.theta1, RoundNearest)
            theta2 = rem2pi(r.theta2, RoundNearest)
            theta3 = rem2pi(r.theta3, RoundNearest)
            return $(rot_type){T}(theta1, theta2, theta3)
        end
    end
end
