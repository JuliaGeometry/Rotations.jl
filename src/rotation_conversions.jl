import Base: convert



# a sgn function
@inline sgn{T}(x::T) = copysign(one(T), x)  # N.B. sign(0.0) != sgn(0.0). sgn(-0.0) == -1 isn't what I want though

"""
wraps theta into the range (-pi, pi]
"""
@inline wrap_angle(theta::Real) = theta + floor((pi-theta) / (2*pi)) * 2*pi # (-pi, pi]


#
# Add all of the transformations from one parameterization to another
#
function build_transform_paths()

    # build transformations from one type to another
    for rT in keys(conversion_path)          # rT is the type to transform from
        for oT in keys(conversion_path)  # oT is the type to transform to
            if (oT != rT) && !any(def_conv -> (def_conv[1] == rT) && (def_conv[2] == oT), defined_conversions)
                inter_conv = conversion_path[rT] == nothing ? conversion_path[oT] : conversion_path[rT]
                qb = quote
                    @inline convert{T <: $(oT)}(::Type{T}, X::$(rT)) = convert(T, convert($(inter_conv), X))  # go via trans_dict[rT]
                end
                #println(qb)
                eval(qb)
            end
        end
    end

    # make sure RotMatrix has call mapped to convert for rotation types
    for rT in keys(conversion_path)
        if rT != RotMatrix
            qb = quote
                @compat @inline (::Type{T}){T <: RotMatrix}(X::$(rT)) = convert(T, X)
            end
            eval(qb)
        end
    end
end
build_transform_paths()


#######################################
# implementations
#######################################

#=
function quat_to_rot(q::Quaternion)  # this version assumes q is normalized

    # get rotation matrix from quaternion
    xx = q.v1 * q.v1
    yy = q.v2 * q.v2
    zz = q.v3 * q.v3
    xy = q.v1 * q.v2
    zw = q.s  * q.v3
    xz = q.v1 * q.v3
    yw = q.v2 * q.s
    yz = q.v2 * q.v3
    xw = q.s  * q.v1

    # initialize rotation part
    return @fsa([1 - 2 * (yy + zz)    2 * (xy - zw)       2 * (xz + yw);
                 2 * (xy + zw)        1 - 2 * (xx + zz)   2 * (yz - xw);
                 2 * (xz - yw)        2 * (yz + xw)       1 - 2 * (xx + yy)])

end
=#


"""
function to convert the input quaternion to a rotation matrix
"""
function quat_to_rot(q::Quaternion)

    ww = (q.s  * q.s)
    xx = (q.v1 * q.v1)
    yy = (q.v2 * q.v2)
    zz = (q.v3 * q.v3)

    mi = 1/sqrt(ww + xx + yy + zz) # norm

    # get rotation matrix from quaternion
    ww *= mi
    xx *= mi
    yy *= mi
    zz *= mi
    xy = (q.v1 * q.v2) * mi
    zw = (q.s  * q.v3) * mi
    xz = (q.v1 * q.v3) * mi
    yw = (q.v2 * q.s)  * mi
    yz = (q.v2 * q.v3) * mi
    xw = (q.s  * q.v1) * mi

    # initialize rotation part
    return @fsa([ww + xx - yy - zz    2 * (xy - zw)       2 * (xz + yw);
                 2 * (xy + zw)        ww - xx + yy - zz   2 * (yz - xw);
                 2 * (xz - yw)        2 * (yz + xw)       ww - xx - yy + zz])

end

"""
 a mutable version filling the 3x3 buffer R using the 4 element vector q

 quat_to_rot!(R, q)

"""
function quat_to_rot!(R, q)

    # get rotation matrix from quaternion
    ww = q[1] * q[1]
    xx = q[2] * q[2]
    yy = q[3] * q[3]
    zz = q[4] * q[4]
    xy = q[2] * q[3]
    zw = q[1] * q[4]
    xz = q[2] * q[4]
    yw = q[3] * q[1]
    yz = q[3] * q[4]
    xw = q[1] * q[2]

    # fill the rotation matrix
    R[1,1], R[1,2], R[1,3] = ww + xx - yy - zz,        2 * (xy - zw),       2 * (xz + yw)
    R[2,1], R[2,2], R[2,3] =     2 * (xy + zw),    ww - xx + yy - zz,       2 * (yz - xw)
    R[3,1], R[3,2], R[3,3] =     2 * (xz - yw),        2 * (yz + xw),   ww - xx - yy + zz


end



"""
function to convert the input rotation matrix into a quaternion.
"""
function rot_to_quat(R::RotMatrix)
    q = Quaternion(sqrt(abs(1  + R[1,1] + R[2,2] + R[3,3])) / 2,
                   copysign(sqrt(abs(1 + R[1,1] - R[2,2] - R[3,3]))/2, R[3,2] - R[2,3]),
                   copysign(sqrt(abs(1 - R[1,1] + R[2,2] - R[3,3]))/2, R[1,3] - R[3,1]),
                   copysign(sqrt(abs(1 - R[1,1] - R[2,2] + R[3,3]))/2, R[2,1] - R[1,2]),
                   true)

    return q

end

"""
 a mutable version filling the 4 element buffer q using the 3x3 matrix R

 rot_to_quat!(q, R)

"""
function rot_to_quat!(q, R)
    q[1] = sqrt(abs(1  + R[1,1] + R[2,2] + R[3,3])) / 2
    q[2] = copysign(sqrt(abs(1 + R[1,1] - R[2,2] - R[3,3]))/2, R[3,2] - R[2,3])
    q[3] = copysign(sqrt(abs(1 - R[1,1] + R[2,2] - R[3,3]))/2, R[1,3] - R[3,1])
    q[4] = copysign(sqrt(abs(1 - R[1,1] - R[2,2] + R[3,3]))/2, R[2,1] - R[1,2])
end

"""
function to find the stereographic projection of a quaternion (see ?SpQuat).
"""
function quat_to_spquat(q::Quaternion)

    # N.B. there's a singularity at q = [-1.0,0,0,0] - avoid it
    # alpha2 = (q.s > -.99999999) ? (1 - q.s) / (1 + q.s) : 0.0
    s = sgn(q.s)  # this way will keep the 2norm of the output spquat <= 1, I think they're better behaved in this range.  Also avoids the singularity
    alpha2 = (1 - s * q.s) / (1 + s * q.s)
    spq = SpQuat(s * q.v1 * (alpha2 + 1)/2,  s * q.v2 * (alpha2 + 1)/2, s * q.v3 * (alpha2 + 1)/2)
end
function quat_to_spquat_naive(q::Quaternion) # makes no attempt to stay away from the singularity
    alpha2 = (1 - q.s) / (1 + q.s)
    spq = SpQuat(q.v1 * (alpha2 + 1)/2,  q.v2 * (alpha2 + 1)/2, q.v3 * (alpha2 + 1)/2)
end

"""
 a mutable version filling the 3 element buffer spq using the 4 element vector q

 quat_to_spquat!(spq, q)

"""
function quat_to_spquat!(spq, q)
    qs = sgn(q[1])  # this way will keep the 2norm of the output spquat <= 1, I think they're better behaved in this range.  Also avoids the singularity
    alpha2 = (1 - qs * q[1]) / (1 + qs * q[1])
    spq[1], spq[2], spq[3]  = qs * q[2] * (alpha2 + 1)/2,  qs * q[3] * (alpha2 + 1)/2, qs * q[4] * (alpha2 + 1)/2
end


"""
function to convert a the stereographic projection of a unit quaternion back into a quaternion (see ?SpQuat).
"""
function spquat_to_quat(spq::SpQuat)
    alpha2 = sum(spq.x .* spq.x + spq.y .* spq.y + spq.z .* spq.z)
    q = Quaternion((1-alpha2) / (alpha2 + 1), 2*spq.x   / (alpha2 + 1),   2*spq.y  / (alpha2 + 1), 2*spq.z / (alpha2 + 1), true)
    q *= sgn(q.s)
    return q
end
function spquat_to_quat_naive(spq::SpQuat) # makes no attempt to stay away from the singularity
    alpha2 = sum(spq.x .* spq.x + spq.y .* spq.y + spq.z .* spq.z)
    q = Quaternion((1-alpha2) / (alpha2 + 1), 2*spq.x   / (alpha2 + 1),   2*spq.y  / (alpha2 + 1), 2*spq.z / (alpha2 + 1), true)
end

"""
 a mutable version filling the 4 element vector q using the 3 element buffer spq

 spquat_to_quat!(q, spq)

"""
function spquat_to_quat!(q, spq)
    alpha2 = spq[1] * spq[1] + spq[2] * spq[2] + spq[3] * spq[3]
    q[1] = (1 - alpha2) / (alpha2 + 1)
    qs = sgn(q[1])
    q[1] *= qs
    q[2], q[3], q[4] = qs * 2 * spq[1] / (alpha2 + 1),   qs * 2 * spq[2]  / (alpha2 + 1),  qs * 2 * spq[3] / (alpha2 + 1)
end


# define an spquat -> rot because its common
"""
 a mutable version filling the 3 x 3 buffer R using the 3 element buffer spq

 sqquat_to_rot!(R, spq)

"""
function spquat_to_rot!(R, spq)

    # use the R buffer to store the quaternion
    spquat_to_quat!(R, spq)
    q = (R[1], R[2], R[3], R[4])  # and grab it

    # now fill R using the quaternion we just grabbed
    quat_to_rot!(R, q)

end




"""
function to convert an arbitrary axis rotation to a quaternion
"""
function angleaxis_to_quat(ar::AngleAxis)
    qtheta = cos(ar.theta / 2)
    s = sin(ar.theta / 2) / sqrt(sum(ar.axis_x * ar.axis_x + ar.axis_y * ar.axis_y + ar.axis_z * ar.axis_z))
    return Quaternion(qtheta, s * ar.axis_x, s * ar.axis_y, s * ar.axis_z)
end

"""
function to convert a quaternion to an arbitrary axis rotation
"""
function quat_to_angleaxis{T}(q::Quaternion{T})                        # TODO: make this differentiable by making a series expansion near s = 0
    s = sqrt(q.v1*q.v1 + q.v2*q.v2 + q.v3*q.v3)
    theta =  2 * atan2(s, q.s)
    return s > 0 ? AngleAxis(theta, q.v1 / s, q.v2 / s, q.v3 / s) : AngleAxis(theta, one(theta), zero(theta), zero(theta))
end

"""
function to convert a rodrigues vector to an angle axis representation
"""
function rodrigues_to_angleaxis{T}(rv::RodriguesVec{T})                 # TODO: make this differentiable by making a series expansion near theta = 0
    theta = norm(rv)
    return theta > 0 ? AngleAxis(theta, rv.sx / theta, rv.sy / theta, rv.sz / theta) : AngleAxis(zero(theta), one(theta), zero(theta), zero(theta))
end

"""
function to convert a rodrigues vector to a quaternion
"""
function rodrigues_to_quat{T}(rv::RodriguesVec{T})
    theta = norm(rv)
    qtheta = cos(theta / 2)
    s = 1/2 * sinc((theta / 2) / pi) # Rework - autodiff can't handle sinc
    return Quaternion(qtheta, s * rv.sx, s * rv.sy, s * rv.sz)
end

"""
function to convert a a quaternion to a rodrigues vector
"""
function quat_to_rodrigues{T}(q::Quaternion{T})
    s2 = q.v1*q.v1 + q.v2*q.v2 + q.v3*q.v3
    if (s2 > 0)
        cos_t2 = sqrt(s2)
        theta = 2 * atan2(cos_t2, q.s)
        sc = theta / cos_t2
    else
        sc = 2                 # N.B. the 2 "should" match the derivitive as cos_t2 -> 0
    end
    return RodriguesVec(sc * q.v1, sc * q.v2, sc * q.v3 )
end

"""
function to convert an angle axis format to a rodrigues vector
"""
function angleaxis_to_rodrigues(aa::AngleAxis)
    s = aa.theta / sqrt(aa.axis_x * aa.axis_x + aa.axis_y * aa.axis_y + aa.axis_z * aa.axis_z)
    return RodriguesVec(s * aa.axis_x, s * aa.axis_y, s * aa.axis_z)
end

"""
Function to test if R is a valid rotation matrix.  This checks I - R * R'
"""
function valid_rotation(R, thresh::Real=1e-9)

    if !isa(R, Mat)
        R = convert_rotation(RotMatrix, R)
    end
    res = R * R'
    et = eltype(res)
    dX = res - Mat{size(res,1), size(res,1), et}(eye(et, size(res,1)))
    return vecnorm(dX) < thresh

end

