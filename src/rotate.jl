
#
# definition for rotating points
#


"""
Rotate the point X using the rotation parameters specified in R (i.e. R * X)

Xout = rotate(R, Xin)

"""
@inline rotate{T}(R::RotMatrix{T}, X::FixedVector{3,T}) =   R * X # rotate via the rotation matrix

# remove this definition when the FixedSizeArrays update makes a tagged version
@inline rotate{T1, T2}(R::RotMatrix{T1}, X::FixedVector{3,T2}) = convert(RotMatrix{promote_type(T1,T2)}, R) * convert(Vec{3, promote_type(T1,T2)}, X)

# Using Quaternions to do the rotation
function rotate(q::Quaternion, X::FixedVector{3})
    # qo = (q*Quaternion(X)*conj(q)) / norm(q)  # expand the above for speed
    s = norm(q)     # SpQuats wouldn't need this step
    qo = (-q.v1 * X[1] - q.v2 * X[2] - q.v3 * X[3],
          q.s   * X[1] + q.v2 * X[3] - q.v3 * X[2],
          q.s   * X[2] - q.v1 * X[3] + q.v3 * X[1],
          q.s   * X[3] + q.v1 * X[2] - q.v2 * X[1])

    Xo = Vec((-qo[1] * q.v1 + qo[2] * q.s  - qo[3] * q.v3 + qo[4] * q.v2) / s,
             (-qo[1] * q.v2 + qo[2] * q.v3 + qo[3] * q.s  - qo[4] * q.v1) / s,
             (-qo[1] * q.v3 - qo[2] * q.v2 + qo[3] * q.v1 + qo[4] * q.s)  / s)

end

# Using Rodrigues formula on an AngleAxis parameterization (assume unit axis length) to do the rotation
# (implementation from: https://ceres-solver.googlesource.com/ceres-solver/+/1.10.0/include/ceres/rotation.h)
function rotate(aa::AngleAxis, X::FixedVector{3})
    w = rotation_axis(aa)
    ct, st = cos(aa.theta), sin(aa.theta)
    w_cross_pt = cross(w, X)
    m = dot(X, w) * (one(w_cross_pt[1]) - ct)
    Xo = Vec(X[1] * ct + w_cross_pt[1] * st + w[1] * m,
             X[2] * ct + w_cross_pt[2] * st + w[2] * m,
             X[3] * ct + w_cross_pt[3] * st + w[3] * m)
    return Xo
end

# Rodrigues formula starting with a Rodrigues vector
function rotate{T}(rv::RodriguesVec{T}, X::FixedVector{3})
    theta = norm(rv)
    if (theta > eps(T)) # use eps here because we have the series expansion defined
        Xo = rotate(AngleAxis(rv), X)
    else
        w_cross_pt = Vec(rv[2] * X[3] - rv[3] * X[2],
                         rv[3] * X[1] - rv[1] * X[3],
                         rv[1] * X[2] - rv[2] * X[1])
        Xo = X + w_cross_pt
    end
end


# and build the rest (having trouble overloading with the RotationTypes union)
for rType in setdiff(RotTypeList, [RotMatrix, Quaternion, AngleAxis, RodriguesVec])
    qb = quote
        rotate{rT <: $(rType), T <: Real}(Rot::rT, X::Vec{3,T}) = rotate($(conversion_path[rType]){promote_type(T, eltype(rT))}(Rot), X)
    end
    eval(qb)
end

#
# in case of mutable input (not the recomended form)
#
@inline rotate{rType <: RotationTypes, T <: Real}(R::rType, X::Vector{T}) =
    Array(rotate(R, Vec{3, promote_type(T, eltype(rType))}(X[1], X[2], X[3])))   # cast to fixed size array

