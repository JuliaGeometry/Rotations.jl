#
# Module to define Jacobian and hessian functions
# all of these should match the implemented conversion
#


#=
# get the partial derivative of the ith input parameter
ith_partial{N}(X::SMatrix{3,N}, i) = SVector(X[1,i], X[2,i], X[3,i])

ith_partial{N}(X::SMatrix{4,N}, i) = SVector(X[1,i], X[2,i], X[3,i], X[4,i])

# reformat to produce the usual 3x3 rotation matrix in this case
ith_partial{N}(X::SMatrix{9,N}, i) = @SMatrix([X[1,i]   X[4,i]   X[7,i];
                                               X[2,i]   X[5,i]   X[8,i];
                                               X[3,i]   X[6,i]   X[9,i]])
=#

#######################################################
# Jacobians for transforming to / from rotation matrices
# (only to rotation matrix is implemented - the other way seems weird)
#######################################################


"""
    jacobian(::Type{output_param}, R::input_param)
Returns the jacobian for transforming from the input rotation parameterization
to the output parameterization, centered at the value of R.

    jacobian(R::rotation_type, X::AbstractVector)
Returns the jacobian for rotating the vector X by R.
"""
function jacobian(::Type{RotMatrix},  q::UnitQuaternion)

    # let q = s * qhat where qhat is a unit quaternion and  s is a scalar,
    # then R = RotMatrix(q) = RotMatrix(s * qhat) = s * RotMatrix(qhat)

    # get R(q)
    # R = q[:] # FIXME: broken with StaticArrays 0.4.0 due to https://github.com/JuliaArrays/StaticArrays.jl/issues/128
    R = SVector(Tuple(q))

    # solve d(s*R)/dQ (because its easy)
    dsRdQ = @SMatrix [ 2*q.w   2*q.x   -2*q.y   -2*q.z ;
                       2*q.z   2*q.y    2*q.x    2*q.w ;
                      -2*q.y   2*q.z   -2*q.w    2*q.x ;

                      -2*q.z   2*q.y    2*q.x   -2*q.w ;
                       2*q.w  -2*q.x    2*q.y   -2*q.z ;
                       2*q.x   2*q.w    2*q.z    2*q.y ;

                       2*q.y   2*q.z    2*q.w    2*q.x ;
                      -2*q.x  -2*q.w    2*q.z    2*q.y ;
                       2*q.w  -2*q.x   -2*q.y    2*q.z ]

    # get s and dsdQ
    # s = 1
    dsdQ = @SVector [2*q.w, 2*q.x, 2*q.y, 2*q.z]

    # now R(q) = (s*R) / s
    # so dR/dQ = (s * d(s*R)/dQ - (s*R) * ds/dQ) / s^2
    #          = (d(s*R)/dQ - R*ds/dQ) / s

    # now R(q) = (R(s*q)) / s   for scalar s, because RotMatrix(s * q) = s * RotMatrix(q)
    #
    # so dR/dQ = (dR(s*q)/dQ*s - R(s*q) * ds/dQ) / s^2
    #          = (dR(s*q)/dQ*s - s*R(q) * ds/dQ) / s^2
    #          = (dR(s*q)/dQ   - R(q) * ds/dQ) / s

    jac = dsRdQ - R * transpose(dsdQ)

    # now reformat for output.  TODO: is the the best expression?
    # return Vec{4,Mat{3,3,T}}(ith_partial(jac, 1), ith_partial(jac, 2), ith_partial(jac, 3), ith_partial(jac, 4))

end


# derivatives of R w.r.t a MRP
function jacobian(::Type{RotMatrix},  X::MRP)

    # get the derivatives of the quaternion w.r.t to the mrp
    dQdX = jacobian(UnitQuaternion,  X)

    # get the derivatives of the rotation matrix w.r.t to the mrp
    dRdQ = jacobian(RotMatrix,  UnitQuaternion(X))

    # and return
    return dRdQ * dQdX

end


#######################################################
# Jacobians for transforming Quaternion <-> MRP
#
#######################################################

function jacobian(::Type{UnitQuaternion},  X::MRP)

    # differentiating
    # q = Quaternion((1-alpha2) / (alpha2 + 1), 2*X.x / (alpha2 + 1),   2*X.y  / (alpha2 + 1), 2*X.z / (alpha2 + 1), true)
    # q *= sgn(q.w)

    # derivatives of alpha2
    vspq = SVector(X.x, X.y, X.z)
    alpha2 = X.x * X.x + X.y * X.y + X.z * X.z
    dA2dX = 2 * vspq

    # f = (1-alpha2) / (alpha2 + 1);
    den2 = (alpha2 + 1) * (alpha2 + 1)
    dQ1dX = (-dA2dX * (alpha2 + 1) - dA2dX * (1-alpha2)) / den2

    # do the on diagonal terms
    # f = 2*x / (alpha2 + 1) => g = 2*x, h = alpha2 + 1
    # df / dx = (dg * h - dh * g) / (h^2)
    dQiDi = 2 * ((alpha2 + 1) .- dA2dX .* vspq) / den2

    # do the off entries
    # f = 2x / (alpha2 + 1)
    dQxDi = -2 * vspq[1] * dA2dX / den2
    dQyDi = -2 * vspq[2] * dA2dX / den2
    dQzDi = -2 * vspq[3] * dA2dX / den2

    # assemble it all
    dQdX = @SMatrix [ dQ1dX[1]  dQ1dX[2]  dQ1dX[3] ;
                      dQiDi[1]  dQxDi[2]  dQxDi[3] ;
                      dQyDi[1]  dQiDi[2]  dQyDi[3] ;
                      dQzDi[1]  dQzDi[2]  dQiDi[3] ]

    return dQdX
end


#
# Jacobian converting from a Quaternion to an MRP
#
function jacobian(::Type{MRP}, q::UnitQuaternion{T}) where T
    den = 1 + q.w
    scale = 1 / den
    dscaledQw = -(scale * scale)
    dSpqdQw = SVector(q.x, q.y, q.z) * dscaledQw
    J0 = @SMatrix [ dSpqdQw[1]  scale   zero(T) zero(T) ;
                    dSpqdQw[2]  zero(T) scale   zero(T) ;
                    dSpqdQw[3]  zero(T) zero(T) scale ]

    # Need to project out norm component of Quat
    dQ = @SVector [q.w, q.x, q.y, q.z]
    return J0 - (J0*dQ)*dQ'
end



#######################################################
# Jacobian for rotating a 3 vectors
#
#######################################################

# Note: this is *not* projected into the orthogonal matrix tangent space.
# can do this by projecting each 3x3 matrix (row of 9) by (jacobian[i] - r * jacabian[i]' * r) / 2   (for i = 1:3)
function jacobian(r::RotMatrix{3}, X::AbstractVector)
    @assert length(X) === 3
    T = promote_type(eltype(r), eltype(X))
    Z = zero(T)

    @inbounds return @SMatrix [ X[1] Z    Z     X[2]  Z     Z     X[3]  Z     Z    ;
                                Z    X[1] Z     Z     X[2]  Z     Z     X[3]  Z    ;
                                Z    Z    X[1]  Z     Z     X[2]  Z     Z     X[3] ]
end

@inline function d_cross(u::AbstractVector{T}) where T
    @assert length(u) === 3
    @inbounds return @SMatrix [ zero(T)  -u[3]      u[2]    ;
                                u[3]      zero(T)  -u[1]    ;
                               -u[2]      u[1]      zero(T) ]
end

# TODO: should this be jacobian(:rotate, q,  X)   # or something?
function jacobian(q::UnitQuaternion, X::AbstractVector)
    @assert length(X) === 3
    T = eltype(q)

    # derivatives ignoring the scaling
    q_im = SVector(2*q.x, 2*q.y, 2*q.z)
    dRdQr  = SVector{3}(2 * q.w * X + cross(q_im, X))
    dRdQim = -X * q_im' + dot(X, q_im) * one(SMatrix{3,3,T}) + q_im * X' - 2*q.w * d_cross(X)

    dRdQs = hcat(dRdQr, dRdQim)

    # include normalization (S, s = norm of quaternion)
    dSdQ = SVector(2*q.w, 2*q.x, 2*q.y, 2*q.z)     # h(x)

    # and finalize with the quotient rule
    Xo = q * X           # N.B. g(x) = s * Xo, with dG/dx = dRdQs
    Xom = Xo * transpose(dSdQ)
    return dRdQs -  Xom
end

function jacobian(spq::MRP, X::AbstractVector)
    dQ = jacobian(UnitQuaternion, spq)
    q = UnitQuaternion(spq)
    return jacobian(q, X) * dQ
end



#=
#######################################################
# Hessians for transforming Quaternion <-> SpQuat
#
#######################################################

#
# 2nd derivative of the SpQuat - > Quaternion transformation
#
"""
1) hessian(::Type{output_param}, R::input_param)
Returns the 2nd order partial derivatives for transforming from the input rotation parameterization to the output parameterization, centered at the value of R.
The output is an N vector of DxD matrices, where N and D are the number of parameters in the output and input parameterizations respectively.
2) hessian(R::rotation_type, X::AbstractVector)
Returns the 2nd order partial derivatives for rotating the vector X by R.
The output is an 3 vector of DxD matrices, where D is the number of parameters of the rotation parameterization.
"""
function hessian(::Type{Quaternion},  X::SpQuat)

    # make it match q = Quaternion(X) which puts the return in the domain with q.w >= 0
    q = spquat_to_quat_naive(X)
    s = sgn(q.w)

    # state with the hessian of the first Quaternion term
    # Q[1] = (1-alpha2) / (alpha2 + 1)
    #
    # let A = 1 + alpha2
    # let B = 1 - alpha2
    #
    # dQ1 / dx = (dB/dx * A - B * dA/dx) / (A^2)
    #
    # let C = (dB/dx * A - B * dA/dx)
    # so  dQ1/dx = C / (A^2)
    #
    # then
    # ddQ1/dxdy = (dC/dy*A^2 - C * 2*A * dA/dy) / (A^4)
    #           = (dC/dy*A   - 2 * C * dA/dy) / (A^3)
    #
    # with
    # dC/dy = dB/dxdy * A + dB/dx*dA/dy - dB/dy*dA/dx - B * dA/dxdy

    # calculate C
    vspq = SVector(X)
    mag = sum(X.x .* X.x + X.y .* X.y + X.z .* X.z)
    A = 1 + mag
    B = 1 - mag
    dAdX = 2 * vspq
    dBdX = -2 * vspq
    C = dBdX * A - B * dAdX


    # dC/dxx terms (middle two terms cancel)
    T = eltype(A)
    dCdxx = one(Mat{3,3,T}) * (-2 * (A + B))

    # dC/dxy terms (dB/dxdy and dA/dxdy) are zero
    # put it all together
    A3 = A*A*A
    d2Q1 = s * (dCdxx * A - 2 * C *  dAdX') / A3

    #
    # now "pure" 2nd derivatives
    #
    # Qx = 2*x / (alpha2 + 1) = 2*x/A
    # dQxdx = (2 * A - 2*x*dA/dx) / (A^2)
    #       = C / A^2
    # with C = 2 * A - 2*x*dA/dX
    #
    # dCdx =  2*dA/dX - (2*dA/dX + 2*x*2) = 2*(dA/dX - dAdX - 2*x)
    #      = -4x
    #
    # d2Qxdxx = (dC/dX*A*A - C*2*A*dA/dX) / (A^4)
    #         = (dC/dX*A - 2*C*dA/dX) / (A^3)
    # d2Qxdxx = (-4*A*vspq - 2*(2*A - 2*vspq.*dAdX) .* dAdX) / A3
    #         =  (-4*A*vspq - 4*dAdX .* (A - vspq.*dAdX)) / A3
    d2Qdxx    = s * -4 * (A*vspq + dAdX .* (A - vspq.*dAdX)) / A3

    #
    # now the mixed 2nd derivatives
    #
    # start at
    # Qx    = 2*x / (alpha2 + 1) = 2*x/A
    # dQxdy = (-2*x*dAdy) / (A^2)
    # d2Qxdxy = (-2*dA/dy*A^2 + 2*x*dA/dy*2*A*dA/dx) / (A^4)
    #         = (-2*dA/dy*A + 2*x*dA/dy*2*dA/dx) / (A^3)
    #         = 2*dA/dy*(2*x*dA/dx - A) / (A^3)
    d2Qxdxy    = s * 2 * dAdX * (2*vspq.*dAdX - A)' / A3

    #
    # now the other 2nd derivatives
    #

    # Qx = 2*x / (alpha2 + 1) = 2*x/A
    # dQxdy =  (-2*x*dAdXy) / (A^2)
    # dQxdyz = -(-2*x*dAdXy*2*A*dAdXz) / (A^4)
    #        =  (4*x*dAdXy.*dAdXz) / (A3)
    d2Qxdyz = s * (4 * vspq[1] * dAdX[2] * dAdX[3]) / A3
    d2Qydxz = s * (4 * vspq[2] * dAdX[1] * dAdX[3]) / A3
    d2Qzdxy = s * (4 * vspq[3] * dAdX[1] * dAdX[2]) / A3

    # Qx = 2*x / (alpha2 + 1) = 2*x/A
    # dQxdy =   (-2*x*dAdXy) / (A^2)
    # let C =   -2 * x * dAdXy
    # dCdy  =   -4 * x

    # dQxdyy =  (dCdy * A^2 + 2*x*dAdXy*2*A*dAdXy) / (A^4)
    #        =  (dCdy * A   + 2*x*dAdXy*2*dAdXy) / (A^3)
    #        =  (-4*x * A   + 2*x*dAdXy*2*dAdXy) / (A^3)
    #        =  4*x*(dAdXy^2 - A) / (A^3)
    d2Qxdyy   = s * 4 * vspq * (dAdX.*dAdX - A)' / A3

    #
    # And form it all
    #
    hess = SVector{4, Mat{3,3,T}}(

      Tuple(d2Q1),

      #d2Qx
      Tuple(@fsa([d2Qdxx[1]       d2Qxdxy[2, 1]   d2Qxdxy[3, 1];
                  d2Qxdxy[2, 1]   d2Qxdyy[1,2]    d2Qxdyz;
                  d2Qxdxy[3, 1]   d2Qxdyz         d2Qxdyy[1, 3]])),

      #d2Qy
      Tuple(@fsa([d2Qxdyy[2, 1]   d2Qxdxy[1, 2]   d2Qydxz;
                  d2Qxdxy[1, 2]   d2Qdxx[2]       d2Qxdxy[3, 2];
                  d2Qydxz         d2Qxdxy[3, 2]   d2Qxdyy[2, 3]])),
      #d2Qz
      Tuple(@fsa([d2Qxdyy[3, 1]   d2Qzdxy         d2Qxdxy[1, 3];
                  d2Qzdxy         d2Qxdyy[3, 2]   d2Qxdxy[2, 3];
                  d2Qxdxy[1, 3]   d2Qxdxy[2, 3]   d2Qdxx[3]]))

           )
end

#
# 2nd derivative of the Quaternion -> SpQuat transformation
#
function hessian(::Type{SpQuat},  X::Quaternion)

    # always apply to the quaternion in the domain with the real part >= 0
    s = sgn(X.s)
    Xim = SVector(s * X.v1, s * X.v2, s * X.v3)

    # A = (1 - qs) / (1 + qs)
    # dAdQ1 = 2 / (qs + 1)^2
    # d2AdQ1 = 4 / (qs +1 )^3
    d = (1 + s * X.s); d2 = d*d; d3 = d2*d
    dAdQs = -2 / d2
    d2AdQs = 4 / d3

    # N.B. only partials of the real part are non-zero (i.e. the output is sparse)
    # partials w.r.t the imaginary parts of the quaternion are all the same
    dSdQxQs = dAdQs / 2

    # now do d/dQ1dQ1 terms
    d2Sx_dQsQs = Xim/2  * d2AdQs

    #
    # And build it
    #
    T = typeof(dSdQxQs)
    z = zero(T)

    hess = SVector{3, Mat{4,4,T}}(

        # d2Sx
        Tuple(@fsa([d2Sx_dQsQs[1]  dSdQxQs  z        z;
                    dSdQxQs        z        z        z;
                    z              z        z        z;
                    z              z        z        z])),

        # d2Sy
        Tuple(@fsa([d2Sx_dQsQs[2]  z        dSdQxQs  z;
                    z              z        z        z;
                    dSdQxQs        z        z        z;
                    z              z        z        z])),

        # d2Sz
        Tuple(@fsa([d2Sx_dQsQs[3]  z        z        dSdQxQs;
                    z              z        z        z;
                    z              z        z        z;
                    dSdQxQs        z        z        z]))

        )

end


#######################################################
# Jacobian for rotating one Quaternion by another
#
#######################################################

#
# get the jacobian of quaternion multiplication w.r.t. the right side quaternion
# each column is the ith partial
#
# TODO: should this be jacobian{T}(:*, const_q::Quaternion{T},  variable_q::Quaternion{T})   # or something?
function jacobian(const_q::Quaternion,  variable_q::Quaternion, right_variable::Type{Val{true}}=Val{true})

    @fsa([const_q.w   -const_q.x     -const_q.y    -const_q.z;
          const_q.x   const_q.w      -const_q.z     const_q.y;
          const_q.y   const_q.z      const_q.w     -const_q.x;
          const_q.z  -const_q.y      const_q.x     const_q.w])

end

function jacobian(variable_q::Quaternion, const_q::Quaternion, right_variable::Type{Val{false}})

    @fsa([ const_q.w   -const_q.x     -const_q.y    -const_q.z;
           const_q.x   const_q.w       const_q.z    -const_q.y;
           const_q.y  -const_q.z      const_q.w      const_q.x;
           const_q.z   const_q.y     -const_q.x     const_q.w])

end

=#


#=

#######################################################
# Hessian for rotating a 3 vector
#
#######################################################

function hessian(q::Quaternion, X::AbstractVector)

    s = norm(q)
    T = typeof(s)
    z = zero(T)

    # f(x)
    Xo = rotate(q, X)

    #
    # first derivative without removing the scaling (g(x) = s * f(x))
    #
    q_im = SVector(2*q.x, 2*q.y, 2*q.z)
    dRdQr  = 2 * q.w * X + cross(q_im, X)
    dRdQim = -X * q_im' + dot(X, q_im) * one(Mat{3,3,T}) + q_im * X' - 2* q.w * d_cross(X)

    #
    # second derivative ignoring the scaling
    #
    X2 = 2*X; Xm2 = -X2
    d2dri = d_cross(Xm2)

    # hessian for the X coord
    d2R1s = @fsa([X2[1]        d2dri[1,1]   d2dri[1,2]   d2dri[1,3];
                  d2dri[1,1]   X2[1]        X2[2]        X2[3];
                  d2dri[1,2]   X2[2]        Xm2[1]       z;
                  d2dri[1,3]   X2[3]        z            Xm2[1]      ])

    # hessian for the Y coord
    d2R2s = @fsa([X2[2]        d2dri[2,1]   d2dri[2,2]   d2dri[2,3];
                  d2dri[2,1]   Xm2[2]       X2[1]        z;
                  d2dri[2,2]   X2[1]        X2[2]        X2[3];
                  d2dri[2,3]   z            X2[3]        Xm2[2]      ])


    # hessian for the Z coord
    d2R3s = @fsa([X2[3]        d2dri[3,1]   d2dri[3,2]   d2dri[3,3];
                  d2dri[3,1]   Xm2[3]       z            X2[1];
                  d2dri[3,2]   z            Xm2[3]       X2[2];
                  d2dri[3,3]   X2[1]        X2[2]        X2[3]      ])

    #
    # Now the scaling part, s(x) = sqrt(x' * x)
    #
    dSdQ = SVector(q) / s
    dSdQt = dSdQ'

    s2 = s*s
    d2SdQ = (one(Mat{4,4,T}) * s -  SVector(q) * dSdQ') / (s2)

    # and combine them
    dSdQ_2 = 2 * dSdQ * dSdQt
    gd1 = SVector(dRdQr[1], dRdQim[1, 1], dRdQim[1, 2], dRdQim[1, 3])
    gd2 = SVector(dRdQr[2], dRdQim[2, 1], dRdQim[2, 2], dRdQim[2, 3])
    gd3 = SVector(dRdQr[3], dRdQim[3, 1], dRdQim[3, 2], dRdQim[3, 3])

    hess = SVector{3, Mat{4,4,T}}(

        # d2X[1]
        Tuple((d2R1s - Xo[1]*d2SdQ)/s + (Xo[1] * dSdQ_2 - gd1 * dSdQt - dSdQ * gd1')/s2),

        # d2X[2]
        Tuple((d2R2s - Xo[2]*d2SdQ)/s + (Xo[2] * dSdQ_2 - gd2 * dSdQt - dSdQ * gd2')/s2),

        # d2X[3]
        Tuple((d2R3s - Xo[3]*d2SdQ)/s + (Xo[3] * dSdQ_2 - gd3 * dSdQt - dSdQ * gd3')/s2)

        )
end

function hessian(spq::SpQuat, X::AbstractVector)

    # converting to a Quaternion
    j1 = jacobian(Quaternion, spq)
    j1t = j1'
    h1 = hessian(Quaternion, spq)

    # and rotating
    q = Quaternion(spq)
    j2 = jacobian(q, X)
    h2 = hessian(q, X)

    # build them
    hess = SVector{3, Mat{3, 3, eltype(j2)}}(

        # d2X[1]
        Tuple((j2[1,1] * h1[1] + j2[1,2] * h1[2] + j2[1,3] * h1[3] + j2[1,4] * h1[4]) + j1t * h2[1] * j1),

        # d2X[2]
        Tuple((j2[2,1] * h1[1] + j2[2,2] * h1[2] + j2[2,3] * h1[3] + j2[2,4] * h1[4]) + j1t * h2[2] * j1),

        # d2X[3]
        Tuple((j2[3,1] * h1[1] + j2[3,2] * h1[2] + j2[3,3] * h1[3] + j2[3,4] * h1[4]) + j1t * h2[3] * j1)

        )

end
=#
