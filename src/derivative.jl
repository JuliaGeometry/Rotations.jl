#
# Module to define Jacobian and hessian functions
# all of these should match the implemented conversion
#


# get the partial derivative of the ith input parameter
ith_partial{N}(X::Mat{3,N}, i) = Vec(X[1,i], X[2,i], X[3,i])

ith_partial{N}(X::Mat{4,N}, i) = Vec(X[1,i], X[2,i], X[3,i], X[4,i])

# reformat to produce the usual 3x3 rotation matrix in this case
ith_partial{N}(X::Mat{9,N}, i) = @fsa([X[1,i]   X[4,i]   X[7,i];
                                       X[2,i]   X[5,i]   X[8,i];
                                       X[3,i]   X[6,i]   X[9,i]])


#######################################################
# Jacobians for transforming to / from rotation matrices
# (only to rotation matrix is implemented - the other way seems weird)
#######################################################

"""

1) jacobian(::Type{output_param}, R::input_param)

Returns the jacobian for transforming from the input rotation parameterization to the output parameterization, centered at the value of R.

2) jacobian(R::rotation_type, X::Union{Vec, Vector})

Returns the jacobian for rotating the vector X by R.

"""
function jacobian{T}(::Type{RotMatrix},  q::Quaternion{T})

    # let q = s * qhat where qhat is a unit quaternion and  s is a scalar,
    # then R = RotMatrix(q) = RotMatrix(s * qhat) = s * RotMatrix(qhat)

    # get R(q)
    t = (RotMatrix(q)...)
    R = Mat(t)  # fsa cant do Vec x [1 x n] matrix

    # solve d(s*R)/dQ (because its easy)
    dsRdQ = @fsa([2*q.s    2*q.v1   -2*q.v2   -2*q.v3;
                  2*q.v3   2*q.v2    2*q.v1    2*q.s;
                 -2*q.v2   2*q.v3   -2*q.s     2*q.v1;

                 -2*q.v3   2*q.v2    2*q.v1   -2*q.s;
                  2*q.s   -2*q.v1    2*q.v2   -2*q.v3;
                  2*q.v1   2*q.s     2*q.v3    2*q.v2;

                  2*q.v2   2*q.v3    2*q.s     2*q.v1;
                 -2*q.v1  -2*q.s     2*q.v3    2*q.v2;
                  2*q.s   -2*q.v1   -2*q.v2    2*q.v3])

    # get s and dsdQ
    s = norm(q)
    dsdQ = @fsa([(q.s/s) (q.v1/s) (q.v2/s) (q.v3/s)])

    # now R(q) = (s*R) / s
    # so dR/dQ = (s * d(s*R)/dQ - (s*R) * ds/dQ) / s^2
    #          = (d(s*R)/dQ - R*ds/dQ) / s

    # now R(q) = (R(s*q)) / s   for scalar s, because RotMatrix(s * q) = s * RotMatrix(q)
    #
    # so dR/dQ = (dR(s*q)/dQ*s - R(s*q) * ds/dQ) / s^2
    #          = (dR(s*q)/dQ*s - s*R(q) * ds/dQ) / s^2
    #          = (dR(s*q)/dQ   - R(q) * ds/dQ) / s

    jac = (dsRdQ - R * dsdQ) / s

    # now reformat for output.  TODO: is the the best expression?
    # return Vec{4,Mat{3,3,T}}(ith_partial(jac, 1), ith_partial(jac, 2), ith_partial(jac, 3), ith_partial(jac, 4))

end


# derivatives of R w.r.t a SpQuat
function jacobian(::Type{RotMatrix},  X::SpQuat)

    # get the derivatives of the quaternion w.r.t to the spquat
    dQdX = jacobian(Quaternion,  X)

    # get the derivatives of the rotation matrix w.r.t to the spquat
    dRdQ = jacobian(RotMatrix,  Quaternion(X))

    # and return
    return dRdQ * dQdX

end


#######################################################
# Jacobians for transforming Quaternion <-> SpQuat
#
#######################################################

function jacobian(::Type{Quaternion},  X::SpQuat)

    # differentiating
    # q = Quaternion((1-alpha2) / (alpha2 + 1), 2*X.x / (alpha2 + 1),   2*X.y  / (alpha2 + 1), 2*X.z / (alpha2 + 1), true)
    # q *= sgn(q.s)

    # derivatives of alpha2
    vspq = Vec(X)
    alpha2 = sum(X.x .* X.x + X.y .* X.y + X.z .* X.z)
    dA2dX = 2 * vspq

    # the function we're differentiating forces q.s +ve
    qs = sgn((1-alpha2) / (alpha2 + 1))

    # f = (1-alpha2) / (alpha2 + 1);
    den2 = (alpha2 + 1) * (alpha2 + 1)
    dQ1dX = qs * (-dA2dX * (alpha2 + 1) - dA2dX * (1-alpha2)) / den2

    # do the on diagonal terms
    # f = 2*x / (alpha2 + 1) => g = 2*x, h = alpha2 + 1
    # df / dx = (dg * h - dh * g) / (h^2)
    dQiDi = qs * 2 * ((alpha2 + 1) - dA2dX .* vspq) / den2

    # do the off entries
    # f = 2x / (alpha2 + 1)
    dQxDi = qs * -2 * vspq[1] * dA2dX / den2
    dQyDi = qs * -2 * vspq[2] * dA2dX / den2
    dQzDi = qs * -2 * vspq[3] * dA2dX / den2

    # assemble it all
    dQdX = @fsa([dQ1dX[1]  dQ1dX[2]  dQ1dX[3];
                 dQiDi[1]  dQxDi[2]  dQxDi[3];
                 dQyDi[1]  dQiDi[2]  dQyDi[3];
                 dQzDi[1]  dQzDi[2]  dQiDi[3]
               ])


    return dQdX
end


#
# Jacobian converting from a Quaternion to an SpQuat
#
function jacobian{T}(::Type{SpQuat},  q::Quaternion{T})

    s = sgn(q.s)
    den = (1 + s * q.s)
    dalpha2dQs = (-s * den - (1 - s*q.s)*s) / (den * den)

    dSpqdQs = s/2 * Vec(q.v1, q.v2, q.v3) * dalpha2dQs

    alpha2 = (1 - s * q.s) / (1 + s * q.s)
    dSpqdQv = s/2 * (alpha2 + 1)

    @fsa([dSpqdQs[1]  dSpqdQv  zero(T)   zero(T);
          dSpqdQs[2]  zero(T)  dSpqdQv   zero(T);
          dSpqdQs[3]  zero(T)  zero(T)   dSpqdQv])
end


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

2) hessian(R::rotation_type, X::Union{Vec, Vector})

Returns the 2nd order partial derivatives for rotating the vector X by R.
The output is an 3 vector of DxD matrices, where D is the number of parameters of the rotation parameterization.

"""
function hessian(::Type{Quaternion},  X::SpQuat)

    # make it match q = Quaternion(X) which puts the return in the domain with q.s >= 0
    q = spquat_to_quat_naive(X)
    s = sgn(q.s)

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
    vspq = Vec(X)
    mag = sum(X.x .* X.x + X.y .* X.y + X.z .* X.z)
    A = 1 + mag
    B = 1 - mag
    dAdX = 2 * vspq
    dBdX = -2 * vspq
    C = dBdX * A - B * dAdX


    # dC/dxx terms (middle two terms cancel)
    T = eltype(A)
    dCdxx = eye(Mat{3,3,T}) * (-2 * (A + B))

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
    hess = Vec{4, Mat{3,3,T}}(

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
    Xim = Vec(s * X.v1, s * X.v2, s * X.v3)

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

    hess = Vec{3, Mat{4,4,T}}(

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

    @fsa([const_q.s   -const_q.v1     -const_q.v2    -const_q.v3;
          const_q.v1   const_q.s      -const_q.v3     const_q.v2;
          const_q.v2   const_q.v3      const_q.s     -const_q.v1;
          const_q.v3  -const_q.v2      const_q.v1     const_q.s])

end

function jacobian(variable_q::Quaternion, const_q::Quaternion, right_variable::Type{Val{false}})

    @fsa([ const_q.s   -const_q.v1     -const_q.v2    -const_q.v3;
           const_q.v1   const_q.s       const_q.v3    -const_q.v2;
           const_q.v2  -const_q.v3      const_q.s      const_q.v1;
           const_q.v3   const_q.v2     -const_q.v1     const_q.s])

end


#######################################################
# Jacobian for rotating a 3 vectors
#
#######################################################

# cross product derivative (w.r.t the right hand side)
d_cross{T}(u::Vec{3,T}) = @fsa([ zero(T) -u[3]      u[2];
                                 u[3]     zero(T)  -u[1];
                                -u[2]     u[1]      zero(T)])

# TODO: should this be jacobian(:rotate, q,  X)   # or something?
function jacobian(q::Quaternion, X::Union{FixedVector{3}, AbstractVector})

    s = norm(q)
    T = eltype(s)

    # derivatives ignoring the scaling
    q_im = Vec(2*q.v1, 2*q.v2, 2*q.v3)
    dRdQr  = Tuple(2 * q.s * X + cross(q_im, X))
    dRdQim = Tuple(-X * q_im' + dot(X, q_im) * eye(Mat{3,3,T}) + q_im * X' - 2*q.s * d_cross(X))
    dRdQs = Mat{3,4, eltype(s)}(dRdQr, dRdQim[1], dRdQim[2], dRdQim[3])

    # include normalization
    dSdQ = Vec(q) / sqrt(s)     # h(x)

    # and finalize with the quotient rule
    Xo = rotate(q, X)           # N.B. g(x) = s * Xo, with dG/dx = dRdQs
    Xom = Mat{3,4,T}(Tuple(dSdQ[1] * Xo), Tuple(dSdQ[2] * Xo), Tuple(dSdQ[3] * Xo), Tuple(dSdQ[4] * Xo)) # Trouble with FSA: outer product of Xo and dSdQ
    dRdQ = (dRdQs -  Xom) / s

end

function jacobian(spq::SpQuat, X::Union{FixedVector{3}, AbstractVector})
    dQ = jacobian(Quaternion, spq)
    q = Quaternion(spq)
    jac = jacobian(q, X) * dQ
end


#######################################################
# Hessian for rotating a 3 vector
#
#######################################################

function hessian(q::Quaternion, X::Union{FixedVector{3}, AbstractVector})

    s = norm(q)
    T = typeof(s)
    z = zero(T)

    # f(x)
    Xo = rotate(q, X)

    #
    # first derivative without removing the scaling (g(x) = s * f(x))
    #
    q_im = Vec(2*q.v1, 2*q.v2, 2*q.v3)
    dRdQr  = 2 * q.s * X + cross(q_im, X)
    dRdQim = -X * q_im' + dot(X, q_im) * eye(Mat{3,3,T}) + q_im * X' - 2* q.s * d_cross(X)

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
    dSdQ = Vec(q) / s
    dSdQt = dSdQ'

    s2 = s*s
    d2SdQ = (eye(Mat{4,4,T}) * s -  Vec(q) * dSdQ') / (s2)

    # and combine them
    dSdQ_2 = 2 * dSdQ * dSdQt
    gd1 = Vec(dRdQr[1], dRdQim[1, 1], dRdQim[1, 2], dRdQim[1, 3])
    gd2 = Vec(dRdQr[2], dRdQim[2, 1], dRdQim[2, 2], dRdQim[2, 3])
    gd3 = Vec(dRdQr[3], dRdQim[3, 1], dRdQim[3, 2], dRdQim[3, 3])

    hess = Vec{3, Mat{4,4,T}}(

        # d2X[1]
        Tuple((d2R1s - Xo[1]*d2SdQ)/s + (Xo[1] * dSdQ_2 - gd1 * dSdQt - dSdQ * gd1')/s2),

        # d2X[2]
        Tuple((d2R2s - Xo[2]*d2SdQ)/s + (Xo[2] * dSdQ_2 - gd2 * dSdQt - dSdQ * gd2')/s2),

        # d2X[3]
        Tuple((d2R3s - Xo[3]*d2SdQ)/s + (Xo[3] * dSdQ_2 - gd3 * dSdQt - dSdQ * gd3')/s2)

        )
end

function hessian(spq::SpQuat, X::Union{FixedVector{3}, AbstractVector})

    # converting to a Quaternion
    j1 = jacobian(Quaternion, spq)
    j1t = j1'
    h1 = hessian(Quaternion, spq)

    # and rotating
    q = Quaternion(spq)
    j2 = jacobian(q, X)
    h2 = hessian(q, X)

    # build them
    hess = Vec{3, Mat{3, 3, eltype(j2)}}(

        # d2X[1]
        Tuple((j2[1,1] * h1[1] + j2[1,2] * h1[2] + j2[1,3] * h1[3] + j2[1,4] * h1[4]) + j1t * h2[1] * j1),

        # d2X[2]
        Tuple((j2[2,1] * h1[1] + j2[2,2] * h1[2] + j2[2,3] * h1[3] + j2[2,4] * h1[4]) + j1t * h2[2] * j1),

        # d2X[3]
        Tuple((j2[3,1] * h1[1] + j2[3,2] * h1[2] + j2[3,3] * h1[3] + j2[3,4] * h1[4]) + j1t * h2[3] * j1)

        )

end

