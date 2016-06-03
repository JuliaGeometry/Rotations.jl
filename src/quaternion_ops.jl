#############################################################
#
#  add extra quaternion operations for Quaternions
#  some should probably be merged back into Quaternions
#
#############################################################


"""
function to find the mean of a set of quaternians

Args:

    Qarray         - Vector of quaternions

    method         - if 0 (default), the returned quaternion minimizes sum_i=1:n (||R(Qi) - R(Qbar)||)  with || the frobius norm, or
                                     equivalently minimizes  sum_i=1:n (sin(phi_i / 2))^2, where phi_i is the angle of Qe_i, with  Qbar = Q_i x Qe_i

                   - if 1, the returned quaternian minimizes sum_i=1:n (phi_i ^2) where phi_i is the angle of Qe_i, with  Qbar = Q_i x Qe_i
                           this should be about the same as method == 0 if phi_i are all small
"""
function mean(Qarray::Vector{Quaternion}, method::Integer = 0)

    if (method == 0)
        M = zeros(4, 4)
        for i = 1:length(Qarray)
            Qi = convert(Vector{Float64}, Qarray[i])      # convert types to ensure we don't get quaternion multiplication
            M[:,:] += Qi * (Qi')
        end
        evec = eigfact(Symmetric(M), 4:4)
        Qbar = Quaternion(evec.vectors[1], evec.vectors[2], evec.vectors[3], evec.vectors[4], true)
    else
        error("I haven't coded this")
    end

    return Qbar

end

"""
function to return the derivitive of the rotation matrix w.r.t the quaternion elements

Args:

    q            - a 4 element quaternion

Returns:

    dRdQ         - an array of four 3 x 3 matrices specifying each partial derivitive, i.e. dRdQ[i] = dR/d(Q[i])

"""
function quattorotderiv{T}(q::Quaternion{T})

    # preallocate
    dRdQ = Vector{Mat{3,3,T}}(4)

    # w.r.t q[1]
    dRdQ[1] = 2*@fsa( [q.s     -q.v3    q.v2;
                       q.v3     q.s    -q.v1;
                      -q.v2     q.v1    q.s])

    # w.r.t q[2]
    dRdQ[2] = 2*@fsa([q.v1     q.v2    q.v3;
                      q.v2    -q.v1   -q.s;
                      q.v3     q.s    -q.v1])

    # w.r.t q[3]
    dRdQ[3] = 2*@fsa([-q.v2     q.v1    q.s;
                       q.v1     q.v2    q.v3;
                      -q.s      q.v3   -q.v2])

    # w.r.t q[4]
    dRdQ[4] = 2*@fsa( [-q.v3    -q.s    q.v1;
                       q.s     -q.v3   -q.v2;
                       q.v1     q.v2    q.v3])


    return dRdQ

end




"""

function to return the 2nd derivitives of the rotation matrix w.r.t the quaternion elements

Args:

    q            - a 4 element quaternion.  The 2nd derivitive is constant for all q so thiis
                     input is just given for comptibility

Returns:

    d2RdQ          - an 4 x 4 array of 3 x 3 matrices specifying each partial derivitive,
                     i.e. d2RdQ[i,j] = d2R/(d(Q[i])d(Q[j]))

Notes:

    Returned matrix is stored compactly (i.e. modifying d2RdQ[i,j] will change d2RdQ[j,i])

"""
function quattorotderiv2{T}(q::Quaternion{T}=Quaternion(1.0))

    # preallocate
    d2RdQ = Array(Mat{3,3,T}, 4, 4)

    ######################################
    # Partials involving Q[1]
    ######################################


    # dR / dQ[1]dQ[1]
    d2RdQ[1,1] = 2*@fsa([1 0 0; 0 1 0; 0 0 1])

    # dR / dQ[1]dQ[2]
    d2RdQ[1,2] = d2RdQ[2,1] = 2*@fsa([0          0         0;
                                      0          0        -1;
                                      0          1         0])

    # dR / dQ[1]dQ[2]
    d2RdQ[1,3] = d2RdQ[3,1] = 2*@fsa([0          0         1;
                                      0          0         0;
                                     -1          0         0])

    # dR / dQ[1]dQ[2]
    d2RdQ[1,4] = d2RdQ[4,1] = 2*@fsa([0         -1        0;
                                      1          0         0;
                                      0          0         0])

    ######################################
    # Partials involving Q[2]
    ######################################


    # dR / dQ[2]dQ[2]
    d2RdQ[2,2] =           2*@fsa([1          0         0;
                                   0         -1         0;
                                   0          0        -1])

    # dR / dQ[2]dQ[3]
    d2RdQ[2,3] = d2RdQ[3,2] = 2*@fsa([0          1         0;
                                      1          0         0;
                                      0          0         0])

    # dR / dQ[2]dQ[4]
    d2RdQ[2,4] = d2RdQ[4,2] = 2*@fsa([0          0         0;
                                      0          0         0;
                                      1          0         0])



    ######################################
    # Partials involving Q[3]
    ######################################


    # dR / dQ[3]dQ[3]
    d2RdQ[3,3] =           2*@fsa([-1          0         0;
                                    0          1         0;
                                    0          0        -1])

    # dR / dQ[3]dQ[4]
    d2RdQ[3,4] = d2RdQ[4,3] = 2*@fsa([0          0         0;
                                      0          0         1;
                                      0          1         0])


    ######################################
    # Partials involving Q[4]
    ######################################

    # dR / dQ[4]dQ[4]
    d2RdQ[4,4] =           2 *@fsa([-1          0         0;
                                     0         -1         0;
                                     0          0         1])

    return d2RdQ

end
