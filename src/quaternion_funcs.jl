#############################################################
#
#  add extra quaternion operations for Quaternions
#  some should probably be merged back into Quaternions
#
#############################################################


"""
function to find the mean of a set of quaternians

Args:

    Qarray         - A vector of quaternions

    method         - if 0 (default), the returned quaternion minimizes sum_i=1:n (||R(Qi) - R(Qbar)||)  with || the frobius norm, or
                                     equivalently minimizes  sum_i=1:n (sin(phi_i / 2))^2, where phi_i is the angle of Qe_i, with  Qbar = Q_i x Qe_i

                   - if 1, the returned quaternian minimizes sum_i=1:n (phi_i ^2) where phi_i is the angle of Qe_i, with  Qbar = Q_i x Qe_i
                           this should be about the same as method == 0 if phi_i are all small
"""
function mean(Qarray::Vector{Quaternion}, method::Integer = 0)

    if (method == 0)
        M = zeros(4, 4)
        for i = 1:length(Qarray)
            Qc = Qarray[i] * Qarray[i]'
            for j = 1:16
                @inbounds M[j] += Qc[j]  # N.B. broadcast seems bad
            end
        end
        evec = eigfact!(Symmetric(M), 4:4)
        Qbar = Quaternion(evec.vectors[1], evec.vectors[2], evec.vectors[3], evec.vectors[4], true)
    else
        error("I haven't coded this")
    end

    return Qbar

end

