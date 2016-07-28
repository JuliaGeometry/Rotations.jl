using Quaternions
using ForwardDiff

thresh = 1e-12

@testset "Derivative checks" begin

    ###################################
    # Test Jacobians
    ###################################

    @testset "Jacobian checks" begin

        # Quaternion to rotation matrix
        @testset "Jacobian (Quaternion -> RotMatrix)" begin

            # transformation function
            function Jr(X)
                R = Tuple(RotMatrix(Quaternion(X...)))
                V = [R[1]..., R[2]..., R[3]...]
            end

            for i = 1:10    # do some repeats

                q = nquatrand()  # a random quaternion

                # test jacobian to a Rotation matrix
                R_jac = Rotations.jacobian(RotMatrix, q)
                FD_jac = ForwardDiff.jacobian(Jr, vec(q))

                # compare
                @test all(abs(Matrix(R_jac) - FD_jac) .< thresh)

            end
        end

        #
        # Quaternion to rotation matrix should also work on non-unit quaternions
        #
        @testset "Jacobian (Non-unit Quaternion -> RotMatrix)" begin

            # transformation function
            function Jr(X)
                R = Tuple(RotMatrix(Quaternion(X...)))
                V = [R[1]..., R[2]..., R[3]...]
            end

            for i = 1:10    # do some repeats

                q = quatrand()  # a random non-unit quaternion

                # test jacobian to a Rotation matrix
                R_jac = Rotations.jacobian(RotMatrix, q)
                FD_jac = ForwardDiff.jacobian(Jr, vec(q))

                # compare
                @test all(abs(Matrix(R_jac) - FD_jac) .< thresh)
            end
        end

        # SpQuat to Quaternion
        @testset "Jacobian (SpQuat -> Quaternion)" begin

            for i = 1:10    # do some repeats

                q = nquatrand()                           # a random quaternion
                spq = Rotations.quat_to_spquat_naive(q)   # I want to test the full domain of SpQuats, not just the one with ||.|| < 1

                # transformation function
                J(X) = vec(Quaternion(SpQuat(X...)))

                # test jacobian to a Rotation matrix
                R_jac = Rotations.jacobian(Quaternion, spq)
                FD_jac = ForwardDiff.jacobian(J, vec(spq))

                # compare
                @test all(abs(Matrix(R_jac) - FD_jac) .< thresh)

            end
        end


        # Quaternion to SpQuat
        @testset "Jacobian (Quaternion -> SpQuat)" begin

            for i = 1:10    # do some repeats

                q = nquatrand()                           # a random quaternion

                # transformation function
                J(X) = vec(SpQuat(Quaternion(X...)))

                # test jacobian to a Rotation matrix
                R_jac = Rotations.jacobian(SpQuat, q)
                FD_jac = ForwardDiff.jacobian(J, vec(q))

                # compare
                @test all(abs(Matrix(R_jac) - FD_jac) .< thresh)

            end
        end

        # SpQuat to rotation matrix
        @testset "Jacobian (SpQuat -> RotMatrix)" begin

            # transformation function
            function Jr(X)
                R = Tuple(RotMatrix(SpQuat(X...)))
                V = [R[1]..., R[2]..., R[3]...]
            end

            for i = 1:10    # do some repeats

                q = nquatrand()  # a random quaternion
                spq = Rotations.quat_to_spquat_naive(q)   # I want to test the full domain of SpQuats, not just the one with ||.|| < 1

                # test jacobian to a Rotation matrix
                R_jac = Rotations.jacobian(RotMatrix, spq)
                FD_jac = ForwardDiff.jacobian(Jr, vec(spq))

                # compare
                @test all(abs(Matrix(R_jac) - FD_jac) .< thresh)

            end
        end


        ##############################
        # Jacobians for rotating stuff
        ##############################

        # Quaternion multiplication
        @testset "Jacobian (Quaternion muliplication w.r.t. the right quaternion)" begin

            for i = 1:10    # do some repeats

                ql = quatrand()    # a random quaternion (should work for non-unit quaternions)
                qr = quatrand()    # a random quaternion (should work for non-unit quaternions)

                # transformation function
                J(X) = vec(ql * Quaternion(X...))

                # test jacobian to a Rotation matrix
                R_jac = Rotations.jacobian(ql, qr)
                FD_jac = ForwardDiff.jacobian(J, vec(qr))

                # compare
                @test all(abs(Matrix(R_jac) - FD_jac) .< thresh)

            end
        end
        @testset "Jacobian (Quaternion muliplication w.r.t. the left quaternion)" begin

            for i = 1:10    # do some repeats

                ql = quatrand()    # a random quaternion (should work for non-unit quaternions)
                qr = quatrand()    # a random quaternion (should work for non-unit quaternions)

                # transformation function
                J(X) = vec(Quaternion(X...) * qr)

                # test jacobian to a Rotation matrix
                R_jac = Rotations.jacobian(ql, qr, Val{false})
                FD_jac = ForwardDiff.jacobian(J, vec(ql))

                # compare
                @test all(abs(Matrix(R_jac) - FD_jac) .< thresh)

            end
        end


        # rotate a point by a quaternion
        @testset "Jacobian (Quaternion rotation)" begin

            for i = 1:10    # do some repeats

                q = nquatrand()    # a random quaternion
                X = randn(Vec{3,Float64})

                # transformation function
                J(q) = Vector(rotate(Quaternion(q...), X))

                # test jacobian to a Rotation matrix
                R_jac = Rotations.jacobian(q, X)
                FD_jac = ForwardDiff.jacobian(J, vec(q))

                # compare
                @test all(abs(Matrix(R_jac) - FD_jac) .< thresh)

            end
        end

        # rotate a point by an SpQuat
        @testset "Jacobian (SpQuat rotation)" begin

            for i = 1:10    # do some repeats

                spq = SpQuat(nquatrand())    # a random quaternion
                X = randn(Vec{3,Float64})

                # transformation function
                J(spq) = Vector(rotate(SpQuat(spq...), X))

                # test jacobian to a Rotation matrix
                R_jac = Rotations.jacobian(spq, X)
                FD_jac = ForwardDiff.jacobian(J, vec(spq))

                # compare
                @test all(abs(Matrix(R_jac) - FD_jac) .< thresh)

            end
        end


    end


    ###################################
    # Test Hessians
    ###################################

    @testset "Hessian checks" begin

        # SpQuat to Quaternion
        @testset "Hessian (SpQuat -> Quaternion)" begin

            for i = 1:10    # do some repeats

                q = nquatrand()                           # a random quaternion
                spq = Rotations.quat_to_spquat_naive(q)   # I want to test the full domain of SpQuats, not just the one with ||.|| < 1

                # test jacobian to a Rotation matrix
                R_hess = Rotations.hessian(Quaternion, spq)
                for d = 1:4
                    H(X) = vec(Quaternion(SpQuat(X...)))[d]   # transformation function
                    FD_hess = ForwardDiff.hessian(H, vec(spq))

                    # compare
                    @test all(abs(Matrix(R_hess[d]) - FD_hess) .< thresh)
                end
            end
        end


        # Quaternion to SpQuat
        @testset "Hessian (Quaternion -> SpQuat)" begin

            for i = 1:10    # do some repeats

                q = nquatrand()                           # a random quaternion

                # test jacobian to a Rotation matrix
                R_hess = Rotations.hessian(SpQuat, q)
                for d = 1:3
                    H(X) = vec(SpQuat(Quaternion(X...)))[d]   # transformation function
                    FD_hess = ForwardDiff.hessian(H, vec(q))

                    # compare
                    @test all(abs(Matrix(R_hess[d]) - FD_hess) .< thresh)
                end
            end
        end

        @testset "Hessian (Quaternion rotation)" begin

            for i = 1:10   # do some repeats

                q = nquatrand()    # a random quaternion
                X = randn(Vec{3,Float64})

                # transformation function
                R_hess = Rotations.hessian(q, X)

                for d = 1:3
                    H(q) = Vector(rotate(Quaternion(q...), X))[d]

                    # test
                    FD_hess = ForwardDiff.hessian(H, vec(q))

                    # compare
                    @test all(abs(Matrix(R_hess[d]) - FD_hess) .< thresh)
                end
            end
        end

        @testset "Hessian (SpQuat rotation)" begin

            for i = 1:10    # do some repeats

                q = nquatrand()    # a random quaternion
                spq = SpQuat(q)
                X = randn(Vec{3,Float64})

                # transformation function
                R_hess = Rotations.hessian(spq, X)

                for d = 1:3
                    H(spq) = Vector(rotate(SpQuat(spq...), X))[d]

                    # test
                    FD_hess = ForwardDiff.hessian(H, vec(spq))

                    # compare
                    @test all(abs(Matrix(R_hess[d]) - FD_hess) .< thresh)
                end
            end
        end
    end
end
