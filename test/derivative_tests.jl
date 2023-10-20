using ForwardDiff

#thresh = 1e-12

@testset "Derivative checks" begin

    ###################################
    # Test Jacobians
    ###################################

    @testset "Jacobian checks" begin

        # Quaternion to rotation matrix
        @testset "Jacobian (QuatRotation -> RotMatrix)" begin
            for i in 1:10    # do some repeats
                q = rand(QuatRotation)  # a random quaternion

                # test jacobian to a Rotation matrix
                R_jac = Rotations.jacobian(RotMatrix, q)
                FD_jac = ForwardDiff.jacobian(x -> SVector{9}(QuatRotation(x[1], x[2], x[3], x[4])),
                                              Rotations.params(q))

                # compare
                @test FD_jac ≈ R_jac
            end
        end

        # MRP to UnitQuternion
        @testset "Jacobian (MRP -> QuatRotation)" begin
            for i in 1:10    # do some repeats
                p = rand(MRP)  # a random MRP

                # test jacobian to a Rotation matrix
                R_jac = Rotations.jacobian(QuatRotation, p)
                FD_jac = ForwardDiff.jacobian(x -> (q = QuatRotation(MRP(x[1],x[2],x[3]));
                                                    Rotations.params(q)),
                                              SVector(p.x, p.y, p.z))

                # compare
                @test FD_jac ≈ R_jac
            end
        end

        @testset "Jacobian (MRP -> QuatRotation) [Corner Cases]" begin
            for p in [MRP(1.0, 0.0, 0.0), MRP(0.0, 1.0, 0.0), MRP(0.0, 0.0, 1.0)]
                # test jacobian to a Rotation matrix
                R_jac = Rotations.jacobian(QuatRotation, p)
                FD_jac = ForwardDiff.jacobian(x -> (q = QuatRotation(MRP(x[1],x[2],x[3]));
                                                    Rotations.params(q)),
                                              SVector(p.x, p.y, p.z))

                # compare
                @test FD_jac ≈ R_jac
            end
        end

        # MRP to QuatRotation
        @testset "Jacobian (QuatRotation -> MRP)" begin
            for i in 1:10    # do some repeats
                q = rand(QuatRotation)  # a random QuatRotation

                # test jacobian to a Rotation matrix
                R_jac = Rotations.jacobian(MRP, q)
                FD_jac = ForwardDiff.jacobian(x -> (p = MRP(QuatRotation(x[1], x[2], x[3], x[4]));
                                                    SVector(p.x, p.y, p.z)),
                                              Rotations.params(q))

                # compare
                @test FD_jac ≈ R_jac
            end
        end

        # MRP to rotation matrix
        @testset "Jacobian (MRP -> RotMatrix)" begin
            for i in 1:10    # do some repeats
                p = rand(MRP)  # a random MRP

                # test jacobian to a Rotation matrix
                R_jac = Rotations.jacobian(RotMatrix, p)
                FD_jac = ForwardDiff.jacobian(x -> SVector{9}(MRP(x[1], x[2], x[3])),
                                              SVector(p.x, p.y, p.z))

                # compare
                @test FD_jac ≈ R_jac
            end
        end


        ##############################
        # Jacobians for rotating stuff
        ##############################
#=
        # Quaternion multiplication
        @testset "Jacobian (Quaternion muliplication w.r.t. the right quaternion)" begin

            for i in 1:10    # do some repeats

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

            for i in 1:10    # do some repeats

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
        end =#


        # rotate a point by a RotMatrix
        @testset "Jacobian (RotMatrix rotation)" begin
            for i in 1:10    # do some repeats
                r = rand(RotMatrix{3,Float64})    # a random quaternion
                v = randn(SVector{3,Float64})

                # test jacobian to a Rotation matrix
                R_jac = Rotations.jacobian(r, v)
                FD_jac = ForwardDiff.jacobian(x -> RotMatrix(x)*v,
                                              SVector(r))

                # compare
                @test FD_jac ≈ R_jac
            end
        end

        # rotate a point by a quaternion
        @testset "Jacobian (QuatRotation rotation)" begin
            for i in 1:10    # do some repeats
                q = rand(QuatRotation)    # a random quaternion
                v = randn(SVector{3,Float64})

                # test jacobian to a Rotation matrix
                R_jac = Rotations.jacobian(q, v)
                FD_jac = ForwardDiff.jacobian(x -> QuatRotation(x[1], x[2], x[3], x[4])*v,
                                              Rotations.params(q))

                # compare
                @test FD_jac ≈ R_jac
            end
        end

        # rotate a point by a MRP
        @testset "Jacobian (MRP rotation)" begin
            for i in 1:10    # do some repeats
                p = rand(MRP)    # a random quaternion
                v = randn(SVector{3,Float64})

                # test jacobian to a Rotation matrix
                R_jac = Rotations.jacobian(p, v)
                FD_jac = ForwardDiff.jacobian(x -> MRP(x[1], x[2], x[3])*v,
                                              SVector(p.x, p.y, p.z))

                # compare
                @test FD_jac ≈ R_jac
            end
        end

        @testset "RotationVec near zero" begin
            rot(x) = RotationVec(x, 0, 0)
            for i = 1:10
                v = randn(SVector{3,Float64})
                rotv(x) = rot(x)*v
                drotv(x) = ForwardDiff.derivative(rotv, x)
                # The following broken tests will be fixed by https://github.com/JuliaDiff/ForwardDiff.jl/pull/669
                @test_broken drotv(0.0) ≈ [0, -v[3], v[2]]
                @test drotv(1e-20) ≈ [0, -v[3], v[2]]
                @test_broken ForwardDiff.derivative(drotv, 0.0) ≈ [0, -v[2], -v[3]]
                @test ForwardDiff.derivative(drotv, 1e-20) ≈ [0, -v[2], -v[3]]
            end
        end
#=
        # rotate a point by an MRP
        @testset "Jacobian (MRP rotation)" begin

            for i in 1:10    # do some repeats

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
=#

    end


    ###################################
    # Test Hessians
    ###################################
    #=
    @testset "Hessian checks" begin

        # SpQuat to Quaternion
        @testset "Hessian (SpQuat -> Quaternion)" begin

            for i in 1:10    # do some repeats

                q = nquatrand()                           # a random quaternion
                spq = Rotations.quat_to_spquat_naive(q)   # I want to test the full domain of SpQuats, not just the one with ||.|| < 1

                # test jacobian to a Rotation matrix
                R_hess = Rotations.hessian(Quaternion, spq)
                for d in 1:4
                    H(X) = vec(Quaternion(SpQuat(X...)))[d]   # transformation function
                    FD_hess = ForwardDiff.hessian(H, vec(spq))

                    # compare
                    @test all(abs(Matrix(R_hess[d]) - FD_hess) .< thresh)
                end
            end
        end


        # Quaternion to SpQuat
        @testset "Hessian (Quaternion -> SpQuat)" begin

            for i in 1:10    # do some repeats

                q = nquatrand()                           # a random quaternion

                # test jacobian to a Rotation matrix
                R_hess = Rotations.hessian(SpQuat, q)
                for d in 1:3
                    H(X) = vec(SpQuat(Quaternion(X...)))[d]   # transformation function
                    FD_hess = ForwardDiff.hessian(H, vec(q))

                    # compare
                    @test all(abs(Matrix(R_hess[d]) - FD_hess) .< thresh)
                end
            end
        end

        @testset "Hessian (Quaternion rotation)" begin

            for i in 1:10   # do some repeats

                q = nquatrand()    # a random quaternion
                X = randn(Vec{3,Float64})

                # transformation function
                R_hess = Rotations.hessian(q, X)

                for d in 1:3
                    H(q) = Vector(rotate(Quaternion(q...), X))[d]

                    # test
                    FD_hess = ForwardDiff.hessian(H, vec(q))

                    # compare
                    @test all(abs(Matrix(R_hess[d]) - FD_hess) .< thresh)
                end
            end
        end

        @testset "Hessian (SpQuat rotation)" begin

            for i in 1:10    # do some repeats

                q = nquatrand()    # a random quaternion
                spq = SpQuat(q)
                X = randn(Vec{3,Float64})

                # transformation function
                R_hess = Rotations.hessian(spq, X)

                for d in 1:3
                    H(spq) = Vector(rotate(SpQuat(spq...), X))[d]

                    # test
                    FD_hess = ForwardDiff.hessian(H, vec(spq))

                    # compare
                    @test all(abs(Matrix(R_hess[d]) - FD_hess) .< thresh)
                end
            end
        end
    end =#
end
