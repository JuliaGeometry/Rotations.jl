# Rotations.jl

[![Build Status](https://travis-ci.org/FugroRoames/Rotations.jl.svg?branch=master)](https://travis-ci.org/FugroRoames/Rotations.jl)

This package implements various 3D rotation parameterizations and defines conversions between them.  Parameterizations encode rotations about the origin (see [CoordinateTransformations.jl](https://github.com/FugroRoames/CoordinateTransformations.jl) or [AffineTransforms.jl](https://github.com/timholy/AffineTransforms.jl) for composing more general transformations).

### Example Usage

```julia

    using Rotations

    # create a rotation matrix
    R = eye(RotMatrix{Float64})

    # create a point
    X = Vec(1.0, 2.0, 3.0)

    # convert to a Quaternion and rotate the point
    q = Quaternion(R)
    Xo = rotate(ea, X)

    # convert to a Stereographic quaternion projection (recommended for applications with differentiation) and rotate
    spq = SpQuat(R)
    Xo = rotate(spq, X)

    # convert to a Rodrigues Vector and rotate
    rv = RodriguesVec(R)
    Xo = rotate(rv, X)

    # convert to Euler Angles (using the default ordering scheme) and rotate
    ea = EulerAngles(R)
    Xo = rotate(ea, X)

    # convert to Euler Angles (specifying an angle order) and rotate
    ea = EulerAngles{Rotations.EulerXYZ}(R)
    Xo = rotate(ea, X)

    # convert to proper Euler Angles (using the default ordering scheme) and rotate
    ea = ProperEulerAngles(R)
    Xo = rotate(ea, X)

    # convert to proper Euler Angles (specifying an angle order) and rotate
    ea = ProperEulerAngles{Rotations.EulerXYX}(R)
    Xo = rotate(ea, X)


```

### Rotation Parameterizations

1. **Rotation Matrix** `RotMatrix{T}`

    A 3 x 3 rotation matrix storing the rotation.  This is a typealias of a [FixedSizeArrays](https://github.com/SimonDanisch/FixedSizeArrays.jl) `Mat{3,3,T}`.  A rotation matrix `R` should have the property `I = R * R'`.


2. **Arbitrary Axis Rotation** `AngleAxis{T}`

    A 4 element immutable array with fields `theta`, `axis_x`, `axis_y`, and `axis_z` to store the rotation angle and axis of the rotation.


3. **EulerAngles** `EulerAngles{Order <: TaitByranOrder, T}`

    A 3 element immutable array which stores the Euler angles with [**Tait Byran**](https://en.wikipedia.org/wiki/Euler_angles#Tait.E2.80.93Bryan_angles) angle ordering (e.g. `EulerXYZ`) encoded by the template parameter `Order`.


4. **ProperEulerAngles** `ProperEulerAngles{Order <: ProperEulerOrder, T}`

    A 3 element immutable array which stores the Euler angles with [**Proper Euler**](https://en.wikipedia.org/wiki/Euler_angles#Conventions) angle ordering (e.g. `EulerXYX`) encoded by the template paramter `Order`.


5. **Quaternions** `Quaternion{T}`

    A 4 element immutable array containing the quaternion representation of the rotation.  This uses the [Quaternions](https://github.com/JuliaGeometry/Quaternions.jl) package.

    A non-unit quaternion ```q``` is treated as a scaled unit quaternion, ```q = s * qhat```, so that ```RotMatrix(q) == s * RotMatrix(qhat)```


6. **Rodrigues Vector** `RodriguesVec{T}`

    A 3 element immutable array encoding an angle axis representation as angle * axis.  This type is used in packages such as [OpenCV](http://docs.opencv.org/2.4/modules/calib3d/doc/camera_calibration_and_3d_reconstruction.html#void%20Rodrigues%28InputArray%20src,%20OutputArray%20dst,%20OutputArray%20jacobian%29).

    Note: If you're differentiating a Rodrigues Vector check the result is what you expect at theta = 0.  The first derivitate of `rotate()` *should* behave, but higher derivitives / parameterization conversions should be tested.  The Stereographic Quaternion Projection is the recommended three parameter format for differentiation.


7. **Stereographic Quaternion Projection** `SpQuat{T}`

    A 3 element immutable array containing the stereographic projection of a unit quaternion.  This projection can be visualized as a pin hole camera, with the pin hole matching the quaternion `[-1,0,0,0]` and the image plane containing the origin and having normal direction `[1,0,0,0]`.  The "no rotation" `Quaternion(1.0,0,0,0)` then maps to the `SpQuat(0,0,0)`

    These are similar to the Rodrigues vector in that the axis direction is stored, and the rotation angle is encoded in the length of the axis.  This type has the nice property that the derivitives of the rotation matrix w.r.t. the `SpQuat` parameters are rational functions.  This makes the `SpQuat` type a good choice for differentiation / optimization.



### Import / Export

All parameterizations can be converted to and from mutable / immutable vectors, e.g.

```julia

    using FixedSizeArrays

    # export
    q = Quaternion(1.0,0,0,0)
    v_mutable = vec(q)
    v_immutable = Vec(q)

    # import
    q2 = Quaternion(v_mutable)
    q2 = Quaternion(v_immutable)

```

### Derivatives

Some derivative calculations are include in this package, e.g.

```julia
    q = Quaternion(1.0,0,0,0)

    # 1st order only for conversion to a rotation matrix
    jac = Rotations.jacobian(RotMatrix, q)

    # 1st and 2nd order for Quaternion <-> SpQuat
    jac = Rotations.jacobian(SpQuat, q)
    hess = Rotations.hessian(SpQuat, q)

    # 1st and 2nd order for rotating points using Quaternions and SpQuats
    using FixedSizeArrays
    X = randn(Vec{3,Float64})
    jac = Rotations.jacobian(q, X)
    hess = Rotations.hessian(q, X)
```

### Notes

This package assumes [active (right handed) rotations](https://en.wikipedia.org/wiki/Active_and_passive_transformation) where applicable.


### Why use immutables / FixedSizeArrays?

They're faster (BLAS isn't great for 3x3 matrices) and don't need preallocating.  A benchmark case is included:

```julia
    cd(Pkg.dir("Rotations") * "/test")
    include("benchmark.jl")
    BenchMarkRotations.benchmark_mutable()
    # Rotating using mutables
    #   0.077123 seconds (3 allocations: 208 bytes)
    #   Rotating using immutables
    #   0.003569 seconds (4 allocations: 160 bytes)
```

```

