# Common Methods for Rotations

## `rotation_angle`, `rotation_axis`
### 3D
A rotation matrix can be expressed with one rotation around an axis and angle.
these function can calculate these values.

**example**
```@setup angle_and_axis
using Rotations
```

```@repl angle_and_axis
R = RotXYZ(2.4, -1.8, 0.5)
θ = rotation_angle(R)
n = rotation_axis(R)
# These matrices are approximately equal.
R ≈ AngleAxis(θ, n...)
```

### 2D
(TBW)

## `Rotations.params`
The parameters of the rotation can be obtained by `Rotations.params`.

**example**
```@setup params
using Rotations
```

```@repl params
R = one(RotMatrix{3})  # Identity matrix
Rotations.params(RotYZY(R))  # Proper Euler angles, (y,z,y)
Rotations.params(RotXYZ(R))  # Tait–Bryan angles, (y,z,y)
Rotations.params(AngleAxis(R))  # Rotation around an axis (theta, axis_x, axis_y, axis_z)
Rotations.params(RotationVec(R))  # Rotation vector (v_x, v_y, v_z)
Rotations.params(UnitQuaternion(R))  # Quaternion (w, x, y, z)
Rotations.params(RodriguesParam(R))  # Rodrigues Parameters (x, y, z)
Rotations.params(MRP(R))  # Modified Rodrigues Parameters (x, y, z)
```

## `isrotation`
Check the given matrix is rotation matrix.

**example**

(TBW)

## `rand`
```@setup rand
using Rotations
```

### `rand` for ``SO(2)``
The following types have the same algebraic structure as ``SO(2)``

* `RotMatrix{2}`
* `Angle2d`
* `RotX`
* `RotY`
* `RotZ`

The random distribution is based on [Haar measure](https://en.wikipedia.org/wiki/Haar_measure).

**example**
```@repl rand
R = rand(Angle2d)
```

### `rand` for ``SO(3)``
The following types have an algebraic structure that is homomorphic to ``SO(3)``.

* `RotMatrix{3}`
* `RotXYZ` (and other Euler angles)
* `AngleAxis`
* `RotationVec`
* `UnitQuaternion`
* `RodriguesParam`
* `MRP`

**example**
```@repl rand
R = rand(RotationVec)
```

The random distribution is based on [Haar measure](https://en.wikipedia.org/wiki/Haar_measure).

### `rand` for `RotXY` and etc.
There also are methods for `rand(::RotXY)` and other 2-axis rotations.

**example**
```@repl rand
R = rand(RotXY)
```

The random distribution is **NOT** based on [Haar measure](https://en.wikipedia.org/wiki/Haar_measure) because the set of `RotXY` doesn't have group structure.

Note that:
* `rand(RotX)` is same as `RotX(2π*rand())`.
* `rand(RotXY)` is same as `RotXY(2π*rand(), 2π*rand())`.
* `rand(RotXYZ)` is **not** same as `RotXYZ(2π*rand(), 2π*rand(), 2π*rand())`.
* But `rand(RotXYZ)` is same as `RotXYZ(rand(UnitQuaternion))`.
