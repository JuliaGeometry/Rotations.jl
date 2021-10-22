# Types in Rotations.jl

## Abstract rotations
A matrix ``R`` is called *rotation matrix* if ``R`` satisfies ``R^\top = R^{-1}`` and ``\det(R)=1``.
In `Rotations.jl`, there's an abstract type for rotations matrix, `Rotation{L}`.
Where `L` is a size of the rotation matrix.

## Type hierarchy

```@setup hierarchy
using InteractiveUtils
```

```@repl hierarchy
using Rotations, StaticArrays
Rotation <: StaticMatrix <: AbstractMatrix
subtypes(Rotation{2})
subtypes(Rotation{3})
```

## 2D rotations
* `RotMatrix2{T}`
    * Rotation matrix in 2 dimensional Euclidean space.
* `Angle2d`
    * Parametrized with rotational angle.

## 3D rotations
* `RotMatrix3{T}`
    * Rotation matrix in 3 dimensional Euclidean space.
* `RotX`, `RotYZ`, `RotXYZ` and etc.
    * Euler angles.
* `AngleAxis`
    * Rotation around given axis and angle.
* `RotationVec`
    * Rotation around given axis. The length of axis vector represents its angle.
* `UnitQuaternion`
    * A 3D rotation parameterized by a unit quaternion.
* `MRP`
    * A 3D rotation encoded by the stereographic projection of a unit quaternion.
