# 3D Rotation Matrix

An ``3 \times 3`` rotation matrix storing the rotation.
This is a simple wrapper for a [StaticArrays](https://github.com/JuliaArrays/StaticArrays.jl) `SMatrix{3,3,T}`.
A rotation matrix ``R`` should have the property ``I = R R^\top``, but this isn't enforced by the constructor.
On the other hand, all the types below are guaranteed to be "proper" rotations for all input parameters (equivalently: parity conserving, in ``SO(3)``, ``\det(R) = 1``, or a rotation without reflection).

## `RotMatrix3`

**example**

(TBW)
