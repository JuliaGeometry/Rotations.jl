# 3D Rotation with Quaternion and Related Parameters

```@setup quaternions
using Rotations
```

## `QuatRotation`
A 3D rotation parametrized by a unit [quaternion](https://en.wikipedia.org/wiki/Quaternion) ([versor](https://en.wikipedia.org/wiki/Versor)).
Note that the constructor will renormalize the quaternion to be a unit quaternion, and that although they follow the same multiplicative *algebra* as quaternions, it is better to think of `QuatRotation` as a ``3 \times 3`` matrix rather than as a quaternion *number*.

```math
\begin{aligned}
    \left\{w_{\mathrm{Q}} + ix_{\mathrm{Q}} + jy_{\mathrm{Q}} + kz_{\mathrm{Q}} \in \mathbb{H} \ | \ x_{\mathrm{Q}}, y_{\mathrm{Q}}, z_{\mathrm{Q}} \in \mathbb{R} \right\}
    \simeq S^3
\end{aligned}
```

**example**
```@repl quaternions
one(QuatRotation)  # identity rotation
α, β, γ = 1.2, -0.8, 0.1;
RotX(α) ≈ QuatRotation(cos(α/2),sin(α/2),0,0)  # These matrices are equal
RotY(β) ≈ QuatRotation(cos(β/2),0,sin(β/2),0)  # These matrices are equal
RotZ(γ) ≈ QuatRotation(cos(γ/2),0,0,sin(γ/2))  # These matrices are equal
```

See also [Quaternions.jl documentation](https://juliageometry.github.io/Quaternions.jl/stable/examples/rotations/).

## `RodriguesParam`
A 3-parameter representation of 3D rotations that has a singularity at ``180^{\circ}``.
They can be interpreted as a projection of the unit quaternion onto the plane tangent to the quaternion identity.
They are computationally efficient and do not have a sign ambiguity of unit quaternion.

```math
\left\{
\begin{aligned}
    \begin{aligned}
    x_{\mathrm{R}} &= x_{\mathrm{Q}}/w_{\mathrm{Q}} \\
    y_{\mathrm{R}} &= y_{\mathrm{Q}}/w_{\mathrm{Q}} \\
    z_{\mathrm{R}} &= z_{\mathrm{Q}}/w_{\mathrm{Q}}
    \end{aligned}
\end{aligned}
\right.
```

**example**
```@repl quaternions
one(RodriguesParam)  # identity rotation
α, β, γ = 1.2, -0.8, 0.1;
RotX(α) ≈ RodriguesParam(tan(α/2),0,0)  # These matrices are equal
RotY(β) ≈ RodriguesParam(0,tan(β/2),0)  # These matrices are equal
RotZ(γ) ≈ RodriguesParam(0,0,tan(γ/2))  # These matrices are equal
```

## `MRP` (Modified Rodrigues Parameters)

A 3D rotation encoded by the stereographic projection of a unit quaternion.
This projection can be visualized as a pin hole camera, with the pin hole matching the quaternion ``-1+0i+0j+0k`` and the image plane containing the origin and having normal direction ``1+0i+0j+0k``.
The "identity rotation" `Quaternion(1.0,0,0,0)` then maps to the `MRP(0,0,0)`

These are similar to the Rodrigues vector in that the axis direction is stored in an unnormalized form, and the rotation angle is encoded in the length of the axis.
This type has the nice property that the derivatives of the rotation matrix w.r.t. the `MRP` parameters are rational functions, making the `MRP` type a good choice for differentiation / optimization.

They are frequently used in aerospace applications since they are a 3-parameter representation whose singularity happens at ``360^\circ``.
In practice, the singularity can be avoided with some switching logic between one of two equivalent MRPs (obtained by projecting the negated quaternion).

```math
\left\{
\begin{aligned}
    \begin{aligned}
    x_{\mathrm{M}} &= \dfrac{x_{\mathrm{Q}}}{w_{\mathrm{Q}}-1} \\
    y_{\mathrm{M}} &= \dfrac{y_{\mathrm{Q}}}{w_{\mathrm{Q}}-1} \\
    z_{\mathrm{M}} &= \dfrac{z_{\mathrm{Q}}}{w_{\mathrm{Q}}-1}
    \end{aligned}
\end{aligned}
\right.
```

**example**
```@repl quaternions
one(MRP)  # identity rotation
α, β, γ = 1.2, -0.8, 0.1;
RotX(α) ≈ MRP(sin(α/2)/(cos(α/2)-1),0,0)  # These matrices are equal
RotY(β) ≈ MRP(0,sin(β/2)/(cos(β/2)-1),0)  # These matrices are equal
RotZ(γ) ≈ MRP(0,0,sin(γ/2)/(cos(γ/2)-1))  # These matrices are equal
```
