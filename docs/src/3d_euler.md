# 3D Rotation with Euler Angles

[Euler angles](https://en.wikipedia.org/wiki/Euler_angles) are a way to represent a rotation matrix with three rotations around cardinal axes.

In Rotations.jl, there are 12 concrete types for Euler angles.

* Proper Euler angles
    * `RotZXZ`, `RotXYX`, `RotYZY`, `RotZYZ`, `RotXZX`, `RotYXY`
* Tait–Bryan angles
    * `RotXYZ`, `RotYZX`, `RotZXY`, `RotXZY`, `RotZYX`, `RotYXZ`

In addition, Rotations.jl provides concrete types that represent rotations in one or two axes.

* one axis
    * `RotX`, `RotY`, `RotZ`
* two axes
    * `RotXY`, `RotYZ`, `RotZX`, `RotXZ`, `RotZY`, `RotYX`

## Rotation around one axis

```math
\begin{aligned}
R_{x}(\alpha)
&= \begin{pmatrix}
1 & 0 & 0 \\
0 & \cos(\alpha) & -\sin(\alpha) \\
0 & \sin(\alpha) & \cos(\alpha) \\
\end{pmatrix}, \\
R_{y}(\alpha)
&= \begin{pmatrix}
\cos(\alpha) & 0 & \sin(\alpha) \\
0 & 1 & 0 \\
-\sin(\alpha) & 0 & \cos(\alpha) \\
\end{pmatrix}, \\
R_{z}(\alpha)
&= \begin{pmatrix}
\cos(\alpha) & -\sin(\alpha) & 0\\
\sin(\alpha) & \cos(\alpha) & 0\\
0 & 0 & 1
\end{pmatrix}
\end{aligned}
```

**example**
```@setup one_axis
using Rotations
```

Here's an example for `RotZ`:
```@repl one_axis
α = 1.2  # Rotation angle
R = RotZ(α)
Q = [cos(α) -sin(α) 0
      sin(α)  cos(α) 0
      0       0      1]
# These matrices are equal
R == Q
```

And more examples for `RotX` and `RotY`:

```@repl one_axis
# These matrices are equal
RotX(α) == [1 0       0
            0 cos(α) -sin(α)
            0 sin(α)  cos(α)]
# These matrices are equal
RotY(α) == [cos(α) 0 sin(α)
            0      1 0
           -sin(α) 0 cos(α)]
```

## Rotation around two axes

```math
\begin{aligned}
R_{xy}(\alpha,\beta)
&= R_{x}(\alpha)R_{y}(\beta),
& R_{yx}(\alpha,\beta)
&= R_{y}(\alpha)R_{x}(\beta), \\
%
R_{yz}(\alpha,\beta)
&= R_{y}(\alpha)R_{z}(\beta),
& R_{zy}(\alpha,\beta)
&= R_{z}(\alpha)R_{y}(\beta), \\
%
R_{zx}(\alpha,\beta)
&= R_{z}(\alpha)R_{x}(\beta),
& R_{xz}(\alpha,\beta)
&= R_{x}(\alpha)R_{z}(\beta)
\end{aligned}
```

**example**
```@setup two_axis
using Rotations
```

```@repl two_axis
α, β = 1.2, 4.7  # Rotation angles
RotX(α) * RotY(β)
RotXY(α, β)
# These matrices are equal
RotX(α) * RotY(β) == RotXY(α, β)
```

## Rotation around three axes (Euler Angles)

**Proper Euler angles**

```math
\begin{aligned}
R_{xyx}(\alpha, \beta, \gamma) &= R_x(\alpha) R_y(\beta) R_x(\gamma),
& R_{yxy}(\alpha, \beta, \gamma) &= R_y(\alpha) R_x(\beta) R_y(\gamma), \\
R_{yzy}(\alpha, \beta, \gamma) &= R_y(\alpha) R_z(\beta) R_y(\gamma),
& R_{zyz}(\alpha, \beta, \gamma) &= R_z(\alpha) R_y(\beta) R_z(\gamma), \\
R_{zxz}(\alpha, \beta, \gamma) &= R_z(\alpha) R_x(\beta) R_z(\gamma),
& R_{xzx}(\alpha, \beta, \gamma) &= R_x(\alpha) R_z(\beta) R_x(\gamma)
\end{aligned}
```

**Tait–Bryan angles**

```math
\begin{aligned}
R_{xyz}(\alpha, \beta, \gamma) &= R_x(\alpha) R_y(\beta) R_z(\gamma),
& R_{yxz}(\alpha, \beta, \gamma) &= R_y(\alpha) R_x(\beta) R_z(\gamma), \\
R_{yzx}(\alpha, \beta, \gamma) &= R_y(\alpha) R_z(\beta) R_x(\gamma),
& R_{zyx}(\alpha, \beta, \gamma) &= R_z(\alpha) R_y(\beta) R_x(\gamma), \\
R_{zxy}(\alpha, \beta, \gamma) &= R_z(\alpha) R_x(\beta) R_y(\gamma),
& R_{xzy}(\alpha, \beta, \gamma) &= R_x(\alpha) R_z(\beta) R_y(\gamma)
\end{aligned}
```

**example**
```@setup three_axis
using Rotations
```

```@repl three_axis
α, β, γ = 1.2, 4.7, -0.4  # Rotation angles
RotYXZ(α, β, γ)
RotY(α)*RotX(β)*RotZ(γ)
# These matrices are equal
RotYXZ(α, β, γ) == RotY(α)*RotX(β)*RotZ(γ)
```
