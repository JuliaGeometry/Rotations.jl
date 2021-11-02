# Rotations

3D rotations made easy in Julia

!!! note "Documentation"
    The documentation is still work in progress.
    For more information, see also
    * [README in the repository](https://github.com/JuliaGeometry/Rotations.jl)
    * [Tests in the repository](https://github.com/JuliaGeometry/Rotations.jl/tree/master/test)
    (TBW) comments are left in areas that have not yet been documented.
    Feel free to [open pull requests](https://github.com/JuliaGeometry/Rotations.jl/pulls) and improve this document!

## Installation
```
pkg> add Rotations
```

## First example

```@repl
using Rotations

# 3D Rotation by Euler Angles
R_euler = RotXYZ(1,2,3)

# Get an angle and an axis of the rotation
rotation_angle(R_euler), rotation_axis(R_euler)

# Convert the rotation to unit quaternion
R_quat = UnitQuaternion(R_euler)

# Get quaternion parameters of the rotation
Rotations.params(R_quat)

# Convert the rotation to MRP (Modified Rodrigues Parameters)
R_mrp = MRP(R_euler)
Rotations.params(R_mrp)

# Get parameters of the MRP
Rotations.params(R_mrp)

# Also supports 2D rotation
R_2d = Angle2d(π/6)

# Also supports some differentiation
Rotations.jacobian(RotMatrix, R_quat)
```
