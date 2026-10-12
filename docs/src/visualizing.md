```@setup ext
using Rotations, Makie
Ext = Base.get_extension(Rotations, :RotationsMakieExt)
```
# Visualizing Rotations

Visualizing rotations is crucial for understanding their behavior and effects. While a general ``3 \times 3`` rotation matrix has 9 parameters, rotations actually live on the 3-dimensional manifold ``SO(3)``, making them suitable for 3D visualization.

This package provides plotting capabilities through [Makie.jl](https://makie.juliaplots.org/) via the `RotationsMakieExt` extension, which is automatically loaded when both `Rotations` and a Makie backend are available.

## Plotting Rotation Matrices

### Basic 2D Rotation Plotting

```julia
using Rotations
using CairoMakie  # or GLMakie for interactive plots

# Create a 2D rotation
rot2d = RotMatrix{2}(π/4)  # 45-degree rotation

# Plot the rotation
fig, ax, plt = plot(rot2d)
```

```@docs
RotationsMakieExt.rotation2plot
RotationsMakieExt.rotation2plot!
RotationsMakieExt.Rotation2Plot
```

The functions are also accessible from the `Visualization` module of the package.

### Basic 3D Rotation Plotting

```julia
# Create a 3D rotation from Euler angles
rot3d = RotXYZ(π/6, π/4, π/3)

# Plot the rotation
fig, ax, plt = plot(rot3d)
```

```@docs
RotationsMakieExt.rotation3plot
RotationsMakieExt.rotation3plot!
RotationsMakieExt.Rotation3Plot
```

### Combining Multiple Rotations

```julia
using GLMakie

# Create several rotations
rotations = [
    RotX(π/6),
    RotY(π/4), 
    RotZ(π/3),
    RotXYZ(π/6, π/4, π/3)
]

# Plot them in a grid
fig = Figure(resolution = (800, 800))

for (i, rot) in enumerate(rotations)
    ax = Axis3(fig[div(i-1, 2) + 1, mod(i-1, 2) + 1], 
               title = "Rotation $i")
    plot!(ax, rot, 
          boxcolor = [:red, :green, :blue, :orange][i],
          axissize = 1.5)
end

fig
```
