module RotationsMakieExt

using Rotations
using Makie
using StaticArrays

Makie.plottype(::Rotation{2}) = Rotation2Plot
Makie.plottype(::Rotation{3}) = Rotation3Plot

@recipe Rotation2Plot (rot::Rotation{2},) begin
    origin = SVector(0, 0)
    boxsize = 0.2
    axissize = 1.0
    axiscolors = [:red, :green]
    boxcolor = :gray
end

function Makie.plot!(p::Rotation2Plot{<:Tuple{<:Rotation{2}}})
    rot = p[1][]
    origin = p.origin[]
    boxsize = p.boxsize[]
    axissize = p.axissize[]
    axiscolors = p.axiscolors[]
    boxcolor = p.boxcolor[]

    l = boxsize
    L = axissize
    e₁ = rot[:, 1]
    e₂ = rot[:, 2]
    ox, oy = origin

    # Draw axes
    lines!(p, [ox, ox + L * e₁[1]], [oy, oy + L * e₁[2]], color=axiscolors[1])
    lines!(p, [ox, ox + L * e₂[1]], [oy, oy + L * e₂[2]], color=axiscolors[2])

    # Draw box
    points = [SVector(ox, oy) + rot * SVector(x, y) for x in (-l, l) for y in (-l, l)]
    poly!(p, points[[1, 2, 4, 3]], color=boxcolor)

    p
end

@recipe Rotation3Plot (rot::Rotation{3},) begin
    origin = SVector(0, 0, 0)
    boxsize = 0.2
    axissize = 1.0
    axiscolors = [:red, :green, :blue]
    boxcolor = :gray
end

function Makie.plot!(p::Rotation3Plot{<:Tuple{<:Rotation{3}}})
    rot = p[1][]
    origin = p.origin[]
    boxsize = p.boxsize[]
    axissize = p.axissize[]
    axiscolors = p.axiscolors[]
    boxcolor = p.boxcolor[]

    l = boxsize
    L = axissize
    e₁ = rot[:, 1]
    e₂ = rot[:, 2]
    e₃ = rot[:, 3]
    ox, oy, oz = origin

    # Draw axes
    lines!(p, [ox, ox + L * e₁[1]], [oy, oy + L * e₁[2]], [oz, oz + L * e₁[3]], color=axiscolors[1])
    lines!(p, [ox, ox + L * e₂[1]], [oy, oy + L * e₂[2]], [oz, oz + L * e₂[3]], color=axiscolors[2])
    lines!(p, [ox, ox + L * e₃[1]], [oy, oy + L * e₃[2]], [oz, oz + L * e₃[3]], color=axiscolors[3])

    # Draw box
    points = [Point3f(ox, oy, oz) + rot * Point3f(x, y, z) for x in (-l, l) for y in (-l, l) for z in (-l, l)]
    faces = [
        1, 2, 4, 4, 3, 1, # bottom
        5, 6, 8, 8, 7, 5, # top
        1, 2, 6, 6, 5, 1, # front
        3, 4, 8, 8, 7, 3, # back
        1, 3, 7, 7, 5, 1, # left
        2, 4, 8, 8, 6, 2  # right
    ]
    mesh!(p, points, faces, color=boxcolor)

    p
end

Makie.preferred_axis_type(::Rotation3Plot) = Makie.Axis3

end
