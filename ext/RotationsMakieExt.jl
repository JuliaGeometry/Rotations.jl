module RotationsMakieExt

using Rotations
using Makie
using StaticArrays

Makie.plottype(::Rotation{2}) = Rotation2Plot
Makie.plottype(::Rotation{3}) = Rotation3Plot

@recipe Rotation2Plot (rot::Rotation{2},) begin
    begin
        arrows_attributes = Makie.documented_attributes(Arrows2D)
        filter!(a -> a.first ∉ [:color], arrows_attributes.d)
        arrows_attributes
    end...
    "Origin point for the 2D rotation visualization"
    origin = SVector(0, 0)
    "Size of the reference square"
    boxsize = 0.2
    "Length of the axis lines"
    axissize = 1.0
    "Colours for the two axes (x, y)"
    axiscolors = [:red, :green]
    "Colour of the reference square"
    boxcolor = :gray
end

function Makie.plot!(p::Rotation2Plot{<:Tuple{<:Rotation{2}}})
    # Use ComputationPipeline for dynamic updates
    map!(p.attributes, [:rot, :origin, :boxsize, :axissize], [:arrows_origin, :arrows1_end, :arrows2_end, :box_points]) do rot, origin, l, L
        e₁ = rot[:, 1]
        e₂ = rot[:, 2]
        ox, oy = origin

        # Compute box corners
        box_points = [Point2f(ox, oy) + rot * Point2f(x, y) for x in (-l, l) for y in (-l, l)][[1, 2, 4, 3]]

        return (Point2f(ox, oy), Point2f(ox + L * e₁[1], oy + L * e₁[2]), Point2f(ox + L * e₂[1], oy + L * e₂[2]), box_points)
    end

    map!(p.attributes, [:axiscolors], [:arrows1_color, :arrows2_color]) do axiscolors
        (axiscolors[1], axiscolors[2])
    end

    arrows2d!(p, p.attributes, p.arrows_origin, p.arrows1_end; color=p.arrows1_color)
    arrows2d!(p, p.attributes, p.arrows_origin, p.arrows2_end; color=p.arrows2_color)
    poly!(p, p.attributes, p.box_points; color=p.boxcolor)

    return p
end

@recipe Rotation3Plot (rot::Rotation{3},) begin
    begin
        arrows_attributes = Makie.documented_attributes(Arrows3D)
        filter!(a -> a.first ∉ [:color], arrows_attributes.d)
        arrows_attributes
    end...
    "Origin point for the 3D rotation visualization"
    origin = SVector(0, 0, 0)
    "Size of the reference cube"
    boxsize = 0.2
    "Length of the axis lines"
    axissize = 1.0
    "Colours for the three axes (x, y, z)"
    axiscolors = [:red, :green, :blue]
    "Colour of the reference cube faces"
    boxcolor = :gray
    "Colour of cube edges"
    edgecolor = :black
end

function Makie.plot!(p::Rotation3Plot{<:Tuple{<:Rotation{3}}})
    # Use ComputationPipeline for dynamic updates
    map!(p.attributes, [:rot, :origin, :boxsize, :axissize], [:arrows_origin, :arrows1_end, :arrows2_end, :arrows3_end, :cube_points, :cube_faces, :cube_edges]) do rot, origin, l, L
        e₁ = rot[:, 1]
        e₂ = rot[:, 2]
        e₃ = rot[:, 3]
        ox, oy, oz = origin

        # Compute cube vertices
        cube_points = [Point3f(ox, oy, oz) + rot * Point3f(x, y, z) for x in (-l, l) for y in (-l, l) for z in (-l, l)]

        # Define cube faces (triangulated)
        cube_faces = [
            1, 2, 4, 4, 3, 1, # bottom
            5, 6, 8, 8, 7, 5, # top
            1, 2, 6, 6, 5, 1, # front
            3, 4, 8, 8, 7, 3, # back
            1, 3, 7, 7, 5, 1, # left
            2, 4, 8, 8, 6, 2  # right
        ]

        # Define cube edges for wireframe
        cube_edges = [
            # Bottom face edges
            cube_points[1], cube_points[2], cube_points[2], cube_points[4],
            cube_points[4], cube_points[3], cube_points[3], cube_points[1],
            # Top face edges
            cube_points[5], cube_points[6], cube_points[6], cube_points[8],
            cube_points[8], cube_points[7], cube_points[7], cube_points[5],
            # Vertical edges
            cube_points[1], cube_points[5], cube_points[2], cube_points[6],
            cube_points[3], cube_points[7], cube_points[4], cube_points[8]
        ]

        return (
            Point3f(ox, oy, oz),
            Point3f(ox + L * e₁[1], oy + L * e₁[2], oz + L * e₁[3]),
            Point3f(ox + L * e₂[1], oy + L * e₂[2], oz + L * e₂[3]),
            Point3f(ox + L * e₃[1], oy + L * e₃[2], oz + L * e₃[3]),
            cube_points, cube_faces, cube_edges
        )
    end

    map!(p.attributes, [:axiscolors], [:arrows1_color, :arrows2_color, :arrows3_color]) do axiscolors
        (axiscolors[1], axiscolors[2], axiscolors[3])
    end

    arrows3d!(p, p.attributes, p.arrows_origin, p.arrows1_end, color=p.arrows1_color)
    arrows3d!(p, p.attributes, p.arrows_origin, p.arrows2_end, color=p.arrows2_color)
    arrows3d!(p, p.attributes, p.arrows_origin, p.arrows3_end, color=p.arrows3_color)
    mesh!(p, p.attributes, p.cube_points, p.cube_faces; color=p.boxcolor)
    linesegments!(p, p.attributes, p.cube_edges; color=p.edgecolor)

    return p
end

Makie.preferred_axis_type(::Rotation3Plot) = Makie.Axis3

end
