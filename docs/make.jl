using Rotations
using Documenter
using Makie

DocMeta.setdocmeta!(Rotations, :DocTestSetup, :(using Rotations); recursive=true)

RotationsMakieExt = Base.get_extension(Rotations, :RotationsMakieExt)
makedocs(;
    modules=[Rotations, RotationsMakieExt],
    repo="https://github.com/JuliaGeometry/Rotations.jl/blob/{commit}{path}#{line}",
    sitename="Rotations.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://JuliaGeometry.github.io/Rotations.jl",
        assets = ["assets/custom.css", "assets/favicon.ico"],
        repolink = "https://github.com/JuliaGeometry/Rotations.jl",
    ),
    pages=[
        "Home" => "index.md",
        "Rotation Types" => "rotation_types.md",
        "2D Rotations" => [
            "Matrix" => "2d_matrix.md",
            "Angle" => "2d_angle.md",
            ],
        "3D Rotations" => [
            "Matrix" => "3d_matrix.md",
            "Euler Angles" => "3d_euler.md",
            "Angle and Axis" => "3d_angleaxis.md",
            "Quaternion and Related Parameters" => "3d_quaternion.md",
            ],
        "Rotation Generators" => [
            "Rotation Generator Types" => "rotation_generator_types.md",
            "2D Rotation Generators" => "2d_rotation_generator.md",
            "3D Rotation Generators" => "3d_rotation_generator.md",
            ],
        "General Dimensional Rotations" => "general_dimensional_rotations.md",
        "Common Methods for Rotations" => "functions.md",
        "Visualizing Rotations" => "visualizing.md",
        "Function Reference" => "functionreference.md",
        "References" => "reference.md",
    ],
)

deploydocs(;
    repo="github.com/JuliaGeometry/Rotations.jl",
)
