using Documenter, SetPropagation_FEM_Examples

DocMeta.setdocmeta!(SetPropagation_FEM_Examples, :DocTestSetup,
                    :(using SetPropagation_FEM_Examples); recursive=true)

# Generate models
include("generate.jl")

makedocs(
    sitename = "SetPropagation_FEM_Examples.jl",
    modules = [SetPropagation_FEM_Examples],
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true",
                             collapselevel = 1,
                             assets = ["assets/juliareach-light.css"]),
    pages = ["Introduction" => "introduction.md",
        #"Harmonic oscillator" => "examples/sdof.md",
        #"Clamped-Free Bar" => "examples/clamped.md",
        #"One-dimensional heat transfer" => "examples/heat1d.md",
        "Concrete casting heat of hydration" => "examples/Heat3D.md"]
)

deploydocs(
    repo = "github.com/JuliaReach/SetPropagation_FEM_Examples.jl.git",
    push_preview = false
)
