# generate examples
import Literate
import SetPropagation_FEM_Examples: @modelpath

MODELS = [
    joinpath(@__DIR__, "..", "examples", "SDOF"),
    #joinpath(@__DIR__, "..", "examples", "Clamped"),
    #joinpath(@__DIR__, "..", "examples", "Heat1D"),
    #joinpath(@__DIR__, "..", "examples", "Heat3D")
]
GENERATEDDIR = joinpath(@__DIR__, "src", "examples")
MODELDIR = joinpath(@__DIR__, "..", "examples")
mkpath(GENERATEDDIR)

macro modelpath(model_path::String, name::String)
    return joinpath(MODELDIR, model_path, name)
end

for model in MODELS
    for file in readdir(model)
        if endswith(file, ".jl")
            input = abspath(joinpath(model, file))
            script = Literate.script(input, GENERATEDDIR, credit=false)
            code = strip(read(script, String))
            mdpost(str) = replace(str, "@__CODE__" => code)
            Literate.markdown(input, GENERATEDDIR, postprocess=mdpost, credit=false)
            Literate.notebook(input, GENERATEDDIR, execute=true, credit=false)
        elseif any(endswith.(file, [".png", ".jpg", ".gif"]))
            cp(joinpath(model, file), joinpath(GENERATEDDIR, file); force=true)
        else
            @warn "ignoring $file in $model"
        end
    end
end
