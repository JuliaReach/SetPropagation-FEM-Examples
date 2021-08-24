module SetPropagation_FEM_Examples

# See docs in NeuralNetworkAnalysis.jl
macro modelpath(model_path::String, name::String)
    __source__.file === nothing && return nothing
    _dirname = dirname(String(__source__.file))
    dir = isempty(_dirname) ? pwd() : abspath(_dirname)
    return joinpath(dir, name)
end

end # module
