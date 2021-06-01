# deactivate plot GUI, which is not available in Docker
ENV["GKSwstype"] = "100"

# instantiate project
import Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

const TARGET_FOLDER = "result"
const RESULTS_FILE = "results.csv"

function main()
    if !isdir(TARGET_FOLDER)
        mkdir(TARGET_FOLDER)
    end
    global io = open(joinpath(TARGET_FOLDER, RESULTS_FILE), "w")

    println("Running RFEM benchmarks...")

    # Production-destruction benchmark
    #println("###\nRunning production-destruction model benchmark\n###")
    #include("models/ProdDestruct/prod_destruct_benchmark.jl")

    print(io, "\n")# deactivate plot GUI, which is not available in Docker
ENV["GKSwstype"] = "100"

# instantiate project
import Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

const TARGET_FOLDER = "result"
const RESULTS_FILE = "results.csv"

function main()
    if !isdir(TARGET_FOLDER)
        mkdir(TARGET_FOLDER)
    end
    global io = open(joinpath(TARGET_FOLDER, RESULTS_FILE), "w")

    println("Running NLN benchmarks...")

    # Production-destruction benchmark
    println("###\nRunning production-destruction model benchmark\n###")
    include("models/ProdDestruct/prod_destruct_benchmark.jl")

    # Coupled Van der Pol benchmark
    println("###\nRunning Van der Pol benchmark\n###")
    include("models/VanDerPol/vanderpol_benchmark.jl")

    # Laub-Loomis benchmark
    println("###\nRunning Laub-Loomis benchmark\n###")
    include("models/LaubLoomis/laubloomis_benchmark.jl")

    # Quadrotor benchmark
    println("###\nRunning Quadrotor benchmark\n###")
    include("models/Quadrotor/quadrotor_benchmark.jl")

    # Lotka-Volterra tangential crossing benchmark
    println("###\nRunning Lotka-Volterra tangential crossing benchmark\n###")
    include("models/LotkaVolterra/lotka_volterra_benchmark.jl")

    # Spacecraft benchmark
    println("###\nSpacecraft benchmark\n###")
    include("models/Spacecraft/spacecraft_benchmark.jl")

    print(io, "\n")
    println("Finished running benchmarks.")
    close(io)
    nothing
end

main()
    println("Finished running benchmarks.")
    close(io)
    nothing
end

main()
