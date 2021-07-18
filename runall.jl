# instantiate project
#import Pkg
#Pkg.activate(@__DIR__)
#Pkg.instantiate()

const OUTPUT_FOLDER = "output"
const RESULTS_FILE = "runtimes.csv"

function main()

    # --------------------------
    # Initialization
    # --------------------------
    if !isdir(OUTPUT_FOLDER)
        mkdir(OUTPUT_FOLDER)
        mkdir(joinpath(OUTPUT_FOLDER, "SDOF"))
        mkdir(joinpath(OUTPUT_FOLDER, "Clamped"))
        mkdir(joinpath(OUTPUT_FOLDER, "Heat1D"))
        mkdir(joinpath(OUTPUT_FOLDER, "Heat3D"))
    end
    global io = open(joinpath(OUTPUT_FOLDER, RESULTS_FILE), "w")
    global TARGET_FOLDER
    println("Running examples...")

    # --------------------------
    # Single-degree-of-freedom
    # --------------------------
    println("\n (1/4) Running single-degree-of-freedom (SDOF) model... \n")
    TARGET_FOLDER = joinpath(OUTPUT_FOLDER, "SDOF")
    include(joinpath("examples", "SDOF", "SDOF_Preliminaries.jl"))
    include(joinpath("examples", "SDOF", "SDOF_Time.jl"))
    include(joinpath("examples", "SDOF", "SDOF_PE_AD.jl"))

    # --------------------------
    # Clamped bar
    # --------------------------
    TARGET_FOLDER = joinpath(OUTPUT_FOLDER, "Clamped")
    println("\n (2/4) Running Clamped bar model\n")
    include(joinpath("examples", "Clamped", "Clamped.jl"))

    # --------------------------
    # Heat 1D
    # --------------------------
    TARGET_FOLDER = joinpath(OUTPUT_FOLDER, "Heat1D")
    println("\n (3/4) Running one-dimensional heat transfer model...\n")
    include(joinpath("examples", "Heat1D", "Heat1D_TempProfile.jl"))
    include(joinpath("examples", "Heat1D", "Heat1D_TempGradient.jl"))

    # --------------------------
    # Heat 3D
    # --------------------------
    TARGET_FOLDER = joinpath(OUTPUT_FOLDER, "Heat3D")
    println("\n (4/4) Running three-dimensional concrete casting model...\n")
    include(joinpath("examples", "Heat3D", "Heat3D.jl"))

    print(io, "\n")
    println("Finished running benchmarks.")
    close(io)
    nothing
end

main()
