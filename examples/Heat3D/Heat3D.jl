# # Concrete casting heat of hydration
#
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/examples/Heat3D.ipynb)
#

using ReachabilityAnalysis, StructuralDynamicsODESolvers, MAT
using LinearAlgebra, SparseArrays
using LazySets.Arrays
using FastExpm

(@isdefined TARGET_FOLDER) ? nothing : TARGET_FOLDER = ""

# load internal (non-exported) functions
using ReachabilityAnalysis: ReachSolution, discretize, normalize, solve, center, step_size

include("/home/mforets/Projects/paper_rfem/examples/krylov_concept.jl")

# reduce default tolerance for vertex enumeration
LazySets.set_ztol(Float64, 1e-12)

using Plots
using Plots.PlotMeasures
using LaTeXStrings

using LazySets: center

function load_heat3d_matrices(; path=joinpath(@__DIR__, "Heat3D.mat"))
    vars = matread(path)

    K = vars["K"]
    C = vars["C"]
    timeIncr = vars["timeIncr"]
    fext = vars["fext"]
    qextTamb = vars["qextTamb"][:] # vector
    QhG =  vars["QhG"][:] # vector
    udots=vars["udots"]
    xs=vars["xs"]

    return K, C, qextTamb, QhG, timeIncr
end

function load_heat3d_params(; path=joinpath(@__DIR__, "params.mat"))
    vars       = matread(path)
    rho        = vars["rho"]
    m          = vars["m"]
    QFH        = vars["QFH"]
    Tambmin    = vars["Tambmin"]
    Tambvar    = vars["Tambvar"]
    Tambvarmin = vars["Tambvarmin"]
    Tambvarmax = vars["Tambvarmax"]
    Tini       = vars["Tini"]
    QFHmin     = vars["QFHmin"]
    QFHmax     = vars["QFHmax"]

    return rho, QFHmin, QFHmax, m, Tini, Tambmin, Tambvarmin, Tambvarmax, Tambvar, QFH
end

function load_heat3d_octaveSols(; path=joinpath(@__DIR__, "solsOctave.mat"))
    vars = matread(path)
    solAOctave = vars["Ts3DA"]
    solBOctave = vars["Ts3DB"]
    return solAOctave, solBOctave
end


# compute the maximum temperature at the given nodes
# by default we compute the max temperature in all the nodes,
# but you can specify a vector of indices to search;
# for example pass nodes=1:100 to make the search only on the first 100 nodes
function maxtemp(sol::ReachSolution; nodes=nothing)
    ## initialize search
    Tmax = -Inf; kmax = -1; imax = -1

    if isnothing(nodes)
        nodes = 1:dim(sol)
    end

    ## for each reach-set we compute the max bounds
    ## in every direction
    for (k, R) in enumerate(sol)
        X = set(R)
        q = center(X) + radius_hyperrectangle(X)
        θi = argmax(view(q, nodes))
        θ = q[θi]
        if θ > Tmax
            Tmax = θ
            kmax = k
            imax = θi
        end
    end

    δ = step_size(sol.alg)
    dtmax = tspan(sol[kmax])
    return (Tmax=Tmax, node=imax, step=kmax, dt=dtmax)
end

# compute the minimum temperature on the given nodes
function mintemp(sol::ReachSolution; nodes=nothing)
    ## initialize search
    Tmin = Inf; kmin = -1; imin = -1

    if isnothing(nodes)
        nodes = 1:dim(sol)
    end

    ## for each reach-set we compute the max bounds
    ## in every direction
    for (k, R) in enumerate(sol)
        X = set(R)
        q = center(X) - radius_hyperrectangle(X)
        θi = argmax(view(q, nodes))
        θ = q[θi]
        if θ < Tmin
            Tmin = θ
            kmin = k
            imin = θi
        end
    end

    δ = step_size(sol.alg)
    dtmin = tspan(sol[kmin])
    return (Tmin=Tmin, node=imin, step=kmin, dt=dtmin)
end

## Control nodes
controlNodeA = 11^2*6+1
controlNodeB = 11^2*9+1
nNodes = 11^3

function AhomCalc(C, K, qextTamb, QhG, rho, m)

    invC = inv( Matrix(C) );
    A = - invC * K |> Matrix
    n = size(A, 1)

    b0   =  invC * qextTamb ;
    bext =  invC * QhG * rho * m ;
    bsin =  invC * qextTamb ;

    omega = pi/12.;

    Ahom = [ A                    b0  bext  bsin      zeros(n) ;
             transpose( zeros(n))  0   0    0         0 ;
             transpose( zeros(n))  0  -m    0         0 ;
             transpose( zeros(n))  0   0    0         1 ;
             transpose( zeros(n))  0   0    -omega^2  0 ] ;

    return sparse(Ahom) # cast to sparse matrix
end

function solve_heat3d_rfem(; model=StepIntersect(Forward(inv=false)), X0=nothing, case=1)

    ## load vectors and params
    K, C, qextTamb, QhG, timeIncr = load_heat3d_matrices()
    rho, QFHmin, QFHmax, m, Tini, Tambmin, Tambvarmin, Tambvarmax, Tambvar, QFH = load_heat3d_params()
    usA, usB = load_heat3d_octaveSols()

    δ = timeIncr
    n = size(K, 1)

    NSTEPS = size(usA,2) - 1

    ## compute "homogeneized" system matrix
    Ahom = AhomCalc(C, K, qextTamb, QhG, rho, m)

    QFH     = (QFHmax+QFHmin)*0.5
    Tambvar = (Tambvarmax+Tambvarmin)*0.5

    if case == 1
        ## opcion singleton
        X0 = Singleton( Tini*ones(n)) ×
            Singleton((Tambmin+Tambvar*0.5)*ones( 1 ) )× Singleton( QFH*ones(1) ) ×
            Singleton((-Tambvar*0.5)*ones( 1 ) )× Singleton( zeros(1) )
        X0 = Singleton(LazySets.center(box_approximation(X0)))

    else
        ## opcion distribuida
        X0 = Singleton( Tini*ones(n)) ×
             Interval(Tambvarmin * 0.5 + Tambmin, Tambvarmax * 0.5 + Tambmin ) ×
             Interval( QFHmin, QFHmax ) ×
             Interval(-Tambvarmax * 0.5, -Tambvarmin*0.5) ×
             Singleton( zeros(1) )
        X0 = box_approximation(X0)
    end

    # set propagation method
    method = :sf
    println("Solving with set prop method (Case $case) $method")
    if method == :box
        prob = @ivp(x' = Ahom*x, x(0) ∈ X0)
        alg = BOX(δ=δ, approx_model=model)
        @time sol = solve(prob, alg=alg, NSTEPS=NSTEPS)

    elseif method == :sf
        @time sol = _solve_heat3d_sf(Ahom, X0, δ, NSTEPS, nodeidx=[controlNodeA, controlNodeB])
    end

    if case == 1
        X0 = Singleton(vcat(Tini*ones(n),
                      (Tambmin+Tambvar*0.5)*ones(1) ,
                      QFHmax*ones(1) ,
                      (-Tambvar*0.5)*ones(1),
                      zeros(1)))

        probOrbitMax = @ivp(x' = Ahom*x, x(0) ∈ X0)

        println("Solving with set prop method (Case 1) ORBIT")
        @time solOrbit = solve( probOrbitMax, NSTEPS=NSTEPS, alg=ORBIT(δ=δ))
        return sol, solOrbit

    elseif case == 2

        ## Orbit min
        X0 = Singleton(
              vcat( Tini*ones(n),
                    (Tambmin+Tambvarmin*0.5)*ones(1) ,
                    QFHmin*ones(1) ,
                    (-Tambvarmin*0.5)*ones(1),
                    zeros(1) ) ) ;

        probOrbit   = @ivp(x' = Ahom*x, x(0) ∈ X0)
        println("Solving with set prop method (Case 2) ORBIT min")
        @time solOrbitmin = solve( probOrbit, NSTEPS=NSTEPS, alg=ORBIT(δ=δ));

        ## Orbit min
        X0 = Singleton(
              vcat( Tini*ones(n),
                    (Tambmin+Tambvarmax*0.5)*ones(1) ,
                    QFHmax*ones(1) ,
                    (-Tambvarmax*0.5)*ones(1),
                    zeros(1) ) ) ;

        probOrbit = @ivp(x' = Ahom*x, x(0) ∈ X0)
        println("Solving with set prop method (Case 2) ORBIT max")
        @time solOrbitmax = solve( probOrbit, NSTEPS=NSTEPS, alg=ORBIT(δ=δ))
        return sol, solOrbitmin, solOrbitmax
   end

end


maxtemp_heat3d(sol; nodes=1:dim(sol), kwargs...) = maxtemp(sol, nodes=nodes)
maxtemp_heat3d(; nodes=nothing, kwargs...) = maxtemp(solve_heat3d(kwargs...), nodes=nodes)
mintemp_heat3d(sol; nodes=1:dim(sol), kwargs...) = mintemp(sol, nodes=nodes)
mintemp_heat3d(; nodes=nothing, kwargs...) = mintemp(solve_heat3d(kwargs...), nodes=nodes)

function solve_heat3d_implicit_euler()

    ## load vectors and params
    K, C, qextTamb, QhG, timeIncr = load_heat3d_matrices()
    rho, QFHmin, QFHmax, m, Tini, Tambmin, Tambvarmin, Tambvarmax, Tambvar, QFH = load_heat3d_params()
    usA, usB = load_heat3d_octaveSols()

    δ = timeIncr;  n = size(K, 1)
    NSTEPS = ( size(usA,2)-1 ) ;

    omega = pi/12.;

    M    = zeros(n, n)
    Tambvar = (Tambvarmin + Tambvarmax ) * 0.5 ;
    QFH = (QFHmin+ QFHmax)*0.5;

    R    = [ QhG*QFH*rho*m*exp(-m*δ*(i-1))+qextTamb*(Tambmin+Tambvar*(0.5 + 0.5*sin( omega*δ*(i-1)-pi/2.0 ))) for i in 1:NSTEPS+1 ] ;

    B   = Diagonal(ones(n));

    sys  =  SecondOrderConstrainedLinearControlContinuousSystem(M, Matrix(C), Matrix(K), B, nothing, R)

    U₀ = Tini*ones(n) ;

    prob = InitialValueProblem(sys, (U₀, U₀))
    alg  = BackwardEuler(Δt=δ)

    println("Solving with Backward Euler (Case 1)")
    @time solBackwardEulerSDCaseA = solve(prob, alg, NSTEPS=NSTEPS)

    R    = [ QhG*QFHmax*rho*m*exp(-m*δ*(i-1))+qextTamb*(Tambmin+Tambvarmax*(0.5 + 0.5*sin( omega*δ*(i-1)-pi/2.0 ))) for i in 1:NSTEPS+1 ] ;

    sys  =  SecondOrderConstrainedLinearControlContinuousSystem(M, Matrix(C), Matrix(K), B, nothing, R)
    prob = InitialValueProblem(sys, (U₀, U₀))

    println("Solving with Backward Euler (Case 2)")
    @time solBackwardEulerSDCaseB = solve(prob, alg, NSTEPS=NSTEPS)


    return solBackwardEulerSDCaseA, solBackwardEulerSDCaseB
end



function _solve_heat3d_sf(A, X0, δ, NSTEPS; nodeidx, mk=70, tol=1e-6)
    num_nodes = isa(nodeidx, Int) ? 1 : length(nodeidx)

    Eplus = _build_Eplus_optim(A, X0, δ; m=mk, tol=tol, cutoff=eps(Float64))

    dirs = Vector{LazySets.Arrays.SingleEntryVector{Float64}}()
    k = size(A, 1)
    for j in 1:num_nodes
        node = nodeidx[j]
        push!(dirs, SingleEntryVector(node, k, 1.0))
        push!(dirs, SingleEntryVector(node, k, -1.0))
    end
    cdirs = CustomDirections(dirs)

    propagacion = :expm_fast
    Ω0, Φᵀ = _discretize(A, δ, propagacion, X0, Eplus, mk, tol)

    Mout = Matrix{Float64}(undef, 2*num_nodes, NSTEPS)


    Aᵀδ = copy(transpose(A .* δ))

    for j in 1:num_nodes
        node = nodeidx[j]
        d = zeros(k)
        d[node] = 1.0
        println("Propagando indice $node ...")
        idx = (2(j-1) + 1):2j

        if isnothing(Φᵀ)
            _reach_homog_krylov_LGG09_modif!(view(Mout, idx, :), Ω0, Aᵀδ, d, NSTEPS; m=mk, tol=tol)
        else
            _reach_homog_dir_LGG09_modif!(view(Mout, idx, :), Ω0, Φᵀ, d, NSTEPS)
        end
    end
    sol = [TemplateReachSet(cdirs, view(Mout, :, k), (0.0 .. δ) + (k-1) * δ) for k in 1:NSTEPS] |> Flowpipe
end

function _discretize(A, δ, propagacion, X0::AbstractSingleton, Eplus, mk, tol)
    ΦX0 = _expv(A, X0.element, 1, δ, m=mk, tol=tol) |> Singleton
    Ω0 = CH(X0, minkowski_sum(ΦX0, Eplus))
    return Ω0, nothing
end

function _discretize(A, δ, propagacion, X0::Hyperrectangle, Eplus, mk, tol)
    if propagacion == :expm_base
        Φ = exp(Matrix(transpose(A) .* δ))
        Φ = sparse(Φ)
        droptol!(Φ, eps(Float64))
    elseif propagacion == :expm_fast
        Φ = real.(fastExpm(A .* δ, threshold=eps(Float64), nonzero_tol=eps(Float64)) )
    end
    Ω0 = CH(X0, Φ*X0 ⊕ Eplus)
    return Ω0, transpose(Φ)
end

function _reach_homog_dir_LGG09_modif!(Mout, Ω₀, Φᵀ, d, NSTEPS)
    rᵢ = copy(d)
    rᵢ₊₁ = similar(rᵢ)

    @inbounds for i in 1:NSTEPS
        Mout[1, i] = ρ(rᵢ, Ω₀)
        Mout[2, i] = ρ(-rᵢ, Ω₀)

        # update cache for the next iteration
        mul!(rᵢ₊₁, Φᵀ, rᵢ)
        copy!(rᵢ, rᵢ₊₁)
    end
    return nothing
end


#-----------------------------------------------------------
# Solutions LSD
solSDCase1, solSDCase2 = solve_heat3d_implicit_euler()

# Numerical solutions
tdom                     = times(solSDCase1)
solnodeSDCase1(node)     = displacements(solSDCase1, node)

tdom                     = times(solSDCase2)
solnodeSDCase2(node)     = displacements(solSDCase2, node)

#-----------------------------------------------------------
# Solutions Octave
solAOctave, solBOctave = load_heat3d_octaveSols()

#-----------------------------------------------------------
# Solutions SetBased

## Case 1
solRFEMCase1, solOrbitCase1 = solve_heat3d_rfem(case=1)

## Case 2
solRFEMCase2, solOrbitCase2min, solOrbitCase2max = solve_heat3d_rfem(case=2);

# Check maximum node temp
# @time maxtemp_heat3d(solRFEMCase2, nodes=1:1331)

#-----------------------------------------------------------
# Plots Case 1

timeMark = 200

indTimeMark = 155

tdom[indTimeMark]
solnodeSDCase1(controlNodeA)[indTimeMark]

fig = plot(tdom, solnodeSDCase1(controlNodeA),
           xlab=L"\textrm{Time (h)}",
           lc=:red, ylab=L"\textrm{Temperature } (^\circ\!C)",
           lab="",
           linestyle=:solid,
           legendfontsize=15, legend=:bottomright,
           tickfont=font(18, "Times"), guidefontsize=18,
           xtick=([0, 50, 100, 150, 200], [L"0", L"50", L"100", L"150", L"200"]),
           ytick=([20, 30, 40, 50, 60, 70, 80], [L"20", L"30", L"40", L"50", L"60", L"70", L"80"]),
           xguidefont=font(18, "Times"),
           yguidefont=font(18, "Times"),
           bottom_margin=6mm, left_margin=8mm, right_margin=4mm, top_margin=3mm,
           size=(900, 600))

plot!( tdom, solnodeSDCase1(controlNodeB), lc=:blue, lab="" )

# markers
plot!( [tdom[indTimeMark]], [solnodeSDCase1(controlNodeA)[indTimeMark]], marker=:utriangle, lc=:red, c=:red, lab=L"\textrm{Backward Euler at A}" )
plot!( [tdom[indTimeMark]], [solnodeSDCase1(controlNodeB)[indTimeMark]], marker=:circle, lc=:blue, c=:blue, lab=L"\textrm{Backward Euler at B}" )

plot!( fig, solRFEMCase1, vars=(0, controlNodeA), alpha=.2, lc=:magenta, c=:magenta, lab=L"\textrm{Flowpipe at A}")
plot!( fig, solRFEMCase1, vars=(0, controlNodeB), alpha=.2, lc=:green, c=:green, lab=L"\textrm{Flowpipe at B}")
#plot!(fig, tdom, solAOctave[controlNodeA,1:end], lc=:blue, lab="Backward Euler point A")

# subida
lens!(fig, [9, 12], [39, 43], inset = (1, bbox(0.15, 0.69, 0.28, 0.24)),
           subplot=5,
           tickfont=font(18, "Times"),
           xticks=([9.0, 10.5, 12.0], [ L"9",L"10.5", L"12"]),
           yticks=([39, 41, 43], [ L"39", L"41", L"43"])
)

lens!(fig, [50, 52], [66.5, 67], inset = (1, bbox(0.2, 0.38, 0.18, 0.18)),
           subplot=6,
           tickfont=font(18, "Times"),
           xticks=([50, 51, 52], [L"50", L"51", L"52"]),
           yticks=([66.5, 66.75, 67], [ L"66.5", "", L"67"])
)

lens!(fig, [50, 52], [83.3, 84.5], inset = (1, bbox(0.7, 0.1, 0.2, 0.2)),
           subplot=7,
           tickfont=font(18, "Times"),
           xticks=([50, 51, 52], [L"50", L"51", L"52"]),
           yticks=([83.5, 84.0, 84.5], [ L"83.5", "", L"84.5"])
)

savefig(fig, joinpath(TARGET_FOLDER, "heat3d_Case1.pdf"))

#-----------------------------------------------------------
# Plots Case 2

fig = plot(legendfontsize=16, legend=:topright,
           tickfont=font(18, "Times"), guidefontsize=18,
           xtick=([0, 50, 100, 150, 200], [L"0", L"50", L"100", L"150", L"200"]),
           ytick=([20, 30, 40, 50, 60, 70, 80], [L"20", L"30", L"40", L"50", L"60", L"70", L"80"]),
           xguidefont=font(18, "Times"),
           yguidefont=font(18, "Times"),
           bottom_margin=6mm, left_margin=8mm, right_margin=4mm, top_margin=3mm, size=(900, 600))

plot!(fig, solRFEMCase2, vars=(0, controlNodeA), alpha=.2, lc=:red, c=:red,
    xlab= L"\textrm{Time (h)}",
    ylab = L"\textrm{Temperature } (^\circ\!C)",
    lab=L"\textrm{Flowpipe at A}")#

#plot!(tdom, solnodeSDCase2(controlNodeA), xlab="Time", lc=:red, ylab="Temperature", lab="Backward Euler")
#plot!(fig, tdom, solBOctave[controlNodeA,1:end], lc=:blue, lab="Backward Euler point A")
#plot!( fig, solOrbitCase2max, vars=(0, controlNodeA), alpha=.2, lc=:red, lab="OrbitMax", setiestype=:path)
#plot!( fig, solOrbitCase2min, vars=(0, controlNodeA), alpha=.2, lc=:red, lab="OrbitMin", setiestype=:path)
plot!(fig,  solRFEMCase2, vars=(0, controlNodeB), alpha=.2, lc=:blue, c=:blue, lab=L"\textrm{Flowpipe at B}")#
plot!(fig, tdom, solnodeSDCase1(controlNodeA), lc=:red, lab=L"\textrm{Backward Euler at A}" )
plot!(fig, tdom, solnodeSDCase1(controlNodeB), lc=:blue, lab=L"\textrm{Backward Euler at B}")

savefig(fig, joinpath(TARGET_FOLDER, "heat3d_Case2.pdf"))

#=

julia>     Ahom = AhomCalc(C, K, qextTamb, QhG, rho, m)
1335×1335 Matrix{Float64}:
 -7.28091  0.77986  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  1.55972  …  0.0  0.0  0.0  0.0  0.0  0.0  4.16147  0.0082213   4.16147    0.0
  ⋮                                ⋮                        ⋮             ⋱  ⋮                        ⋮

  0.0      0.0      0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0         0.0  0.0  0.0  0.0  0.0  0.0  0.0      0.0        -0.0685389  0.0

julia> sparse(Ahom)
1335×1335 SparseMatrixCSC{Float64, Int64} with 10787 stored entries:
⠻⣦⡀⠀⠀⢸
⠀⠈⠻⣦⡀⢸
⠀⠀⠀⠈⠻⣾

julia> 1335 * 1335
1782225

julia> Ahom * Ahom;

julia> Ahom * Ahom;

julia> @time Ahom * Ahom;
  0.100113 seconds (2 allocations: 13.597 MiB)

julia> @time Ahom * Ahom;
  0.109501 seconds (2 allocations: 13.597 MiB)

julia> Ahom_sp = sparse(Ahom);

julia> @time Ahom_sp * Ahom_sp;
  0.133171 seconds (126.53 k allocations: 8.742 MiB, 99.41% compilation time)

julia> @time Ahom_sp * Ahom_sp;
  0.000643 seconds (7 allocations: 1.460 MiB)

=#
