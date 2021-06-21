# Load package and utility scripts

using ReachabilityAnalysis
using StructuralDynamicsODESolvers
const SD = StructuralDynamicsODESolvers
using ReachabilityAnalysis: solve

using Plots
using Plots.PlotMeasures
using LaTeXStrings

include("utils.jl")
include("data.jl")

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

    return Ahom
end

function solve_heat3d_rfem(; model=StepIntersect(Forward(inv=false)),
                             X0=nothing, case=1)

    # load vectors and params
    K, C, qextTamb, QhG, timeIncr = load_heat3d_matrices()
    rho, QFHmin, QFHmax, m, Tini, Tambmin, Tambvarmin, Tambvarmax, Tambvar, QFH = load_heat3d_params()
    usA, usB = load_heat3d_octaveSols()

    δ = timeIncr;  n = size(K, 1);

    NSTEPS = ( size(usA,2)-1 ) ;

    # compute "homogeneized" system matrix
    Ahom = AhomCalc( C, K, qextTamb, QhG, rho, m ) ;

    QFH     = (QFHmax+QFHmin)*0.5 ;
    Tambvar = (Tambvarmax+Tambvarmin)*0.5 ;

    if case == 1
        # opcion singleton
        X0 = Singleton( Tini*ones(n)) ×
            Singleton((Tambmin+Tambvar*0.5)*ones( 1 ) )× Singleton( QFH*ones(1) ) ×
            Singleton((-Tambvar*0.5)*ones( 1 ) )× Singleton( zeros(1) ) ;

        prob = @ivp(x' = Ahom*x, x(0) ∈ X0)
        alg = BOX(δ=δ, approx_model=model)
        sol = solve(prob, alg=alg, NSTEPS=NSTEPS)

        X0 = Singleton( vcat( Tini*ones(n),
                    (Tambmin+Tambvar*0.5)*ones(1) ,
                    QFHmax*ones(1) ,
                    (-Tambvar*0.5)*ones(1),
                    zeros(1) ) ) ;

        probOrbitMax = @ivp(x' = Ahom*x, x(0) ∈ X0)
        solOrbit = solve( probOrbitMax, NSTEPS=NSTEPS, alg=ORBIT(δ=δ));
        # --------------------------------
        return sol, solOrbit

    elseif case == 2

        X0 = Singleton( Tini*ones(n)) ×
             Interval(Tambvarmin * 0.5 + Tambmin, Tambvarmax * 0.5 + Tambmin ) ×
             Interval( QFHmin, QFHmax ) ×
             Interval(-Tambvarmax * 0.5, -Tambvarmin*0.5) ×
             Singleton( zeros(1) ) ;

        prob = @ivp(x' = Ahom*x, x(0) ∈ X0)
        alg = BOX(δ=δ, approx_model=model)
        sol = solve(prob, alg=alg, NSTEPS=NSTEPS)
        # --------------------------------

        # --- Orbit min ---
        X0 = Singleton(
              vcat( Tini*ones(n),
                    (Tambmin+Tambvarmin*0.5)*ones(1) ,
                    QFHmin*ones(1) ,
                    (-Tambvarmin*0.5)*ones(1),
                    zeros(1) ) ) ;

        probOrbit   = @ivp(x' = Ahom*x, x(0) ∈ X0)
        solOrbitmin = solve( probOrbit, NSTEPS=NSTEPS, alg=ORBIT(δ=δ));
        # --------------------------------

            # --- Orbit min ---
        X0 = Singleton(
              vcat( Tini*ones(n),
                    (Tambmin+Tambvarmax*0.5)*ones(1) ,
                    QFHmax*ones(1) ,
                    (-Tambvarmax*0.5)*ones(1),
                    zeros(1) ) ) ;

        probOrbit   = @ivp(x' = Ahom*x, x(0) ∈ X0)
        solOrbitmax = solve( probOrbit, NSTEPS=NSTEPS, alg=ORBIT(δ=δ));
        # --------------------------------

        return sol, solOrbitmin, solOrbitmax
   end

end

# -------------------------------------------------

maxtemp_heat3d(sol; nodes=1:dim(sol), kwargs...) = maxtemp(sol, nodes=nodes)
maxtemp_heat3d(; nodes=nothing, kwargs...) = maxtemp(solve_heat3d(kwargs...), nodes=nodes)
mintemp_heat3d(sol; nodes=1:dim(sol), kwargs...) = mintemp(sol, nodes=nodes)
mintemp_heat3d(; nodes=nothing, kwargs...) = mintemp(solve_heat3d(kwargs...), nodes=nodes)

# -------------------------------------------------
function solve_heat3d_implicit_euler()

    # load vectors and params
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

    solBackwardEulerSDCaseA = SD.solve(prob, alg, NSTEPS=NSTEPS)

    R    = [ QhG*QFHmax*rho*m*exp(-m*δ*(i-1))+qextTamb*(Tambmin+Tambvarmax*(0.5 + 0.5*sin( omega*δ*(i-1)-pi/2.0 ))) for i in 1:NSTEPS+1 ] ;

    sys  =  SecondOrderConstrainedLinearControlContinuousSystem(M, Matrix(C), Matrix(K), B, nothing, R)
    prob = InitialValueProblem(sys, (U₀, U₀))

    solBackwardEulerSDCaseB = SD.solve(prob, alg, NSTEPS=NSTEPS)


    return solBackwardEulerSDCaseA, solBackwardEulerSDCaseB
end

#-----------------------------------------------------------
# Solutions LSD
@time solSDCase1, solSDCase2 = solve_heat3d_implicit_euler();

#-----------------------------------------------------------
# Numerical solutions
tdom                     = times(solSDCase1);
solnumSDCase1            = displacements(solSDCase1)[1:end];
solnodeSDCase1(node)     = [ x[node] for x in solnumSDCase1 ] ;

tdom                     = times(solSDCase2);
solnumSDCase2            = displacements(solSDCase2)[1:end];
solnodeSDCase2(node)     = [ x[node] for x in solnumSDCase2 ] ;

#-----------------------------------------------------------
# Solutions Octave
solAOctave, solBOctave = load_heat3d_octaveSols();

#-----------------------------------------------------------
# Solutions SetBased

## Case 1
@time solRFEMCase1, solOrbitCase1 = solve_heat3d_rfem( case = 1 ) ;

## Case 2
@time solRFEMCase2, solOrbitCase2min, solOrbitCase2max = solve_heat3d_rfem( case = 2 ) ;
# Check maximum node temp
# @time maxtemp_heat3d(solRFEMCase2, nodes=1:1331)

#-----------------------------------------------------------
# Plots Case 1

timeMark = 200;

indTimeMark = 155 ;

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

savefig(fig, "heat3d_Case1.pdf")

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

savefig(fig,"heat3d_Case2.pdf")
