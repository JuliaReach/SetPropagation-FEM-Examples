# # Single degree-of-freedom

# ### Equations of motion

using ReachabilityAnalysis
using StructuralDynamicsODESolvers
using ReachabilityAnalysis: discretize
using Plots, Plots.PlotMeasures, LaTeXStrings

# Struct that holds a problem describing an harmonic oscillator with frequency ω:
# x''(t) + ω^2 x(t) = 0
#
# solution x(t) = Acos(ωt + B), v(t) = -ωAsin(ωt + B)
# x(0) = Acos(B)
# v(0) = -ωAsin(B)
#
# special case: if x(0) = A, v(0) = 0 => x(t) = A cos ωt
struct SDOF{T, ST}
    ivp::T
    X0::ST
    ω::Float64
    Ampl::Float64
    T::Float64
end

amplitude(p::SDOF) = p.Ampl
frequency(p::SDOF) = p.ω
period(p::SDOF) = p.T
MathematicalSystems.initial_state(p::SDOF) = p.X0
MathematicalSystems.InitialValueProblem(p::SDOF) = p.ivp
MathematicalSystems.state_matrix(p::SDOF) = state_matrix(p.ivp)

function analytic_solution(p::SDOF; A=p.Ampl, B=0.0, x0=nothing, v0=nothing)
    @assert p.X0 isa Singleton "the analytic solution requires a singleton initial condition"

    ω = p.ω
    if !isnothing(x0) && !isnothing(v0)
        A = sqrt(x0^2 + v0^2 / ω^2)
        B = atan(-v0 / (ω*x0))
    end
    return t -> A * cos(ω * t + B)
end

function analytic_derivative(p::SDOF; A=p.Ampl, B=0.0, x0=nothing, v0=nothing)
    @assert p.X0 isa Singleton "the analytic solution requires a singleton initial condition"

    ω = p.ω
    if !isnothing(x0) && !isnothing(v0)
        A = sqrt(x0^2 + v0^2 / ω^2)
        B = atan(-v0 / (ω*x0))
    end
     return t -> -ω * A * sin(ω * t + B)
end

function InitialValueProblem_quad(p::SDOF)
    M = hcat(1.0)
    ω = p.ω
    K = hcat(ω^2)
    C = zeros(1, 1)
    R = zeros(1)
    S = SecondOrderAffineContinuousSystem(M, C, K, R)
    @assert p.X0 isa Singleton
    x0 = element(p.X0)
    U0 = [x0[1]]
    U0dot = [x0[2]]
    return IVP(S, (U0, U0dot))
end

function sdof(; T=0.5,     # period
                Ampl=1.0,  # amplitude
                X0 = Singleton([Ampl, 0.0])) # initial condition
    ## frequency
    ω = 2π / T

    ## cast as a first-order system:
    ## x' = v
    ## v' = -ω^2 * x
    A = [ 0.0     1.0;
         -ω^2     0.0]

    prob = @ivp(X' = AX, X(0) ∈ X0)
    return SDOF(prob, X0, ω, Ampl, T)
end

## Funciones para calcular PE y AD

function consecutivos(indis)

    diffdx1 = diff(indis);
    lista_indx = [];

    j = 1
    while j <= length(indis)
        inds = [indis[j]]
        while (j < length(indis)) && (diffdx1[j]==1)
            j+=1
            push!(inds, indis[j])
        end
        j+=1
        push!(lista_indx, inds)
    end

    return lista_indx

end

function period_evolution2(T, tf, α; X0 = Singleton([1.0, 0.0]), polardirs = PolarDirections(200))

    ω = 2π / T

    ## x' = v
    ## v' = -ω^2 * x
    A = [0.0     1.0;
         -ω^2    0.0]

    prob = @ivp(X' = AX, X(0) ∈ X0)
    guard = HPolyhedron([HalfSpace([ 0.0,  1.0], 0.0),    # v ≤ 0
                         HalfSpace([ 0.0, -1.0], 0.0),    # v ≥ 0
                         HalfSpace([-1.0,  0.0], 0.0)])   # x ≥ 0

    alg(αi) = VREP(δ=αi*T, approx_model=StepIntersect(setops=:concrete))
    ## alg(αi) = GLGM06(δ=αi*T))
    ## alg(αi) = LGG09(δ=αi*T, template=polardirs))

    conjuntos = []
    for i in 1:length(α)
        sol = solve(prob, tspan=(0.0, tf), alg=alg(α[i]))
        ind_i = findall(R -> !isdisjoint(set(R), guard), array(sol)) #TODO sol[2:end]
        push!(conjuntos, [sol[x] for x in consecutivos(ind_i[2:end])])
    end

    return conjuntos
end

function AD_numerica(xb, N, npers)
    ts_barrido = Int(N); #ultimos N periodos
    return (1-(maximum(xb[end-ts_barrido:end]))^(1/(Int(npers)-N)))*100;
end

function _compute_AD(T, N, numpers, α)
    tf = numpers*T

    @time conjs = period_evolution2(T, tf, α)

    Ts = [[tspan(R) for R in c] for c in conjs] #intervalos de tiempo en los que se da una vuelta
    periodos = [[R/i for (i, R) in enumerate(c)] for c in Ts] #divido entre el numero de vueltas hechas para tener los periodos

    amplitudes = [[overapproximate(project(overapproximate(ConvexHullArray(set.(R)), Hyperrectangle), [1]), Interval).dat for R in c] for c in conjs]

    amps_min = [inf.(a1) for a1 in amplitudes]
    ## amps_max = [sup.(a1) for a1 in amplitudes]

    amps_min_num = AD_numerica.(amps_min, N, numpers)

    return amps_min_num
end

function _compute_PE(T, N, numpers, α)

    tf = numpers*T

    conjs = period_evolution2(T, tf, α)

    Ts = [[tspan(R) for R in c] for c in conjs] #intervalos de tiempo en los que se da una vuelta
    periodos = [[R/i for (i, R) in enumerate(c)] for c in Ts]

    intersecciones = [];
    for i in 1:length(periodos)
        z = [];
        for j in 1:length(periodos[i])
            push!(z, reduce(∩, periodos[i][1:j]))
        end
        push!(intersecciones, z)
    end

    ## ¿Todos los periodos limite contienen al T del sistema?
    @assert all([T ∈ interss[end] for interss in intersecciones])

    x = [Singleton([alphai]) × Interval((Ai[end]-T)/T*100) for (alphai, Ai) in zip(α, intersecciones)];
    y1 = UnionSetArray([ConvexHull(x[k], x[k+1]) for k in 1:length(x)-1]);

    PE_RFEM_sup = ([sup(Ai[end]) for Ai in intersecciones].-T)./T*100;
    PE_RFEM_inf = ([inf(Ai[end]) for Ai in intersecciones].-T)./T*100;

    return PE_RFEM_inf, PE_RFEM_sup
end

## Setup

p = sdof()

prob = IVP(p)
q = InitialValueProblem_quad(p)

@show x0 = initial_state(prob)

@show T = period(p)

tmax = 4.0
dom = 0:0.001:tmax

# some plotting options
COLOR_initial = :green
COLOR_transformed = :yellow
COLOR_convexhull = :blue
COLOR_Ω0 = :orange1
COLOR_Ω0_box = :lightblue
COLOR_analytic = :magenta
COLOR_numeric = :red
COLOR_zono = :chartreuse;

#-

# ## Singleton initial conditions

# ### Discretization

model = StepIntersect(setops=:concrete)

α = 0.05
@show δ = α*T

@show x0
probn = normalize(prob)
probd = discretize(probn, δ, model)

Φ = state_matrix(probd)
@show Ω0 = initial_state(probd);


function plot_singleton_discretization()
    fig = plot(legend=:topleft, xlab=L"u(t)", ylab=L"v(t)",
               legendfontsize=15,
               tickfont=font(20, "Times"),
               guidefontsize=20,
               xguidefont=font(30, "Times"),
               yguidefont=font(30, "Times"),
               xtick=([0.95, 0.96, 0.97, 0.98, 0.99, 1.0], [L"0.95", L"0.96", L"0.97", L"0.98", L"0.99", L"1.0"]),
               ytick=([0.0, -1.0, -2.0, -3.0, -4.0], [L"0.0", L"-1.0", L"-2.0", L"-3.0", L"-4.0"]),
               xlims=(0.94, 1.0),
               ylims=(-4.0, 1.0),
               bottom_margin=6mm, left_margin=2mm, right_margin=4mm, top_margin=3mm, size=(900, 600))
    plot!(fig, box_approximation(Ω0), lab=L"\textrm{Box}(\Omega_0)", c=COLOR_Ω0_box)
    plot!(fig, Ω0, lab=L"\Omega_0", c=COLOR_Ω0, alpha=0.6)

    ## analytic solution of the ODE
    tdom1 = 0:0.0001:α*T

    asol_x = analytic_solution(p)
    asol_v = analytic_derivative(p)

    ## refine the domain of interest
    plot!(fig, asol_x.(tdom1), asol_v.(tdom1), seriestype=:path, c=COLOR_analytic, lw=1.5, lab=L"\textrm{Analytic} = \Omega_0^e")

    ## initial state: singleton
    plot!(fig, x0, c=COLOR_initial, markersize=12, alpha=0.5, markershape=:circle, lab=L"\textrm{Initial } \mathcal{X}_0")

    ## transformed state: linear map of the singleton
    plot!(fig, Φ * x0, c=COLOR_transformed, markersize=12, alpha=0.5, markershape=:circle, lab=L"\textrm{Transformation }\Phi \mathcal{X}_0")

    ## initial state: singleton
    plot!(fig, CH(x0, Φ * x0), c=COLOR_convexhull, markersize=0, linestyle=:solid, lw=2.0, markershape=:circle, lab=L"\textrm{Convex hull CH}(\mathcal{X}_0, \Phi \mathcal{X}_0)")

    ## numerical solution with smaller step size
    algn = Newmark(Δt=α*T/5, α=0.25, δ=0.5)
    soln = solve(q, algn, NSTEPS=round(Int, tmax / (α*T)));
    plot!(fig, soln.U[1:10], soln.U′[1:10], markershape=:utriangle, markersize=8, color=COLOR_numeric, lab="")
    plot!(fig, soln.U[11:11], soln.U′[11:11], markershape=:utriangle, markersize=8, color=COLOR_numeric, lab=L"\textrm{Numeric}")

    return fig
end

fig = plot_singleton_discretization()
savefig(fig, "sdof_singleton_discretization.pdf")

# ### Flowpipe

model = StepIntersect(setops=:concrete)
α = 0.05

alg = VREP(δ=α*T, approx_model=model)
solvrep = solve(prob, tspan=(0.0, tmax), alg=alg);

alg = GLGM06(δ=α*T, approx_model=model)
solglg = solve(prob, tspan=(0.0, tmax), alg=alg);

alg = BOX(δ=α*T, approx_model=model)
solbox = solve(prob, tspan=(0.0, tmax), alg=alg);

function plot_singleton_flowpipe()
    fig = plot(legend=:topleft, xlab=L"u(t)", ylab=L"v(t)",
               legendfontsize=15,
               tickfont=font(20, "Times"),
               guidefontsize=20,
               xguidefont=font(30, "Times"),
               yguidefont=font(30, "Times"),
               xtick=([-1.0, -0.5, 0.0, 0.5, 1.0], [L"-1.0", L"-0.5", L"0.0", L"0.5", L"1.0"]),
               ytick=([0.0, -3.0, -6.0, -9.0, -12.0], [L"0.0", L"-3.0", L"-6.0", L"-9.0", L"-12.0"]),
               xlims=(-1.1, 1.2),
               ylims=(-14, 1.4),
               bottom_margin=6mm, left_margin=2mm, right_margin=4mm, top_margin=8mm, size=(900, 600))

    NMAX = 9

    ## propagation of Ω0 without overapproximation
    plot!(fig, solbox[1:NMAX], vars=(1, 2), c=COLOR_Ω0_box, lab=L"\textrm{Flowpipe (hyperrectangles)}")

    ## propagation of Ω0 using zonotopes
    plot!(fig, solglg[1:NMAX], vars=(1, 2), c=COLOR_zono, lw=1.0, alpha=.6, lab=L"\textrm{Flowpipe (zonotopes)}")

    ## propagation of Ω0 withot overapproximation
    plot!(fig, solvrep[1:NMAX], vars=(1, 2), c=COLOR_Ω0, alpha=0.8, lw=1.0, lab=L"\textrm{Flowpipe (polygons)}")

    ## analytic solution of the ODE
    tdom1 = 0:0.0001:NMAX*α*T
    asol_x = analytic_solution(p)
    asol_v = analytic_derivative(p)
    ## refine the domain of interest
    plot!(fig, asol_x.(tdom1), asol_v.(tdom1), seriestype=:path, c=COLOR_analytic, lw=2.0, lab=L"\textrm{Analytic} = F^e")

    ## numerical solution with smaller step size
    algn = Newmark(Δt=α*T/5, α=0.25, δ=0.5)
    soln = solve(q, algn, NSTEPS=round(Int, tmax / (α*T)));
    plot!(fig, soln.U[1:1:NMAX*5], soln.U′[1:1:NMAX*5], markersize=3.0, markershape=:utriangle, color=COLOR_numeric, lab="")
    plot!(fig, soln.U[NMAX*10+1:NMAX*10+1], soln.U′[NMAX*10+1:NMAX*10+1], markershape=:utriangle, color=COLOR_numeric, lab=L"\textrm{Numeric}")
    ##  \Delta t = \delta / 5

    fig
end

fig = plot_singleton_flowpipe()
savefig(fig, "sdof_singleton_flowpipe.pdf")

# ## Distributed initial conditions

# ### Discretization

p = sdof()

prob = IVP(p)
q = InitialValueProblem_quad(p)

X0 = Hyperrectangle(low=[0.9, -0.1], high=[1.1, 0.1])

probset = IVP(prob.s, X0);

@show T = period(p)

tmax = 4.0
dom = 0:0.001:tmax

model = StepIntersect(setops=:concrete)

α = 0.05
δ = α*T

probset_d = discretize(normalize(probset), δ, model)

Φ = state_matrix(probset_d)
Ω0 = initial_state(probset_d);

#-

import Random

Random.seed!(3333)

function plot_distributed_discretization()
    fig = plot(legend=:left, xlab=L"u(t)", ylab=L"v(t)",
               legendfontsize=15,
               tickfont=font(20, "Times"),
               guidefontsize=20,
               xguidefont=font(30, "Times"),
               yguidefont=font(30, "Times"),
               xtick=([0.85, 0.90, 0.95, 1.00, 1.05, 1.10], [L"0.85", L"0.90", L"0.95", L"1.00", L"1.05", L"1.10"]),
               ytick=([0.0, -1.0, -2.0, -3.0, -4.0], [L"0.0", L"-1.0", L"-2.0", L"-3.0", L"-4.0"]),
               xlims=(0.84, 1.12),
               ylims=(-4.5, 0.5),
               bottom_margin=6mm, left_margin=2mm, right_margin=4mm, top_margin=3mm, size=(900, 600))

    plot!(fig, box_approximation(Ω0), lab=L"\textrm{Box}(\Omega_0)", c=COLOR_Ω0_box)
    plot!(fig, Ω0, lab=L"\Omega_0", c=COLOR_Ω0, alpha=0.6)

    ## convex hull of the initial state and the transformed state
    plot!(fig, CH(X0, Φ * X0), c=COLOR_convexhull, lab=L"\textrm{Convex hull CH}(\mathcal{X}_0, \Phi \mathcal{X}_0)")

    ## initial state: singleton
    plot!(fig, X0, c=COLOR_initial, alpha=1., lab=L"\textrm{Initial } \mathcal{X}_0")


    ## transformed state: linear map of the singleton
    plot!(fig, Φ * X0, c=COLOR_transformed, alpha=1., lab=L"\textrm{Transformation }\Phi \mathcal{X}_0")

    ## -------- Random sample of analytic solutions ------
    tdom1 = 0:0.0001:α*T

    nsamples = 50
    X0sampled = sample(X0, nsamples);
    ω = p.ω

    s = X0sampled[1]
    asol_x = analytic_solution(p, x0=s[1], v0=s[2])
    asol_v = analytic_derivative(p, x0=s[1], v0=s[2])
    plot!(fig, asol_x.(tdom1), asol_v.(tdom1), seriestype=:path, lw=1.5, c=COLOR_analytic, lab=L"\textrm{Analytic (random)}")

    for k in 2:nsamples
        s = X0sampled[k]
        asol_x = analytic_solution(p, x0=s[1], v0=s[2])
        asol_v = analytic_derivative(p, x0=s[1], v0=s[2])
        plot!(fig, asol_x.(tdom1), asol_v.(tdom1), seriestype=:path, c=COLOR_analytic, lw=1.5, lab="")
    end

    return fig
end

fig = plot_distributed_discretization()
savefig(fig, "sdof_distributed_discretization.pdf")

# ### Flowpipe

model = StepIntersect(setops=:concrete)
α = 0.05

alg = VREP(δ=α*T, approx_model=model)
solvrepset = solve(probset, tspan=(0.0, tmax), alg=alg);

alg = GLGM06(δ=α*T, approx_model=model)
solglgset = solve(probset, tspan=(0.0, tmax), alg=alg);

alg = BOX(δ=α*T, approx_model=model)
solboxset = solve(probset, tspan=(0.0, tmax), alg=alg);

import Random

Random.seed!(3333)

function plot_distributed_flowpipe()
    fig = plot(legend=:top, xlab=L"u(t)", ylab=L"v(t)",
               legendfontsize=15,
               tickfont=font(20, "Times"),
               guidefontsize=20,
               xguidefont=font(30, "Times"),
               yguidefont=font(30, "Times"),
               xtick=([-1.0, -0.5, 0.0, 0.5, 1.0], [L"-1.0", L"-0.5", L"0.0", L"0.5", L"1.0"]),
               ytick=([0.0, -3.0, -6.0, -9.0, -12.0], [L"0.0", L"-3.0", L"-6.0", L"-9.0", L"-12.0"]),
               xlims=(-1.15, 1.2),
               ylims=(-15, 1.1),
               bottom_margin=6mm, left_margin=2mm, right_margin=4mm, top_margin=3mm, size=(900, 600))
    NMAX = 9

    ## propagation of Ω0 without overapproximation
    plot!(fig, solboxset[1:NMAX], vars=(1, 2), c=COLOR_Ω0_box, lab=L"\textrm{Flowpipe (hyperrectangles)}")

    ## propagation of Ω0 using zonotopes
    plot!(fig, solglgset[1:NMAX], vars=(1, 2), c=COLOR_zono, lw=1.0, alpha=.6, lab=L"\textrm{Flowpipe (zonotopes)}")

    ## propagation of Ω0 withot overapproximation
    plot!(fig, solvrepset[1:NMAX], vars=(1, 2), c=COLOR_Ω0, alpha=0.8, lw=1.0, lab=L"\textrm{Flowpipe (polygons)}")

    ## -------- Random sample of analytic solutions ------
    tdom1 = 0:0.0001:NMAX*α*T

    nsamples = 50
    X0sampled = sample(X0, nsamples);
    ω = p.ω

    s = X0sampled[1]
    asol_x = analytic_solution(p, x0=s[1], v0=s[2])
    asol_v = analytic_derivative(p, x0=s[1], v0=s[2])
    plot!(fig, asol_x.(tdom1), asol_v.(tdom1), seriestype=:path, lw=1.0, c=COLOR_analytic, lab=L"\textrm{Analytic (random)}")

    for k in 2:nsamples
        s = X0sampled[k]
        asol_x = analytic_solution(p, x0=s[1], v0=s[2])
        asol_v = analytic_derivative(p, x0=s[1], v0=s[2])
        plot!(fig, asol_x.(tdom1), asol_v.(tdom1), seriestype=:path, c=COLOR_analytic, lw=1.0, lab="")
    end

    return fig
end

fig = plot_distributed_flowpipe()
savefig(fig, "sdof_distributed_flowpipe.pdf")
