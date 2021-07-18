include("SDOF_Model.jl")

using Plots, Plots.PlotMeasures, LaTeXStrings

# ## Setup

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
