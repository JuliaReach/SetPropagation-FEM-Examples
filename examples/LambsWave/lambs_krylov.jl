using MAT
using ReachabilityAnalysis
using SparseArrays
using LinearAlgebra
using StructuralDynamicsODESolvers
using LazySets.Arrays

using Plots
using LaTeXStrings
using Plots.PlotMeasures

# utils for plots
__toL(x, digits) = L"%$(round(x, digits=digits))"
xticklatex(vec, digits) = (vec, __toL.(vec, Ref(digits)))

LazySets.set_rtol(Float64, eps(Float64))
LazySets.set_ztol(Float64, eps(Float64))
LazySets.set_atol(Float64, eps(Float64))

include("../krylov_concept.jl")

include("lambs_modelo.jl")

# cases:
#    - "constant" (Heaviside)
#    - "rickler" (step-load approximation of the rickler wavelet)
#
function sol_numerica_lamb(p::FEMProblem; case="heaviside", δt=nothing, δtfac=1.0,
                           alg=Trapezoidal, X0=nothing, kwargs...)

    if isnothing(δt)
        δt = step_size_lamb(p, δtfac)
    end

    Tmax = 1.0 # seconds
    NSTEPS = round(Int, Tmax/δt)
    tsp = range(0.0, Tmax, length=NSTEPS+1)

    m = size(p.s.K, 1)
    B = Matrix(1.0I, m, m)
    X = Universe(m)

    H = Heaviside
    if case == "constant" || case == "heaviside"
        Fval = -1e6
        Fval = [p.b * Fval for ti in tsp]

    elseif case == "ricker"
        ancho_escalon = 0.05 # 50 ms de ancho
        # por cuestion de simetria se usa 2e6 / 2 respecto del articulo de Accurate Wave Prop
        fr(t) = -1e6 * (H(ancho_escalon*3 - t) - 3*H(ancho_escalon*2 - t) + 3*H(ancho_escalon - t))
        Fval = [p.b * fr(ti) for ti in tsp] # "step-loading" en el articulo

    else
        throw(ArgumentError("case $case not defined"))
    end
    sys = SecondOrderConstrainedLinearControlContinuousSystem(sparse(p.s.M), p.s.C, p.s.K, B, Universe(m), Fval)

    if isnothing(X0)
        X0 = (zeros(m), zeros(m))
    end
    prob = @ivp(sys, u(0) ∈ X0)
    alg = alg(Δt=δt)
    sol = solve(prob, alg, NSTEPS=NSTEPS)

    # al final no invertimos
    # invertimos el signo de todos los desplazamientos para comparar con el articulo Accurate Wave prop
    #sol.U .= -1 .* sol.U
    # lo mismo con las velocidades
    #sol.U′ .= -1 .* sol.U′

    return sol
end

function sol_numerica_lamb_statistics(p::FEMProblem; case="ricker", δt=nothing, δtfac=1.0,
                                      alg=Trapezoidal, X0=nothing,
                                      sampler, idx, vecsamples,
                                      kwargs...)
    if isnothing(δt)
        δt = step_size_lamb(p, δtfac)
    end

    Tmax = 1.0 # seconds
    NSTEPS = round(Int, Tmax/δt)
    tsp = range(0.0, Tmax, length=NSTEPS+1)

    m = size(p.s.K, 1)
    B = Matrix(1.0I, m, m)
    X = Universe(m)

    H = Heaviside
    if case == "constant" || case == "heaviside"
        Fval = -1e6
        Fval = [p.b * Fval for ti in tsp]

    elseif case == "ricker"
        ancho_escalon = 0.05 # 50 ms de ancho
        # por cuestion de simetria se usa 2e6 / 2 respecto del articulo de Accurate Wave Prop
        fr(t) = -1e6 * (H(ancho_escalon*3 - t) - 3*H(ancho_escalon*2 - t) + 3*H(ancho_escalon - t))
        Fval = [p.F * fr(ti) for ti in tsp] # "step-loading" en el articulo

    else
        throw(ArgumentError("case $case not defined"))
    end
    sys = SecondOrderConstrainedLinearControlContinuousSystem(sparse(p.s.M), p.s.C, p.s.K, B, Universe(m), Fval)

    if isnothing(X0)
        X0 = (zeros(m), zeros(m))
    end
    prob = @ivp(sys, u(0) ∈ X0)
    alg = alg(Δt=δt)

    sol_statistics = StructuralDynamicsODESolvers._solve_statistics(alg, prob, NSTEPS, sampler, idx, vecsamples)

    # invertimos el signo de todos los desplazamientos para comparar con el articulo Accurate Wave prop
    #sol.U .= -1 .* sol.U

    return sol_statistics
end



function plot_nodo(sol; nodeidx=3, direction="Horizontal", method="Newmark")

    fig = plot(xlab=L"\textrm{Time (s)}",
               ylab=L"\textrm{%$direction displacement} (\times 10^{-5} \textrm{m})",
               legend=:topright,
               legendfontsize=25,
               tickfont=font(20, "Times"),
               guidefontsize=10,
               xguidefont=font(20, "Times"),
               yguidefont=font(20, "Times"),
               xtick = xticklatex([0, 0.2, 0.4, 0.6, 0.8, 1.0], 2),
               ytick = ([-4e-5, -2e-5, 0, 2e-5, 4e-5, 6e-5], [L"-4", L"-2", L"0", L"2", L"4", L"6"]),
               xlims=(0.0, 1.0),
               ylims=(-5e-5, 8e-5),
               bottom_margin=10mm, left_margin=6mm, right_margin=12mm, top_margin=3mm, size=(900, 600))

    plot!(fig, sol, vars=(0, nodeidx), seriestype=:path, markershape=:none, lw=2.0, lab=L"\textrm{%$method}")
    fig
end

function plot_lambs(p::FEMProblem; nodeidx=3, direction="Horizontal", case="heaviside", kwargs...)

    fig = plot(xlab=L"\textrm{Time (s)}",
               ylab=L"\textrm{%$direction displacement} (\times 10^{-5} \textrm{m})",
               legend=:topright,
               legendfontsize=25,
               tickfont=font(20, "Times"),
               guidefontsize=10,
               xguidefont=font(20, "Times"),
               yguidefont=font(20, "Times"),
               xtick = xticklatex([0, 0.2, 0.4, 0.6, 0.8, 1.0], 2),
               ytick = ([-4e-5, -2e-5, 0, 2e-5, 4e-5, 6e-5], [L"-4", L"-2", L"0", L"2", L"4", L"6"]),
               xlims=(0.0, 1.0),
               ylims=(-5e-5, 8e-5),
               bottom_margin=10mm, left_margin=6mm, right_margin=12mm, top_margin=3mm, size=(900, 600))

    # resolver numericas
    sol_newmark = sol_numerica_lamb(p; nodeidx=nodeidx, case=case, alg=Trapezoidal, kwargs...)
    sol_bathe = sol_numerica_lamb(p; nodeidx=nodeidx, case=case, alg=Bathe, kwargs...)

    # resolver flowpipe
    out = _sol_reachability_lamb(p; case=case, nodeidx=nodeidx, kwargs...)

    plot!(fig, out, vars=(0, nodeidx), lc=:blue, c=:blue, lab=L"\textrm{Flowpipe}")
    plot!(fig, sol_newmark, vars=(0, nodeidx), seriestype=:path, markershape=:none, lw=2.0, lab=L"\textrm{Newmark}")
    plot!(fig, sol_bathe, vars=(0, nodeidx), seriestype=:path, markershape=:none, lw=2.0, lab=L"\textrm{Bathe}")
    fig
end

function generar_figura_malla2_singleton_position()

    # =================
    # setear parametros
    # =================
    p = load_lambs(mesh=2);
    @show m = size(p.s.M, 1)
    nodeidx = 3
    δtfac = 1.0
    δt = step_size_lamb(p, δtfac) # 0.001
    case = "heaviside"
    condicion_inicial = "singleton"
    mk = 140
    tol = 1e-7

    # ====================
    # resolver numericas
    # ====================
    println("Resolviendo con Newmark... ")
    @time sol_newmark = sol_numerica_lamb(p; nodeidx=nodeidx, case=case, alg=Trapezoidal, δt=δt, condicion_inicial=condicion_inicial)

    println("Resolviendo con Bathe... ")
    @time sol_bathe = sol_numerica_lamb(p; nodeidx=nodeidx, case=case, alg=Bathe, δt=δt, condicion_inicial=condicion_inicial)

    # =================
    # resolver flowpipe
    # =================
    println("Resolviendo Flowpipe... ")
    @time out = _sol_reachability_lamb(p; case=case, nodeidx=nodeidx, δt=δt, condicion_inicial=condicion_inicial, mk=mk, tol=tol)

    # =================
    # graficar
    # =================
    fig = plot(xlab=L"\textrm{Time (s)}",
               ylab=L"\textrm{Displacement} (\times 10^{-5} \textrm{m})",
               legend=:bottomright,
               legendfontsize=15,
               tickfont=font(25, "Times"),
               guidefontsize=15,
               xguidefont=font(25, "Times"),
               yguidefont=font(25, "Times"),
               xtick = xticklatex([0, 0.2, 0.4, 0.6, 0.8, 1.0], 2),
               ytick = ([-1e-5, 0, 1e-5, 2e-5, 3e-5, 4e-5], [L"-1", L"0", L"1", L"2", L"3", L"4"]),
               xlims=(0.0, 1.0),
               ylims=(-1.2e-5, 4.2e-5),
               bottom_margin=10mm, left_margin=6mm, right_margin=12mm, top_margin=3mm, size=(900, 600))

    plot!(fig, out, vars=(0, nodeidx), c=:lightblue, lc=:black, lw=0.5, alpha=1., lab=L"\textrm{Flowpipe}")

    # grafico continuo
    plot!(fig, sol_newmark, vars=(0, nodeidx), seriestype=:path, marker=:none, lw=2.0, color=:red, lab="")
    plot!(fig, sol_bathe, vars=(0, nodeidx), seriestype=:path, marker=:none, lw=2.0, color=:green, lab="")

    # con markers salteados
    paso = 15
    npoints = length(sol_newmark.U)
    plot!(fig, sol_newmark.t[1:paso:npoints], [sol_newmark.U[i][nodeidx] for i in 1:paso:npoints], vars=(0, nodeidx), linetype=:scatter, marker=:circle, markersize=8, lw=2.0, color=:red, lab=L"\textrm{Newmark}")
    plot!(fig, sol_bathe.t[1:paso:npoints], [sol_bathe.U[i][nodeidx] for i in 1:paso:npoints], vars=(0, nodeidx), linetype=:scatter, marker=:utriangle, markersize=8, lw=2.0, color=:green, lab=L"\textrm{Bathe}")

    # grafico bounding box donde hago zoom
    tzoom_left = 0.40
    tzoom_right = 0.50
    yzoom_bottom = 2e-5
    yzoom_top = 4e-5

    PBotLeft  = [tzoom_left, yzoom_bottom]
    PBotRight = [tzoom_right, yzoom_bottom]
    PTopLeft  = [tzoom_left, yzoom_top]
    PTopRight = [tzoom_right, yzoom_top]

    L1 = LineSegment(PBotLeft, PBotRight)
    L2 = LineSegment(PBotRight, PTopRight)
    L3 = LineSegment(PTopRight, PTopLeft)
    L4 = LineSegment(PTopLeft, PBotLeft)

    plot!(fig, L1, color=:black, style=:dash, lw=3.0, seriestype=:solid, marker=:none, alpha=1.)
    plot!(fig, L2, color=:black, style=:dash, lw=3.0, seriestype=:solid, marker=:none, alpha=1.)
    plot!(fig, L3, color=:black, style=:dash, lw=3.0, seriestype=:solid, marker=:none, alpha=1.)
    plot!(fig, L4, color=:black, style=:dash, lw=3.0, seriestype=:solid, marker=:none, alpha=1.)

    savefig(fig, "malla2_singleton_pos_nodeidx_3.pdf")

    # =============================
    # graficar parte zoomeada
    # =============================

    # como todos los solvers usan el mismo step size, puedo usar Newmark para encontrar
    # los indices donde hago zoom
    a = findfirst(>(0), sol_newmark.t .- tzoom_left)
    b = findfirst(>(0), sol_newmark.t .- tzoom_right)

    fig = plot(xlab=L"\textrm{Time (s)}",
               ylab=L"\textrm{Displacement} (\times 10^{-5} \textrm{m})",
               legend=:bottomright,
               legendfontsize=15,
               tickfont=font(25, "Times"),
               guidefontsize=15,
               xguidefont=font(25, "Times"),
               yguidefont=font(25, "Times"),
               xtick = xticklatex([0.4, 0.425, 0.45, 0.475, 0.5], 3),
               ytick = ([2.5e-5, 2.75e-5, 3.0e-5, 3.25e-5, 3.5e-5], [L"2.5", L"2.75", L"3.0", L"3.25", L"3.5"]),
               xlims=(0.395, 0.505),
               ylims=(2.4e-5, 3.7e-5),
               bottom_margin=10mm, left_margin=6mm, right_margin=12mm, top_margin=3mm, size=(900, 600))

    plot!(fig, out[a:b], vars=(0, nodeidx), c=:lightblue, lc=:black, lw=0.5, alpha=1., lab=L"\textrm{Flowpipe}")

    plot!(fig, sol_newmark.t[a:b], [sol_newmark.U[i][nodeidx] for i in a:b], vars=(0, nodeidx), seriestype=:path, markershape=:none, lw=2.0, color=:red, lab="", alpha=1.)
    plot!(fig, sol_bathe.t[a:b], [sol_bathe.U[i][nodeidx] for i in a:b], seriestype=:path, markershape=:none, lw=2.0, color=:green, lab="", alpha=1.)

    # con markers salteados
    paso = 2
    plot!(fig, sol_newmark.t[a:paso:b], [sol_newmark.U[i][nodeidx] for i in a:paso:b], vars=(0, nodeidx), linetype=:scatter, marker=:circle, markersize=10, lw=2.0, color=:red, lab=L"\textrm{Newmark}")
    plot!(fig, sol_bathe.t[a:paso:b], [sol_bathe.U[i][nodeidx] for i in a:paso:b], vars=(0, nodeidx), linetype=:scatter, marker=:utriangle, markersize=10, lw=2.0, color=:green, lab=L"\textrm{Bathe}")


    savefig(fig, "malla2_singleton_pos_nodeidx_3_zoom.pdf")

    return nothing
end

function generar_figura_malla2_singleton_velocity()

    # =================
    # setear parametros
    # =================
    p = load_lambs(mesh=2);
    @show m = size(p.s.M, 1)
    nodeidx = 3
    δtfac = 1.0
    δt = step_size_lamb(p, δtfac) # 0.001
    case = "heaviside"
    condicion_inicial = "singleton"
    mk = 140
    tol = 1e-7

    # ====================
    # resolver numericas
    # ====================
    println("Resolviendo con Newmark... ")
    @time sol_newmark = sol_numerica_lamb(p; nodeidx=nodeidx, case=case, alg=Trapezoidal, δt=δt, condicion_inicial=condicion_inicial)

    println("Resolviendo con Bathe... ")
    @time sol_bathe = sol_numerica_lamb(p; nodeidx=nodeidx, case=case, alg=Bathe, δt=δt, condicion_inicial=condicion_inicial)

    # =================
    # resolver flowpipe
    # =================
    println("Resolviendo Flowpipe... ")
    @time out = _sol_reachability_lamb(p; case=case, nodeidx=nodeidx+m, δt=δt, condicion_inicial=condicion_inicial, mk=mk, tol=tol)

    # =================
    # graficar
    # =================
    fig = plot(xlab=L"\textrm{Time (s)}",
               ylab=L"\textrm{Velocity} (\times 10^{-4} \textrm{m/s})",
               legend=:topright,
               legendfontsize=15,
               tickfont=font(25, "Times"),
               guidefontsize=15,
               xguidefont=font(25, "Times"),
               yguidefont=font(25, "Times"),
               xtick = xticklatex([0, 0.2, 0.4, 0.6, 0.8, 1.0], 2),
               ytick = ([-1e-4, 0, 1e-4, 2e-4, 3e-4, 4e-4, 5e-4], [L"-1", L"0", L"1", L"2", L"3", L"4", L"5"]),
               xlims=(-0.02, 1.0),
               ylims=(-2.1e-4, 5.5e-4),
               bottom_margin=10mm, left_margin=6mm, right_margin=12mm, top_margin=3mm, size=(900, 600))

    plot!(fig, out, vars=(0, nodeidx+m), c=:lightblue, lc=:black, lab=L"\textrm{Flowpipe}", alpha=1.)

    # grafico continuo
    plot!(fig, sol_newmark.t, [sol_newmark.U′[i][nodeidx] for i in 1:length(sol_newmark.U′)], vars=(0, nodeidx), seriestype=:path, markershape=:none, lw=2.0, lab="")
    plot!(fig, sol_bathe.t, [sol_bathe.U′[i][nodeidx] for i in 1:length(sol_bathe.U′)], seriestype=:path, markershape=:none, lw=2.0, lab="")

    # con markers salteados
    paso = 15
    npoints = length(sol_newmark.U)
    plot!(fig, sol_newmark.t[1:paso:npoints], [sol_newmark.U′[i][nodeidx] for i in 1:paso:npoints], vars=(0, nodeidx), linetype=:scatter, marker=:circle, markersize=10, lw=2.0, c=:red, lab=L"\textrm{Newmark}")
    plot!(fig, sol_bathe.t[1:paso:npoints], [sol_bathe.U′[i][nodeidx] for i in 1:paso:npoints], linetype=:scatter, marker=:utriangle, markersize=10, lw=2.0, c=:green, lab=L"\textrm{Bathe}")

    # grafico bounding box donde hago zoom
    tzoom_left = 0.8
    tzoom_right = 0.99
    yzoom_bottom = -1e-4
    yzoom_top = 1e-4

    PBotLeft  = [tzoom_left, yzoom_bottom]
    PBotRight = [tzoom_right, yzoom_bottom]
    PTopLeft  = [tzoom_left, yzoom_top]
    PTopRight = [tzoom_right, yzoom_top]

    L1 = LineSegment(PBotLeft, PBotRight)
    L2 = LineSegment(PBotRight, PTopRight)
    L3 = LineSegment(PTopRight, PTopLeft)
    L4 = LineSegment(PTopLeft, PBotLeft)

    plot!(fig, L1, color=:black, style=:dash, lw=3.0, seriestype=:solid, marker=:none, alpha=1.)
    plot!(fig, L2, color=:black, style=:dash, lw=3.0, seriestype=:solid, marker=:none, alpha=1.)
    plot!(fig, L3, color=:black, style=:dash, lw=3.0, seriestype=:solid, marker=:none, alpha=1.)
    plot!(fig, L4, color=:black, style=:dash, lw=3.0, seriestype=:solid, marker=:none, alpha=1.)

    savefig(fig, "malla2_singleton_vel_nodeidx_3.pdf")

    # =============================
    # graficar parte zoomeada
    # =============================

    # como todos los solvers usan el mismo step size, puedo usar Newmark para encontrar
    # los indices donde hago zoom
    a = findfirst(>(0), sol_newmark.t .- tzoom_left)
    b = findfirst(>(0), sol_newmark.t .- tzoom_right)

    fig = plot(xlab=L"\textrm{Time (s)}",
               ylab=L"\textrm{Velocity} (\times 10^{-4} \textrm{m/s})",
               legend=:topright,
               legendfontsize=15,
               tickfont=font(25, "Times"),
               guidefontsize=15,
               xguidefont=font(25, "Times"),
               yguidefont=font(25, "Times"),
               xtick = xticklatex([0.8, 0.9, 1.0], 3),
               ytick = ([-1.0e-4, -0.5e-4, 0, 0.5e-4, 1e-4], [L"-1.0", L"-0.5", L"0.0", L"0.5", L"1.0"]),
               xlims=(0.79, 1.01),
               ylims=(-1.05e-4, 1.05e-4),
               bottom_margin=10mm, left_margin=6mm, right_margin=12mm, top_margin=3mm, size=(900, 600))

    plot!(fig, out[a:b-1], vars=(0, nodeidx+m), lc=:black, c=:lightblue, lw=0.5, alpha=1., lab=L"\textrm{Flowpipe}")
    plot!(fig, sol_bathe.t[a:b], [sol_bathe.U′[i][nodeidx] for i in a:b], seriestype=:path, markershape=:none, c=:red, lw=2.0, lab="", alpha=1.)
    plot!(fig, sol_newmark.t[a:b], [sol_newmark.U′[i][nodeidx] for i in a:b], vars=(0, nodeidx), seriestype=:path, markershape=:none, c=:green, lw=2.0, lab="", alpha=1.)

    # con markers salteados
    paso = 2
    plot!(fig, sol_newmark.t[a:paso:b], [sol_newmark.U′[i][nodeidx] for i in a:paso:b], vars=(0, nodeidx), linetype=:scatter, marker=:circle, markersize=10, lw=2.0, color=:red, alpha=1., lab=L"\textrm{Newmark}")
    plot!(fig, sol_bathe.t[a:paso:b], [sol_bathe.U′[i][nodeidx] for i in a:paso:b], vars=(0, nodeidx), linetype=:scatter, marker=:utriangle, markersize=10, lw=2.0, color=:green, alpha=1., lab=L"\textrm{Bathe}")

    savefig(fig, "malla2_singleton_vel_nodeidx_3_zoom.pdf")

    return nothing
end

# ========================================
# Metodo de reachability
# ========================================

function _sol_reachability_lamb(p::FEMProblem; δt=nothing, nodeidx=[3, 4],
                                case="ricker", δtfac=1.0, alg=Trapezoidal,
                                mk=150, tol=1e-12,
                                Eplus_optim=true,
                                kwargs...)

    if case == "heaviside"
        return _sol_reachability_lamb_heaviside(p; δt=δt, nodeidx=nodeidx, δtfac=δtfac, alg=alg, mk=mk, tol=tol, Eplus_optim=Eplus_optim, kwargs...)
    else
        error("ver lambs_krylov_old.jl")
    end
end

function _sol_reachability_lamb_heaviside(p::FEMProblem;
                                          δtfac=1.0,
                                          δt=step_size_lamb(p, δtfac),
                                          nodeidx=[3, 4],
                                          alg=Trapezoidal,
                                          mk=150, tol=1e-12,
                                          Eplus_optim=true,
                                          condicion_inicial="singleton",

                                          # para condicion inicial "distribuida"
                                          radio_desplazamientos = nothing, # 3.5e-7
                                          radio_velocidades = nothing, # 5e-6

                                          kwargs...)

    # ===================
    # Definir la carga
    # ===================

    # @assert δt <= δtnom

    Tmax = 1.0 # seconds
    NSTEPS = round(Int, Tmax/δt)

    H = Heaviside

    ancho_escalon = 0.05 # 50 ms de ancho
    # por cuestion de simetria se usa 2e6 / 2 respecto del articulo de Accurate Wave Prop
    Δ = ancho_escalon
    α = -1e6
    fload(t) = α * H(t)

    # ==============================
    # Homogeneizar
    # ==============================
    m = size(p.s.K, 1)
    Shom = _homogeneize(p)
    Aext = Shom.A
    @assert size(Aext, 1) == 2m+1

    # =====================================
    # Preparacion de la condicion inicial
    # =====================================

    if condicion_inicial == "singleton"
        X0_init = Singleton(zeros(2m+1))
        X0_init.element[end] = α

    elseif condicion_inicial == "distribuida"
        X0_init = Hyperrectangle(zeros(2m+1), zeros(2m+1))
        if isnothing(radio_desplazamientos) || isnothing(radio_velocidades)
            error("debe especificar `radio_desplazamientos` y `radio_velocidades`")
        end
        X0_init.center[end] = α # fuerza
        X0_init.radius[1:m] .= radio_desplazamientos # desplazamientos iniciales
        X0_init.radius[m+1:2m] .= radio_velocidades # velocidades iniciales

    else
        error("condicion_inicial = $condicion_inicial no existe")
    end

    # ==============================
    # Propagacion a tiempo continuo
    # ==============================
    Aᵀδ = copy(transpose(Aext .* δt))
    num_nodes = isa(nodeidx, Int) ? 1 : length(nodeidx)

    # alg = LGG09(δ=δt, n=2m+1, vars=nodeidx)
    dirs = Vector{LazySets.Arrays.SingleEntryVector{Float64}}()
    for j in 1:num_nodes
        node = nodeidx[j]
        push!(dirs, SingleEntryVector(node, 2m+1, 1.0))  # 0, 0, 1, ... , 0, 0, 0, 0, 0, 0
        push!(dirs, SingleEntryVector(node, 2m+1, -1.0))
    end
    cdirs = CustomDirections(dirs)

    println("Calculando el factor de bloating...")

    Eplus = _build_Eplus_optim(Aext, X0_init, δt; m=mk, tol=tol, cutoff=eps(Float64))

    println("Propagando a tiempo continuo...")

    if condicion_inicial == "singleton"
        ΦX0_init = _expv(Aext, X0_init.element, 1, δt, m=mk, tol=tol) |> Singleton
        Ω0 = CH(X0_init, ΦX0_init ⊕ Eplus)

        Mout = Matrix{Float64}(undef, 2*num_nodes, NSTEPS)
        for j in 1:num_nodes
            node = nodeidx[j]
            d = zeros(2m+1)
            d[node] = 1.0
            println("Propagando indice $node ...")
            idx = (2(j-1) + 1):2j
            _reach_homog_krylov_LGG09_modif!(view(Mout, idx, :), Ω0, Aᵀδ, d, NSTEPS; m=mk, tol=tol)
        end
        F = [TemplateReachSet(cdirs, view(Mout, :, k), (0.0 .. δt) + (k-1) * δt) for k in 1:NSTEPS]

    elseif condicion_inicial == "distribuida"
        # Ω0 = CH(X0_init, ΦX0_init * X0_init ⊕ Eplus)

        Mout = Matrix{Float64}(undef, 2*num_nodes, NSTEPS)
        for j in 1:num_nodes
            node = nodeidx[j]

            d = zeros(2m+1)
            d[node] = 1.0

            println("Propagando indice $node ...")
            idx = (2(j-1) + 1):2j
            _reach_homog_krylov_LGG09_modif!(view(Mout, idx, :), X0_init, Aext, δt, Eplus, d, NSTEPS; m=mk, tol=tol)
        end
        F = [TemplateReachSet(cdirs, view(Mout, :, k), (0.0 .. δt) + (k-1) * δt) for k in 1:NSTEPS]

    else
        error("condicion_inicial = $condicion_inicial no existe")
    end

    # niego para ajustar el sistema de coordenadas con el articulo
    proyectar = get(kwargs, :proyectar, false)
    if proyectar
        # OJO si uso esto, el flowpipe tiene una sola variable 1
        @assert length(nodeidx) == 1
        out = [ReachSet(Interval(-Mout[2, k], Mout[1, k]), (0.0 .. δt) + (k-1) * δt) for k in 1:size(Mout, 2)]
    else
        #α = -1.0
        α = 1.0
        out = [ReachSet(α * set(R), tspan(R)) for R in F]
    end

    return out
end

# computes the support functions of Ω0 = CH(X0, Φ*X0 ⊕ Eplus)
# along +d and -d for NSTEPS
#
# for k = 0, 1, ..., NSTEPS-1
#
# ρ(d, Φ^k Ω0) = ρ((Φ^T)^k d, CH(X0, Φ * X0 ⊕ Eplus))
#              = max(ρ((Φ^T)^k d, X0), ρ((Φ^T)^k d, Φ * X0) + ρ((Φ^T)^k d, Eplus))
#              = max(ρ((Φ^T)^k d, X0), ρ((Φ^T)^(k+1) d, X0) + ρ((Φ^T)^k d, Eplus))
#              = max(ρ(vk, X0), ρ((Φ^T)^k * vk, X0) + ρ(vk, Eplus)), vk := (Φ^T)^k d
#
function _reach_homog_krylov_LGG09_modif!(out, X0, A, δ, Eplus,
                                          d::AbstractVector, NSTEPS;
                                          hermitian=false, m=min(30, size(Aᵀδ, 1)), tol=1e-7)

    @assert size(out, 1) == 2
    @assert size(out, 2) == NSTEPS

    Aᵀδ = transpose(A) .* δ

    # initialization of the krylov subspace
    # esta parte genera una base de
    # {Aᵀδ * d, (Aᵀδ)^2 * d, (Aᵀδ)^3 * d, ..., (Aᵀδ)^m * d}
    TA, Tb = eltype(Aᵀδ), eltype(d)
    T = promote_type(TA, Tb)
    Ks = KrylovSubspace{T, real(T)}(length(d), m)
    arnoldi!(Ks, Aᵀδ, d; m=m, ishermitian=hermitian, tol=tol)

    # rᵢ stores is the cache for each vector: (Φᵀ)^i d and (Φᵀ)^{i+1} d
    rᵢ = deepcopy(d)
    rᵢ₊₁ = similar(d)

    #
    # TODO se puede optimizar el loop ya que el ρ(rᵢ₊₁, X0) se repite en ρ(rᵢ, X0)
    #
    @inbounds for i in 1:NSTEPS
        expv!(rᵢ₊₁, i*1.0, Ks) # rᵢ₊₁ <- exp(A^T δ * i) * d
        out[1, i] = max(ρ(rᵢ, X0), ρ(rᵢ₊₁, X0) + ρ(rᵢ, Eplus))
        out[2, i] = max(ρ(-rᵢ, X0), ρ(-rᵢ₊₁, X0) + ρ(-rᵢ, Eplus))
        copyto!(rᵢ, rᵢ₊₁)
    end
    return out
end


# ===========================
# estadisticas con numericas
# ===========================

using Sobol
using StructuralDynamicsODESolvers: SolutionExtrema

sol_numericas_random_statistics(p; m, vecsamples::Int, idx, radio_desplazamientos, radio_velocidades) = sol_numericas_random_statistics(p; m=m,
    vecsamples=[vecsamples], idx=idx, radio_desplazamientos=radio_desplazamientos, radio_velocidades=radio_velocidades)

# ejemplo: vecsamples = [10, 20, 100]
function sol_numericas_random_statistics(p; m, vecsamples::Vector{Int}, idx, δt, radio_desplazamientos, radio_velocidades)
    #X0rand_desplazamientos = [[2(rand() - 0.5) * 5e-7 for _ in 1:825] for i in 1:NTRAJ]
    #X0rand_velocidades = [[2(rand() - 0.5) * 5e-7 for _ in 1:825] for i in 1:NTRAJ]
    #remake = [(X0d, X0v) for (X0d, X0v) in zip(X0rand_desplazamientos, X0rand_velocidades)]

    H0 = Hyperrectangle(zeros(2m), vcat(fill(radio_desplazamientos, m), fill(radio_velocidades, m)))

    sampleado = :vertices

    if sampleado == :sobol
        # para samplear dentro de H0
        sampler = SobolSeq(H0)
    elseif sampleado == :vertices
        # para samplear en los vertices
        sampler = VertexSampler(H0)
    end

    @time out = sol_numerica_lamb_statistics(p; case="heaviside", alg=Trapezoidal, δt=δt, X0=nothing, sampler=sampler, idx=idx, vecsamples=vecsamples)
end

Sobol.SobolSeq(H::Hyperrectangle) = SobolSeq(low(H), high(H))

struct VertexSampler{T}
    H::T
end

function StructuralDynamicsODESolvers._next!(S::VertexSampler, unew)
    c = S.H.center
    r = S.H.radius
    for i in 1:length(unew)
        a = ifelse(rand() > 0.5, 1, -1)
        unew[i] = c[i] + a*r[i]
    end
end

StructuralDynamicsODESolvers._next!(s::ScaledSobolSeq, unew) = next!(s, unew)

function generar_figura_malla2_envolventes_position(; vecsamples = nothing)

    # =================
    # setear parametros
    # =================
    p = load_lambs(mesh=2);
    @show m = size(p.s.M, 1)
    nodeidx = 3
    δt = step_size_lamb(p, 1.0) # 0.001
    case = "heaviside"
    condicion_inicial = "distribuida"
    mk = 140
    tol = 1e-7

    radio_desplazamientos = 3.5e-7
    radio_velocidades = 5e-6

    # ====================
    # correr simulaciones
    # ====================
    out_extrema_malla2 = nothing
    if !isnothing(vecsamples)
        out_extrema_malla2 = sol_numericas_random_statistics(p, m=m, idx=nodeidx, vecsamples=vecsamples, δt=δt, radio_desplazamientos=radio_desplazamientos, radio_velocidades=radio_velocidades)
    end

    # ====================
    # calcular flowpipe
    # ====================
    @time out_flowpipe = _sol_reachability_lamb(p, case=case, condicion_inicial=condicion_inicial, nodeidx=nodeidx, mk=mk, tol=tol, δt=δt, proyectar=true, radio_desplazamientos=radio_desplazamientos, radio_velocidades=radio_velocidades)

    # ====================
    # graficar flowpipe
    # ====================
    # q = length(vecsamples)
    # c = distinguishable_colors(q);
    fig = plot()

    plot!(fig, out_flowpipe, vars=(0, 1), lw=0.0, lc=:lightblue, c=:lightblue, lab=L"\textrm{Flowpipe}")

    if !isnothing(vecsamples)
        [plot!(fig, out_extrema_malla2.t, out_extrema_malla2.hist[k][:U_min], lw=2.0, lab="$k") for (i, k) in enumerate(vecsamples)]
        [plot!(fig, out_extrema_malla2.t, out_extrema_malla2.hist[k][:U_max], lw=2.0, lab="") for (i, k) in enumerate(vecsamples)]
    end
    savefig(fig, "malla2_distribuida_pos_nodeidx_3.pdf")

    out_flowpipe, out_extrema_malla2, fig
end

function generar_figura_malla2_envolventes_velocity(; vecsamples = nothing)

    # =================
    # setear parametros
    # =================
    p = load_lambs(mesh=2);
    @show m = size(p.s.M, 1)
    nodeidx = 3
    δt = step_size_lamb(p, 1.0)
    case = "heaviside"
    condicion_inicial = "distribuida"
    mk = 140
    tol = 1e-7

    radio_desplazamientos = 3.5e-7
    radio_velocidades = 5e-6

    # ====================
    # correr simulaciones
    # ====================
    println("Corriendo simulaciones...")
    out_extrema_malla2 = nothing
    if !isnothing(vecsamples)
        out_extrema_malla2 = sol_numericas_random_statistics(p, m=m, idx=nodeidx, vecsamples=vecsamples, δt=δt, radio_desplazamientos=radio_desplazamientos, radio_velocidades=radio_velocidades)
    end

    # ====================
    # calcular flowpipe
    # ====================
    println("Calculando el flowpipe...")
    @time out_flowpipe = _sol_reachability_lamb(p, case=case, condicion_inicial=condicion_inicial, nodeidx=nodeidx+m, mk=mk, tol=tol, δt=δt, proyectar=true, radio_desplazamientos=radio_desplazamientos, radio_velocidades=radio_velocidades)

    # ====================
    # graficar flowpipe
    # ====================
    # q = length(vecsamples)
    # c = distinguishable_colors(q);
    fig = plot()

    plot!(fig, out_flowpipe, vars=(0, 1), lw=0.0, lc=:lightblue, c=:lightblue, lab=L"\textrm{Flowpipe}")
    [plot!(fig, out_extrema_malla2.t, out_extrema_malla2.hist[k][:Udot_min], lw=2.0, lab="$k") for (i, k) in enumerate(vecsamples)]
    [plot!(fig, out_extrema_malla2.t, out_extrema_malla2.hist[k][:Udot_max], lw=2.0, lab="") for (i, k) in enumerate(vecsamples)]

    savefig(fig, "malla2_distribuida_vel_nodeidx_3.pdf")

    out_flowpipe, out_extrema_malla2, fig
end

function _norma_inf_envolvente_numericas(sol, vecsamples)
    sol_norma_inf = Vector{Float64}()
    for k in vecsamples
        Udotmin_k = sol.hist[k][:Udot_min] .|> abs
        Udotmax_k = sol.hist[k][:Udot_max] .|> abs

        # norma infinito de las envolventes numericas
        sol_norma_inf_k = maximum(max.(Udotmin_k, Udotmax_k))
        push!(sol_norma_inf, sol_norma_inf_k)
    end
    return sol_norma_inf
end

function _norma_1_envolvente_numericas(sol, vecsamples)
    sol_norma_1 = Vector{Float64}()
    for k in vecsamples
        Udotmin_k = sol.hist[k][:Udot_min] .|> abs
        Udotmax_k = sol.hist[k][:Udot_max] .|> abs

        # norma infinito de las envolventes numericas
        sol_norma_1_k = sum(max(Udotmin_k[i], Udotmax_k[i]) for i in 1:length(Udotmin_k))
        push!(sol_norma_1, sol_norma_1_k)
    end
    δt = out_extrema_malla2_vel.alg.Δt
    return sol_norma_1 * δt
end

function generar_figura_malla2_envolvente_velocidad()

    results = matread("results.mat")

    vecsamples = "N" .* string.(results["vecsamples"])

    times = results["out_extrema_malla2_vel_t"]
    vmin = [results["out_extrema_malla2_vel_hist"][k]["NUdot_min"] for k in vecsamples]
    vmax = [results["out_extrema_malla2_vel_hist"][k]["NUdot_max"] for k in vecsamples];

    xx = [Interval(a, b) for (a, b) in zip(results["out_flowpipe_vel_lo"], results["out_flowpipe_vel_hi"])];
    tt = [IntervalArithmetic.Interval(a, b) for (a, b) in zip(results["out_flowpipe_vel_tstart"], results["out_flowpipe_vel_tend"])];
    F = [ReachSet(a, b) for (a, b) in zip(xx, tt)];

    cols = distinguishable_colors(length(vecsamples))

    fig = plot(xlab=L"\textrm{Time (s)}",
               ylab=L"\textrm{Velocity} (\times 10^{-4} \textrm{m/s})",
               legend=:outerright,
               legendfontsize=20,
               tickfont=font(20, "Times"),
               guidefontsize=10,
               xguidefont=font(20, "Times"),
               yguidefont=font(20, "Times"),
               xtick = xticklatex([0, 0.2, 0.4, 0.6, 0.8, 1.0], 2),
               ytick = ([-10e-4, -5e-4, 0, 5e-4, 10e-4], [L"-10", L"-5", L"0", L"5", L"10"]),
               xlims=(0.0, 1.0),
               ylims=(-1.3e-3, 1.2e-3),
               bottom_margin=10mm, left_margin=6mm, right_margin=12mm, top_margin=3mm, size=(900, 600))

    labs = [L"10^0", L"10^1", L"10^2", L"10^3", L"10^4", L"10^5"]

    plot!(fig, F, vars=(0, 1), lw=1., lc=:lightblue, c=:lightblue, lab=L"\textrm{Flowpipe}")
    [plot!(fig, times, vmin[i], lw=2.0, c=cols[i], lab=labs[i], ls=:solid) for i in 1:length(vmin)]
    [plot!(fig, times, vmax[i], lw=2.0, c=cols[i], lab="", ls=:solid) for i in 1:length(vmin)]

    savefig(fig, "fig_vel.pdf")

    fig
end
