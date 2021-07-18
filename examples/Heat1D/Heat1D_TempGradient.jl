using Plots, LaTeXStrings, Plots.Measures
using LinearAlgebra, MAT

# load the model
include("Heat1D_Model.jl")

import IntervalArithmetic
const IA = IntervalArithmetic

file = joinpath(@__DIR__, "heat1d_99.mat")

δ = 1e-5
k = 1_000
tol = 0.1

@time prob, alg, NSTEPS = _solve_heat1d(NSTEPS=1, δ=δ , file=file, model=StepIntersect(Forward()))
n = statedim(prob);

LGG_solu_dirs = LGG_dirs_solucion(n)
LGG_grad_dirs = LGG_dirs_gradiente(n)

L = 1.0
Δx = L/(n+2)

# warm up
solve(prob, alg=LGG09(δ=δ, template=LGG_grad_dirs, approx_model=Forward(setops=:box)), NSTEPS=1)
varias_sol(file=file, δ=δ, tol=tol, NSTEPS=1, NCURVAS=1);

@time grad_lgg = solve(prob, alg=LGG09(δ=δ, template=LGG_grad_dirs, approx_model=Forward(setops=:box)), NSTEPS=k);

NCURVAS = 1_000
@time Temps, gradsT = varias_sol(file=file, δ=δ, tol=tol, NSTEPS=k, NCURVAS=NCURVAS);

# spatial gradient
function plot_gradiente_espacial(; paso = 1_000, δ = 1E-5, NCURVAS=length(gradsT), extremos_y = 6.0)

    fig = plot(xlab=L"\textrm{Position}",
               ylab=L"\textrm{Temperature gradient}",
               legendfontsize=18, legend=:top,
               tickfont=font(25, "Times"), guidefontsize=20,
               xguidefont=font(30, "Times"),
               yguidefont=font(30, "Times"),
               ylims=(-extremos_y-0.5, extremos_y),
               xlims=(0.0, 1.0),
               xtick=([0.0, 0.25, 0.50, 0.75, 1.0], [L"0.0", L"0.25", L"0.50", L"0.75", L"1.0"]),
               ytick=([-10.0, -5.0, 0.0, 5.0, 10.0], [L"-10", L"-5", L"0", L"5", L"10"]),
               bottom_margin=6mm, left_margin=6mm, right_margin=4mm, top_margin=6mm, size=(900, 600))

    fp_lgg_grad = gradiente_LGG_graddirs(grad_lgg = grad_lgg, paso_temporal = paso)
    plot!(fig, fp_lgg_grad, vars=(0,1), lw=.5, c=:lightblue)

    xsnumer = range(0, L, length=length(gradsT[1][1])+1)
    xsnumer = 0.5*(xsnumer[1:end-1]+xsnumer[2:end])
    for j = 1:NCURVAS
        plot!(fig, xsnumer, gradsT[j][paso], c=:magenta, lw=1.0, lab="")
    end
    fig

end

num_paso = 100
fig = plot_gradiente_espacial(paso = num_paso, extremos_y = 10.0)

savefig(fig, joinpath(TARGET_FOLDER, "grad_temp_barra_paso$num_paso.pdf"))

# Temporal gradient
function plot_gradiente_temporal(; nodo_grad = 66, δ = 1E-5, kmax = length(grad_lgg), NCURVAS=length(gradsT), extremos_y = [-20, 25])

    fig = plot(xlab=L"\textrm{Time } (\times 10^{-4})",
               ylab=L"\textrm{Temperature gradient}",
               legendfontsize=18, legend=:top,
               tickfont=font(25, "Times"), guidefontsize=20,
               xguidefont=font(30, "Times"),
               yguidefont=font(30, "Times"),
               ylims=(extremos_y[1], extremos_y[2]),
               xlims=(0.0, kmax*δ + 0.1e-4),
               xticks=([0.0, 0.00025, 0.0005, 0.00075, 0.001], [L"0", L"2.5", L"5.0", L"7.5", L"10.0"]),
               yticks=([-20, -10, 0, 10, 20], [L"-20", L"-10", L"0", L"10", L"20"]),
               bottom_margin=6mm, left_margin=6mm, right_margin=4mm, top_margin=3mm, size=(900, 600))

    kmax = length(grad_lgg)
    xsnumer = range(0, L, length=length(gradsT[1][1])+1)
    xsnumer = 0.5*(xsnumer[1:end-1]+xsnumer[2:end])
    x_nodo_grad = xsnumer[nodo_grad]

    fp_lgg_grad = gradiente_LGG_graddirs(grad_lgg = grad_lgg, paso_temporal = 1)

    fp_grad_nodo = [ReachSet(set(fp_lgg_grad(x_nodo_grad)), tspan(grad_lgg[1]))]
    for j in 2:kmax
        fp_lgg_grad = gradiente_LGG_graddirs(grad_lgg = grad_lgg, paso_temporal = j)
        push!(fp_grad_nodo, ReachSet(set(fp_lgg_grad(x_nodo_grad)), tspan(grad_lgg[j])))
    end
    fff = Flowpipe(fp_grad_nodo)

    plot!(fig, fff, vars=(0,1), lw=0.5, lc=:black, c=:lightblue)

    tiempo = range(0., step=δ, length=kmax)
    long_t = length(gradsT[1])
    for j in 1:NCURVAS
        y = [gradsT[j][i][nodo_grad] for i in 1:long_t]
        plot!(fig, tiempo, y[1:end-1], c=:magenta, lab="", lw=0.5)
    end

    fig
end

nnodo = 66
fig = plot_gradiente_temporal(nodo_grad=nnodo, kmax=100)
savefig(fig, joinpath(TARGET_FOLDER, "grad_temp_tiempo_paso$nnodo.pdf"))
