using Plots, Plots.PlotMeasures, LaTeXStrings

include("SDOF_Model.jl")

(@isdefined TARGET_FOLDER) ? nothing : TARGET_FOLDER = ""

p = sdof()
prob = IVP(p)
asol = analytic_solution(p)
T = period(p)

tmax = 4.0
dom = 0:0.001:tmax

f(α; alg=VREP(δ=α*T, approx_model=StepIntersect(setops=:concrete))) = solve(prob, tspan=(0.0, tmax), alg=alg)

# solution using numerical scheme
q = InitialValueProblem_quad(p)

g(α; alg=Newmark(Δt=α*T, α=0.25, δ=0.5)) = solve(q, alg, NSTEPS=round(Int, tmax / (α*T)));

function plot_x_vs_t(α)
    fig = plot(xlab=L"t", ylab=L"u(t)",
               legendfontsize=12, legend=:none,
               tickfont=font(20, "Times"),
               guidefontsize=20,
               xguidefont=font(30, "Times"),
               yguidefont=font(30, "Times"),
               xtick=([0.0, 1.0, 2.0, 3.0, 4.0], [L"0.0", L"1.0", L"2.0", L"3.0", L"4.0"]),
               ytick=([-1.0, -0.5, 0.0, 0.5, 1.0], [L"-1.0", L"-0.5", L"0.0", L"0.5", L"1.0"]),
               xlims=(0.0, 4.5), ylims=(-1.2, 1.3),
               bottom_margin=6mm, left_margin=2mm, right_margin=4mm, top_margin=3mm, size=(900, 600))

    L1 = LineSegment([3.3, -1.2], [3.7, -1.2])
    L2 = LineSegment([3.3, 1.2], [3.7, 1.2])
    L3 = LineSegment([3.3, -1.2], [3.3, 1.2])
    L4 = LineSegment([3.7, -1.2], [3.7, 1.2])

    #plot!(fig, UnionSetArray([L1, L2, L3, L4]), color=:black, style=:dash, lw=1.5, seriestype=:solid, marker=:none, alpha=1.)
    plot!(fig, L1, color=:black, style=:dash, lw=2.5, seriestype=:solid, marker=:none, alpha=1.)
    plot!(fig, L2, color=:black, style=:dash, lw=2.5, seriestype=:solid, marker=:none, alpha=1.)
    plot!(fig, L3, color=:black, style=:dash, lw=2.5, seriestype=:solid, marker=:none, alpha=1.)
    plot!(fig, L4, color=:black, style=:dash, lw=2.5, seriestype=:solid, marker=:none, alpha=1.)

    # Flowpipe
    println("Flowpipe solution (VREP)")
    fp = f(α) # warm-up
    @time fp = f(α)
    plot!(fig, fp, vars=(0, 1), lab=L"\alpha = %$α~\textrm{Flowpipe}", color=:lightblue, lw=0.5)

    # Numerical solution
    println("Numerical solution (Newmark)")
    y = g(α, alg=Newmark(Δt=α*T, α=0.25, δ=0.5)) # warm-up
    @time y = g(α, alg=Newmark(Δt=α*T, α=0.25, δ=0.5))
    tdom = 0:(α*T):tmax
    plot!(fig, tdom, displacements(y, 1), linestyle=:solid, marker=:circle, lab=L"\alpha = %$α~\textrm{Newmark}", color=:red)

    println("Numerical solution (Bathe)")
    y = g(α, alg=Bathe(Δt=α*T)) # warm-up
    @time y = g(α, alg=Bathe(Δt=α*T))
    plot!(fig, tdom, displacements(y, 1), linestyle=:solid, marker=:utriangle, lab=L"\alpha = %$α~\textrm{Bathe}", color=:green)

    # Analytic solution
    plot!(fig, dom, asol.(dom), lab=L"\textrm{Analytic}", color=:magenta, lw=1.2)

    fig
end

function plot_x_vs_t_zoom(α)
    fig = plot(xlab=L"t", ylab=L"u(t)", legend=:topleft,
               legendfont=font(15, "Times"),
               tickfont=font(20, "Times"),
               guidefontsize=20,
               xguidefont=font(30, "Times"),
               yguidefont=font(30, "Times"),
               xtick=([3.3, 3.5, 3.7], [L"3.3", L"3.5", L"3.7"]),
               ytick=([-1.0, -0.5, 0.0, 0.5, 1.0], [L"-1.0", L"-0.5", L"0.0", L"0.5", L"1.0"]),
               xlims=(3.28, 3.72), ylims=(-1.2, 1.2),
               bottom_margin=6mm, left_margin=2mm, right_margin=4mm, top_margin=3mm, size=(900, 600))

    a = 3.3; b = 3.7

    # flowpipe solutions
    aux = f(α)(a .. b)
    plot!(fig, aux[1:end-1], vars=(0, 1), lab=L"\textrm{Flowpipe}", color=:lightblue)

    # numerical solutions
    tdom = a:(α*T):b
    ai = round(Int, a / (α * T))
    bi = round(Int, b / (α * T))

    y = g(α, alg=Newmark(Δt=α*T, α=0.25, δ=0.5)) |> displacements
    plot!(fig, tdom, [yi[1] for yi in y[ai:bi]], linestyle=:solid, lw=2.5, marker=:circle, markersize=10.0, lab=L"\textrm{Newmark}", color=:red)

    y = g(α, alg=Bathe(Δt=α*T)) |> displacements
    plot!(fig, tdom, [yi[1] for yi in y[ai:bi]], linestyle=:solid, lw=2.5, marker=:utriangle, markersize=10.0, lab=L"\textrm{Bathe}", color=:green)

    # analytic solution
    domw = a:0.001:b
    plot!(fig, domw, asol.(domw), lab=L"\textrm{Analytic}", color=:magenta, lw=2.0)

    fig
end

# ==========================
# Figures for α = 0.05
# ==========================
println("Solving for α = 0.05")
fig = plot_x_vs_t(0.05)
savefig(fig, joinpath(TARGET_FOLDER, "sdof5a.pdf"))

fig = plot_x_vs_t_zoom(0.05)
savefig(fig, joinpath(TARGET_FOLDER, "sdof5b.pdf"))

# ==========================
# Figures for α = 0.1
# ==========================

println("Solving for α = 0.1")
fig = plot_x_vs_t(0.1)
savefig(fig, joinpath(TARGET_FOLDER, "sdof6a.pdf"))

fig = plot_x_vs_t_zoom(0.1)
savefig(fig, joinpath(TARGET_FOLDER, "sdof6b.pdf"))
