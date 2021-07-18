using Plots, LaTeXStrings, Plots.Measures
using ReachabilityAnalysis: dim
using MAT

(@isdefined TARGET_FOLDER) ? nothing : TARGET_FOLDER = ""

# load the model
include("Heat1D_Model.jl")

# plot the center
# flowpipe
δ = 1e-5
file = joinpath(@__DIR__, "heat1d_99.mat")
@time sol = solve_heat1d(NSTEPS=30_000, δ=δ, file=file, model=StepIntersect(Forward()));

# numerical solutions
@time sol₊ = solve_heat1d_implicit_euler(file=file, δ=δ, tol=0.1, NSTEPS=30_000);
@time sol₋ = solve_heat1d_implicit_euler(file=file, δ=δ, tol=-0.1, NSTEPS=30_000);

node = 51 # 1 a 99 (nodo medio = 51)

sold₊ = displacements(sol₊)
sold₋ = displacements(sol₋)

tdom = times(sol₊)

sol₊_num(node) = [x[node] for x in sold₊[1:end-1]];
sol₋_num(node) = [x[node] for x in sold₋[1:end-1]];

sol₊_node = sol₊_num(node)
sol₋_node = sol₋_num(node);

# temperature evolution at the central node
function plot_node_vs_t()

    fig = plot(xlab=L"\textrm{Time (s)}", ylab=L"\textrm{Temperature at the center}",
               legendfontsize=18, legend=:topright,
               tickfont=font(20, "Times"), guidefontsize=20,
               xguidefont=font(30, "Times"),
               yguidefont=font(30, "Times"),
               #title=L"\textrm{SDOF with } \alpha = %$α",
               xtick=[0.0, 0.1, 0.2, 0.3], # , ytick=[-1.0, 0.0, 1.0],
               xlims=(0.0, 0.33), ylims=(0.0, 0.2, 0.4, 0.6, 0.8),
               bottom_margin=6mm, left_margin=6mm, right_margin=4mm, top_margin=3mm, size=(900, 600))

    #    fig = plot(xlab="Time (s)", ylab="Temperature at the center")

    # \delta = 10^{-5}
    plot!(fig, sol[1:5:end], vars=(0, node), alpha=.1, lc=:blue, lab=L"\textrm{Flowpipe}")

    plot!(fig, tdom[1:350:end], sol₊_node[1:350:end], linestyle=:solid, marker=:utriangle, lab=L"\textrm{Backward Euler}", color=:red)
    plot!(fig, tdom[1:350:end], sol₋_node[1:350:end], lab="", marker=:utriangle, linestyle=:solid, color=:red)

    # Analytic solution
    sol₊_analytic = evolucion_temp_teo(node/100., 1.1, tdom)
    sol₋_analytic = evolucion_temp_teo(node/100., 0.9, tdom)

    plot!(fig, tdom, sol₊_analytic, c=:magenta, linestyle = :solid, lab=L"\textrm{Analytic (PDE)}", lw=1.2)
    plot!(fig, tdom, sol₋_analytic, c=:magenta, linestyle = :solid, lab="")

    #xtick=([0.010, 0.015, 0.020, 0.025, 0.030], [L"0.010", L"0.015", L"0.020", L"0.025", L"0.030"])
    lens!(fig, [0.01, 0.03], [0.6, 0.85], inset = (1, bbox(0.58, 0.33, 0.4, 0.4)),
           subplot=2,
           tickfont=font(18, "Times"),
           xticks=([0.01, 0.020, 0.030], [L"0.01", L"0.02", L"0.03"]),
           yticks=([0.60, 0.65, 0.70, 0.75, 0.80, 0.85], ["", L"0.65", "", L"0.75", "", L"0.85"]))
    fig
end

fig = plot_node_vs_t()
savefig(fig, joinpath(TARGET_FOLDER, "temp_center.pdf"))


# temperature profile for different instants of time
function plot_espacial()

    fig = plot(xlab=L"\textrm{Position}", ylab=L"\textrm{Temperature}",
               legendfontsize=18, legend=:top,
               tickfont=font(20, "Times"), guidefontsize=20,
               xguidefont=font(30, "Times"),
               yguidefont=font(30, "Times"),
               ylims=(-0.05, 1.25),
               xlims=(-0.02, 1.02),
               #title=L"\textrm{SDOF with } \alpha = %$α",
               xtick=[0.0, 0.25, 0.50, 0.75, 1.0], # , ytick=[-1.0, 0.0, 1.0],
               ytick=[0.0, 0.25, 0.50, 0.75, 1.0, 1.25],
               bottom_margin=6mm, left_margin=6mm, right_margin=4mm, top_margin=3mm, size=(900, 600))

    paso = 5
    ancho_linea = 0.75
    col_linea = :blue
    marker_col = :blue

    tiempo_perfil = 0.0
    col = :blue

    fp = flowpipe_snapshot(sol, tiempo_perfil)
    xs, Tteo₋ = perfil_temp_teo(101, 0.9, tiempo_perfil)
    xs, Tteo₊ = perfil_temp_teo(101, 1.1, tiempo_perfil)
    plot!(fig, fp, lw=.1, c=col)
    plot!(fig, xs, Tteo₊, lc=col_linea, lw=ancho_linea, linestyle=:solid, lab="")
    plot!(fig, xs[1:paso:end], Tteo₊[1:paso:end], seriestype=:scatter, marker=:utriangle, c=marker_col, lab=L"t=0")
    plot!(fig, xs, Tteo₋, lc=col_linea, lw=ancho_linea, linestyle=:solid, lab="")
    plot!(fig, xs[1:paso:end], Tteo₋[1:paso:end], seriestype=:scatter, marker=:utriangle, c=marker_col, lab="")

    tiempo_perfil = 0.03
    col = :green
    col_linea = :green
    marker_col = :green

    fp = flowpipe_snapshot(sol, tiempo_perfil)
    xs, Tteo₋ = perfil_temp_teo(101, 0.9, tiempo_perfil)
    xs, Tteo₊ = perfil_temp_teo(101, 1.1, tiempo_perfil)
    plot!(fig, fp, lw=.1, c=col)
    plot!(fig, xs, Tteo₊, lc=col_linea, lw=ancho_linea, linestyle=:solid, lab="")
    plot!(fig, xs[1:paso:end], Tteo₊[1:paso:end], marker=:utriangle, c=marker_col, seriestype=:scatter, lab=L"t=0.03")
    plot!(fig, xs, Tteo₋, lc=col_linea, lw=ancho_linea, linestyle=:solid, lab="")
    plot!(fig, xs[1:paso:end], Tteo₋[1:paso:end], marker=:utriangle, c=marker_col, seriestype=:scatter, lab="")

    tiempo_perfil = 0.1
    col = :yellow
    col_linea = :orange
    marker_col = :orange

    fp = flowpipe_snapshot(sol, tiempo_perfil)
    xs, Tteo₋ = perfil_temp_teo(101, 0.9, tiempo_perfil)
    xs, Tteo₊ = perfil_temp_teo(101, 1.1, tiempo_perfil)
    plot!(fig, fp, lw=.1, c=col)
    plot!(fig, xs, Tteo₊, lc=col_linea, lw=ancho_linea, linestyle=:solid, lab="")
    plot!(fig, xs[1:paso:end], Tteo₊[1:paso:end], marker=:utriangle, c=marker_col, seriestype=:scatter, lab=L"t=0.1")
    plot!(fig, xs, Tteo₋, lc=col_linea, lw=ancho_linea, linestyle=:solid, lab="")
    plot!(fig, xs[1:paso:end], Tteo₋[1:paso:end], marker=:utriangle, c=marker_col, seriestype=:scatter, lw=0.0, lab="")

    return fig

end

fig = plot_espacial()
savefig(fig, joinpath(TARGET_FOLDER, "temp_barra.pdf"))
