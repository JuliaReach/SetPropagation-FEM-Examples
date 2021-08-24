using LazySets
using LazySets.Arrays
using Plots, LaTeXStrings, Plots.PlotMeasures
using LinearAlgebra

(@isdefined TARGET_FOLDER) ? nothing : TARGET_FOLDER = ""
include("Clamped_Model.jl")

NSTEPS = 12000
nMasses = 1000

# computed parameters
finalTime = NSTEPS * dtn
tdom = range(0, finalTime, length = NSTEPS + 1)

PosTwoThirds = 700
PosMiddle    = round(Int, nMasses / 2)
VelTwoThirds = PosTwoThirds + nMasses
VelMiddle    = PosMiddle + nMasses

U₀ = zeros(nMasses)
V₀ = zeros(nMasses)
p = InitialValueProblem(s1000, (U₀, V₀))

println("Solving with Bathe method")
@time solBathe = solve(p, Bathe(dtn), NSTEPS=NSTEPS)

println("Solving with Newmark method")
@time solNewmark = solve(p, Trapezoidal(dtn), NSTEPS=NSTEPS)

uBathe_PosTwoThirds = displacements(solBathe, PosTwoThirds)
uBathe_VelTwoThirds = velocities(solBathe, PosTwoThirds)

uNewmark_PosTwoThirds = displacements(solNewmark, PosTwoThirds)
uNewmark_VelTwoThirds = velocities(solNewmark, PosTwoThirds)

# =====================
# Reachability results
# =====================

p = InitialValueProblem(s1000, Singleton(zeros(2*nMasses)))
println("Solving with ORBIT method")
@time sol_nobloating = solve(p, NSTEPS=NSTEPS, alg=ORBIT(δ=dtn))

println("Solving with LGG09 method")
@time Mout = _sol_reachability_clamped(s1000);
sol_PosTwoThirds = [ReachSet(Interval(-Mout[2, k], Mout[1, k]), (0.0 .. dtn) + (k-1) * dtn) for k in 1:size(Mout, 2)] |> Flowpipe
sol_VelTwoThirds = [ReachSet(Interval(-Mout[4, k], Mout[3, k]), (0.0 .. dtn) + (k-1) * dtn) for k in 1:size(Mout, 2)] |> Flowpipe


# =========================
# graficos desplazamientos
# =========================

# time window where the analytic solution is computed
tdomasol = range(0.0, finalTime, length=10_000)

asol_PosTwoThirds, asol_VelTwoThirds = _sol_analitica_clamped(PosTwoThirds, tdomasol, 500)

function plot_displacements()

    fig = plot(xlab=L"\textrm{Time}", ylab=L"\textrm{Displacement of node } %$PosTwoThirds",
               legend=:bottomright,
               legendfontsize=15,
               tickfont=font(20, "Times"),
               guidefontsize=20,
               xguidefont=font(20, "Times"),
               yguidefont=font(20, "Times"),
               xtick=([0, 0.0025, 0.005, 0.0075, 0.01],[L"0.000", L"0.0025", L"0.0050", L"0.0075", L"0.01"]),
               ytick=([0, 0.02, 0.04, 0.06, 0.08, 0.1], [L"0.00", L"0.02", L"0.04", L"0.06", L"0.08", L"0.1"]),
               xlims=(0.0, 0.012), ylims=(-0.002, 0.1),
               bottom_margin=6mm, left_margin=20mm, right_margin=4mm, top_margin=3mm, size=(900, 600))

    # flowpipe
    plot!(fig, sol_PosTwoThirds, vars=(0, 1), c=:lightblue, lc=:lightblue, lw=4.0, alpha=1., lab="") # lab=L"\textrm{Flowpipe}")

    # plot!(fig, tdom, uNewmark, seriestype=:path, marker=:none, color=:red, lab=L"\textrm{Newmark}", lw=1.0)
    # plot!(fig, tdom, uBathe, seriestype=:path, marker=:none, color=:green, lab=L"\textrm{Bathe}", lw=1.0)
    # plot!(fig, sol_nobloating, vars=(0, PosTwoThirds), seriestype=:path, marker=:none, lw = 2.0, c=:blue, lab=L"\textrm{Analytic (ODE)}")

    tZoomDispLeft  = 0.0055
    tZoomDispRight = 0.0058

    yZoomDispBott  = 0.08
    yZoomDispTop = 0.098

    # == Analytic ODE ==
    plot!(fig, sol_nobloating, vars=(0, PosTwoThirds), seriestype=:path, marker=:none, lw = 1.5, alpha=1, c=:blue, lab="")
    #plot!(fig, sol_nobloating[ai:paso:bi], vars=(0, VelTwoThirds), linetype=:scatter, marker=:square, markersize=5,  c=:blue, lab=L"\textrm{Analytic (ODE)}")


    plot!(fig, tdomasol, asol_PosTwoThirds, c=:magenta, lw=1., lab="") # , lab=L"\textrm{Analytic (PDE)}")

    PBotLeft  = [tZoomDispLeft,yZoomDispBott]
    PBotRight = [tZoomDispRight,yZoomDispBott]
    PTopLeft  = [tZoomDispLeft,yZoomDispTop]
    PTopRight = [tZoomDispRight,yZoomDispTop]

    L1 = LineSegment( PBotLeft, PBotRight )
    L2 = LineSegment( PBotRight, PTopRight )
    L3 = LineSegment( PTopRight, PTopLeft )
    L4 = LineSegment( PTopLeft, PBotLeft )

    plot!(fig, L1, color=:black, style=:dash, lw=3.0, seriestype=:solid, marker=:none, alpha=1.)
    plot!(fig, L2, color=:black, style=:dash, lw=3.0, seriestype=:solid, marker=:none, alpha=1.)
    plot!(fig, L3, color=:black, style=:dash, lw=3.0, seriestype=:solid, marker=:none, alpha=1.)
    plot!(fig, L4, color=:black, style=:dash, lw=3.0, seriestype=:solid, marker=:none, alpha=1.)
end

fig = plot_displacements()
savefig(fig, joinpath(TARGET_FOLDER, "u_node700.pdf"))

# ================================
# graficos desplazamientos (ZOOM)
# ================================

function plot_displacements_zoom()

    tZoomDispLeft  = 0.00559
    tZoomDispRight = 0.00572

    idxtleft = findfirst(t -> t > tZoomDispLeft, tdom)
    idxtright = findfirst(t -> t > tZoomDispRight, tdom)

    yZoomDispBott  = 0.0925
    yZoomDispTop = 0.0937

    fig = plot(xlab=L"\textrm{Time } (\times 10^{-3})", ylab=L"\textrm{Displacement of node } %$PosTwoThirds",
               legend=:bottomright,
               legendfontsize=15,
               tickfont=font(20, "Times"),
               guidefontsize=20,
               xguidefont=font(20, "Times"),
               yguidefont=font(20, "Times"),
               xtick=([0.00560, 0.00562, 0.00564, 0.00566, 0.00568, 0.00570, 0.00572], [L"5.60", L"5.62", L"5.64", L"5.66", L"5.68", L"5.70", L"5.72"]),
               ytick=([0.09250, 0.09275, 0.09300, 0.09325, 0.09350], [L"0.09250", L"0.09275", L"0.09300", L"0.09325", L"0.09350"]),
               ylims=(yZoomDispBott-0.00005, 0.09355), xlims=(5.6e-3, tZoomDispRight),
               bottom_margin=6mm, left_margin=30mm, right_margin=10mm, top_margin=3mm, size=(900, 600))

    # flowpipe
    plot!(fig, sol_PosTwoThirds(tZoomDispLeft .. tZoomDispRight), vars=(0, 1), c=:lightblue, lc=:black, lw=0.5, alpha=1., lab=L"\textrm{Flowpipe}")

    # solucion analitica
    tdomasol_zoom = range(tZoomDispLeft, tZoomDispRight, length=15_000)

    asol_PosTwoThirds_zoom, asol_VelTwoThirds_zoom = _sol_analitica_clamped(PosTwoThirds, tdomasol_zoom, 2000);

    plot!(fig, tdomasol_zoom, asol_PosTwoThirds_zoom, c=:magenta, lw=1.5, lab=L"\textrm{Analytic (PDE)}")

    paso = 8

    ai = idxtleft
    bi = idxtright

    # == Newmark ==
    # continua sin markers
    plot!(fig, tdom[ai:bi], uNewmark_PosTwoThirds[ai:bi], seriestype=:path, marker=:none, lab="", color=:red, lw=2.0)

    # con markers pero salteados
    plot!(fig, tdom[ai:paso:bi], uNewmark_PosTwoThirds[ai:paso:bi], linetype=:scatter, marker=:circle, markersize=5, lab=L"\textrm{Newmark}", color=:red)

    # == Bathe ==
    # continua sin markers
    plot!(fig, tdom[ai:bi], uBathe_PosTwoThirds[ai:bi], seriestype=:path, marker=:none, lab="", color=:green, lw=2.0)

    # con markers pero salteados
    plot!(fig, tdom[ai:paso:bi], uBathe_PosTwoThirds[ai:paso:bi], linetype=:scatter, marker=:utriangle, markersize=5, lab=L"\textrm{Bathe}", color=:green)

    # == Analytic ODE ==
    plot!(fig, sol_nobloating, vars=(0, PosTwoThirds), seriestype=:path, marker=:none, lw = 1.5, alpha=1, c=:blue, lab=L"\textrm{Analytic (ODE)}")
    plot!(fig, sol_nobloating[ai:paso:bi], vars=(0, PosTwoThirds), linetype=:scatter, marker=:square, markersize=5,  c=:blue, lab="")

end

fig = plot_displacements_zoom()
savefig(fig, joinpath(TARGET_FOLDER, "u_node700_zoom.pdf"))

# =========================
# velocidad
# =========================

function plot_velocity()

    fig = plot(xlab=L"\textrm{Time}", ylab=L"\textrm{Velocity of node } %$PosTwoThirds",
               legend=:bottomright,
               legendfontsize=15,
               tickfont=font(20, "Times"),
               guidefontsize=20,
               xguidefont=font(20, "Times"),
               yguidefont=font(20, "Times"),
               xtick=([0, 0.0025, 0.005, 0.0075, 0.01],[L"0.000", L"0.0025", L"0.0050", L"0.0075", L"0.01"]),
               ytick=([-75, -50, -25, 0, 25, 50, 75], [L"-75", L"-50", L"0.0", L"25", L"50", L"75"]),
               xlims=(0.0, 0.012), # , ylims=(-90.0, 90.0),
               bottom_margin=6mm, left_margin=20mm, right_margin=4mm, top_margin=3mm, size=(900, 600))

    # flowpipe
    plot!(fig, sol_VelTwoThirds, vars=(0, 1), c=:lightblue, lc=:lightblue, lw=4.0, alpha=1., lab="") # lab=L"\textrm{Flowpipe}")

    # plot!(fig, tdom, uNewmark, seriestype=:path, marker=:none, color=:red, lab=L"\textrm{Newmark}", lw=1.0)

    # plot!(fig, tdom, uBathe, seriestype=:path, marker=:none, color=:green, lab=L"\textrm{Bathe}", lw=1.0)
    # plot!(fig, sol_nobloating, vars=(0, PosTwoThirds), seriestype=:path, marker=:none, lw = 2.0, c=:blue, lab=L"\textrm{Analytic (ODE)}")

    tZoomDispLeft  = 0.008175
    tZoomDispRight = 0.008375

    yZoomDispBott  = 40
    yZoomDispTop = 92

    # == Analytic ODE ==
    plot!(fig, sol_nobloating, vars=(0, VelTwoThirds), seriestype=:path, marker=:none, lw = 1.5, alpha=1, c=:blue, lab="")
    #plot!(fig, sol_nobloating[ai:paso:bi], vars=(0, VelTwoThirds), linetype=:scatter, marker=:square, markersize=5,  c=:blue, lab=L"\textrm{Analytic (ODE)}")

    # == Analytic PDE ==
    tdomasol = range(0, finalTime, length=1000)

    _, asol_VelTwoThirds = _sol_analitica_clamped(PosTwoThirds, tdomasol, 1000);

    plot!(fig, tdomasol, asol_VelTwoThirds, c=:magenta, lw=1.5, lab="") # , lab=L"\textrm{Analytic (PDE)}")

    # ZOOM BOX
    PBotLeft  = [tZoomDispLeft,yZoomDispBott]
    PBotRight = [tZoomDispRight,yZoomDispBott]
    PTopLeft  = [tZoomDispLeft,yZoomDispTop]
    PTopRight = [tZoomDispRight,yZoomDispTop]

    L1 = LineSegment( PBotLeft, PBotRight )
    L2 = LineSegment( PBotRight, PTopRight )
    L3 = LineSegment( PTopRight, PTopLeft )
    L4 = LineSegment( PTopLeft, PBotLeft )

    plot!(fig, L1, color=:black, style=:dash, lw=3.0, seriestype=:solid, marker=:none, alpha=1.)
    plot!(fig, L2, color=:black, style=:dash, lw=3.0, seriestype=:solid, marker=:none, alpha=1.)
    plot!(fig, L3, color=:black, style=:dash, lw=3.0, seriestype=:solid, marker=:none, alpha=1.)
    plot!(fig, L4, color=:black, style=:dash, lw=3.0, seriestype=:solid, marker=:none, alpha=1.)

end

fig = plot_velocity()
savefig(fig, joinpath(TARGET_FOLDER, "v_node700.pdf"))

# ==============================
# graficos velocidad (ZOOM)
# ==============================

function plot_velocity_zoom()


    tZoomDispLeft  = 0.008175
    tZoomDispRight = 0.008375

    yZoomDispBott  = 30
    yZoomDispTop = 92

    idxtleft = findfirst(t -> t > tZoomDispLeft, tdom)
    idxtright = findfirst(t -> t > tZoomDispRight, tdom)


    fig = plot(xlab=L"\textrm{Time } (\times 10^{-3})", ylab=L"\textrm{Velocity of node } %$PosTwoThirds",
               legend=:bottomright,
               legendfontsize=15,
               tickfont=font(20, "Times"),
               guidefontsize=20,
               xguidefont=font(20, "Times"),
               yguidefont=font(20, "Times"),
               xtick=([0.00820, 0.00825, 0.00830, 0.00835], [L"8.20", L"8.25", L"8.30", L"8.35"]),
               ytick=([30, 40, 50, 60, 70, 80, 90], [L"30", L"40", L"50", L"60", L"70", L"80", L"90"]),
               ylims=(yZoomDispBott, yZoomDispTop),
               xlims=(tZoomDispLeft, tZoomDispRight),
               bottom_margin=6mm, left_margin=30mm, right_margin=10mm, top_margin=3mm, size=(900, 600))

    # flowpipe
    plot!(fig, sol_VelTwoThirds(tZoomDispLeft .. tZoomDispRight), vars=(0, 1), c=:lightblue, lc=:black, lw=0.5, alpha=1., lab=L"\textrm{Flowpipe}")

    # solucion analitica
    tdomasol_zoom = range(tZoomDispLeft, tZoomDispRight, length=103)

    asol_VelTwoThirds_zoom, asol_VelTwoThirds_zoom = _sol_analitica_clamped(PosTwoThirds, tdomasol_zoom, 1000);

    plot!(fig, tdomasol_zoom, asol_VelTwoThirds_zoom, c=:magenta, lw=1.5, lab=L"\textrm{Analytic (PDE)}")

    paso = 8

    ai = idxtleft
    bi = idxtright

    # == Newmark ==
    # continua sin markers
    plot!(fig, tdom[ai:bi], uNewmark_VelTwoThirds[ai:bi], seriestype=:path, marker=:none, lab="", color=:red, lw=2.0)

    # con markers pero salteados
    plot!(fig, tdom[ai:paso:bi], uNewmark_VelTwoThirds[ai:paso:bi], linetype=:scatter, marker=:circle, markersize=5, lab=L"\textrm{Newmark}", color=:red)

    # == Bathe ==
    # continua sin markers
    plot!(fig, tdom[ai:bi], uBathe_VelTwoThirds[ai:bi], seriestype=:path, marker=:none, lab="", color=:green, lw=2.0)

    # con markers pero salteados
    plot!(fig, tdom[ai:paso:bi], uBathe_VelTwoThirds[ai:paso:bi], linetype=:scatter, marker=:utriangle, markersize=5, lab=L"\textrm{Bathe}", color=:green)

    # == Analytic ODE ==
    plot!(fig, sol_nobloating, vars=(0, VelTwoThirds), seriestype=:path, marker=:none, lw = 1.5, alpha=1, c=:blue, lab=L"\textrm{Analytic (ODE)}")
    plot!(fig, sol_nobloating[ai:paso:bi], vars=(0, VelTwoThirds), linetype=:scatter, marker=:square, markersize=5,  c=:blue, lab="")
end

fig = plot_velocity_zoom()
savefig(fig, joinpath(TARGET_FOLDER, "v_node700_zoom.pdf"))
