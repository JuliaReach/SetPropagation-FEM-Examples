using ReachabilityAnalysis, Plots, Statistics
using LaTeXStrings
using Plots.PlotMeasures
using MAT

include("SDOF.jl")

vars = matread("calculos_octave_PE_AD.mat");
α       = vars["alphas"]';
N       = vars["N"]'; # numero de periodos usados para el calculo
numpers = vars["npers"]'; # numero de periodos de la senal
AD_Bathe_teo_octave   = vars["AD_Bathe_teo"];
AD_Bathe_num_octave   = vars["AD_Bathe_num"];
AD_Newmark_teo_octave = vars["AD_Newmark_teo"];
AD_Newmark_num_octave = vars["AD_Newmark_num"];
PE_Bathe_octave       = vars["vectorPETsOctave_Bathe"];
PE_Newmark_octave     = vars["vectorPETsOctave_Newmark"];
PE_Bathe_teo_octave   = vars["PE_Bathe_teo"];
PE_Newmark_teo_octave = vars["PE_Newmark_teo"];

T = 1/2;
PE_RFEM_inf, PE_RFEM_sup = _compute_PE(T, N, numpers, α);
AD_RFEM_inf = _compute_AD(T, N, numpers, α);

function _plot_PE()
    fig = plot(xlab=L"\alpha", ylab=L"\% PE", legend=:topleft,
               legendfontsize=15,
               tickfont=font(20, "Times"),
               guidefontsize=20,
               xguidefont=font(30, "Times"),
               yguidefont=font(30, "Times"),
               xtick=([0.0, 0.05, 0.10, 0.15, 0.20], [L"0.0", L"0.05", L"0.10", L"0.15", L"0.20"]),
               ytick=([0.0, 2.5, 5.0, 7.5, 10.0], [L"0.0", L"2.5", L"5.0", L"7.5", L"10.0"]),
               xlims=(0.0, 0.22), ylims=(-0.3, 11.0),
               bottom_margin=6mm, left_margin=20mm, right_margin=4mm, top_margin=3mm, size=(900, 600))
    plot!(fig, α[1:2:end], PE_Bathe_octave[1:2:end], c=:green, linestyle=:solid, lw=2.5, markersize=6, lab=L"\textrm{Bathe}", marker=:utriangle)
    plot!(fig, α[1:2:end], PE_Newmark_octave[1:2:end], c=:red, linestyle=:solid, lw=2.5, markersize=6, lab=L"\textrm{Newmark}", marker=:circle)
    plot!(fig, α[1:2:end], PE_RFEM_sup[1:2:end], c=:blue, linestyle=:solid, lw=2.5, markersize=6, lab=L"\textrm{Flowpipe}", marker=:square)
    fig
end

fig = _plot_PE()
savefig(fig, "PE_all.pdf")

function _plot_AD()
    fig = plot(xlab=L"\alpha", ylab=L"\% AD", legend=:topleft,
               legendfontsize=15,
               tickfont=font(20, "Times"),
               guidefontsize=20,
               xguidefont=font(30, "Times"),
               yguidefont=font(30, "Times"),
               xtick=([0.0, 0.05, 0.10, 0.15, 0.20], [L"0.0", L"0.05", L"0.10", L"0.15", L"0.20"]),
               ytick=([0.0, 1.0, 2.0, 3.0], [L"0.0", L"1.0", L"2.0", L"3.0"]),
               xlims=(0.0, 0.22), ylims=(-0.3, 3.2),
               bottom_margin=6mm, left_margin=20mm, right_margin=4mm, top_margin=3mm, size=(900, 600))
    plot!(fig, α[1:2:end], AD_Bathe_num_octave[1:2:end], c=:green, linestyle=:solid, lw=2.5, markersize=6, lab=L"\textrm{Bathe}", marker=:utriangle)
    plot!(fig, α[1:2:end], AD_Newmark_num_octave[1:2:end], c=:red, linestyle=:solid, lw=2.5, markersize=6, lab=L"\textrm{Newmark}", marker=:circle)
    plot!(fig, α[1:2:end], AD_RFEM_inf[1:2:end], c=:blue, linestyle=:solid, lw=2.5, markersize=6, lab=L"\textrm{Flowpipe}", marker=:square)
    fig
end
fig = _plot_AD()
savefig(fig, "AD_all.pdf")
