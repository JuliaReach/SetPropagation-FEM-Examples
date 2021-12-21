using Plots
using Colors
using Statistics
using MAT
import Random

include("lambs_krylov.jl")

δt = step_size_lamb(load_lambs(mesh=2), 1.0);

# ---------------------------------------------
# Generar figuras condicion inicial singleton
# ---------------------------------------------
generar_figura_malla2_singleton_position()
generar_figura_malla2_singleton_velocity()

# ---------------------------------------------
# Generar datos envolvente
# ---------------------------------------------

vecsamples = [1, 10, 100, 1000, 10_000, 100_000]

# corridas de warm-up
out_flowpipe_pos, out_extrema_malla2_pos, fig_pos = generar_figura_malla2_envolventes_position(vecsamples=[1, 10]);
out_flowpipe_vel, out_extrema_malla2_vel, fig_vel = generar_figura_malla2_envolventes_velocity(vecsamples=[1, 10]);

# corridas largas
Random.seed!(1234)
out_flowpipe_pos, out_extrema_malla2_pos, fig_pos = generar_figura_malla2_envolventes_position(vecsamples=vecsamples);
out_flowpipe_vel, out_extrema_malla2_vel, fig_vel = generar_figura_malla2_envolventes_velocity(vecsamples=vecsamples);

# ---------------------------------------------
# Calculo de normas del flowpipe
# ---------------------------------------------

# dt siempre el mismo para cada tramo, y la ventana de tiempo va hasta 1 segundo,
# entonces el dt se saca de factor comun

# max(max_t |v_max(t)|, max_t |v_min(t)|)
out_flowpipe_vel_norma_inf_envolvente = max(abs(ρ([1.0], out_flowpipe_vel)), abs(ρ([-1.0], out_flowpipe_vel)))
# salida: 0.0012224968252669921

# δt * sum_{i=1}^N max(|v_max(t)|, |v_min(t)|)
δt = step_size_lamb(load_lambs(mesh=2), 1.0)
@assert tspan(out_flowpipe_vel[1]).hi ≈ δt
out_flowpipe_vel_norma_1_envolvente = sum(max(abs(ρ([1.0], set(X))), abs(ρ([-1.0], set(X)))) for X in out_flowpipe_vel) * δt
# salida: 0.0008133221497265329

# ---------------------------------------------
# Calculo de normas de numericas
# ---------------------------------------------

out_extrema_malla2_vel_norma_inf_envolvente = _norma_inf_envolvente_numericas(out_extrema_malla2_vel, vecsamples)
#=
6-element Vector{Float64}:
 0.00056983336205989
 0.0005752808846025316
 0.0005758510647681003
 0.0005821999057203168
 0.0006118261295338235
 0.0006220764784910799
 =#

out_extrema_malla2_vel_norma_1_envolvente = _norma_1_envolvente_numericas(out_extrema_malla2_vel, vecsamples)
#=
julia> out_extrema_malla2_vel_norma_1_envolvente = _norma_1_envolvente_numericas(out_extrema_malla2_vel, vecsamples)
6-element Vector{Float64}:
 9.268765778974471e-5
 0.00013515950042666908
 0.0001661339910196396
 0.00018520059379205858
 0.0001997982127229552
 0.00021423496225308065
 =#

# ------------------------------
# Guardar datos para las tablas
# ------------------------------

# convert the keys of a dictionary into strings; the values might themselves
# be dictionaries, in which case the operation is applied recursively
# moreover, all keys are appended a capital letter N, e.g. "1000", it is converted to "N1000"
function _to_str_keys(d::Dict)
    strkeys = "N" .* string.(keys(d))

    if first(d)[2] isa Dict
        return Dict(strkeys .=> _to_str_keys.(values(d)))
    else
        return Dict(strkeys .=> values(d))
    end
end

using MAT

results  = Dict{String, Any}("vecsamples" => vecsamples,
                             "dt" => δt,
                             "tiempo_flowpipe" => "8.504446 seconds (9.82 k allocations: 292.719 MiB, 1.10% gc time)",

                             "out_flowpipe_pos_lo" => [set(X).dat.lo for X in out_flowpipe_pos],
                             "out_flowpipe_pos_hi" => [set(X).dat.hi for X in out_flowpipe_pos],
                             "out_flowpipe_pos_tstart" => [tstart(X) for X in out_flowpipe_pos],
                             "out_flowpipe_pos_tend" => [tend(X) for X in out_flowpipe_pos],

                             "out_flowpipe_vel_lo" => [set(X).dat.lo for X in out_flowpipe_vel],
                             "out_flowpipe_vel_hi" => [set(X).dat.hi for X in out_flowpipe_vel],
                             "out_flowpipe_vel_tstart" => [tstart(X) for X in out_flowpipe_vel],
                             "out_flowpipe_vel_tend" => [tend(X) for X in out_flowpipe_vel],

                             "out_flowpipe_vel_norma_inf_envolvente" => out_flowpipe_vel_norma_inf_envolvente,
                             "out_flowpipe_vel_norma_1_envolvente" => out_flowpipe_vel_norma_1_envolvente,
                             "out_extrema_malla2_vel_norma_inf_envolvente" => out_extrema_malla2_vel_norma_inf_envolvente,
                             "out_extrema_malla2_vel_norma_1_envolvente" => out_extrema_malla2_vel_norma_1_envolvente,

                             "out_extrema_malla2_pos_times" => out_extrema_malla2_pos.idx,
                             "out_extrema_malla2_pos_U_min" => out_extrema_malla2_pos.U_min,
                             "out_extrema_malla2_pos_U_max" => out_extrema_malla2_pos.U_max,
                             "out_extrema_malla2_pos_Udot_min" => out_extrema_malla2_pos.U′_min,
                             "out_extrema_malla2_pos_Udot_max" => out_extrema_malla2_pos.U′_max,
                             "out_extrema_malla2_pos_hist" => _to_str_keys(out_extrema_malla2_pos.hist),
                             "out_extrema_malla2_pos_t" => collect(out_extrema_malla2_pos.t),
                             "out_extrema_malla2_pos_times" => _to_str_keys(out_extrema_malla2_pos.times),

                             "out_extrema_malla2_vel_times" => out_extrema_malla2_vel.idx,
                             "out_extrema_malla2_vel_U_min" => out_extrema_malla2_vel.U_min,
                             "out_extrema_malla2_vel_U_max" => out_extrema_malla2_vel.U_max,
                             "out_extrema_malla2_vel_Udot_min" => out_extrema_malla2_vel.U′_min,
                             "out_extrema_malla2_vel_Udot_max" => out_extrema_malla2_vel.U′_max,
                             "out_extrema_malla2_vel_hist" => _to_str_keys(out_extrema_malla2_vel.hist),
                             "out_extrema_malla2_vel_t" => collect(out_extrema_malla2_vel.t),
                             "out_extrema_malla2_vel_times" => _to_str_keys(out_extrema_malla2_vel.times)
                             );

matwrite("results.mat", results)

# ---------------
# Generar figura
# ---------------

generar_figura_malla2_envolvente_velocidad()
