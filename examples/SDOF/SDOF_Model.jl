# # Single degree-of-freedom

using ReachabilityAnalysis, StructuralDynamicsODESolvers
using ReachabilityAnalysis:discretize

# ### Equations of motion

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

# ================================
# Funciones para calcular PE y AD
# ================================

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

function compute_amplitude_decay(T, N, numpers, α)
    tf = numpers*T
    conjs = period_evolution2(T, tf, α)

    ## intervalos de tiempo en los que se da una vuelta
    Ts = [[tspan(R) for R in c] for c in conjs]

    ## divido entre el numero de vueltas hechas para tener los periodos
    periodos = [[R/i for (i, R) in enumerate(c)] for c in Ts]

    amplitudes = [[overapproximate(project(overapproximate(ConvexHullArray(set.(R)), Hyperrectangle), [1]), Interval).dat for R in c] for c in conjs]
    amps_min = [inf.(a1) for a1 in amplitudes]
    amps_min_num = AD_numerica.(amps_min, N, numpers)

    return amps_min_num
end

function compute_period_elongation(T, N, numpers, α)

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
    ##@. PE_RFEM_sup = (PE_RFEM_sup/T - 1) * 100

    PE_RFEM_inf = ([inf(Ai[end]) for Ai in intersecciones].-T)./T*100;

    return PE_RFEM_inf, PE_RFEM_sup
end
