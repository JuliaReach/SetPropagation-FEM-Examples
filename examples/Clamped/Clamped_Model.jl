#=
This model is taken from [1]. The theoretical derivations of the analytic solution can be found in [2].

[1] Malakiyeh, Mohammad Mahdi, Saeed Shojaee, and Klaus-Jürgen Bathe. "The Bathe time integration method revisited for prescribing desired numerical dissipation." Computers & Structures 212 (2019): 289-298.

[2] Mechanical Vibrations, Gerardin et al, page 250-251.
=#

using ReachabilityAnalysis, LinearAlgebra, LazySets, SparseArrays
using StructuralDynamicsODESolvers
using LazySets.Arrays

include("../krylov_concept.jl")
using FastExpm

LazySets.set_ztol(Float64, 1e-15)
LazySets.set_atol(Float64, 1e-15)
LazySets.set_rtol(Float64, 1e-15)

function clamped_bar(; N=1000, E=30e6, ρ=7.3e-4, A=1, L=200)
    ℓ = L / N
    K = (E * A / ℓ) * SymTridiagonal(fill(2.0, N), fill(-1.0, N))
    K[end, end] = E * A / ℓ

    M = (ρ * A * ℓ / 2) * Diagonal(vcat(fill(2.0, N-1), 1.0))
    M[end, end] = ρ * A * ℓ / 2

    return M, K
end

"""
    clamped_free(; N)

Free (no forcing term) instance of the clamped bar.
"""
function clamped_free(; N)
    M, K = clamped_bar(N=N)
    C = spzeros(N, N) # no damping
    sys = SecondOrderLinearContinuousSystem(M, C, K)
end

"""
    clamped_forced(; N, F=10e3, E=30e6, A=1)

Forced (constant force `F`) instance of the clamped bar.
"""
function clamped_forced(; N, F=10e3, E=30e6, A=1)
    M, K = clamped_bar(N=N)
    C = spzeros(N, N) # no damping
    F = vcat(zeros(N-1), F) # the right-most node is excited
    sys = SecondOrderAffineContinuousSystem(M, C, K, F)
end

# "nominal" step size
const dtn = 9.88e-7

# problem instances
s5 = clamped_forced(N=5);
s10 = clamped_forced(N=10);
s20 = clamped_forced(N=20);
s30 = clamped_forced(N=30);
s50 = clamped_forced(N=50);
s100 = clamped_forced(N=100);
s200 = clamped_forced(N=200);
s500 = clamped_forced(N=500);
s1000 = clamped_forced(N=1000);

# box initial condition
function prob_box(s, e)
    two_n = 2statedim(s)
    return InitialValueProblem(s, X0(two_n) + BallInf(two_n, e))
end

# solve with GLGM06
function solve_glg(s; δ=dtn, NSTEPS=10, model=Forward())
    # singleton initial condition
    prob = InitialValueProblem(s, Singleton(zeros(2statedim(s))))
    alg = GLGM06(δ=δ, approx_model=model)
    sol = solve(ivp, NSTEPS=NSTEPS, alg=alg, homogeneize=true)
end

# solve with LGG09 for nodeidx
function solve_lgg(s; δ=dtn, NSTEPS=100, nodeidx, model=Forward())
    # singleton initial condition
    prob = InitialValueProblem(s, Singleton(zeros(2statedim(s))))

    n = statedim(s)
    eplus = SingleEntryVector(nodeidx, 2n+1, 1.0)
    eminus = SingleEntryVector(nodeidx, 2n+1, -1.0)

    dirs = [eminus, eplus]

    alg = LGG09(δ=δ, template=dirs, approx_model=model)
    sol = solve(ivp, NSTEPS=NSTEPS, alg=alg, homogeneize=true)
end

# solve with BOX
# usage: sol = solve_box(s10, NSTEPS=100)
function solve_box(s; δ=dtn, NSTEPS=100, model=Forward())

    # ph = homogeneize(normalize(p))
    # @time pd = discretize(ph, dt, model) # ~35 sec for 1000 nodes

    alg = BOX(δ=δ, approx_model=model)

    # singleton initial condition
    prob = InitialValueProblem(s, Singleton(zeros(2statedim(s))))

    sbox = solve(p, alg=alg, NSTEPS=NSTEPS, homogeneize=true)
    return sbox
end

# Ec. (4.181) libro Geradin & Rixen (2015)
function _sol_analitica_clamped(x, t, Ncorte)

    F = 10E3
    E = 30E6
    dens = 7.3E-4
    L = 200
    A = 1

    # reescalar
    nMasses = 1000
    x = x*L/nMasses

    m = A*dens
    c = sqrt(E*A/m)

    α = 8*F*L/π^2/E/A

    β = π/2/L
    u = sin.( β * x ) * ( 1.0 .- cos.( β * c * t ) )
    v = sin.( β * x ) * β * c * sin.( β * c * t )

    for s in 2:Ncorte
        u .+= (-1)^(s-1)/(2*s-1)^2 * sin.( (2*s-1)* β * x ) * ( 1.0 .- cos.( (2*s-1)* β * c * t ) )
        v .+= (-1)^(s-1)/(2*s-1)^2 * sin.( (2*s-1)* β * x ) * ( (2*s-1)* β * c * sin.( (2*s-1)* β * c * t ) )
    end

    return u*α, v*α

end

function _sol_reachability_clamped(p; δt=dtn, NSTEPS=NSTEPS,
                                      nodeidx=[PosTwoThirds, VelTwoThirds],
                                      mk=150, tol=1e-12)

    # ==============================
    # Homogeneize
    # ==============================
    m = size(p.K, 1)
    Shom = _homogeneize(p)
    Aext = Shom.A
    @assert size(Aext, 1) == 2m+1

    # la condicion inicial es nula y se agrega la fuerza al final
    X0 = Vector{Float64}(undef, 2m+1)
    X0 = zeros(2m+1)
    X0[end] = 1.0

    # ==============================
    # Propagacion a tiempo continuo
    # ==============================
    #Aᵀδ = copy(transpose(Aext .* δt))
    num_nodes = isa(nodeidx, Int) ? 1 : length(nodeidx)

    #println("Calculando factor de bloating...")

    Eplus = _build_Eplus_optim(Aext, X0, δt; m=mk, tol=tol, cutoff=eps(Float64))

    #println("Propagando a tiempo continuo...")
    _resolver_clamped_expm(Aext, X0, m, δt, mk, tol, Eplus, num_nodes, NSTEPS, nodeidx)
end

function _resolver_clamped_expm(Aext, X0, m, δt, mk, tol, Eplus, num_nodes, NSTEPS, nodeidx)
    fastexpm = true

    if fastexpm
        Φᵀ = real.( fastExpm(transpose(Aext) .* δt, threshold=eps(Float64), nonzero_tol=eps(Float64)) )
    else
        Φᵀ = exp(Matrix(transpose(Aext) .* δt))
        Φᵀ = sparse(Φᵀ)
        droptol!(Φᵀ, eps(Float64))
    end

    q0 = transpose(transpose(X0) * Φᵀ) # para no calcular Φ y Φ^T

    # build Ω0 plus
    q0_msum_Eplus = Hyperrectangle(q0, Eplus.radius) # q0 ⊕ Eplus
    Ω0 = CH(Singleton(X0), q0_msum_Eplus)

    step_intersect = true
    if step_intersect
        Eplusstep = _build_Eplus_optim(Aext, q0, δt; m=mk, tol=tol, cutoff=eps(Float64))
        Ω0minus = CH(Singleton(q0), Singleton(X0) ⊕ Eplusstep)
        Ω0 = Ω0 ∩ Ω0minus
    end

    # Ω0box = box_approximation(Ω0, Hyperrectangle)

    Mout = Matrix{Float64}(undef, 2*num_nodes, NSTEPS)

    d = zeros(2m+1)
    @inbounds for j in 1:num_nodes
        node = nodeidx[j]
        d[node] = 1.0
    #    println("Propagando indice $node ...")
        idx = (2(j-1) + 1):2j
        reach_homog_dir_LGG09_modif!(view(Mout, idx, :), Ω0, Φᵀ, d, NSTEPS)
        d[node] = 0.0
    end
    δt0 = interval(0, δt)
    return Mout
end

function reach_homog_dir_LGG09_modif!(Mout, Ω₀, Φᵀ, d, NSTEPS)
    rᵢ = copy(d)
    rᵢ₊₁ = similar(rᵢ)

    @inbounds for i in 1:NSTEPS
        Mout[1, i] = ρ(rᵢ, Ω₀)
        Mout[2, i] = ρ(-rᵢ, Ω₀)

        # update cache for the next iteration
        mul!(rᵢ₊₁, Φᵀ, rᵢ)
        copy!(rᵢ, rᵢ₊₁)
    end
    return nothing
end
