using ReachabilityAnalysis
using StructuralDynamicsODESolvers

# returns A = -inv(C) * K
function load_heat1d(; file="heat1d_99.mat")
    path = joinpath("mats", file)
    vars = matread(path)

    K = vars["K"]
    C = vars["C"]

    return K, C
end

function initial_distribution(n=100, tol=0.1)
    xs = range(0, 1, length=n+2)[2:end-1]
    c0 = @. sin(pi*xs) + 0.5*sin(3*pi*xs)
    r0 = tol * c0
    X0 = Hyperrectangle(c0, r0)
    return X0
end

function solve_heat1d(; kwargs...)
    prob, alg, NSTEPS = _solve_heat1d(; kwargs...)
    return solve(prob, alg=alg, NSTEPS=NSTEPS)
end

function _solve_heat1d(; file="heat1d_99.mat",
                         δ=1e-5,
                         NSTEPS=10_000,
                         model=Forward(inv=true),
                         X0=nothing,
                         alg=BOX(δ=δ, approx_model=model))

    ## load model, we convert to a dense matrix to take the inverse in the discretization
    K, C = load_heat1d(file=file)
    A = -inv(Matrix(C)) * K |> Matrix
    n = size(A, 1)

    if isnothing(X0)
        X0 = initial_distribution(n)
    end

    prob = @ivp(x' = A*x, x(0) ∈ X0)

    return prob, alg, NSTEPS
end


maxtemp_heat1d(sol; nodes=1:dim(sol), kwargs...) = maxtemp(sol, nodes=nodes)
maxtemp_heat1d(; nodes=nothing, kwargs...) = maxtemp(solve_heat1d(kwargs...), nodes=nodes)

mintemp_heat1d(sol; nodes=1:dim(sol), kwargs...) = mintemp(sol, nodes=nodes)
mintemp_heat1d(; nodes=nothing, kwargs...) = mintemp(solve_heat1d(kwargs...), nodes=nodes)

# dado el numero de puntos, discretiza el dominio [0, 1] devuelve el perfil de temperaturas analitico
# en los tiempos dados en un tiempo fijo
function perfil_temp_teo(npuntos, Tampl, t::AbstractFloat)
    rho       = 1. ;
    cSpHe     = 1. ;
    kCond     = 1. ;

    xs = range(0, 1, length=npuntos)
    alpha = kCond / ( rho * cSpHe ) ;
    Tout = Tampl * ( exp(-(  pi )^2 * alpha * t ) * sin.(pi * xs) + exp(-(3*pi)^2 * alpha * t ) * 0.5 * sin.( 3 * pi * xs ) );
    return xs, Tout
end

# para un nodo dado, me da la temperatura en funcion del tiempo en los instantes
# de tiempo dados por el vector t
# xs : posicion
# Tampl = 1 : perfil de temperaturas con suma de senos
# Tampl = 1.1 : lo anterior pero con un 10 % mas
function evolucion_temp_teo(xs::AbstractFloat, Tampl, t::AbstractVector)
    rho       = 1. ;
    cSpHe     = 1. ;
    kCond     = 1. ;

    alpha = kCond / ( rho * cSpHe ) ;
    Tout = Tampl * ( exp.(-(  pi )^2 * alpha * t ) * sin.(pi * xs) + exp.(-(3*pi)^2 * alpha * t ) * 0.5 * sin.( 3 * pi * xs ) );
    return Tout
end

# Solucion numerica con Backward Euler
# use tol = 0.1 for 10% increase in the initial uncertainty
# use tol = -0.1 for 10% decrease in the initial uncertainty
function solve_heat1d_implicit_euler(; file="heat1d_99.mat", δ=1e-5, tol=0.1, NSTEPS=10_000, U₀=nothing)

    K, C = load_heat1d(file=file)

    n = size(K, 1)
    M = zeros(n, n)
    R = zeros(n)
    sys = SecondOrderAffineContinuousSystem(M, C, K, R)
    X0 = initial_distribution(n, abs(tol))

    ## extremal value
    x₀ = initial_state_singleton(X0, tol, U₀)

    prob = InitialValueProblem(sys, (x₀, x₀))

    alg = BackwardEuler(Δt=δ)
    return solve(prob, alg, NSTEPS=NSTEPS)
end

function initial_state_singleton(X0, tol, U₀)
    if isnothing(U₀)
        x₀ = X0.center + sign(tol) * X0.radius
    elseif U₀ == :sample
        x₀ = hcat(sample(X0, 1)...) # sample(X0, 1)
    elseif U₀ == :zigzag₊
        x₀ = [X0.center[i] + (-1)^i*X0.radius[i] for i in 1:n]
    elseif U₀ == :zigzag₋
        x₀ = [X0.center[i] + (-1)^(i-1)*X0.radius[i] for i in 1:n]
    end
    return x₀
end

## no-op
ReachabilityAnalysis.convexify(R::ReachabilityAnalysis.AbstractReachSet) = R

function _gradient_cpa(sol, t::AbstractFloat)
    X = sol(t) |> convexify
    n = dim(X)  # for 100 nodes, n = 99

    ST = Interval{Float64, IA.Interval{Float64}}
    grad = Vector{ST}(undef, n+1)

    dx = 1/n

    ## condiciones de borde (cero en los bordes)
    R = overapproximate(Projection(X, 1:1), Interval)
    grad[1] = set(R).dat * (1/dx) |> Interval
    L = overapproximate(Projection(X, n:n), Interval)
    grad[n+1] = -set(L).dat * (1/dx) |> Interval

    for i in 1:n-1
        L = overapproximate(Projection(X, i:i), Interval)
        R = overapproximate(Projection(X, i+1:i+1), Interval)
        grad[i+1] = Interval((set(R).dat - set(L).dat) * (1/dx))
    end
    dt = tspan(X)

    return ReachSet(CartesianProductArray(grad), dt)
end

function _gradient(sol, t::AbstractFloat)
    Tg = _gradient_cpa(sol, t)

    arr = array(set(Tg))
    n = dim(sol) + 1
    xs = range(0, 1, length=n+1)
    dx = [Interval(xs[i], xs[i+1]) for i in 1:n]
    u = UnionSetArray([dxi × ai for (dxi, ai) in zip(dx, arr)])
end

## function para calcular el gradiente de una sol numerica
## TODO improve using single loop
function gradiente_numerica(sol, dx)
    U = sol.U
    N = length(U) # numero de pasos temporales
    out = [vcat(x[1], diff(x, dims=1)[:], -x[end]) .* (1/dx) for x in U]
    return out
end

struct SolBE{S, G}
    displacement::S
    gradient::G
end

# example uso _varias_sol_be:
#
#   num_paso = 100
#   n_curvas = 100
#   tol = 0.1
#   δ = 1E-5
#   xs, out = _varias_sol_be(δ=δ, tol=tol, num_paso = num_paso, n_curvas=n_curvas)
#   figg = plot()
#   for i in 1:length(out)-2
#       plot!(figg, xs, out[i], lab="")
#   end
#   plot!(figg, xs, out[end-1], c=:magenta, lab="")
#   plot!(figg, xs, out[end], c=:magenta, lab="")
#   fig
function varias_sol(; file="heat1d_99.mat", δ, tol, NSTEPS, NCURVAS)
    ## SM = Vector{Matrix{Float64}}
    ## SV = Vector{Vector{Float64}}
    gradientes = [] # Vector{SolBE{SM, SV}}()
    temperaturas = []

    ## dx = 1/(99 + 1)

    ## calcular sol random
    ## separo el primer paso para calcular dx
    sol = solve_heat1d_implicit_euler(file=file, δ=δ, tol=tol, NSTEPS=NSTEPS, U₀=:sample)
    dx = 1/(length(sol.U[1])+1)
    grad = gradiente_numerica(sol, dx)
    push!(gradientes, grad) # SolBE(sol.U, gr))
    push!(temperaturas, sol)

    ## pasos siguientes
    for k in 2:NCURVAS
        sol = solve_heat1d_implicit_euler(file=file, δ=δ, tol=tol, NSTEPS=NSTEPS, U₀=:sample)
        grad = gradiente_numerica(sol, dx)
        push!(gradientes, grad) # SolBE(sol.U, gr))
        push!(temperaturas, sol)
    end

    ## calcular sol extremas
    for t in [tol, -tol]
        sol = solve_heat1d_implicit_euler(file=file, δ=δ, tol=t, NSTEPS=NSTEPS, U₀=nothing)
        grad = gradiente_numerica(sol, dx)
        push!(gradientes, grad) # SolBE(sol.U, gr))
    end

    for Uopt in [:zigzag₊, :zigzag₋]
        sol = solve_heat1d_implicit_euler(file=file, δ=δ, tol=tol, NSTEPS=NSTEPS, U₀=Uopt)
        grad = gradiente_numerica(sol, dx)
        push!(gradientes, grad) # SolBE(sol.U, gr))
        push!(temperaturas, sol)
    end

    return temperaturas, gradientes
end

## Codigo Para analisis de gradientes
function LGG_dirs_solucion(n)
    ## n = statedim(prob)

    M = Diagonal(ones(n))
    M = [Vector(M[:, i]) for i in 1:size(M, 1)]

    return vcat(M, .-M)

end

function LGG_dirs_gradiente(n)
    ## n = statedim(prob)

    B = Bidiagonal(ones(n), -ones(n-1), :U)
    B = [Vector(B[:, i]) for i in 2:size(B, 2)]
    vini = zeros(n); vini[1] = 1.0
    vend = zeros(n); vend[end] = -1.0
    insert!(B, 1, vini);
    push!(B, vend);

    return vcat(B, .-B);

end

## grad_lgg = solucion con LGG usando LGG_grad_dirs
function gradiente_LGG_graddirs(; grad_lgg, L=1.0, paso_temporal)
    n = Int(0.5*length(grad_lgg[1].sf))
    Δx = L/n
    xs = range(0, L, step=Δx)
    α = 1/Δx
    intervs_xs = [IA.Interval(xs[i], xs[i+1]) for i in 1:n]
    Flowpipe([ReachSet(
                  Interval(-grad_lgg[paso_temporal].sf[i+n] * α, grad_lgg[paso_temporal].sf[i] * α),
                  intervs_xs[i])
              for i in 1:n])
end

## solu_lgg = solucion con LGG usando LGG_solu_dirs (direcciones cartesianas)
function gradiente_LGG_soludirs(; solu_lgg, LGG_grad_dirs, L=1.0, paso_temporal)
    n = Int(0.5*length(LGG_grad_dirs))
    Δx = L/n
    xs = range(0, L, step=Δx)
    α = 1/Δx
    intervs_xs = [IA.Interval(xs[i], xs[i+1]) for i in 1:n]
    Flowpipe([ReachSet(
            Interval(-α*ρ(LGG_grad_dirs[i+n], solu_lgg[paso_temporal]), α*ρ(LGG_grad_dirs[i], solu_lgg[paso_temporal])),
            intervs_xs[i])
              for i in 1:n])
end
