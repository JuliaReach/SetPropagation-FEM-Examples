using ReachabilityAnalysis
const RA = ReachabilityAnalysis
const IA = IntervalArithmetic
using StructuralDynamicsODESolvers # ODE solvers
using MAT
using LinearAlgebra
using SparseArrays

using ReachabilityAnalysis: ReachSolution, discretize, normalize, solve, center, step_size

# reduce default tolerance for vertex enumeration
LazySets.set_ztol(Float64, 1e-12)

# compute the maximum temperature at the given nodes
# by default we compute the max temperature in all the nodes,
# but you can specify a vector of indices to search;
# for example pass nodes=1:100 to make the search only on the first 100 nodes
function maxtemp(sol::ReachSolution; nodes=nothing)
    # initialize search
    Tmax = -Inf; kmax = -1; imax = -1

    if isnothing(nodes)
        nodes = 1:dim(sol)
    end

    # for each reach-set we compute the max bounds
    # in every direction
    for (k, R) in enumerate(sol)
        X = set(R)
        q = center(X) + radius_hyperrectangle(X)
        θi = argmax(view(q, nodes))
        θ = q[θi]
        if θ > Tmax
            Tmax = θ
            kmax = k
            imax = θi
        end
    end

    δ = step_size(sol.alg)
    dtmax = tspan(sol[kmax])
    return (Tmax=Tmax, node=imax, step=kmax, dt=dtmax)
end

# compute the minimum temperature on the given nodes
function mintemp(sol::ReachSolution; nodes=nothing)
    # initialize search
    Tmin = Inf; kmin = -1; imin = -1

    if isnothing(nodes)
        nodes = 1:dim(sol)
    end

    # for each reach-set we compute the max bounds
    # in every direction
    for (k, R) in enumerate(sol)
        X = set(R)
        q = center(X) - radius_hyperrectangle(X)
        θi = argmax(view(q, nodes))
        θ = q[θi]
        if θ < Tmin
            Tmin = θ
            kmin = k
            imin = θi
        end
    end

    δ = step_size(sol.alg)
    dtmin = tspan(sol[kmin])
    return (Tmin=Tmin, node=imin, step=kmin, dt=dtmin)
end

# ===============================================================
# Methods to plot flowpipe "snapshots"
# ===============================================================

# no-op
ReachabilityAnalysis.convexify(R::ReachabilityAnalysis.AbstractReachSet) = R

function flowpipe_snapshot(sol, t)
    n = dim(sol)
    solt = sol(t) |> convexify
    X = [set(overapproximate(Projection(solt, i:i), Interval)) for i in 1:n]
    R = ReachSet(CartesianProductArray(X), tspan(solt))
    _flowpipe_espacial(R)
end

# cartesian product array of intervals
function _flowpipe_snapshot(X::ReachSet{N, CartesianProductArray{N, Interval{N, IA.Interval{N}}}}) where {N}
    n = dim(X)
    xs = range(0, 1, length=n+2)

    fp_espacial = array(set(X))
    n_x = length(fp_espacial)
    fpx1 = [Interval(xs[k+1], xs[k+1]) × fp_espacial[k] for k in 1:n_x]

    T0 = 0.0 # temperatura en las puntas
    aux = [Interval(xs[1], xs[1]) × Interval(T0[1], T0[1])]
    for k = 1:n_x
        push!(aux, fpx1[k])
    end
    push!(aux, Interval(xs[end], xs[end]) × Interval(T0[end], T0[end]))
    fpx2 = UnionSetArray([ConvexHull(aux[k], aux[k+1]) for k in 1:length(aux)-1])

    return fpx2
end
