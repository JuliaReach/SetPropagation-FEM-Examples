using ReachabilityAnalysis
using ExponentialUtilities
using Plots
using SparseArrays
using ReachabilityAnalysis: dim

_homogeneize(p::SecondOrderAffineContinuousSystem) = _homogeneize(p.M, p.K, p.b)

function _homogeneize(M, K, F)
    m = size(K, 1)
    Zm = spzeros(m, m)
    Im = sparse(1.0I, m, m)
    M = Diagonal(M)
    invM = inv(M)

    u0 = invM * F
    Aext = [Zm           Im      spzeros(m) ;
           -invM*K       Zm         u0      ;
            spzeros(1, 2m+1)                ]

    S = @system(x' = Aext*x) # sistema homogeneizado
end


function _reach_homog_krylov_LGG09_modif!(out, Ω₀::LazySet, Aᵀδ::AbstractMatrix,
                                          ℓ::AbstractVector, NSTEPS;
                                          hermitian=false, m=min(30, size(Aᵀδ, 1)), tol=1e-7)

    @assert size(out, 1) == 2
    @assert size(out, 2) == NSTEPS

    # initialization of the krylov subspace
    # esta parte genera una base de
    # \{Aᵀδ * ℓ, (Aᵀδ)^2 * ℓ, (Aᵀδ)^3 * ℓ, ..., (Aᵀδ)^m * ℓ}
    TA, Tb = eltype(Aᵀδ), eltype(ℓ)
    T = promote_type(TA, Tb)
    Ks = KrylovSubspace{T, real(T)}(length(ℓ), m)
    arnoldi!(Ks, Aᵀδ, ℓ; m=m, ishermitian=hermitian, tol=tol)

    # rᵢ stores is the cache for each vector: (Φᵀ)^i ℓ
    rᵢ = deepcopy(ℓ)

    @inbounds for i in 1:NSTEPS
        out[1, i] = ρ(rᵢ, Ω₀)
        out[2, i] = ρ(-rᵢ, Ω₀)
        expv!(rᵢ, i*1.0, Ks) # rᵢ <- exp(A^T δ * i) * ℓ
    end
    return out
end

# ------------------------------------------------------------------------
# Auxiliary methods that are not yet available in ReachabilityAnalysis.jl
# ------------------------------------------------------------------------

# calcula: out <- exp(A * NSTEPS * dt) * b
function _expv(A, b, NSTEPS, dt; hermitian=false, m=min(30, size(A, 1)), tol=1e-7)

    # initialization of the krylov subspace
    TA, Tb = eltype(A), eltype(b)
    T = promote_type(TA, Tb)
    Ks = KrylovSubspace{T, real(T)}(length(b), m)
    arnoldi!(Ks, A, b; m=m, ishermitian=hermitian, tol=tol)

    out = similar(b)
    expv!(out, NSTEPS*dt, Ks)

    return out
end

# calcula: out <- exp(A * NSTEPS * dt) * b
function _phiv(A, b, NSTEPS, dt; hermitian=false, m=min(30, size(A, 1)), tol=1e-7)

    # initialization of the krylov subspace
    TA, Tb = eltype(A), eltype(b)
    T = promote_type(TA, Tb)
    Ks = KrylovSubspace{T, real(T)}(length(b), m)
    arnoldi!(Ks, A, b; m=m, ishermitian=hermitian, tol=tol)

    out = Matrix{Float64}(undef, size(A, 1), 3)
    phiv!(out, NSTEPS*dt, Ks, 2)

    return view(out, :, 3) .* dt^2
end

#
# E₊ = ⊡(Φ₂ ⊡(A^2 X₀))
#
function _build_Eplus(A, X0, δt; m=min(30, size(A, 1)), tol=1e-7, Φ₂can=nothing)
    if isa(X0, Vector)
        X0 = Singleton(X0)
    end
    n = dim(X0)

    # build box E+
    H = symmetric_interval_hull(A^2 * X0)
    E₊_high = Vector{Float64}(undef, n)

    if isnothing(Φ₂can)
        Φ₂can = Φ₂_canonical(A, δt, m=m, tol=tol)
    end

    for i in 1:n
        Φ₂col = view(Φ₂can, :, i)
        E₊_high[i] = max(ρ(Φ₂col, H), abs(ρ(-Φ₂col, H)))
    end

    E₊_low = -E₊_high # since it is a SIH
    E₊ = Hyperrectangle(low=E₊_low, high=E₊_high)
end

function Φ₂_canonical(A, δt; m=min(size(A, 1)), tol=1e-7)
    n = size(A, 1)
    P3n = ReachabilityAnalysis._P_3n(abs.(A), δt, n)
    P3n_tr = copy(transpose(P3n))

    Φ₂can = Matrix{Float64}(undef, n, n)
    b = zeros(3n)
    for i in 1:n
        b[i] = 1.0
        out = _expv(P3n_tr, b, 1, 1.0, m=m, tol=tol)
        Φ₂can[:, i] .= view(out, 2n+1:3n)
        b[i] = 0.0
    end
    return Φ₂can
end

function _support_function_Eplus_baseline(d, A, δ, V)
    P = ReachabilityAnalysis.Φ₂(abs.(A), δ)
    E = symmetric_interval_hull(P * V)
    return ρ(d, E)
end

function _support_function_Eplus_krylov(d, A::AbstractMatrix{N}, δ,
                           V::AbstractHyperrectangle;
                           cutoff=eps(N), mk=min(30, size(A, 1)), tol=1e-7) where {N}

    # use only the indices in which V is not flat
    idx = findall(ri -> abs(ri) > cutoff, radius_hyperrectangle(V))
    Aabsᵀ = copy(transpose(abs.(A))) # sparse

    out = zero(N)
    n = size(A, 1)
    ei = zeros(N, n)
    @inbounds for i in idx
        ei[i] = one(N)

        R1 = phiv(δ, Aabsᵀ, ei, 2, m=mk, tol=tol) .* δ^2
        α = ρ(view(R1, :, 3), V)

        R2 = phiv(δ, Aabsᵀ, -ei, 2, m=mk, tol=tol) .* δ^2
        β = ρ(view(R2, :, 3), V)

        out += abs(d[i]) * max(α, abs(β))
        ei[i] = zero(N)
    end
    return out
end

function _support_function_Eplus_krylov_optim(d::AbstractVector{N}, A::AbstractMatrix{N}, δ,
                                              V::AbstractHyperrectangle;
                                              cutoff=eps(N), mk=min(30, size(A, 1)), tol=1e-7) where {N}

    # use only the indices in which V is not flat
    idx = findall(ri -> abs(ri) > cutoff, radius_hyperrectangle(V))
    Aabs = copy(abs.(A)) # sparse

    v = V.radius[idx]
    Maux = phiv(δ, Aabs, v, 2, m=mk, tol=tol)
    Pv = view(Maux, :, 3) .* δ^2 # rescale

    return dot(abs.(d), Pv)
end

function _build_Eplus_optim(A::SparseMatrixCSC{N, D}, X0::AbstractVector{N}, δt; m=min(30, size(A, 1)), tol=1e-7, cutoff=eps(N)) where {N, D}
    _build_Eplus_optim(A, Singleton(X0), δt; m=m, tol=tol, cutoff=cutoff)
end

# optimized version of the constructor of the Eplus set using a single Krylov subspace evaluation
# uses the idea in _support_function_Eplus_krylov_optim
function _build_Eplus_optim(A::SparseMatrixCSC{N, D}, X0::AbstractSingleton{N}, δt; m=min(30, size(A, 1)), tol=1e-7, cutoff=eps(N)) where {N, D}
    n = dim(X0)

    # compute V
    AX0 = linear_map(A, X0)
    V = symmetric_interval_hull(A * AX0)
    v = V.radius

    # compute Pv using Krylov
    Aabs = copy(abs.(A))
    Pv = _phiv(Aabs, v, 1, δt; m=m, tol=tol)

    # return the symmetric interval hull whose radius is Pv
    E⁺ = Hyperrectangle(zeros(n), Pv)

    return E⁺
end

function _build_Eplus_optim(A::SparseMatrixCSC{N, D}, X0::Hyperrectangle{N}, δt; m=min(30, size(A, 1)), tol=1e-7, cutoff=eps(N)) where {N, D}
    n = dim(X0)
    A2 = A * A; # fast if A sparse
    V = symmetric_interval_hull(A2 * X0)
    v = V.radius
    Aabs = copy(abs.(A))
    Pv = _phiv(Aabs, v, 1, δt; m=m, tol=tol)
    E⁺ = Hyperrectangle(zeros(n), Pv)
end
