struct FEMProblem
    s # system
    coords
    connec
    neumdofs
    fixeddofs
    b # Forcing term
    Fnotzerocoord
end

_homogeneize(p::FEMProblem) = _homogeneize(p.s.M, p.s.K, p.b)

# nodes of interest:
#    nodeidx_hor = 3
#    nodeidx_ver = 4
function load_lambs(; mesh=1,
                      file_matrices=joinpath(@__DIR__, "mesh$mesh", "matrices.mat"),
                      file_mesh=joinpath(@__DIR__, "mesh$mesh", "mesh.mat"),
                      file_loads=joinpath(@__DIR__, "mesh$mesh", "loads.mat"))

    # load matrices
    d = matread(file_matrices)
    K = d["KTred"]
    M = d["massMatred"]
    @assert isdiag(M)
    M = Diagonal(M) # it is diagonal
    n = size(K, 1)
    C = spzeros(n, n)
    sys = SecondOrderLinearContinuousSystem(M, C, K)

    neumdofs = d["neumdofs"]
    neumdofs = Int.(neumdofs)[:]

    # load mesh
    mesh = matread(file_mesh)
    coords = mesh["mesh"]["nodesCoords"]
    connec = mesh["mesh"]["conecCell"]
    connec = Int.(connec)

    # fixed degrees of freedom
    fixeddofs = setdiff(1:6*size(coords, 1), neumdofs)

    # load external force definition
    df = matread(file_loads)
    F = df["factorLoadsFextCell"][1, 3][:]
    F = F[neumdofs]
    Fnotzerocoord = findall(!iszero, F)

    return FEMProblem(sys, coords, connec, neumdofs, fixeddofs, F, Fnotzerocoord)
end

# ======================
# Funciones de carga
# ======================

Heaviside(t) = t >= 0 ? 1 : 0

# funcion auxiliar que devuelve el factor de carga aplicado en la coordenada coord
function factor_de_carga(m, NSTEPS; T=1.0, coord=5, func=Heaviside)
    tiempos = range(0.0, T, length=NSTEPS)

    F = Vector{Vector{Float64}}(undef, NSTEPS)
    for i in 1:NSTEPS
        aux = zeros(m)
        aux[coord] .= func(tiempos[i])

        F[i] = aux
    end
    return F
end

function step_size_lamb(p, δtfac)
    nElem = sum(p.connec[:, 2] .== 3) # number of triangular elements
    CFL = 0.125
    cP  = 3200
    h   = sqrt(2*3200^2 / nElem) * δtfac
    δt = CFL * h / cP
end
