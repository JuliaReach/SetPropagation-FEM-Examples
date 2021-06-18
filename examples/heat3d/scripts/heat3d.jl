

# -------------------------------------------------
function AhomCalc( C, K, qextTamb, QhG, rho, m)
    
    invC = inv( Matrix(C) );
    A = - invC * K |> Matrix
    n = size(A, 1)
    
    b0   =  invC * qextTamb ;
    bext =  invC * QhG * rho * m ;
    bsin =  invC * qextTamb ;

    omega = pi/12.;

    Ahom = [ A                    b0  bext  bsin      zeros(n) ;
             transpose( zeros(n))  0   0    0         0 ;
             transpose( zeros(n))  0  -m    0         0 ;
             transpose( zeros(n))  0   0    0         1 ;
             transpose( zeros(n))  0   0    -omega^2  0 ] ;

    return Ahom
end
# -------------------------------------------------
    
    
# -------------------------------------------------
function solve_heat3d_rfem(; model=StepIntersect(Forward(inv=false)),
                             X0=nothing, case=1)

    # load vectors and params
    K, C, qextTamb, QhG, timeIncr = load_heat3d_matrices()
    rho, QFHmin, QFHmax, m, Tini, Tambmin, Tambvarmin, Tambvarmax, Tambvar, QFH = load_heat3d_params()
    usA, usB = load_heat3d_octaveSols()
    
    δ = timeIncr;  n = size(K, 1);
    
    NSTEPS = ( size(usA,2)-1 ) ;
        
    # compute "homogeneized" system matrix
    Ahom = AhomCalc( C, K, qextTamb, QhG, rho, m ) ;

    QFH     = (QFHmax+QFHmin)*0.5 ;
    Tambvar = (Tambvarmax+Tambvarmin)*0.5 ;
    
    if case == 1
        # opcion singleton
        X0 = Singleton( Tini*ones(n)) ×
            Singleton((Tambmin+Tambvar*0.5)*ones( 1 ) )× Singleton( QFH*ones(1) ) ×
            Singleton((-Tambvar*0.5)*ones( 1 ) )× Singleton( zeros(1) ) ;

        prob = @ivp(x' = Ahom*x, x(0) ∈ X0)
        alg = BOX(δ=δ, approx_model=model)
        sol = solve(prob, alg=alg, NSTEPS=NSTEPS)
        
        X0 = Singleton( vcat( Tini*ones(n), 
                    (Tambmin+Tambvar*0.5)*ones(1) , 
                    QFHmax*ones(1) ,
                    (-Tambvar*0.5)*ones(1), 
                    zeros(1) ) ) ;

        probOrbitMax = @ivp(x' = Ahom*x, x(0) ∈ X0)
        solOrbit = solve( probOrbitMax, NSTEPS=NSTEPS, alg=ORBIT(δ=δ));
        # --------------------------------
        return sol, solOrbit

    elseif case == 2

        X0 = Singleton( Tini*ones(n)) ×
             Interval(Tambvarmin * 0.5 + Tambmin, Tambvarmax * 0.5 + Tambmin ) × 
             Interval( QFHmin, QFHmax ) ×
             Interval(-Tambvarmax * 0.5, -Tambvarmin*0.5) ×
             Singleton( zeros(1) ) ;
        
        prob = @ivp(x' = Ahom*x, x(0) ∈ X0)
        alg = BOX(δ=δ, approx_model=model)
        sol = solve(prob, alg=alg, NSTEPS=NSTEPS)
        # --------------------------------
        
        # --- Orbit min ---
        X0 = Singleton( 
              vcat( Tini*ones(n), 
                    (Tambmin+Tambvarmin*0.5)*ones(1) , 
                    QFHmin*ones(1) ,
                    (-Tambvarmin*0.5)*ones(1), 
                    zeros(1) ) ) ;

        probOrbit   = @ivp(x' = Ahom*x, x(0) ∈ X0)
        solOrbitmin = solve( probOrbit, NSTEPS=NSTEPS, alg=ORBIT(δ=δ));
        # --------------------------------

            # --- Orbit min ---
        X0 = Singleton( 
              vcat( Tini*ones(n), 
                    (Tambmin+Tambvarmax*0.5)*ones(1) , 
                    QFHmax*ones(1) ,
                    (-Tambvarmax*0.5)*ones(1), 
                    zeros(1) ) ) ;

        probOrbit   = @ivp(x' = Ahom*x, x(0) ∈ X0)
        solOrbitmax = solve( probOrbit, NSTEPS=NSTEPS, alg=ORBIT(δ=δ));
        # --------------------------------
     
        return sol, solOrbitmin, solOrbitmax
   end
    
end
# -------------------------------------------------



# -------------------------------------------------

maxtemp_heat3d(sol; nodes=1:dim(sol), kwargs...) = maxtemp(sol, nodes=nodes)
maxtemp_heat3d(; nodes=nothing, kwargs...) = maxtemp(solve_heat3d(kwargs...), nodes=nodes)
mintemp_heat3d(sol; nodes=1:dim(sol), kwargs...) = mintemp(sol, nodes=nodes)
mintemp_heat3d(; nodes=nothing, kwargs...) = mintemp(solve_heat3d(kwargs...), nodes=nodes)
# -------------------------------------------------

using StructuralDynamicsODESolvers
const SD = StructuralDynamicsODESolvers
using ReachabilityAnalysis: solve
# -------------------------------------------------



# -------------------------------------------------
function solve_heat3d_implicit_euler()
    
    # load vectors and params
    K, C, qextTamb, QhG, timeIncr = load_heat3d_matrices()
    rho, QFHmin, QFHmax, m, Tini, Tambmin, Tambvarmin, Tambvarmax, Tambvar, QFH = load_heat3d_params()
    usA, usB = load_heat3d_octaveSols()
    
    δ = timeIncr;  n = size(K, 1)
    NSTEPS = ( size(usA,2)-1 ) ;
    
    omega = pi/12.;
    
    M    = zeros(n, n)
    @show Tambvar = (Tambvarmin + Tambvarmax ) * 0.5 ;
    @show QFH = (QFHmin+ QFHmax)*0.5;
    
    R    = [ QhG*QFH*rho*m*exp(-m*δ*(i-1))+qextTamb*(Tambmin+Tambvar*(0.5 + 0.5*sin( omega*δ*(i-1)-pi/2.0 ))) for i in 1:NSTEPS+1 ] ;

    B   = Diagonal(ones(n));
    
    sys  =  SecondOrderConstrainedLinearControlContinuousSystem(M, Matrix(C), Matrix(K), B, nothing, R)

    U₀ = Tini*ones(n) ;
    
    prob = InitialValueProblem(sys, (U₀, U₀))
    alg  = BackwardEuler(Δt=δ)
    
    solBackwardEulerSDCaseA = SD.solve(prob, alg, NSTEPS=NSTEPS)

    R    = [ QhG*QFHmax*rho*m*exp(-m*δ*(i-1))+qextTamb*(Tambmin+Tambvarmax*(0.5 + 0.5*sin( omega*δ*(i-1)-pi/2.0 ))) for i in 1:NSTEPS+1 ] ;

    sys  =  SecondOrderConstrainedLinearControlContinuousSystem(M, Matrix(C), Matrix(K), B, nothing, R)
    prob = InitialValueProblem(sys, (U₀, U₀))

    solBackwardEulerSDCaseB = SD.solve(prob, alg, NSTEPS=NSTEPS)

    
    return solBackwardEulerSDCaseA, solBackwardEulerSDCaseB
end

