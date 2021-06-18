# -------------------------------------------------
function load_heat3d_matrices()

    path = joinpath("mats", "Heat3D.mat")
    vars = matread(path)

    K = vars["K"];
    C = vars["C"];
    timeIncr = vars["timeIncr"];
    fext = vars["fext"];
    qextTamb = vars["qextTamb"][:]; # vector
    QhG =  vars["QhG"][:]; # vector
    udots=vars["udots"];
    xs=vars["xs"];

    return K, C, qextTamb, QhG, timeIncr
end
# -------------------------------------------------


# -------------------------------------------------
function load_heat3d_params()

    # load mat
    path    = joinpath("mats", "params.mat")
    vars    = matread(path)

    # save variables
    rho        = vars["rho"] ;
    m          = vars["m"] ;
    QFH        = vars["QFH"] ;
    Tambmin    = vars["Tambmin"] ;
    Tambvar    = vars["Tambvar"] ;
    Tambvarmin = vars["Tambvarmin"] ;
    Tambvarmax = vars["Tambvarmax"] ;
    Tini       = vars["Tini"]   ;
    QFHmin     = vars["QFHmin"] ;
    QFHmax     = vars["QFHmax"] ;

    return rho, QFHmin, QFHmax, m, Tini, Tambmin, Tambvarmin, Tambvarmax, Tambvar, QFH
end
# -------------------------------------------------


# -------------------------------------------------
function load_heat3d_octaveSols()
    vars             = matread( joinpath("mats", "solsOctave.mat") );
    solAOctave       = vars["Ts3DA"]   ;
    solBOctave       = vars["Ts3DB"]   ;
    return solAOctave, solBOctave
end
# -------------------------------------------------
