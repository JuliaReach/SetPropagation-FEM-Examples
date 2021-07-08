% Example meshing InclusionWave using ONSAS

% command to generate msh file from Octave:
%   system('gmsh -2 inclusionCirc.geo')

clear all, close all

global spitMatrices
spitMatrices = true

addpath(genpath(getenv('ONSAS_PATH')))

% scalar parameters (parametros promedio del tejido mamario)
E = 20e6     % 20 MPa
nu = 0.49 ;  %
p = 1.0 ;    %
thickness = 1 ;
density = 950 ; % parametro promedio del tejido mamario

materials.hyperElasModel  = {'linearElastic'; 'linearElastic'} ;
materials.hyperElasParams = { [ E nu ]; [ 10*E nu ] }      ;
materials.density         = { density ; density }      ;

elements.elemType = { 'node', 'edge', 'triangle' } ;

elements.elemTypeParams = { []; [] ; 2  } ;
elements.elemTypeGeometry = { []; thickness ; thickness } ;

boundaryConds.loadsCoordSys = {[]; []; 'global'  } ;
boundaryConds.loadsTimeFact = { []; []; @(t) 1e4*( 3*(t<5e-5) - 3*(t<10e-5) + 1*(t<15e-5) ) } ;
boundaryConds.loadsBaseVals = { []; []; [ 0 0 -p 0  0 0 ]  } ;
boundaryConds.imposDispDofs = { [1 3] ; [3] ; []  } ;
boundaryConds.imposDispVals = { [0 0] ; [0] ; []  } ;
%md

initialConds = struct();
%md

[ mesh.nodesCoords, mesh.conecCell ] = meshFileReader( 'inclusionCirc.msh' ) ;

analysisSettings.methodName    = 'newmark' ;
analysisSettings.stopTolIts    = 30      ;
analysisSettings.stopTolDeltau = 1.0e-12 ;
analysisSettings.stopTolForces = 1.0e-12 ;
analysisSettings.finalTime      = 150e-5    ;
analysisSettings.alphaNM       = 0.25    ;
analysisSettings.deltaT        = .25e-5    ;
analysisSettings.deltaNM       = 0.5     ;
%md
%md
%md### Output parameters
otherParams.problemName = 'inclusionWave' ;
otherParams.plotsFormat = 'vtk' ;
%md
[matUs, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;
