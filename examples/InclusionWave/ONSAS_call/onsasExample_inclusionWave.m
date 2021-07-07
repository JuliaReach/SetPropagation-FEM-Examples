% Example meshing InclusionWave using ONSAS

% command to generate msh file from Octave:
%   system('gmsh -2 inclusionCirc.geo')

clear all, close all

addpath(genpath(getenv('ONSAS_PATH')))

% scalar parameters (parametros promedio del tejido mamario)
E = 40e6     % 40 MPa
nu = 0.49 ;  %
p = 1.0 ;    %
thickness = 1 ;
# parametro promedio del tejido mamario

materials.hyperElasModel  = {'linearElastic'; 'linearElastic'} ;
materials.hyperElasParams = { [ E nu ]; [ E nu ] }      ;
materials.density         = { 1.0 ; 1.0 }      ;

elements.elemType = { 'node', 'edge', 'triangle' } ;

elements.elemTypeParams = { []; [] ; 2  } ;
elements.elemTypeGeometry = { []; thickness ; thickness } ;

boundaryConds.loadsCoordSys = {[]; []; 'global'  } ;
boundaryConds.loadsTimeFact = { []; []; @(t) t  } ;
boundaryConds.loadsBaseVals = { []; []; [ 0 0 -p 0  0 0 ]  } ;
boundaryConds.imposDispDofs = { [1 3] ; [3] ; []  } ;
boundaryConds.imposDispVals = { [0 0] ; [0] ; []  } ;
%md

initialConds = struct();
%md

[ mesh.nodesCoords, mesh.conecCell ] = meshFileReader( 'inclusionCirc.msh' )

analysisSettings.methodName    = 'alphaHHT' ;
analysisSettings.stopTolIts    = 30      ;
analysisSettings.stopTolDeltau = 1.0e-12 ;
analysisSettings.stopTolForces = 1.0e-12 ;
analysisSettings.finalTime      = 2       ;
analysisSettings.alphaHHT      = 0      ;
analysisSettings.deltaT        = 1      ;
%md
%md
%md### Output parameters
otherParams.problemName = 'inclusionWave' ;
otherParams.plotsFormat = 'vtk' ;
otherParams.spitMatrices = true ;
%md
[matUs, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;
