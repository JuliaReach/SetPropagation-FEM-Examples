% Example meshing InclusionWave using ONSAS

clear all, close all
addpath( genpath( '../../../ONSAS' ) )
% scalar parameters
E = 1e2 ; nu = 0.25 ; p = .5e-4 ; thickness = 1 ;

materials.hyperElasModel  = {'linearElastic'} ;
materials.hyperElasParams = { [ E nu ] }      ;

elements.elemType = { 'node', 'edge', 'triangle' } ;

elements.elemTypeParams = { []; [] ; 2  } ;
elements.elemTypeGeometry = { []; thickness ; thickness } ;

boundaryConds.loadsCoordSys = {[]; []; 'local'  } ;
boundaryConds.loadsTimeFact = { []; []; @(t) t  } ;
boundaryConds.loadsBaseVals = { []; []; [ p 0  0 0  0 0 ]  } ;
boundaryConds.imposDispDofs = { [1] ; [3] ; []  } ;
boundaryConds.imposDispVals = { [0] ; [0] ; []  } ;
%md

initialConds = struct();
%md

[ mesh.nodesCoords, mesh.conecCell ] = meshFileReader( 'inclusionCirc.msh' )

analysisSettings.methodName    = 'newtonRaphson' ;
analysisSettings.stopTolIts    = 30      ;
analysisSettings.stopTolDeltau = 1.0e-12 ;
analysisSettings.stopTolForces = 1.0e-12 ;
analysisSettings.finalTime      = 2       ;
analysisSettings.deltaT        = 1      ;
%md
%md
%md### Output parameters
otherParams.problemName = 'linearPlaneStrain' ;
otherParams.plotsFormat = 'vtk' ;
%md
[matUs, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;
