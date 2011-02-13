function [disp fit optim misc] = decompose(opt)
%DCCOMPOSE Method of DCOptions separates the options into categories. 
%
%   [DISP FIT OPTIM MISC] = DECOMPOSE(OPT) takes the DCOptions object and breaks it into four
%   structures, each containing the options/settings indicated by their name.
%
%   This method is used by PolyDataset to determine display and fitting settings.


%options for display
disp.display = opt.display;
disp.guiHandle = opt.guiHandle;
disp.plotFitProgress = opt.plotFitProgress;

%options for fitting
%fit.analysisMode = opt.analysisMode;
fit.surrogateType = opt.surrogateType;
fit.surfaceType = opt.surfaceType;
fit.surfaceFittingMode = opt.surfaceFittingMode;
fit.surfaceTransformation = opt.surfaceTransformation;
fit.fitNorm = opt.fitNorm;
fit.minFitIter = opt.minFitIter;
fit.maxFitIter = opt.maxFitIter;
fit.nSuccessfulFitIter = opt.nSuccessfulFitIter;
fit.maxPnts4Fit = opt.maxPnts4Fit;
fit.useAllPnts4Fit = opt.useAllPnts4Fit;
fit.fitConvergenceTol = opt.fitConvergenceTol;
fit.nPntsPerCoeff4OneShot = opt.nPntsPerCoeff4OneShot;
fit.nPntsPerParam4ActiveParamSel = opt.nPntsPerParam4ActiveParamSel;
fit.activeParamSelCutOff = opt.activeParamSelCutOff;
fit.nPntsPerCoeff4Validation = opt.nPntsPerCoeff4Validation;
% fit.nPntsPerCoeff4Metamodel = opt.nPntsPerCoeff4Metamodel;
% fit.nPntsPerCoeff4MetamodelValidate = opt.nPntsPerCoeff4MetamodelValidation;
fit.nLocalValidationSearches = opt.nLocalValidationSearches;
% fit.nLocalValidationSearches4Metamodel = opt.nLocalValidationSearches4Metamodel;
fit.maxDenSwing = opt.maxDenSwing;
fit.nComputer = opt.nComputer;

%options for optimization
optim.omitInnerBound = opt.omitInnerBound;
optim.omitOuterBound = opt.omitOuterBound;
optim.maxBranchBoundIter = opt.maxBranchBoundIter;
optim.branchBoundTermTol = opt.branchBoundTermTol;
optim.branchBoundQGapTol = opt.branchBoundQGapTol;
optim.nRestart = opt.nRestart;
optim.tolFun = opt.tolFun;
optim.tolCon = opt.tolCon;
optim.sedumiParEps = opt.sedumiParEps;
optim.fullBlockHRelaxation = opt.fullBlockHRelaxation;
optim.constraints = opt.constraints;
optim.paramOptimObjectiveFctnNorm = opt.paramOptimObjectiveFctnNorm;

%misc options
misc.fileName2Save = opt.fileName2Save;
    
    
   
   
    
