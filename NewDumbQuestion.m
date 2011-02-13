%% Consistency doesn't seem to match up

%% Import class
import DClab.*

%% Load GRI polys
D = createGRI;
nexp = length(D.ModelAndObservationPair)

%% Copy DataSet, and eliminate F5 from copy
D.ModelAndObservationPair(37).name
E = D;
E.ModelAndObservationPair(37) = [];


load GRIDataStruct
trg = GRIDataStruct; %this is essentially Michael's GRI target object
allnums = unique([trg.var]);

%% Extract Model for F5 from original
Mf5 = D.ModelAndObservationPair(37).ResponseModel;

%%
opts = DCOptions;
opts.maxBranchBoundIter = 1;
C = ConsistencyTest(E,opts);

%% Prediction of F5 on E
opts = DCOptions;
opts.maxBranchBoundIter = 5;
P = ResponsePrediction(Mf5,E,opts);
LowerBndRange = 10.^[P.LBo P.LBi]
UpperBndRange = 10.^[P.UBi P.UBo]

%% WarmStart, getting more B&B
% 32.836 34.135; 35.609 36.052, exit on 8
opts.maxBranchBoundIter = 10;
P10 = ResponsePrediction(P,opts);
LowerBndRange = 10.^[P10.LBo P10.LBi]
UpperBndRange = 10.^[P10.UBi P10.UBo]

%% Skip inner bounds, tighter tolerance
% still exited after 8, seemed to do IB anyway,
% and got the same outer bound as above, 32.836, 36.052
opt.maxBranchBoundIter = 10;
opt.branchBoundTermTol = 0.001; %need it small, otherwise it meets the
opt.omitInnerBound = 1; %This is for speed (and obviously optional).  It
Ptrent = ResponsePrediction(Mf5,E,opts);
LowerBndRange = 10.^[Ptrent.LBo Ptrent.LBi]
UpperBndRange = 10.^[Ptrent.UBi Ptrent.UBo]

%% Compare with pg 117 in thesis
% Note that outer bounds are significantly worse than what is described in
% the thesis.