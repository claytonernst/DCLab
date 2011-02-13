%% Script for Trent/Andy/Michael paper
% August 14, 2010

% addpath('c:\unzipped\SeDuMi_1_3\SeDuMi_1_3')
% addpath('c:\unzipped\SeDuMi_1_3\SeDuMi_1_3\conversion')
% addpath('c:\unzipped\SeDuMi_1_3\SeDuMi_1_3\examples')
% rmpath('C:\Documents and Settings\musyn\InstalledMatlabTools\DClabV2\SeDuMi_1_21')

%% Import class
import DClab.*

%% Load GRI polys
D = createGRI;

%% Check number of parameters and number of exp/data pairs
D.nParameters
%%
nexp = D.nPairs

%% Create H simply by eliminating all ModelsObservationPairs
H = D;
H.ModelAndObservationPair(1:nexp) = [];

%% Copy, and eliminate F5 (prior analysis indicates it is a problem)
D.ModelAndObservationPair(37).name
E = D;
E.ModelAndObservationPair(37) = [];

%% Check Consistency of E
C = ConsistencyTest(E);
[C.LB C.UB]

%% Extract Model for F5, look at some properties
Mf5 = D.ModelAndObservationPair(37).ResponseModel;
D.ModelAndObservationPair(37).constraintLowerLimit
D.ModelAndObservationPair(37).constraintUpperLimit
D.ModelAndObservationPair(37).name
D.ModelAndObservationPair(37).observationUncertainty
D.ModelAndObservationPair(37).observedValue
D.ModelAndObservationPair(37).ResponseModel
D.ModelAndObservationPair(37).ResponseObservation

%% Prediction of F5 on H
PonH = ResponsePrediction(Mf5,H);
LowerBndRange = 10.^[PonH.LBo PonH.LBi]
UpperBndRange = 10.^[PonH.UBi PonH.UBo]

%% Prediction of F5 on E
P = ResponsePrediction(Mf5,E);
LowerBndRange = 10.^[P.LBo P.LBi]
UpperBndRange = 10.^[P.UBi P.UBo]

%% Look at Prediction Sensitivities
% LowerBound prediction to lower bounds of ExperimentalObservations
[sLL,sLLidx] = max(abs(P.outerBndSens.lower.expl));
sLL
sLLidx
bar(abs(P.outerBndSens.lower.expl))

%%
% LowerBound prediction to upper bounds of ExperimentalObservations
[sLU,sLUidx] = max(abs(P.outerBndSens.lower.expu));
sLU
sLUidx

%%
% UpperBound prediction to lower bounds of ExperimentalObservations
[sUL,sULidx] = max(abs(P.outerBndSens.upper.expl));
sUL
sULidx

%%
% UpperBound prediction to upper bounds of ExperimentalObservations
[sUU,sUUidx] = max(abs(P.outerBndSens.upper.expu));
sUU
sUUidx
bar(abs(P.outerBndSens.upper.expu))

%% Reduced Prediction, based on Sensitivities
% Prediction bounds are most sensitive to ModelObservation pairs 36 and 37
% within DataSet E.  Keep these as reduced DataSet, and redo prediction
R = E;
R.ModelAndObservationPair([1:35 38:end]) = [];
RedP = ResponsePrediction(Mf5,R);
LowerBndRange = 10.^[RedP.LBo RedP.LBi]
UpperBndRange = 10.^[RedP.UBi RedP.UBo]

%% Increase all experimental uncertainties by 20%
% The function |createGRIb| is just a copy of |createGRI|, which loads the
% polynomials, and the uncertainties that Xiaoqing derived.  The new
% version, |createGRIb| simply scales the uncertainty by the 4th argument,
% as a factor, in this case increasing by 20%.
Di = createGRIb([],[],[],1.10);

%% Copy, and eliminate F5
Ei = Di;
Ei.ModelAndObservationPair(37) = [];

%% Prediction of F5 on Ei
% Note that there is not much difference between these at the original
% prediction interval.  This illustrates that what is shown in the paper is not
% simply that F5 is extremely sensitive to small increases in the feasible
% set.   Indeed, the paper shows that the approximations put forth, while
% seemingly harmless, incur very large changes in prediction, unlike modest
% percentage increases in all of the experimental constraints that define
% the feasible set.
Pi = ResponsePrediction(Mf5,Ei);
LowerBndRange = 10.^[Pi.LBo Pi.LBi]
UpperBndRange = 10.^[Pi.UBi Pi.UBo]

%% Comments
% I didn't compute it, but we can back out the global bounds that the
% sensitivities imply, and couple those with the absolute changes made to
% the experimental uncertainties, and convince ourselves that the
% predictions for Pi (with relaxed uncertainties) are within the provable
% bounds implied by the sensitivity measures.


