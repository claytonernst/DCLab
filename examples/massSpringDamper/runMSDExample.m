%The first time running this script it will run fairly slowly, as it has
%to evaluate the model at many points and make fits.  However, a second
%running of this script should be significantly faster, as the
%ConsistencyTest and ResponsePrediction will load saved points.

DCsetup;

RO1 = ResponseObservation(1.82,0.2); 
RO2 = ResponseObservation(9,1);
RO3 = ResponseObservation(1.12,0.2); 

RM1 = ResponseModel(@msdLinearModel,'peakDis',1); % (function,feature,extra input (force))
RM2 = ResponseModel(@msdLinearModel,'rt95',1); 
RM3 = ResponseModel(@msdLinearModel,'ssMean',1); 

FP1 = FreeParameter('m',7.5,2.5); %[5 10]
FP2 = FreeParameter('b',1.125,0.875); %[0.25 2]
FP3 = FreeParameter('k',1.25,0.75); %[0.5 2]

Pair1 = ModelAndObservationPair(RO1,RM1,'peakDis');
Pair2 = ModelAndObservationPair(RO2,RM2,'rt95');
Pair3 = ModelAndObservationPair(RO3,RM3,'ssMean');

D = DCDataset([Pair1;Pair2;Pair3],[FP1;FP2;FP3]);

opt = DCOptions('surfaceTransformation',{'linXlogY'},'display','final','useAllPnts4Fit',0);

% Initial consistency analysis
ctestObj1 = ConsistencyTest(D,opt);

% Prediction for the peak displacement with a 2N force
RM_peakDis2N = ResponseModel(@msdLinearModel,'peakDis',2);
D_onlyFP = DCDataset([],[FP1;FP2;FP3]);

%These two predictions can use all the points generated in the consistency
%test to quickly regenerate the surfaces
predObj1 = ResponsePrediction(RM_peakDis2N,D_onlyFP,opt); %prediction without experiment constraints
predObj2 = ResponsePrediction(RM_peakDis2N,D,opt); %with experiment constraints - notice improvement

