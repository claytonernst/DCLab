%% create GRI-Mech 3.0 Dataset

dset = createGRI

%This dataset has 102 FreeParameters and 77 ModelAndObservationPairs
n = dset.nParameters
m = dset.nPairs


%% Consistency

ctest = ConsistencyTest(dset)

%Dataset is inconsistant.  Let's find the parameter or experiment that
%consistency measure is most sensitive to.

%% Sensitivity

sens = ctest.upperBndSens;

%create a 2x2 plot with the sensitivities
M = max(abs([sens.paraml; sens.paramu; sens.expl; sens.expu]));
subplot(2,2,1)
bar(sens.paramu)
ylim([0,M])
ylabel('Sens. to up. bnd.')
subplot(2,2,2)
bar(sens.expu)
ylim([0,M])
subplot(2,2,3)
bar(sens.paraml)
ylim([-M,0])
ylabel('Sens. to low. bnd.')
xlabel('Parameter index')
subplot(2,2,4)
bar(sens.expl)
ylim([-M,0])
xlabel('Experiment index')

%Examining the senstivities shows that the consistency measure is highly
%sensitive to the upper bound in the 36th experiment and the lower bound in
%the 37th experiment

%There are several things we could do here.  We could adjust the
%uncertainty in these experiments and see if consistency is achieved.  But
%instead we see if removal of either of these two experiments changes the
%consistency.

%% Experiment removall

%start with a copy of the DCDataset
dsetCopy = dset;

%First let's remove the 36th experiment
dsetCopy.ModelAndObservationPair(36) = [];

%And check the consistency again
ctest2 = ConsistencyTest(dsetCopy)

%Looks like the consistency is undetermined.  Increasing the number of
%branch and bound iterations could help
opt = DCOptions('maxBranchBoundIter',2);
ctest3 = ConsistencyTest(ctest2,opt) %warm start

%So the dataset is still inconsistant even with experiment 36 removed.
%Therefore, experiment 36 does not seem to be the problem (at least not by
%itself).  Let's try removing experiment 37 instead
dsetCopy = dset;
dsetCopy.ModelAndObservationPair(37) = [];
ctest4 = ConsistencyTest(dsetCopy)

%The dataset is now consistent.  This implies that the model and
%experimental data from experiment 37 should be examined more closely.

%% Prediction

%Since we've removed experiment 37, it is actually a good candidate for a
%response prediction.  We will take only its ResponseModel and assume there
%is no data for a moment.

rm37 = dset.ModelAndObservationPair(37).ResponseModel
pred = ResponsePrediction(rm37,dsetCopy)

%The inner and outer bounds on the minimum are a little far apart, and more
%branch and bound iterations could tighten that gap, but we'll leave that
%for now.

%The data and uncertainty assert the following range for experiment 37
range37 = dset.ModelAndObservationPair(37).ResponseObservation.range

%Notice that the experimental range the range of values the model predicts
%over the feasible set constrained by other experiments do not over lap.
%This was the source of the inconsistency we saw earlier.  It is now up
%domain experts to determine why this might be.

%NOTE: The experiment uncertainties used in this example were not determined by the
%experimenters, but rather were ascertained as an educated guess by a
%domain expert.  Therefore, the results presented here should be taken as
%an example of the technique only, not as any assertion about the work they
%represent.