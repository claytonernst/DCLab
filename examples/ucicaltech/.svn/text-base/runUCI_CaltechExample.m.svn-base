%ucicaltech.m

disp('This example has not been updated, and does not run properly.')
disp('Please try one of the other examples.')
return

%Free Parameters
paramNames = {'k1';'k2';'k3';'k4';'k5';'k6';'k7';'k8';'Gt'};
paramRanges = [...
  0.8e6     1.2e6; ...
  0.8e-2   1.2e-2; ...
  3.2e-4   4.8e-4; ...
  3.2e-3   4.8e-3; ...
   3.8      4.2; ...
  7e-6     13e-6; ...
  0.07      0.13; ...
   0.8      1.2; ...
  0.8e4     1.2e4;]; 

n = 9;
P = cell(n,1);

for i1 = 1:n
  P{i1} = FreeParameter(paramNames{i1},mean(paramRanges(i1,:)),paramRanges(i1,2)-mean(paramRanges(i1,:)));
end
FreeParam = vertcat(P{:});

inputLevel = {1;2;5;10;20;50;100;1000};
observations = [0.083 0.122 0.240 0.352 0.384 0.397 0.400 0.397];
unc = 0.1; %use 10% relative uncertainty

%define the ResponseModel and ResponseObservation objects
Pairs = cell(8,1);
for i1 = 1:8
  RM = ResponseModel(@G_protein,inputLevel{i1}*1e-9);
  RO = ResponseObservation(observations(i1),observations(i1)*unc);
  Pairs{i1,1} = ModelAndObservationPair(RO,RM,[num2str(inputLevel{i1}) 'nM']);
end
Pairs = vertcat(Pairs{:});

%create dataset
D = DCDataset(Pairs,FreeParam);
opt2 = DCOptions('display','all','fitConvergenceTol',0.1,'maxBranchBoundIter',3,'analysisMode','metamodelBasedA');

% Create a metamodel
dom(9,1) = struct;
for i1 = 1:9
  dom(i1).name = paramNames{i1};
  dom(i1).range = paramRanges(i1,:);
end

meta = DCMetamodel(D.ModelAndObservationPair(5).ResponseModel,dom,[0.3 0.46],[0.2 0.6]);

tic,
PD = PolyDataset(D,opt2);
toc

%Small
scrambledFreeParam = vertcat(P{3},P{5},P{1},P{2},P{9},P{7},P{8},P{4},P{6});
Dsmall = Dataset(ModelAndObservationPair(RO,RM,[num2str(inputLevel{8}) 'nM']),scrambledFreeParam);
PDsmall = PolyDataset(Dsmall,opt2);

Dsmall2 = Dataset(ModelAndObservationPair(RO,RM,[num2str(inputLevel{8}) 'nM']),FreeParam);
PDsmall2 = PolyDataset(Dsmall2,opt2);

%%Stopped here



%create options structure
opt = DCOptions('display','all','maxBranchBoundIter',1,'trans',{'logXlogY'});

%test the consistency
[lb ub obj] = ConsistTest(D,opt);

%view the getable test object properties
get(obj)
