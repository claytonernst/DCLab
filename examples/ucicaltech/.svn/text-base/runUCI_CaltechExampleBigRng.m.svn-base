%ucicaltech.m

disp('This example has not been updated, and does not run properly.')
disp('Please try one of the other examples.')
return

%Free Parameters
paramNames = {'k1';'k2';'k3';'k4';'k5';'k6';'k7';'k8';'Gt'};
nom = [1e6; 1e-2; 4e-4; 4e-3; 4; 1e-5; 0.1; 1; 1e4];
%unc = 0.5;
%lb = 10.^(log10(nom)-unc);
%ub = 10.^(log10(nom)+unc);
%
%paramRanges = [...
%  0.8e6     1.2e6; ...
%  0.8e-2   1.2e-2; ...
%  3.2e-4   4.8e-4; ...
%  3.2e-3   4.8e-3; ...
%   3.8      4.2; ...
%  7e-6     13e-6; ...
%  0.07      0.13; ...
%   0.8      1.2; ...
%  0.8e4     1.2e4;]; 

n = 9;
P = cell(n,1);

unc.value = 0.5;
unc.type = 'absolute';
unc.transformation = 'log10';

for i1 = 1:n
  P{i1} = FreeParameter(paramNames{i1},nom(i1),unc);
end

FreeParam = vertcat(P{:});

inputLevel = {1;2;5;10;20;50;100;1000};
observations = [0.083 0.122 0.240 0.352 0.384 0.397 0.400 0.397];
unc = 0.1; %use 10% relative uncertainty

%define the ResponseModel and ResponseObservation objects
Pairs = cell(8,1);
for i1 = 1:8
  critRng = observations(i1)*[0.85 1.15];
  trnRng = [-inf inf];
  RM = ResponseModel(@G_protein,inputLevel{i1}*1e-9);
  RM = setPrivate(RM,'name',['BigH:_' num2str(inputLevel{i1}) 'nM']);
  RO = ResponseObservation(observations(i1),observations(i1)*unc);
  Pairs{i1,1} = ModelAndObservationPair(RO,RM,[num2str(inputLevel{i1}) 'nM'],critRng,trnRng);
end
Pairs = vertcat(Pairs{:});

%create dataset
D = Dataset(Pairs,FreeParam);
opt2 = DCOptions('display','all','fitConvergenceTol',0.1,'maxBranchBoundIter',3,'analysisMode','metamodelBasedA');
tic,
PD = PolyDataset(D,opt2);
toc

% % Create a metamodel
% dom(9,1) = struct;
% rng = getPrivate(FreeParam,'range');
% for i1 = 1:9
%   dom(i1).name = paramNames{i1};
%   dom(i1).range = rng(i1,:);
% end
% 
% meta = DCMetamodel(D.ModelAndObservationPair(5).ResponseModel,dom,[0.3 0.46],[0.2 0.6]);




