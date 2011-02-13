% This is a brief template to check the consistency of a dataset

% Create the model assertion who's range you want to predict. The
% typical constructor syntax is 
%   RM = ResponseModel(@dcModel,addlInput1,...,)
%RM0 = ...


% Define the remaining model assertions, one for each experiment
% assertion.
%RM1 = ...
%RMm = ...

% Define the experiment assertions. The typical constructor syntax
% is    
%       RO = ResponseObservation(data,unc)  
% where unc is a positive scalar. See help for asymmetric
% uncertainty bounds.
%
%RO1 = ...
%ROm = ...

% Define the parameter assertion. The typical constructor syntax is
%    FP = FreeParameter('name',nominal,uncertainty)
% where range is 1x2
%
%FP1 = ...
%FPn = ...

%Define the dataset. First the dataset units, with the typical
%constructor syntax
%    Pair = ModelAndObservationPair(RO,RM,'name')
%
%Pair1 = ...
%Pairm = ...
%
%then the dataset.
%Dset = Dataset([Pair1;...;Pairm],[FP1;...;FPn])

%Define any DCOptions. It is often useful to set the maximum number
%of branch and bound iterations to 1 initially to gauge how long the total
%execution will be. Additionally you can set the display level
%here. Valid values are ['off'|'notify'|'iter'|'all'|'ALL']
opt = DCOptions('maxBranchBoundIter',1,'display','iter')

%Call the Prediction constructor
PredObj = ResponsePrediction(MA0,Dset,opt);

%This object has lots of properties available to explore.
