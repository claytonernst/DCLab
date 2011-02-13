function out = badDCModel6(flag,paramMatrix,varargin)
% badDCModel6 defines a responseSimulateSegment that is inconsistent
% with the dimensions of responseList

if ~ischar(flag)
  error('Inputs: The first input to a dcModel must be a char array');
end

switch flag
 case 'simulate'    % Simulate the model to get a prediction for Y
  % Determine how many horizontally concatenated parameter vectors were
  % supplied
  N = size(paramMatrix,2);
  out = zeros(1,N);
  for i1 = 1:N
    % Define the parameter values
    paramVect = paramMatrix(:,i1);

    % Add code to simulate your model at paramVect to product an output y.
    % Any necessary subfunctions may be added to the end of this file. 
    out(i1) = sum(paramVect);
  end
  
 case 'getModelDomain'    
  n = 2; % The number of model parameters
  % Initialize the nx1 model domain structure array with fields .name and
  % .range
  modelDomain = struct('name',cell(n,1),'range',cell(n,1));
  modelDomain(1).name =  'p1'; % Name of first parameter, e.g., 'param1'
  modelDomain(1).range = [-1 1]; % Range of first parameter, e.g., [1 10] or [-inf inf]
  modelDomain(2).name = 'p2';
  modelDomain(2).range = [-inf inf];
  % ... fill in commands to define the name and range of all n parameters
  out = modelDomain;

 case 'getOutputUncertainty'   
  out = [-4 2]; % Define a bound on the modeling error. Be sure
%           % to account for all cases presented by varargin.

 case 'getOutputUncertaintyType'   
  out = 'absolute'; 

 case 'getOutputUncertaintyTransformation'   
  out = 'log10'; 

 case 'isMultipleResponsesEnabled'
  out = true; % Enter a value of 0 or 1. If 1, the output in the 'simulate'
%           % case should be Nresp x size(paramMatrix,2). Additionally a
%           % responseList must be defined in this file and the first element
%           % of varargin must be a member of responseList.


 case 'getResponseList'
  out = {'resp1';'resp2'}; % Enter a column cell array of single line chars. The first
%           % element of varargin must be an element of this cell array.

%---getResponseSimulationSegment an optional condition that is only
%   relevant (but not required) when isMultipleResponsesEnabled == 1 ---
 case 'getResponseSimulationSegment'
  out = 1; % Enter a column vector of scalars. For cases when the
%           % simulation can be broken up into multiple segments. Each
%           % element of this vector indicates the segment which computes
%           % the corresponding element of responseList. This feature is
%           used to get as much information as possible from a single
%           simulation run, but no more than needed. Only responses
%           % with segment number less than or equal to that of varargin{1}
%           % should be produced in the output of the 'simulate' case. At
%           % that point the simulation can be terminated.

 otherwise
  error(['Behavior for flag value ' flag ' is not defined']);
end

%===Add subfunctions here===
