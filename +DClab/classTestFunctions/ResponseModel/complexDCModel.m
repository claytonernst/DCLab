function out = complexDCModel(flag,paramMatrix,varargin)
% A fairly complecated dcModel mfile.
%
% Inputs:
%   flag: A char array indicating the task to be performed. Valid
%         values are 'simulate', 'getModelDomain', 'getOutputUncertainty',
%         'getOutputUncertaintyType', 'getOutputUncertaintyTransformation', 
%         'getDispCallback', 'isSaveEnabled', 'isMultipleResponsesEnabled',
%         'getResponseList', 'getResponseSimulationSegment'. This function
%         MUST support 'simulate' and 'getModelDomain'. The output returned
%         for each of these flag values is described below.
%   paramMatrix: a nxNevals matrix, where each column contains a value
%         of the n dimensional parameter vector. 
%   varargin: any additional inputs given to the ResponseModel constructor
%
% Outputs: for possible flag values
%   'simulate', out is an 1xNevals vector containing the model predictions
%       for the response of interest at each supplied value of the parameter
%       vector. In the special case where multiple responses can be
%       obtained from a single simulation, out can be nRespxNevals. Further
%       details are provided below.
%   'getModelDomain', out is an nx1 structure array with fields .name and .range.
%   'getOutputUncertainty', out is a nonnegative scalar indicating an upper bound on
%       the modeling error. Alternatively, if the modeling error is
%       asymmetric,  out can be a 1x2 vector with a negative and positive
%       component.
%   'getOutputUncertaintyType', out is either 'absolute' or 'relative'.
%       If relative, the uncertainty should be a relative fraction, i.e.,
%       0.1 rather than 10(%).
%   'getOutputUncertaintyTransformation', out is either 'none' or 'log10'.
%       The letting y(x) be the output of the true model at x, The
%       expressions for the bounds on outputUncertainty in each possible
%       case are given below:
%         symmetric case (uncertainty value = a scalar)
%           absolute, no trans: |y(x) - M(x)| <= outputUnc
%           absolute, log10: |log10(y(x)) - log10(M(x))| <= outputUnc
%           relative, no trans: |y(x) - M(x)| <= outputUnc*M(x)
%           relative, log10: |log10(y(x)) - log10(M(x))| <= outputUnc*log10(M(x))
%         asymmetric case (uncertainty value = a 1x2, 1st component negative)
%           absolute, no trans: outputUnc(1) <= y(x) - M(x) <= outputUnc(2)
%           absolute, log10: outputUnc(1) <= log10(y(x)) - log10(M(x)) <= outputUnc(2)
%           relative, no trans: (1+outputUnc(1))*M(x) <= y(x) <= (1+outputUnc(2))*M(x)
%           relative, log10: (1+outputUnc(1))*log10(M(x)) <= log10(y(x)) <= (1+outputUnc(2))*log10(M(x))
%   'getGUIDisplayCallback', out is a function_handle referencing a
%       function on the MATLAB path. This function is executed when a GUI
%       application requests more information about a ResponseModel
%       created from this dcModel file. 
%   'isSaveEnabled', out is either 0 or 1, with 1 indicating the DC
%       software should save evaluations of the model for possible future
%       use. This functionality requires that any arguments passed in
%       varargin be single row character arrays
%   'isMultipleResponsesEnabled', out is a boolean. When true, the output for the
%       'simulate' case will be multidimensional (fat) matrix containing an output for
%       all responses produced by the simulate run. Additionally this file
%       must define a list of these responses, and the first additional
%       input to the file (varargin{1}) must be the name of a response on
%       this list. It is possible (likely?) that the code contained in this
%       file will not use this first additional input.
%   'getResponseList', out will be a column cell array of strings. Calling
%       syntax with this flag must be supported if isMultRespEnabled ==true.
%   'getResponseSimulationSegment', out is a vector of positive integers of
%       the same size as responseList. This is used when multiple responses
%       are enabled and the simulation can be broken up into segments,
%       computing some but not all of the responses.

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

    % Define the output for each response, taking into account the values
    % of responseSimulationSegment
    switch varargin{1}
      case {'resp1','resp2'} %both are in segment 1
        out(1,i1) = exp(sum(paramVect));
        out(2,i1) = 5*exp(sum(paramVect));
      case 'resp3'; %this is in segment 2
        out(1,i1) = exp(sum(paramVect));
        out(2,i1) = 5*exp(sum(paramVect));
        out(3,i1) = 10^sum(paramVect);
      otherwise
        error('First element of varargin not a member of responseList')
    end
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
  out = 0.1; % Define a bound on the modeling error y(x) - M(x) where M(x)
           % is the model defined in this file and y(x) is the output of the actual
           % model. Be sure to account for all cases presented by varargin.

 case 'getOutputUncertaintyType'   
  out = 'relative'; % Define whether the provided value for output uncertainty is
           % 'absolute' or 'relative'.

 case 'getOutputUncertaintyTransformation'   
  out = 'none'; % Define a transformation into the coordinates in which the
           % output uncertainty is expresses. Either 'none' or 'log10'

 case 'getGUIDisplayCallback' 
  out = @disp; % Enter a handle to a function that accepts a
           % ResponseModel object as its input. This function will
           % be called if a GUI application requests info about a 
           % ResponseModel object based on this file.

 case 'isSaveEnabled'
  out = true; % Enter a value of 0 or 1. If 1, model evaluations invoked
           % by the 'simulate' flag will be saved and any inputs in 
           % varargin must be single row character arrays or scalars.

 case 'isMultipleResponsesEnabled'
  out = true'; % Enter a value of 0 or 1. If 1, the output in the 'simulate'
           % case should be Nresp x size(paramMatrix,2). Additionally a
           % responseList must be defined in this file and the first element
           % of varargin must be a member of responseList.

 case 'getResponseList'
  out = {'resp1';'resp2';'resp3'}; % Enter a column cell array of single line chars. The first
           % element of varargin must be an element of this cell array.

 case 'getResponseSimulationSegment'
  out = [1;1;2]; % Enter a column vector of scalars. For cases when the
           % simulation can be broken up into multiple segments. Each
           % element of this vector indicates the segment which computes
           % the corresponding element of responseList. This feature is
           % used to get as much information as possible from a single
           % simulation run, but no more than needed. Only responses
           % with segment number less than or equal to that of varargin{1}
           % should be produced in the output of the 'simulate' case. At
           % that point the simulation can be terminated.

 otherwise
  error(['Behavior for flag value ' flag ' is not defined']);
end

%===Add subfunctions here===
