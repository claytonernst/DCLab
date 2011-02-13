function out = dcModel(flag,paramMatrix,varargin)
% A prototype dcModel mfile. The file can be found at
% ../DClabV1p1/fileTemplates/dcModel.m
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
%           relative, log10: (1+outputUnc(1))*log10(M(x)) <= log10(y(x)) <=
%           (1+outputUnc(2))*log10(M(x))
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

      % Add code to simulate your model at paramVect to product an output y.
      % Any necessary subfunctions may be added to the end of this file. 
      out(i1) = y;
    end
  
  case 'getModelDomain'    
    n = ?? % The number of model parameters
    % Initialize the nx1 model domain structure array with fields .name and
    % .range
    modelDomain = struct('name',cell(n,1),'range',cell(n,1));
    modelDomain(1).name =  ?? % Name of first parameter, e.g., 'param1'
    modelDomain(1).range = ?? % Range of first parameter, e.g., [1 10] or [-inf inf]
    modelDomain(2).name = ??
    modelDomain(2).range = ??
    % ... fill in commands to define the name and range of all n parameters
    out = modelDomain;
    % the modelDomain may be dependent on the additional inputs as they are
    % provided when this flag is called.

    
%---getOutputUncertainty is an optional condition, uncomment if required---
%  case 'getOutputUncertainty'   
%    out = ?? % Define a bound on the modeling error y(x) - M(x) where M(x)
%             % is the model defined in this file and y(x) is the output of
%             % the actual model. Be sure to account for all cases presented
%             % by varargin.

%---getOutputUncertaintyType is an optional condition (default is
%    'absolute') that is only relevant (but not required) when an output
%    uncertainty is provided, uncomment if required---
%  case 'getOutputUncertaintyType'   
%    out = ?? % Define whether the provided value for output uncertainty is
%             % 'absolute' or 'relative'.

%---getOutputUncertaintyTransformation is an optional condition (default is
%    'none') that is only relevant (but not required) when an output
%    uncertainty is provided, uncomment if required---
%  case 'getOutputUncertaintyTransformation'   
%    out = ?? % Define a transformation into the coordinates in which the
%             % output uncertainty is expresses. Either 'none' or 'log10'

%---isSaveEnabled is an optional condition, uncomment if required---
%  case 'isSaveEnabled'
%    out = ?? % Enter a value of 0 or 1. If 1, model evaluations invoked
%             % by the 'simulate' flag will be saved and any inputs in 
%             % varargin must be single row character arrays or scalars.

%---isMultipleResponsesEnabled is an optional condition, uncomment if required---
%  case 'isMultipleResponsesEnabled'
%    out = ?? % Enter a value of 0 or 1. If 1, the output in the 'simulate'
%             % case should be Nresp x size(paramMatrix,2). Additionally a
%             % responseList must be defined in this file and the first element
%             % of varargin must be a member of responseList.

%---getResponseList is a required condition when isMultipleResponsesEnabled == 1 ---
%  case 'getResponseList'
%    out = ?? % Enter a column cell array of single line chars. The first
%             % element of varargin must be an element of this cell array.

%---getResponseSimulationSegment an optional condition that is only
%   relevant (but not required) when isMultipleResponsesEnabled == 1 ---
%  case 'getResponseSimulationSegment'
%    out = ?? % Enter a column vector of scalars. For cases when the
%             % simulation can be broken up into multiple segments. Each
%             % element of this vector indicates the segment which computes
%             % the corresponding element of responseList. This feature is
%             % used to get as much information as possible from a single
%             % simulation run, but no more than needed. Only responses
%             % with segment number less than or equal to that of varargin{1}
%             % should be produced in the output of the 'simulate' case. At
%             % that point the simulation can be terminated.

%---getName is an optional condition, uncomment if require---
%     case 'getName'
%         out = ?? % If the model can be different depending on the additional
%                  % inputs then this let's you define a character array
%                  % (string) that will be used to name files when saving
%                  % points.  Note: if isMultipleResponseEnabled is true (see
%                  % below) this should NOT reflect a specific feature.
%                  % Feature names are handled automatically.  This name can
%                  % depend on the other input arguments to this function as
%                  % they are passed in when getting this string.

  otherwise
    error(['Behavior for flag value ' flag ' is not defined']);
end

%===Add subfunctions here===
