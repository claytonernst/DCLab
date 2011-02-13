function out = lessSimpleDCModel(flag,paramMatrix,varargin)
% Less Simple example/test dcModel mfile having the required
% functionality and a few extras

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
  out = 0.1; % Define a bound on the modeling error. Be sure
             % to account for all cases presented by varargin.

 case 'getOutputUncertaintyType'   
  out = 'relative';

 case 'getGUIDispayCallback' 
  out = @disp; % Enter the name of a function on the path that accepts a
           % ResponseModel object as its input. This function will
           % be called if a GUI application requests info about a 
           % ResponseModel object based on this file.

 case 'isSaveEnabled'
  out = false; % Enter a value of 0 or 1. If 1, model evaluations invoked
               % by the 'simulate' flag will be saved and any inputs in 
               % varargin must be single row character arrays or scalars.

 otherwise
  error(['Behavior for flag value ' flag ' is not defined']);
end

%===Add subfunctions here===
