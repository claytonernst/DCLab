function out = model2(flag,paramMatrix)
% Inputs:
%   flag: A char array indicating the task to be performed. Valid
%         values are 'simulate', 'getModelDomain', 'isSaveEnabled', 
%         'getModelErr', or 'getDispCallback'. The output returned for
%         each of these flag values is described below.
%   paramMatrix: a nxNevals matrix, where each column contains a value
%                of the n dimensional parameter vector. 
%
% Outputs:
%   if flag=='simulate', out is an 1xNevals vector containing the
%             model predictions for the quantity of interest at each 
%             supplied value of the parameter vector.
%   if flag=='getModelDomain', out is an nx1 structure array with fields
%             .name and .range.

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

    % implement z1_dot = x1*exp(-x2*z1), z1(0) = 0
    odefun = @(t,z) paramVect(1)-paramVect(2)*z;
    [tout zout] = ode45(odefun,[0 5],0);
    y = mean(zout);

    out(i1) = y;
  end
  
 case 'getModelDomain'    
  n = 2; % The number of model parameters
  % Initialize the nx1 model domain structure array with fields .name and
  % .range
  modelDomain = struct('name',cell(n,1),'range',cell(n,1));
  modelDomain(1).name =  'X1'; % Name of first parameter, e.g., 'param1'
  modelDomain(1).range = [-inf inf]; % Range of first parameter, e.g., [1 10] or [-inf inf]
  modelDomain(2).name = 'X3';
  modelDomain(2).range = [-inf inf];
  
  out = modelDomain;

 otherwise
  error(['Behavior for flag value ' flag ' is not defined']);
end

