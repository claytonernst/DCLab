function out = badDCModel1(flag,paramMatrix,varargin)
% badDCModel1 does not define the model's domain

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

 otherwise
  error(['Behavior for flag value ' flag ' is not defined']);
end

%===Add subfunctions here===
