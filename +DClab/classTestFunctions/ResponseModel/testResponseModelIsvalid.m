function passed = testResponseModelIsvalid(dispBool)
% This file thrashes ResponseModel/isvalid if dispBool == true, messages
% about the progree will be displayed to the screen. If dispBool = false,
% this function will be "silent".

if nargin == 0
  dispBool = true;
end

passed = true;

%Get current warning state so we can revert to it at the end of this file.
warnState = warning('query','DCLAB:ResponseModel:isvalid');
if dispBool
  warning('on','DCLAB:ResponseModel:isvalid')
else
  warning('off','DCLAB:ResponseModel:isvalid')
end

%Create some good objects to break. Then test that isvalid detects that
%they are broken.
coeffMatrix = [3 4 5]';
domain = struct('name',{'p1';'p2'},'range',{[-1 1];[-inf inf]});
linObj = ResponseModel(coeffMatrix,domain); %create good object
quadObj = ResponseModel(coeffMatrix*coeffMatrix',domain);
dcObj = ResponseModel(@simpleDCModel);

if dispBool
  disp('===Testing ResponseModel isvalid method===')
end

if dispBool
  disp('  == Linear algebraic model case ==')
end

badObj = linObj;
badObj.model = coeffMatrix';
if isvalid(badObj)
  if dispBool
    disp('   ERROR: should be warned about wrong model vector orientation')
  end
  passed = false;
else
  if dispBool
    disp('   OK')
  end
end

badObj = linObj;
badObj.variableTransformations = {'log10';'log10'};
if isvalid(badObj)
  if dispBool
    disp('   ERROR: since lower limit of bounds is neg., should be warned about illegal log10 transformation')
  end
  passed = false;
else
  if dispBool
    disp('   OK')
  end
end

badObj = linObj;
badObj.domain = struct('name',{'p1';'p2';'p3'},'range',{[-1 1];[-inf inf];[-2 2]});
if isvalid(badObj)
  if dispBool
    disp('   ERROR: should be warned about too many components in domain')
  end
  passed = false;
else
  if dispBool
    disp('   OK')
  end
end

badObj = linObj;
badObj.domain = struct('name',{'p1';'p2'},'range',{[-1 1];1});
if isvalid(badObj)
  if dispBool
    disp('   ERROR: should be warned about lack of 1x2 vectors')
  end
  passed = false;
else
  if dispBool
    disp('   OK')
  end
end

badObj = linObj;
badObj.domain = struct('name',{'p1';'p2'},'range',{[-1 1];[2 1]});
if isvalid(badObj)
  if dispBool
    disp('   ERROR: should be warned about empty domain')
  end
  passed = false;
else
  if dispBool
    disp('   OK')
  end
end

badObj = linObj;
badObj.model = rand(3,2);
if isvalid(badObj)
  if dispBool
    disp('   ERROR: should be warned about incorrect dimensions for algebraic model')
  end
  passed = false;
else
  if dispBool
    disp('   OK')
  end
end

badObj = linObj;
badObj.domain = {5};
if isvalid(badObj)
  if dispBool
    disp('   ERROR: should be warned about nonstructure domain')
  end
  passed = false;
else
  if dispBool
    disp('   OK')
  end
end

badObj = linObj;
badObj.outputUncertainty = -1;
if isvalid(badObj)
  if dispBool
    disp('   ERROR: should be warned about negative outputUncertainty')
  end
  passed = false;
else
  if dispBool
    disp('   OK')
  end
end

badObj = setPrivate(linObj,'outputUncertainty',[1 2]);
if isvalid(badObj)
  if dispBool
    disp('   ERROR: should be warned about interval not containing zero')
  end
  passed = false;
else
  if dispBool
    disp('   OK')
  end
end

% It is not possible to set ouputUncertaintyType or
% outputUncertaintyTransformation to illegal values.

badObj = setPrivate(linObj,'guiDisplayCallback',3);
if isvalid(badObj)
  if dispBool
    disp('   ERROR: should be warned improper guiDisplayCallback')
  end
  passed = false;
else
  if dispBool
    disp('   OK')
  end
end

if dispBool
  disp('  == Quadratic algebraic model case ==')
end

badObj = setPrivate(quadObj,'model',rand(3));
if isvalid(badObj)
  if dispBool
    disp('   ERROR: should be warned about nonsymmetric quadratic matrix')
  end
  passed = false;
else
  if dispBool
    disp('   OK')
  end
end

badObj = setPrivate(quadObj,'domain',struct('name',{'p1';'p2';'p3'},'range',{[-1 1];[-inf inf];[-2 2]}));
if isvalid(badObj)
  if dispBool
    disp('   ERROR: should be warned about too many components in domain')
  end
  passed = false;
else
  if dispBool
    disp('   OK')
  end
end

if dispBool
  disp('  == DCModel case ==')
end

badObj = setPrivate(dcObj,'model',4);
if isvalid(badObj)
  if dispBool
    disp('   ERROR: should be warned that model is not a function_handle')
  end
  passed = false;
else
  if dispBool
    disp('   OK')
  end
end

% badObj = setPrivate(dcObj,'model',@badDCModel1);
% if isvalid(badObj)
%   if dispBool
%     disp('   ERROR: should be warned that model domain is not returned')
%   end
%   passed = false;
% else
%   if dispBool
%     disp('   OK')
%   end
% end
% 

try
  badObj = ResponseModel(@badDCModel2); %#ok
  if dispBool
    disp('   ERROR: call should have failed, domain in DCModel is not 1-by-2 vectors')
  end
catch
  if dispBool
    disp('   OK')
  end
end

%badDCModel3 has a negative scalar outputUncertainty
try
  badObj = ResponseModel(@badDCModel3); %#ok
  if dispBool
    disp('   ERROR: call should have failed, DCModel file has a negative scalar outputUncertainty')
  end
catch
  if dispBool
    disp('   OK')
  end
end

badObj = setPrivate(dcObj,'model',@complexDCModel);
badObj = setPrivate(badObj,'additionalInputs',{'resp1',{'a';'b'}});
if isvalid(badObj)
  if dispBool
    disp('   ERROR: should be warned about additional inputs not being convertible to one-line strings')
  end
  passed = false;
else
  if dispBool
    disp('   OK')
  end
end

badObj = setPrivate(dcObj,'model',@complexDCModel);
badObj = setPrivate(badObj,'additionalInputs',{'asdf'});
if isvalid(badObj)
  if dispBool
    disp('   ERROR: should be warned about first additional inputs not among response list')
  end
  passed = false;
else
  if dispBool
    disp('   OK')
  end
end

%badDCModel5 has multipleResponsesEnabled, but doesn't contain a feature
%list.
try 
  badObj = ResponseModel(@badDCModel5); %#ok
  if dispBool
    disp('   ERROR: call should have failed, DCModel file lacks response list')
  end
catch
  if dispBool
    disp('   OK')
  end
end

badObj = setPrivate(dcObj,'model',@badDCModel6);
badObj = setPrivate(badObj,'additionalInputs',{'resp1'});
if isvalid(badObj)
  if dispBool
    disp('   ERROR: should be warned about improper responseSimulation segment dimension')
  end
  passed = false;
else
  if dispBool
    disp('   OK')
  end
end

badObj = setPrivate(dcObj,'model',@badDCModelLocal);
badObj = setPrivate(badObj,'saveEnabled',true);

if isvalid(badObj)
  if dispBool
    disp('   ERROR: should be warned that model is not a function_handle to an mfile')
  end
  passed = false;
else
  if dispBool
    disp('   OK')
  end
end

warning(warnState); %revert warning state

function out = badDCModelLocal(flag,paramMatrix,varargin)

% this function works, but is used to test enabling save with a local
% function

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
 otherwise
  error(['Behavior for flag value ' flag ' is not defined']);
end
