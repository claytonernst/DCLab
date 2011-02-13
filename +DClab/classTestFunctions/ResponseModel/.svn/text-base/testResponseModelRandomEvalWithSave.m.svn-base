function passed = testResponseModelRandomEvalWithSave(dispBool)
%TESTRESPONSEMODELRANDOMEVALWITHSAVE returns true if it appears OK
%
%   PASSED = TESTRESPONSEMODELRANDOMEVALWITHSAVE returns true or false.
%   Several messages will be displayed to the screen.
%
%   PASSED = TESTRESPONSEMODELRANDOMEVALWITHSAVE(DISPBOOL) will suppress
%   all screen displays if DISPBOOL==false.

if nargin == 0
  dispBool = true;
end

passed = true;

if dispBool
  disp('===Testing ResponseModel randomEvalWithSave method===')
end

coeffMatrix = [0 1]';
domain = struct('name',{'p1'},'range',{[-1 1]});
linRM = ResponseModel(coeffMatrix,domain);

if dispBool
  disp('    1 var linear model, no transformaion on x')
end
[x y] = randomEvalWithSave(linRM,[-1 1],{'none'},10);
if ~isequal(x,y)
  passed = false;
  if dispBool
    disp('   ERROR: computation incorrect')
  end
end

if dispBool
  disp('    1 var linear model, log10 transformaion on x')
end
[x y] = randomEvalWithSave(linRM,[0.1 1],{'log10'},100);
if ~isequal(x,y) || min(x) < 0.1 || max(x) > 1
  passed = false;
  if dispBool
    disp('   ERROR: computation incorrect')
  end
end

if dispBool
  disp('    1 var linear model, too large of range on x')
end
try %#ok
  [x y] = randomEvalWithSave(linRM,[-2 1],{'none'},100); %#ok
  passed = false;
  if dispBool
    disp('   ERROR: previous call should have bombed')
  end
end

if dispBool
  disp('    1 var linear model, illegal log10 transformaion on x')
end
try %#ok
  [x y] = randomEvalWithSave(linRM,[-1 1],{'log10'},100); %#ok
  passed = false;
  if dispBool
    disp('   ERROR: previous call should have bombed')
  end
end

if dispBool
  disp('    2 var linear model, log10 transformaion on x, one singleton dimension')
end
coeffMatrix = [0 1 2]';
domain = struct('name',{'p1';'p2'},'range',{[-1 1];[0 1]});
linRM = ResponseModel(coeffMatrix,domain);
[x y] = randomEvalWithSave(linRM,[0 0; 0.1 1],{'none';'log10'},100); 
if max(x(1,:)) > 0 || min(x(1,:)) < 0 || max(x(2,:)) > 1 || min(x(2,:)) < 0.1 || ~isequal(2*x(2,:),y)
  passed = false;
  if dispBool
    disp('   ERROR: computation incorrect')
  end
end

if dispBool
  disp('    simple dcmodel')
end
simpleRM = ResponseModel(@simpleDCModel);
[x y] = randomEvalWithSave(simpleRM,[0 0; 0.1 1],{'none';'log10'},100); 
if max(x(1,:)) > 0 || min(x(1,:)) < 0 || max(x(2,:)) > 1 || min(x(2,:)) < 0.1 || ~isequal(10.^sum(x),y)
  passed = false;
  if dispBool
    disp('   ERROR: computation incorrect')
  end
end

if dispBool
  disp('    simple dcmodel with save enabled')
end
simpleRM = set(simpleRM,'saveEnabled',true');
% Delete any saved evaluations
deleteSavedEvaluations(simpleRM);
[x y] = randomEvalWithSave(simpleRM,[0 0; 0.1 1],{'none';'log10'},100); %#ok

% Check for new saved evaluations
file = which('savedEvaluationsDir');
evalPath = fileparts(file);
dirName = fullfile(evalPath,'simpleDCModel');
if ~exist(dirName,'dir') || ~exist([dirName filesep 'simpleDCModel.mat'],'file') || ~exist([dirName filesep '1.mat'],'file')
  passed = false;
  if dispBool
    disp('   ERROR: computation incorrect')
  end
else
  %load them up for fun
  s = load([dirName filesep '1.mat']);
  x = s.a.paramValues;
  y = s.a.responses;
  xmin = min(x,[],2);
  xmax = max(x,[],2);
  if xmax(1) > 0 || xmin(1) < 0 || xmax(2) > 1 || xmin(2) < 0.1 || ~isequal(10.^sum(x),y)
    passed = false;
    if dispBool
      disp('   ERROR: computation incorrect')
    end
  end
end
% Delete the just created evaluations
deleteSavedEvaluations(simpleRM);

%==============
if dispBool
  disp('    Not testing multiple computers, uncomment to do so')
end
% 
% %Note: for these next few lines to work, you must have
% %makeComputerDCSlave(1) running on another machine that shares the same
% %file server as this machine. Alternatively you can call
% %makeComputerDCSlave(1) using another instance of MATLAB on this machine.
% %Additionally, ../classTestFunctions/ResponseModel/simpleDCModel must be on
% %the path of the machine running makeComputerDCSlave.
% if dispBool
%   disp('    simple dcmodel with multiple computers')
% end
% 
% simpleRM = ResponseModel(@simpleDCModel);
% [x y] = randomEvalWithSave(simpleRM,[0 0; 0.1 1],{'none';'log10'},100,1); 
% if max(x(1,:)) > 0 || min(x(1,:)) < 0 || max(x(2,:)) > 1 || min(x(2,:)) < 0.1 || ~isequal(10.^sum(x),y)
%   passed = false;
%   if dispBool
%     disp('   ERROR: computation incorrect')
%   end
% end
% %Clean up by removing the jobQueue file
% file = [fileparts(which('multiComputerTempDir')) filesep 'jobQueueObj.mat'];
% if exist(file,'file')
%   delete(file,'s');
% end
% 

%=============

if dispBool
  disp('    complex dcmodel with save enabled')
end
cplxRM = ResponseModel(@complexDCModel,'resp2');
cplxRM = set(cplxRM,'saveEnabled',true');
[x y] = randomEvalWithSave(cplxRM,[0 1; 0.1 1],{'none';'log10'},100); 
if max(x(1,:)) > 1 || min(x(1,:)) < 0 || max(x(2,:)) > 1 || min(x(2,:)) < 0.1 || ~all(abs(exp(sum(x))-y) < 1e-10)
  passed = false;
  if dispBool
    disp('   ERROR: computation incorrect')
  end
end

% Delete the just created evaluations
deleteSavedEvaluations(cplxRM);
cplxRM = ResponseModel(@complexDCModel,'resp1');
deleteSavedEvaluations(cplxRM);

