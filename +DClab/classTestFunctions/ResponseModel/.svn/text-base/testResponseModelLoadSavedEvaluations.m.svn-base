function passed = testResponseModelLoadSavedEvaluations(dispBool)
%TESTRESPONSEMODELLOADSAVEDEVALUATIONS returns true if it appears OK
%
%   PASSED = TESTRESPONSEMODELLOADSAVEDEVALUTATIONS returns true or false.
%   Several messages will be displayed to the screen.
%
%   PASSED = TESTRESPONSEMODELLOADSAVEDEVALUATIONS(DISPBOOL) will suppress
%   all screen displays if DISPBOOL==false.

if nargin == 0
  dispBool = true;
end

passed = true;

if dispBool
  disp('===Testing ResponseModel loadSavedEvaluations method===')
end

coeffMatrix = [0 1]';
domain = struct('name',{'p1'},'range',{[-1 1]});
linRM = ResponseModel(coeffMatrix,domain);

if dispBool
  disp('    1 var linear model')
end
[x y] = loadSavedEvaluations(linRM,[-1 1]);
if ~isempty(y)
  passed = false;
  if dispBool
    disp('   ERROR: computation incorrect, no (x,y) pairs should be found')
  end
end

%Create some saved points to load
simpleRM = ResponseModel(@simpleDCModel);
simpleRM = set(simpleRM,'saveEnabled',true');

%Delete any existing points
deleteSavedEvaluations(simpleRM);

%Domain is simpleRM is [-1 1;-inf inf]
randomEvalWithSave(simpleRM,[-1 1; 0.1 1],{'none';'log10'},100); 

if dispBool
  disp('    Loading points from simpleDCModel')
end
%Try to load them
[x y] = loadSavedEvaluations(simpleRM,[-1 1; 0.1 1],{'none';'log10'},100);
if length(y)~=100
  passed = false;
  if dispBool
    disp('   ERROR: computation incorrect, 100 (x,y) pairs should be found')
  end
end

%Load points from a slightly smaller domain
[x y] = loadSavedEvaluations(simpleRM,[-0.5 0.5; 0.1 1],{'none';'log10'},100);
if isempty(y)
  passed = false;
  if dispBool
    disp('   ERROR: computation incorrect, some (x,y) pairs should be found')
  end
end

%Create some more points with different trans/domain to load
randomEvalWithSave(simpleRM,[0.5 0.75; 5 10],{'log10';'log10'},100);
%Load points with transformations that should be unavailable
[x y] = loadSavedEvaluations(simpleRM,[0.5 0.75; 5 10],{'none';'log10'},100);
if ~isempty(y)
  passed = false;
  if dispBool
    disp('   ERROR: computation incorrect, no (x,y) pairs should be found')
  end
end

%Load points disregaring the transformations.
[x y] = loadSavedEvaluations(simpleRM,[0.5 0.75; 5 10],[],100);
if length(y)~=100
  passed = false;
  if dispBool
    disp('   ERROR: computation incorrect, 100 (x,y) pairs should be found')
  end
end

%Load points for which the proper transformations should be available.
[x y] = loadSavedEvaluations(simpleRM,[0.5 0.75; 5 10],{'log10';'log10'},100);
if length(y)~=100
  passed = false;
  if dispBool
    disp('   ERROR: computation incorrect, 100 (x,y) pairs should be found')
  end
end

%Create some more points with one singleton dimensions
randomEvalWithSave(simpleRM,[0.1 0.1; 5 10],{'log10';'none'},100);
%Try to load a subset of these
[x y] = loadSavedEvaluations(simpleRM,[0.1 0.1; 6 9],{'log10';'none'});
if isempty(y)
  passed = false;
  if dispBool
    disp('   ERROR: computation incorrect, some (x,y) pairs should be found')
  end
end

%Delete all created points
deleteSavedEvaluations(simpleRM);

%Create some saved points to load
cplxRM = ResponseModel(@complexDCModel,'resp3');
cplxRM = set(cplxRM,'saveEnabled',true');

%Delete any existing points
deleteSavedEvaluations(cplxRM);

%Domain of cplxRM is [-1 1;-inf inf]
% This will create points for resp1, resp2, and resp3
randomEvalWithSave(cplxRM,[-1 1; 0.1 1],{'none';'log10'},100); 

if dispBool
  disp('    Loading points from cplxDCModel_resp1')
end
%See if we can load points for the first response.
cplxRM = ResponseModel(@complexDCModel,'resp1');
[x y] = loadSavedEvaluations(cplxRM,[-0.9 1; 0.1 1],{'none';'log10'});
if isempty(y) || max(x(1,:)) > 1 || min(x(1,:)) < -0.9 || max(x(2,:)) > 1 || min(x(2,:)) < 0.1 || ~all(abs(sum(x)-y) < 1e-10)
  passed = false;
  if dispBool
    disp('   ERROR: computation incorrect, some (x,y) pairs should be found')
  end
end
deleteSavedEvaluations(cplxRM);

if dispBool
  disp('    Loading points from cplxDCModel_resp2')
end
%See if we can load points for the second response.
cplxRM = ResponseModel(@complexDCModel,'resp2');
[x y] = loadSavedEvaluations(cplxRM,[-0.9 1; 0.1 1],{'none';'log10'});
if isempty(y) || max(x(1,:)) > 1 || min(x(1,:)) < -0.9 || max(x(2,:)) > 1 || min(x(2,:)) < 0.1 || ~all(abs(exp(sum(x))-y) < 1e-10)
  passed = false;
  if dispBool
    disp('   ERROR: computation incorrect, some (x,y) pairs should be found')
  end
end
deleteSavedEvaluations(cplxRM);

if dispBool
  disp('    Loading points from cplxDCModel_resp3')
end
%See if we can load points for the third response.
cplxRM = ResponseModel(@complexDCModel,'resp3');
[x y] = loadSavedEvaluations(cplxRM,[-0.9 1; 0.1 1],{'none';'log10'},10); %skip file 10
if isempty(y) || max(x(1,:)) > 1 || min(x(1,:)) < -0.9 || max(x(2,:)) > 1 || min(x(2,:)) < 0.1 || ~all(abs(10.^(sum(x))-y) < 1e-10)
  passed = false;
  if dispBool
    disp('   ERROR: computation incorrect, some (x,y) pairs should be found')
  end
end
deleteSavedEvaluations(cplxRM);




