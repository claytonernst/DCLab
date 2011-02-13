function MM = reevalerr(MM,RM,domain,criticalRange,opt)

domrng = vertcat(domain.range);
na = length(MM.variableTransformations~=0);

Nv = 50*(na+1)*(na+2)/2;

n = length(MM.variableTransformations);
varTransChar = cell(n,1);
for i1 = 1:n
  if MM.variableTransformations(i1)==0
    varTransChar{i1} = 'log10';
  else
    varTransChar{i1} = 'none';
  end
end
respTrans = MM.responseTransformation;
varTransNum = MM.variableTransformations;

[xV yV] = loadSavedEvaluations(RM,domrng,varTransChar,Nv);

Nsamples = length(yV);

if Nsamples < Nv
  [newXV, newYV] = randomEvalWithSave(RM,domrng,varTransChar,Nv-Nsamples,opt.nComputer);
  xV = [xV newXV];
  yV = [yV newYV];
end

% Reject data that didn't lie in the critical range.
if ~any(isinf(criticalRange)) 
  junkDataV = isnan(yV) | yV > criticalRange(2) | yV < criticalRange(1);
  yV = yV(~junkDataV);
  xV = xV(:,~junkDataV);
elseif isinf(criticalRange(1))
  junkDataV = isnan(yV) | yV > criticalRange(2);
  yV = yV(~junkDataV);
  xV = xV(:,~junkDataV);
elseif isinf(criticalRange(2))
  junkDataV = isnan(yV) | yV < criticalRange(1);
  yV = yV(~junkDataV);
  xV = xV(:,~junkDataV);
else
  junkDataV = isnan(yV);
  yV = yV(~junkDataV);
  xV = xV(:,~junkDataV);
end

if length(yV) < 0.05*Nv
  disp('More than 95 %% of the validation data was rejected.')
  keyboard
end

% Determine the error on the validation points.
yHatV = evalPrivate(MM.model,xV,respTrans,varTransNum);
errV = yHatV-yV;

fitInfo.nSamplePoints4Validation = Nv;
fitInfo.averageErrorOnSamplePoints4Validation = mean(abs(errV));
fitInfo.peakErrorOnSamplePoints4Validation = [min(errV) max(errV)];

% Optimizes the fitting error.

xAll = [xV];
errAll = [errV];

%if any(isnan(yV)) && ( isempty(criticalRange) || isempty(trainingRange) )
%  error('Response Model give NaN outputs, but either criticalRange or trainingRange was undefined.')
%end

% Reject data that didn't lie in the critical range.
%if ~isempty(criticalRange)
%  junkData = isnan(yV) | yV > criticalRange(2) | yV < criticalRange(1);
%  yV = yV(~junkData);
%  xV = xV(:,~junkData);
%  if isempty(yV)
%    disp('No validation data in the critical range')
%    keyboard
%  end
%end
%   
% yHatV = evalPrivate(model.model,xV,respTrans,varTransNum);
% errV = yHatV-yV;
% 
% fitInfo.nSamplePoints4Validation = Nv;
% fitInfo.averageErrorOnSamplePoints4Validation = mean(abs(errV));
% fitInfo.peakErrorOnSamplePoints4Validation = [min(errV) max(errV)];
% 
% % Optimize fitting error:
% 
% % Reject data that didn't lie in the critical range.
% if ~isempty(criticalRange)
%   junkData = isnan(yT) | yT > criticalRange(2) | yT < criticalRange(1);
%   errT = errT(~junkData);
%   xT = xT(:,~junkData);
% end
% 
% xAll = [xT xV];
% errAll = [errT errV];

% Create a handle to the ResponseModel in coded variables.
xCoded2xSTR = '';
for i1 = 1:n
  if varTransNum(i1) == 2
    xCoded2xSTR = [xCoded2xSTR '10.^( 0.5*log10(domrng(' num2str(i1) ',2)/domrng(' num2str(i1) ',1))*(x(' num2str(i1) ',:) + 1) + log10(domrng(' num2str(i1) ',1)) );']; %#ok
  else
    xCoded2xSTR = [xCoded2xSTR '0.5*(diff(domrng(' num2str(i1) ',:)))*(x(' num2str(i1) ',:) + 1) + domrng(' num2str(i1) ',1);']; %#ok
  end
end
xCoded2xSTR(end) = ']';
xCoded2xSTR = ['[' xCoded2xSTR];

eval(['M_A = @(x) rapideval(RM,' xCoded2xSTR ');']);

% Handle to a modified model that is everywhere defined that can be used
% with fmincon to assess the fitting error without difficulty. If the
% ResponseModel evaluates to NaN, that value will be replaced with
% criticalRange +/- 50;

% TODO: we need to verify that the user's criticalRange(1) is posivite if
% we decide upon a log10 transformation for Y.
if ~isinf(criticalRange(2)) 
  M_A_clipped = @(x) clipy(x,M_A,criticalRange(2)+50);
elseif ~isinf(criticalRange(1))
  % Manipulate to avoid problems with log10 transform.
  if strcmp(respTrans,'log10')
    M_A_clipped = @(x) clipy(x,M_A,0.5*criticalRange(1));
  else
    M_A_clipped = @(x) clipy(x,M_A,criticalRange(1)-50);
  end
else
  M_A_clipped = M_A; %assume the model is everywhere defined if no crit range given.
  if any(junkDataT) || any(junkDataV)
    disp('NaN model outputs detected, but the critical range is infinite. Error optimization may encounter difficulties.')
  end
end

% Create handle to metamodel in coded variables.
eval(['meta = @(x) evalPrivate(MM.model,' xCoded2xSTR ',respTrans,varTransNum);']);

% Number of local searches for metamodel. This should be an option and can
% be zero.
n2Opt = 3;

fmopt = optimset('fmincon');
fmopt = optimset(fmopt,'display','off','RelLineSrchBnd',0.05,'RelLineSrchBndDuration',20,'DiffMinChange',1e-4,'maxfunevals',600,'largescale','off');
fmopt2 = optimset(fmopt,'RelLineSrchBnd',0.005); %We reduce the linesearch step size if we end up infeasible. Since we start feasible, this should help...
  
eFlagH = zeros(1,n2Opt);
eFlagL = zeros(1,n2Opt);
ERRH = zeros(1,n2Opt);
ERRL = zeros(1,n2Opt);
XH = zeros(n,n2Opt);
XL = zeros(n,n2Opt);

% Initial seeds:
xHseed = x2Coded(DClab.findExtremeX(xAll,errAll,n2Opt,'max'),domrng,varTransNum);
xLseed = x2Coded(DClab.findExtremeX(xAll,errAll,n2Opt,'min'),domrng,varTransNum);

% Objective functions
if strcmp(respTrans,'log10')
  objfunH = @(x) log10(M_A_clipped(x))-log10(meta(x));
  objfunL = @(x) log10(meta(x))-log10(M_A_clipped(x));
else
  objfunH = @(x) M_A_clipped(x)-meta(x);
  objfunL = @(x) meta(x)-M_A_clipped(x);
end

% Constrant function
if all(isinf(criticalRange))
  confun = [];
elseif isinf(criticalRange(1))
  confun = @(x) criticalRangeUConstFun(x,M_A_clipped,criticalRange(2));
elseif isinf(criticalRange(2))
  confun = @(x) criticalRangeLConstFun(x,M_A_clipped,criticalRange(1));
else
  confun = @(x) criticalRangeConstFun(x,M_A_clipped,criticalRange);
end
 
t1 = cputime;
for i1 = 1:n2Opt
  [XH(:,i1) ERRH(i1) eFlagH(i1)] = fmincon(objfunH,xHseed(:,i1),[],[],[],[],-ones(n,1),ones(n,1),confun,fmopt);
  [XL(:,i1) ERRL(i1) eFlagL(i1)] = fmincon(objfunL,xLseed(:,i1),[],[],[],[],-ones(n,1),ones(n,1),confun,fmopt);

  % If the optimization exited with an infeasible solution, try again
  % with a smaller line search stepsize.
  if eFlagH(i1) == -2
    [XH(:,i1) ERRH(i1) eFlagH(i1)] = fmincon(objfunH,xHseed(:,i1),[],[],[],[],-ones(n,1),ones(n,1),confun,fmopt2);
  end
  if eFlagL(i1) == -2
    [XL(:,i1) ERRL(i1) eFlagL(i1)] = fmincon(objfunL,xLseed(:,i1),[],[],[],[],-ones(n,1),ones(n,1),confun,fmopt2);
  end

  % If maxfunevals occurred but the final value is (essentially) feasible, keep it.
  
  if eFlagH(i1) == 0
    if isempty(confun)
      eFlagH(i1) = 1;
    else
      if  max(confun(XH(:,i1))) <= 1e-4
        eFlagH(i1) = 1;
      end
    end
  end
  if eFlagL(i1) == 0 
    if isempty(confun)
      eFlagL(i1) = 1;
    else
      if  max(confun(XL(:,i1))) <= 1e-4
        eFlagL(i1) = 1;
      end
    end
  end
end
fitInfo.optimErrTime = cputime-t1;

ERRH = -ERRH;
ERRH = ERRH(eFlagH>0);
ERRL = ERRL(eFlagL>0);
XH = XH(:,eFlagH>0);
XL = XL(:,eFlagL>0);

if isempty(ERRL)
  fitInfo.optimFailed(1,1) = true;
  fitInfo.peakErrorFromOptimization(1,1) = min(errAll);
  fitInfo.xGivingErrMin = xLseed(:,1);
else
  fitInfo.optimFailed(1,1) = false;
  [fitInfo.peakErrorFromOptimization(1,1) idx] = min(ERRL);
  fitInfo.xGivingErrMin = XL(:,idx);
end

if isempty(ERRH)
  fitInfo.optimFailed(1,2) = true;
  fitInfo.peakErrorFromOptimization(1,2) = max(errAll);
  fitInfo.xGivingErrMax = xHseed(:,1);
else
  fitInfo.optimFailed(1,2) = false;
  [fitInfo.peakErrorFromOptimization(1,2) idx] = max(ERRH);
  fitInfo.xGivingErrMax = XH(:,idx);
end
fitInfo.xGivingErrMin = coded2X(fitInfo.xGivingErrMin,domrng,varTransNum);
fitInfo.xGivingErrMax = coded2X(fitInfo.xGivingErrMax,domrng,varTransNum);

MM.fitInfo.peakErrorFromOptimization = fitInfo.peakErrorFromOptimization;
MM.fitInfo.xGivingErrMin = fitInfo.xGivingErrMin;
MM.fitInfo.xGivingErrMax = fitInfo.xGivingErrMax;


% Local functions

function xC = x2Coded(x,domrng,varTransNum)

xC = zeros(size(x));
for i1 = 1:size(domrng,1)
  if varTransNum(i1) == 2;
    xC(i1,:) = 2/log10(domrng(i1,2)/domrng(i1,1))*(log10(x(i1,:)) - log10(domrng(i1,1))) - 1;
  else
    xC(i1,:) = 2*(1/diff(domrng(i1,:)))*(x(i1,:) - domrng(i1,1)) - 1;
  end
end

function x = coded2X(xC,domrng,varTransNum)

x = zeros(size(xC));
for i1 = 1:size(domrng,1)
  if varTransNum(i1) == 2;
    x(i1,:) = 10.^( 0.5*log10(domrng(i1,2)/domrng(i1,1))*(xC(i1,:) + 1) + log10(domrng(i1,1)) );
  else
    x(i1,:) = 0.5*(diff(domrng(i1,:)))*(xC(i1,:) + 1) + domrng(i1,1); 
  end
end

function [cval eval] = criticalRangeConstFun(x,fx,critRng)
y = fx(x);
cval = [y-critRng(2); critRng(1)-y];
eval = 0;

function [cval eval] = criticalRangeUConstFun(x,fx,critRngU)
y = fx(x);
cval = y-critRngU;
eval = 0;

function [cval eval] = criticalRangeLConstFun(x,fx,critRngL)
y = fx(x);
cval = critRngL-y;
eval = 0;
