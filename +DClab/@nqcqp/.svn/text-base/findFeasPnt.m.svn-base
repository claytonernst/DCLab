function obj = findFeasPnt(obj,opt)
%findFeasPnt locates 5 points in the feasible set
%
%   NQCQPOBJ = findFeasPnt(NQCQPOBJ) searches for a points in the "center"
%   of the feasible set that gives "small" values of the objective
%   function. The function is used for predictions over crude
%   approximations of the feasible set. The hope is that x's toward the
%   middle of the crude approximation will fall in the actual feasible set.
%   Among such points, we'd like to chose ones that give small oosts. To
%   accomplish this, this method populates the xfeas field of NQCQPOBJ with
%   a 1x5 cell array of feasible points (or empties if none are found).
%   From left to right, these give successively smaller costs as the degree
%   to which the points are forced to the center of the feasible set is
%   relaxed. The method simply returns NQCQPOBJ if this field is already
%   populated.
%
%   NQCQPOBJ = checkSProcInfeas(NQCQPOBJ,OPT) uses options taken from the
%   DCOptions object OPT. This method uses OPT.display, OPT.guihandle,
%   OPT.nRestart, OPT.tolCon, and OPT.tolFun.

%TODO, using a slack variable to find the "center" of the feasible set is
%only a good idea if the constraints are normalized. Currently they are
%not, so this function needs some additional work to be truly effective.

error(nargchk(1,2,nargin));

% Special cases that warrant a speedy exit:

% Provably infeasible
if ~isempty(obj.xfeas)
  %nothing to do, quick exit
  return
end

% Only box constraints
if isempty(obj.Zquad) && isempty(obj.ZrelaxedEq)
  xfeas = cell(1,5);
  LB = obj.variableBnds(:,1); 
  UB = obj.variableBnds(:,2);
  bnded = ~isinf(LB);
  for i1= 1:5
    x = zeros(obj.nVars,1);
    x(bnded) = (UB(bnded)-LB(bnded).*rand(sum(bnded),1)) + LB(bnded);
    xfeas{i1} = x;
  end
  obj.xfeas = xfeas;
  return
end

% Optional input initialization:
if nargin==1
  opt = DClab.DCOptions;
end

% Variable initialization:
xfeas = cell(1,5);

% Useful constants:
guihand = opt.guiHandle;
n = obj.nVars;
nQuad = size(obj.Zquad,1);
nRelaxedEq = size(obj.ZrelaxedEq,1);

% Check if the S-procedure can prove the feasible set is empty. We want
% to do this here, rather than as part of the upperBnd method so that value
% of .sProcInfeas is available in this function.
obj = checkSProcInfeas(obj);

% Attempt to find a feasible point. We will change each quadratic
% constraint from Z(x) <= 0 to Z(x,t) <= t and then minimize t. To
% accomplish this, we will make a new nqcqp object and use the upperBnd 
% method.

Znot = zeros(n+2);
Znot(1,end) = 0.5;
Znot(end,1) = 0.5;

Zquad = obj.Zquad;
for i1 = 1:nQuad
  Zquad{i1} = blkdiag(Zquad{i1},0);
  Zquad{i1}(1,end) = -0.5;
  Zquad{i1}(end,1) = -0.5;
end

ZrelaxedEq = obj.ZrelaxedEq;
for i1 = 1:nRelaxedEq
  ZrelaxedEq{i1} = blkdiag(ZrelaxedEq{i1},0);
  ZrelaxedEq{i1}(1,end) = -0.5;
  ZrelaxedEq{i1}(end,1) = -0.5;
end

tmpObj = DClab.nqcqp;
tmpObj.Znot = Znot;
tmpObj.Zquad = Zquad;
tmpObj.ZrelaxedEq = ZrelaxedEq;
tmpObj.nVars = n+1;
tmpObj.variableBnds = [obj.variableBnds; -Inf Inf];
tmpObj.linXlogX = obj.linXlogX;
tmpObj.knownFeas = obj.knownFeas;
tmpObj.sProcInfeas = obj.sProcInfeas;

% Find an initial point. We will randomly generate an initial x, and then
% add on the max constraint violation as the initial t.

% Use 0 for unbounded variables. By the requirements of the object they
% will be unbounded above and below
LB = obj.variableBnds(:,1); 
UB = obj.variableBnds(:,2);
bnded = ~isinf(LB);
xinit = zeros(n,1);
xinit(bnded) = (UB(bnded)-LB(bnded).*rand(sum(bnded),1)) + LB(bnded);

xx = [1; xinit; 0];
vio1 = zeros(nQuad,1);
vio2 = zeros(nRelaxedEq,1);
for i1 = 1:nQuad
  vio1 = xx'*Zquad{i1}*xx;
end
for i1 = 1:nRelaxedEq
  vio2 = xx'*ZrelaxedEq{i1}*xx;
end

xinit = [xinit; max([vio1; vio2])];
[ub xopt] = upperBnd(tmpObj,xinit,opt);

%If no strictly feasible point (t <-1e-4) was found, try again up to
%nRestarts more times  
if ub < -1e-4
  noneFound = false;
else 
  noneFound = true;
end

nTried = 1;
while noneFound && nTried <= opt.nRestart+1
  if strcmpi(opt.display,'all')
    str = '    fmincon failed to find a feasible point in nqcqp\findFeas, retrying';
    DClab.dcdispstr(str,guihand,false)
  else
    drawnow
  end
  [ub xopt] = upperBnd(tmpObj,[],opt);
  nTried = nTried+1;
  if ub < -1e-4
    noneFound = false;
  end
end %end while

if noneFound
  if ismember(opt.display,{'notify','all','ALL'})
    str = '    Although the S-procedure could not prove infeasibility, fmincon failed to find a feasible point in nqcqp\findFeasPnt';
    DClab.dcdispstr(str,guihand,false)
  end
else
  %We have a feasible point and t is strictly less then zero so the
  %feasible set has volume. Now optimize original objective function
  %subject to 
    %a) Si(x) <= 0.95*t
    %b) Si(x) <= 0.8*t
    %c) Si(x) <= 0.6*t
    %d) Si(x) <= 0.4*t
    %e) Si(x) <= 0.2*t
  %
  %Each of these should be feasible, and we'll start them with a
  %random point. If they turn out to be infeas (it's local search)
  %we'll restart them once with goodx. %We start with a random point
  %because while goodx is feasible, it probably doesn't give a good cost.
  %If we fail to find a feasible point with the random seed, then we revert
  %to starting from goodx.
  
  goodx = xopt(1:end-1); %eliminate t variable.
  t = xopt(end); %should be equal to ub.
  
  % Eval the real objective function at goodx so we can track whether the
  % random initial point outperforms the known-to-be-feasible goodx.
  objVal = [1 goodx']*obj.Znot*[1;goodx];
  scales = [0.95 0.8 0.6 0.4 0.2];
  
  for i1 = 1:5
    % Adjust the quadratic constraints
    tmpObj = obj;
    for i2 = 1:nQuad
      tmpObj.Zquad{i2}(1,1) = tmpObj.Zquad{i2}(1,1) - t*scales(i1);
    end
    % Solve optimization with a randomly generated initial point.
    [ub1 xopt1] = upperBnd(tmpObj,[],opt);
    
    %If we did worse than objVal, or worse yet failed to find a feasible
    %solution, try once more starting from goodx
    if ub1 > objVal
      if strcmpi(opt.display,'all')
        str = '    Random point failed to beat goodx in nqcqp\findFeasPnt, retrying with goodx';
        DClab.dcdispstr(str,guihand,false)
      else
        drawnow
      end
      [ub2 xopt2] = upperBnd(tmpObj,goodx,opt);
      if isinf(ub2)
        if ismember(opt.display,{'notify','all','ALL'})
          str = '    goodx was infeasible in nqcqp\findFeasPnt, this is weird';
          DClab.dcdispstr(str,guihand,false)
        end
        %for some weird reason we couldn't find a feasible point, return
        %empty
        xfeas{1,i1} = [];
      else
        if ub1 < ub2
          if strcmpi(opt.display,'all')
            str = '    random point beat goodx in nqcqp\findFeasPnt';
            DClab.dcdispstr(str,guihand,false)
          end
          xfeas{1,i1} = xopt1;
        else
          if strcmpi(opt.display,'all')
            str = '    goodx beat random point in nqcqp\findFeasPnt';
            DClab.dcdispstr(str,guihand,false)
          end
          xfeas{1,i1} = xopt2;
        end
        %Save the best current objective value to use as our bench mark for
        %deciding whether to retry. As i1 increases the quadratic
        %constraints are relaxed, so objVal should decrease. If it doesn't,
        %a restart is warranted, and we may as well use goodx. Perhaps
        %restarting from random would be better, but I lack the test cases
        %to test this.
        objVal = min([objVal ub1 ub2]); 
      end
    else
      xfeas{1,i1} = xopt1;
    end %if ub1 > objVal
  end  %for i1 = 1:5
end

obj.xfeas = xfeas;
