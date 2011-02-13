function [ub,xopt,mult] = qpub_FMINCON(obj,xinit,opt,relaxEq)
%QPUB_FMINCON produces an upper bound on the optimal value.
%
%   See the help on nqcqp/upperBnd for a description of the input
%   arguments.
%
%   UB = LOWERBND(NQCQPOBJ) is the most simple syntax.
%
%   [UB XOPT MULT] = LOWERBND(NQCQPOBJ,OPT,RELAXEQ) is the most complicated
%   syntax.


% Check if the constraint set is provably empty using the S-procedure. If
% so, make a quick exit. This in the only "return" in this file.
obj = checkSProcInfeas(obj);
if obj.sProcInfeas
  ub = Inf;
  xopt = [];
  mult.lower = zeros(obj.nVars,1);
  mult.upper = zeros(obj.nVars,1);
  mult.quadineq = zeros(size(obj.Zquad,1),1);
  mult.relaxedeq = zeros(size(obj.ZrelaxedEq,1),1);
  return
end

% Initialize inputs
ni = nargin;
if ni==1
  xinit = [];
  opt = [];
  relaxEq = [];
elseif ni==2
  opt = [];
  relaxEq = [];
elseif ni==3
  relaxEq = [];
else
  error(nargchk(1,4,ni));
end

% Problem Data
n = obj.nVars;      % Number of variables
nQuad = size(obj.Zquad,1);
nRelaxEq = size(obj.ZrelaxedEq,1);
Znot = obj.Znot;
Zquad = obj.Zquad;
LB = obj.variableBnds(:,1); 
UB = obj.variableBnds(:,2); 
linXlogX = obj.linXlogX;

% Initialize empty inputs
if isempty(xinit)
  % Create a random initial point. Use 0 for unbounded variables. By
  % the requirements of the object they will be unbounded above and
  % below
  bnded = ~isinf(LB);
  xinit = zeros(n,1);
  xinit(bnded) = (UB(bnded)-LB(bnded).*rand(sum(bnded),1)) + LB(bnded);
  
  % Make sure xinit at least satisfies the equality constraints.
  % Let [i j c1 c2 c3] denote a row of LINXLOGX.
  % Now tweak it to ensure at least the equality constraints are satisfied.
  %           UB(i)-LB(i)
  %    X(i) = ---------- c2 ( c1^[(X(j)-LB(j))/(UB(j)-LB(j))] - 1 ) + LB(i)
  %              c3
  for i1 = 1:size(linXlogX,1)
    ii = linXlogX(i1,1);
    jj = linXlogX(i1,2);
    c1 = linXlogX(i1,3);
    c2 = linXlogX(i1,4);
    c3 = linXlogX(i1,5);
    xinit(ii) = [(UB(ii)-LB(ii))/c3]*c2*(c1^[(xinit(jj)-LB(jj))/(UB(jj)-LB(jj))]-1 ) + LB(ii); %#ok
  end
  
end
if isempty(opt)
  opt = DClab.DCOptions;
end
if isempty(relaxEq)
  relaxEq = false;
end

% would we ever want the notify or final messages to splash to the screen?
if strcmp(opt.display,'ALL')
  dispMode = 'iter';
else
  dispMode = 'off';
end

%Determine which optimization case we're dealing with
% case1: there are no linXlogX equality constraints  and all quadratic
%   constraints and the objective are actually liner. use linprog.
% case2: there are no linXlogX equality constraints, but either the
%   objective or a constraint is quadratic. use fmincon.
% case3: we are supposed to use the quadratic relaxation of the linXlogX
%   equality constraints. use fmincon
% case4: we are supposed to enforce the linXlogX equality constraints. use
%   fmincon.

% Define solveCase, A,B,Aeq,Beq.
A = [];
B = [];
Aeq = [];
Beq = [];
if isempty(obj.ZrelaxedEq)
  linIdx = false(nQuad,1);
  for i2 = 1:length(Zquad)
    if isequal(Zquad{i2}(2:n+1,2:n+1),zeros(n))
      A = [A; 2*Zquad{i2}(1,2:end)];
      B = [B; -Zquad{i2}(1,1)];
      linIdx(i2) = true;
    end
  end
  if sum(linIdx) == length(Zquad) && isequal(Znot(2:n+1,2:n+1),zeros(n))
    solveCase = 1;
  else
    solveCase = 2;
  end
  %remove from the quadratic contraint list those that were actually linear
  Zquad(linIdx) = [];
else
  if relaxEq
    solveCase = 3;
    linIdx = false(nQuad,1);
    for i2 = 1:length(Zquad)
      if isequal(Zquad{i2}(2:n+1,2:n+1),zeros(n))
        A = [A; 2*Zquad{i2}(1,2:end)];
        B = [B; -Zquad{i2}(1,1)];
        linIdx(i2) = true;
      end
    end
    %remove from the quadratic contraint list those that were actually linear
    Zquad(linIdx) = [];
  else
    solveCase = 4;    
  end
end
switch solveCase
  case 1
    lpoptions = optimset('linprog');
    lpoptions = optimset(lpoptions,'Display',dispMode,'TolFun',opt.tolFun,'TolX',opt.tolCon,'TolCon',opt.tolCon);

    % First attempt to solve
    [x,fval,exitflg,output,lambda] = ...
        linprog(2*Znot(1,2:end),A,B,[],[],LB,UB,[],lpoptions);

    % Second attempt to solve. If this also fails, we give up.
    if exitflg <= 0
      if ismember(opt.display,{'notify','all','ALL'})
        str = '  LINPROG failed to find a feasible point, retrying';
        DClab.dcdispstr(str,opt.guiHandle,false)
      end

      %Don't use largescale method this time.
      lpoptions = optimset(lpoptions,'largescale','off');
      % Create a random initial point. Use uniform on [-500, 500] for 
      % unbounded variables.
      idx = ~isinf(obj.variableBnds(:,1));
      xinit = 1000*rand(size(xinit)) - 500;
      xinit(idx) = (UB(idx)-LB(idx)).*rand(sum(idx),1) + LB(idx); 
      [x,fval,exitflg,output,lambda] = ...
        linprog(2*Znot(1,2:end),A,B,[],[],LB,UB,xinit,lpoptions);
    end
    fval = fval+Znot(1,1); %add on affine portion of Znot. 

  case {2,3}

    confun = @(x) DClab.nqcqp.fmconfunRelax(x,[Zquad; obj.ZrelaxedEq]);
    objfun = @(x) DClab.nqcqp.fmobjfunRelax(x,Znot);

    % FMINCON options (large scale doesn't work with nonlin constraints)
    fmoptions = optimset('fmincon');
    try
        fmoptions = optimset(fmoptions,'GradConstr','on','GradObj','on','TolX',opt.tolCon,...
            'LargeScale','off','MaxIter',500,'MaxFunEval',12500,'algorithm','active-set',...
            'TolFun',opt.tolFun,'MaxSQPIter',10000,'Display',dispMode,'TolCon',opt.tolFun);
    catch
        fmoptions = optimset(fmoptions,'GradConstr','on','GradObj','on','TolX',opt.tolCon,...
            'LargeScale','off','MaxIter',500,'MaxFunEval',12500,'algorithm','active-set',...
            'TolFun',opt.tolFun,'MaxSQPIter',10000,'Display',dispMode,'TolCon',opt.tolFun);
        disp('Had issues with optimset in @nqcqp/private/qpub_FMINCON')
    end

    % First attempt to solve
    state = warning('off','optim:fmincon:NLPAlgLargeScaleConflict');
    [x,fval,exitflg,output,lambda]=...
        fmincon(objfun,xinit,full(A),full(B),full(Aeq),full(Beq),full(LB),full(UB),confun,fmoptions);
    warning(state);

    % Second attempt to solve. If this also fails, we give up.
    if exitflg <= 0
      if ismember(opt.display,{'notify','all','ALL'})
        str = '  FMINCON failed to find a feasible point, retrying';
        DClab.dcdispstr(str,opt.guiHandle,false)
      end

      % Create a random initial point. Use uniform on [-500, 500] for
      % unbounded variables.
      idx = ~isinf(obj.variableBnds(:,1));
      xinit = 1000*rand(size(xinit)) - 500;
      xinit(idx) = (UB(idx)-LB(idx)).*rand(sum(idx),1) + LB(idx); 
      [x,fval,exitflg,output,lambda] = ...
          fmincon(objfun,xinit,A,B,Aeq,Beq,LB,UB,confun,fmoptions);
    end

  case 4
    %Here the optimization solver is only aware of the smaller decision
    %vector that results from using the special purpose equality
    %constraints to reduce the number of optimization variables. However,
    %in the objective and constraint functions, we recover the entire
    %decision vector in order to evaluate these. 

    % Define some useful indexing vectors
    x2XElim = obj.linXlogX(:,2);         %xElim = x(x2XElim)
    x2XEqFcnXElim = obj.linXlogX(:,1);   %xEqFcnXElim = x(x2XEqFcnXElim)
    x2XRem = setdiff((1:n)',x2XElim);       %xRem = x(x2XRem)

    % Determine where in xRem it is related by xElim by the equality
    % constraint.

    [trash a2sort b2sort] = intersect(x2XRem,x2XEqFcnXElim);
    xRem2XEqFcnXElim(b2sort,1) = x2XRem(a2sort);

    nElim = length(x2XElim);
    c = zeros(nElim,7);
    c(:,[1 2 3]) = obj.linXlogX(:,[3 4 5]);
    c(:,4) = LB(x2XEqFcnXElim);
    c(:,5) = UB(x2XEqFcnXElim) - LB(x2XEqFcnXElim);
    c(:,6) = LB(x2XElim);
    c(:,7) = UB(x2XElim) - LB(x2XElim);

    % The equality constraint is
    %
    % f_c^{-1}(xl) = c7/log10(c1) * log10( c3/(c2*c5)*(xl-c4) + 1 ) + c6
    %
    % To reduce the computational time, simplify.
    % Let t1 = c7/log10(c1),
    %     t2 = c3/(c2*c5),
    %     t3 = c4;
    %     t4 = c6;
    %
    % Thus the equality constraint is
    %
    % f_c^{-1}(xl) = t1 * log10(t2*(xl-t3) + 1) + t4;

    t = zeros(nElim,4);
    t(:,1) = c(:,7)./log10(c(:,1));
    t(:,2) = c(:,3)./(c(:,2).*c(:,5));
    t(:,3) = c(:,4);
    t(:,4) = c(:,6);

    confun = @(xRem) DClab.nqcqp.fmconfun(xRem,Zquad,x2XRem,x2XElim,xRem2XEqFcnXElim,t);
    objfun = @(xRem) DClab.nqcqp.fmobjfun(xRem,Znot,x2XRem,x2XElim,xRem2XEqFcnXElim,t);

    A=[]; B=[]; %No generalized LP constraints
    Aeq=[]; Beq=[];  %No equality constraints

    % FMINCON options (large scale doesn't work with nonlin constraints)
    fmoptions = optimset('fmincon');
    fmoptions = optimset(fmoptions,'GradConstr','on','GradObj','on','TolX',opt.tolCon,...
        'LargeScale','off','MaxIter',500,'MaxFunEval',12500, ...
        'TolFun',opt.tolFun,'MaxSQPIter',10000,'Display',dispMode,'TolCon',opt.tolFun);

    %First attempt to solve.  
    state = warning('off','optim:fmincon:NLPAlgLargeScaleConflict');
    [xRem,fval,exitflg,output,lambda]=...
        fmincon(objfun,xinit(x2XRem),A,B,Aeq,Beq,LB(x2XRem),UB(x2XRem),confun,fmoptions);
    warning(state);

    %Second attempt to solve
    if exitflg <= 0
      if ismember(opt.display,{'notify','all','ALL'})
        str = '  fmincon failed to find a feasible point, retrying';
        DClab.dcdispstr(str,opt.guiHandle,false)
      end

      % Create a random initial point. Use uniform on [-500, 500] for 
      % unbounded variables.
      idx = ~isinf(obj.variableBnds(:,1));
      xinit = 1000*rand(size(xinit)) - 500;
      xinit(idx) = (UB(idx)-LB(idx)).*rand(sum(idx),1) + LB(idx);
      [xRem,fval,exitflg,output,lambda]=...
          fmincon(objfun,xinit(x2XRem),A,B,Aeq,Beq,LB(x2XRem),UB(x2XRem),confun,fmoptions);
    end

    %Recover full decision vector.
    x = zeros(n,1);
    x(x2XRem) = xRem;
    x(x2XElim) = t(:,1).*log10( t(:,2).*(xRem(xRem2XEqFcnXElim)-t(:,3)) + 1 ) + t(:,4);

end %switch solveCase

%if still infeas, return Inf. Set lambda=[] so we can catch this later
if exitflg <= 0
  if ismember(opt.display,{'notify','all','ALL'})
    str = '  FMINCON failed to find a feasible point, returning infinite bound';
    DClab.dcdispstr(str,opt.guiHandle,false)
  end
  ub = inf;
  xopt = x; %still return x because fmincon tries to minimize 
            %constraint violation on infeas problems. this x "may" be in the true feasible set 
  xopt = []; %no this creates an error!
  lambda = [];
else
  % Output information
  ub = fval;
  xopt = x;
end

%parse lambda to create the desired multiplier structure.

if isempty(lambda)
  %solution bombed return all zeros
  mult.lower = zeros(n,1);
  mult.upper = zeros(n,1);
  mult.quadineq = zeros(nQuad,1);
  mult.relaxedeq = zeros(nRelaxEq,1);

else
  switch solveCase
    case 1
      mult.lower = lambda.lower;
      mult.upper = lambda.upper;
      mult.quadineq = lambda.ineqlin;
      mult.relaxedeq = zeros(0,1);
    case 2
      mult.lower = lambda.lower;
      mult.upper = lambda.upper;
      mult.quadineq = zeros(nQuad,1);
      mult.quadineq(linIdx) = lambda.ineqlin;
      mult.quadineq(~linIdx) = lambda.ineqnonlin;
      mult.relaxedeq = zeros(0,1);
    case 3
      mult.lower = lambda.lower;
      mult.upper = lambda.upper;
      mult.quadineq = zeros(nQuad,1);
      mult.quadineq(linIdx) = lambda.ineqlin;
      nQuadAsQuad = sum(~linIdx);
      mult.quadineq(~linIdx) = lambda.ineqnonlin(1:nQuadAsQuad);
      mult.relaxedeq = lambda.ineqnonlin(nQuadAsQuad+1:end);
    case 4
      mult.lower = zeros(n,1);
      mult.lower(x2XRem) = lambda.lower;
      mult.upper = zeros(n,1);
      mult.upper(x2XRem) = lambda.upper;
      mult.quadineq = lambda.ineqnonlin;
      mult.relaxedeq = zeros(nRelaxEq,1);
  end
end

%TODO, for times when exitflg==0, max iterations have occured, but the
%returned point might be feasible. If so, we should return non inf values.


%       %===Solver could not find a solution. Get user assistance===
% 
%       if exitflg <= 0 && strcmp(solver,'fmincon')
%         disp(['fmincon failed to find a soln, adjust options and run the following to ' ...
%               'retry']);
%         disp(['[x,fval,exitflg,output,multipliers]= fmincon' ...
%           '(@fmobjfun,xinit{i1},A,B,Aeq,Beq,LB,UB,@fmconfun,options,Z0,Zi)']);
%         disp('Alternatively, type return to terminate with no solution found')
%         x_ = x; fval_ = fval; lambda_ = lambda; 
%         fval = -inf; x = []; lambda = [];
%         keyboard
%       end
% 
%       if exitflg <= 0 && strcmp(solver,'linprog')
%         disp(['linprog failed to find a soln, adjust options and run the following to ' ...
%               'retry']);
%         disp(['[x,fval,exitflg,output,lambda]=' ...
%           'linprog(2*Z0(1,2:end),A,B,[],[],LB,UB,[],options)']);
%         disp('Alternatively, type return to terminate with no solution found')
%         x_ = x; fval_ = fval; lambda_ = lambda; 
%         fval = -inf; x = []; lambda = [];
%         keyboard
%       end
%       %===End user assistance===














