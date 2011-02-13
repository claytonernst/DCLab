function nqcqpObj = checkSProcInfeas(nqcqpObj,opt)
%checkSProcInfeas determines if the constraint set is provably empty using the S-procedure.
%
%   NQCQPOBJ = checkSProcInfeas(NQCQPOBJ) populates the .sProcInfeas field of
%   NQCQP with 1 if the S-procedure proves the constraint set is empty, and
%   with 0 otherwise. The method simply returns NQCQPOBJ if this field is
%   already populated.
%
%   NQCQPOBJ = checkSProcInfeas(NQCQPOBJ,OPT) uses options taken from the
%   DCOptions object OPT. This method uses OPT.display and
%   OPT.sedumiParEps.
%
%
% TODO: old help below
%
% Function to check the feasibility of the constraint set using the s-proc
% relaxation. If type is consist or optimXXXnorm, the function quickly
% exists, because these problems are always feasible. It type=pred, then we
% call sedumi to determine if the constraints are feasible. This function
% is called before initiating a local search, to see if it is a waste of
% time.


error(nargchk(1,2,nargin))
error(nargoutchk(0,1,nargout))
if nargin==1
  opt = DClab.DCOptions;
end

% Only do stuff if the .sProcInfeas field is empty
if isempty(nqcqpObj.sProcInfeas)
  
  % Problem Data
  n = nqcqpObj.nVars;

  %convert variableBnds to quadratic constraints
  xbnds = nqcqpObj.variableBnds;

  %==eliminate any variables with infinite bounds
  temp = isinf(xbnds(:,1));
  xbnds = xbnds(~temp,:);
  %==end eliminate inf bounds

  % Create quadratic representations of the box constraints
  nn = size(xbnds,1);
  Zbox = cell(nn,1);
  [Zbox{:}] = deal(spalloc(n+2,n+2,4));
  affineTerm = xbnds(:,1).* xbnds(:,2);
  linTerm = -xbnds(:,1) - xbnds(:,2);
  for i2 = 1:nn
    Zbox{i2}(1,1) = affineTerm(i2);
    Zbox{i2}(i2+1,1) = linTerm(i2)/2;
    Zbox{i2}(1,i2+1) = linTerm(i2)/2;
    Zbox{i2}(i2+1,i2+1) = 1;
  end
  
  % We know the box constraints by themselves are feasible. Add a slack
  % variable t to the Zquad and ZrelaxedEq contraints. This gives us a
  % feasible constraint set over which we can minimize t. If the lower
  % bound on the optimal t is positive, the s-procedure has shown the
  % original constraint set is empty. 
  %
  % TODO: since we just need feasibility, we might gain some speed-up by
  % setting t = 0 and just searching for multipliers, rather than
  % optimizing t subject to such multipliers existing.
  
  Zquad = nqcqpObj.Zquad;
  ZrelaxedEq = nqcqpObj.ZrelaxedEq;
  
  % Add slack variable t to constraint matrices (t is the last variable)
  for i2 = 1:size(Zquad,1)
    Zquad{i2} = blkdiag(Zquad{i2},0);
    Zquad{i2}(1,end) = -0.5;
    Zquad{i2}(end,1) = -0.5;
  end
  for i2 = 1:size(ZrelaxedEq,1)
    ZrelaxedEq{i2} = blkdiag(ZrelaxedEq{i2},0);
    ZrelaxedEq{i2}(1,end) = -0.5;
    ZrelaxedEq{i2}(end,1) = -0.5;
  end
  
  % The objective function is t
  obj = spalloc(n+2,n+2,2);
  obj(1,end) = 0.5;
  obj(end,1) = 0.5;

  Zi = [Zbox; Zquad; ZrelaxedEq];
  N = length(Zi);
  
  % Sedumi (dual-form) notation:
  %     max b'y  s.t. c-A'*y >= 0 
  % Define y:=[lambda_1 ... lambda_N gamma]
  b = spalloc(N+1,1,1); b(end) = 1;

  % N linear constraints:
  %  lambda >=0
  K.l = N;
  cl = spalloc(N,1,0);
  Al = [speye(N) spalloc(N,1,0)];

  % (n+1)x(n+1) LMI constraint:
  %    Z0 - [gamma 0; 0 0] + sum_k lambda_k Z_k >= 0
  K.s = n+2;
  clmi = sparse(obj(:));
  Almi = [];
  if N > 0
    Almi = zeros(numel(Zi{1}),N);
    for i2 = 1:N
      Almi(:,i2) = Zi{i2}(:);
    end
  end
  temp = spalloc(n+2,n+2,1);
  temp(1,1) = 1;
  Almi = sparse([Almi -temp(:)]);

  % Concatenate the constraints
  c = [cl; clmi];
  At = -[Al; Almi];

  % Call Sedumi
  if strcmp(opt.display,'ALL')
    pars.fid = 1;
  else
    pars.fid=0;    % fid=0: no output; fid=1: display output
  end
  pars.eps=opt.sedumiParEps; % Desired accuracy, sedumi default = 1e-8;

  [x,y,info]=sedumi(At,b,c,K,pars); %#ok

  % Output information
  lb = y(end);

  %this function is not super critical, check only for severe numerical
  %problem and display nothing otherwise
  if info.numerr == 2 && ismember(opt.display,{'notify','all','ALL'})
    str = '  SeDuMi numerical error (noncritical) in nqcqp/checkSProcFeas. Problem may be badly scaled';
    DClab.dcdispstr(str,opt.guiHandle,false)
    lb(i1,1) = inf;
  end

  if lb > 0
    nqcqpObj.sProcInfeas = true; %provably infeasible
  else
    nqcqpObj.sProcInfeas = false; %feasibility inconclusive, but s-proc relaxation is feas
  end
end


% old function below

% 
% 
% 
% constr = nqcqpObj(i1).constraints;
% vBnds = nqcqpObj(i1).variableBnds;
% type = nqcqpObj(i1).type;
% 
% if ~isempty(nqcqpObj(i1).sProcFeas)
%   type = 'unnec'; %set type to unnec for a quick exit
% end
% 
% switch type
%  case 'unnec'
%   bool = nqcqpObj(i1).sProcFeas; %do nothing, just leave existing value
% 
%  case {'consist';'optim1norm';'optim2norm';'optimInfnorm'}
%   bool = true; %these are always feasible
%  case 'pred'
% 
%   nvar = size(vBnds,1);
%   %==eliminate any variables with infinite bounds
%   temp = isinf(vBnds);
%   vBnds = vBnds(sum(temp,2)==0,:);
%   nn = size(vBnds,1);
%   %==end eliminate variables with infinite bounds
% 
%   % Create quadratic matrices representing the variable bounds
%   Zx = cell(nn,1);
%   [Zx{:}] = deal(spalloc(nvar+2,nvar+2,4));
%   affineTerm = vBnds(:,1).* vBnds(:,2);
%   linTerm = ( -vBnds(:,1) - vBnds(:,2) )/2;
%   for i2 = 1:nn
%     Zx{i2}(1,1) = affineTerm(i2);
%     Zx{i2}(i2+1,1) = linTerm(i2);
%     Zx{i2}(1,i2+1) = linTerm(i2);
%     Zx{i2}(i2+1,i2+1) = 1;
%   end
% 
%   % Add slack variable to constraint matrices (t is the last variable)
%   for i2 = 1:length(constr)
%     constr{i2} = blkdiag(constr{i2},0);
%     constr{i2}(1,end) = -0.5;
%     constr{i2}(end,1) = -0.5;
%   end
% 
%   % The objective function is t
%   obj = spalloc(nvar+2,nvar+2,2);
%   obj(1,end) = 0.5;
%   obj(end,1) = 0.5;
% 
%   Zi = [constr; Zx];
%   N = length(Zi);
% 
%   % Sedumi (dual-form) notation:
%   %     max b'y  s.t. c-A'*y >= 0 
%   % Define y:=[lambda_1 ... lambda_N gamma]
%   b = spalloc(N+1,1,1); b(end) = 1;
% 
%   % N linear constraints:
%   %  lambda >=0
%   K.l = N;
%   cl = spalloc(N,1,0);
%   Al = [speye(N) spalloc(N,1,0)];
% 
%   % (n+1)x(n+1) LMI constraint:
%   %    Z0 - [gamma 0; 0 0] + sum_k lambda_k Z_k >= 0
%   K.s = nvar+2;
%   clmi = sparse(obj(:));
%   Almi = [];
%   if N > 0
%     Almi = zeros(numel(Zi{1}),N);
%     for i2 = 1:N
%       Almi(:,i2) = Zi{i2}(:);
%     end
%   end
%   temp = spalloc(nvar+2,nvar+2,1);
%   temp(1,1) = 1;
%   Almi = sparse([Almi -temp(:)]);
% 
%   % Concatenate the constraints
%   c = [cl; clmi];
%   At = -[Al; Almi];
% 
%   % Call Sedumi
%   if strcmp(opt.display,'ALL')
%     pars.fid = 1;
%   else
%     pars.fid=0;    % fid=0: no output; fid=1: display output
%   end
%   pars.eps=opt.sedumiParEps; % Desired accuracy, sedumi default = 1e-8;
% 
%   [x,y,info]=sedumi(At,b,c,K,pars); %#ok
% 
%   % Output information
%   lb = y(end);
% 
%   %this function is not super critical, check only for severe numerical
%   %problem and display nothing otherwise
%   if info.numerr == 2 && ismember(opt.display,{'notify','all','ALL'})
%     str = '  SeDuMi numerical error (noncritical) in nqcqp/checkSProcFeas. Problem may be badly scaled';
%     DClab.dcdispstr(str,opt.guiHandle,false)
%     lb(i1,1) = inf;
%   end
% 
%   if lb > 0
%     bool = false; %provably infeasible
%   else
%     bool = true; %feasibility inconclusive, but s-proc relaxation is feas
%   end
% 
%  otherwise
%   error('Invalid type')
% end
% 
% nqcqpObj(i1).sProcFeas = bool;

