function xinit = findStartPt(PD,node,activeIdx,startPt,Q)
%FINDSTARTPT determines a seed for nqcqp/upperBnd
%
%   XINIT = FINDSTARTPT(PD,NODE,ACTIVEIDX) generates a length(ACTIVEIDX)
%   column vector XINIT by sampling from a uniform distribution over the
%   NODE subdomain of the PolyDataset PD.
%
%   XINIT = FINDSTARTPT(PD,NODE,ACTIVEIDX,STARTPT) massages the n-by-1
%   vector STARTPT to produce the corresponding length length(ACTIVEIDX)
%   column vector XINIT suitable to initialized nqcqp/upperBnd. Note that
%   XINIT may not lie in the NODE subdomain if STARTPT does not.
%
%   XINIT = FINDSTARTPT(PD,ACTIVEIDX,STARTPT,Q) has the same behaviour as
%   the previous unless STARTPT=[]. In this later case, XINIT will be the
%   "mean" vector taken from the optimal Q matrix from the rank relaxation
%   performed by nqcqp/lowerBnd.


% n = PD.nParameters;
% %Determine if we need to include multiplier for logX variables.
% if max(activeIdx) > n
%   xinit = zeros(2*n,2);
% else
%   xinit = zeros(n,2);
% end
ni = nargin;
switch ni
  case 3
    startPt = [];
    Q = [];
  case 4
    Q = [];
  otherwise
   error(nargchk(3,5,ni)) 
end

if ~isempty(startPt)
  n = PD.nParameters;
  bnds = PD.PiecewiseSurrogateModelTree(node).domainRange;
  
  % Make sure supplied start point lies in the subdomain.
  if DClab.issubset([startPt startPt],bnds)

    %Determine if we need to include initial points for logX variables.
    if max(activeIdx) > n
      normStartPtLinX = 2*diag(1./diff(bnds,[],2))*(startPt-bnds(:,1)) - 1;
      normStartPtLogX = 2*diag(1./log10(bnds(:,2)./bnds(:,1)))*log10(startPt./bnds(:,1)) - 1;
      normStartPt = [normStartPtLinX; normStartPtLogX];
      xinit = normStartPt(activeIdx);
    else
      normStartPtLinX = 2*diag(1./diff(bnds,[],2))*(startPt-bnds(:,1)) - 1;
      xinit = normStartPtLinX(activeIdx);
    end
  else
    % Warn, and call recursively to generate a random point.
    warning('Supplied start point does not lie in the corresponding subdomain, a random point will be used') %#ok  
    xinit = findStartPt(PD,node,activeIdx);
  end
  
  % This xinit will be feasible for the special purpose equality
  % constraints.
  
elseif ~isempty(Q)
 xinit = Q(2:end,1);
 % This xinit may not be feasible for the special purpose equality
 % constraints.

else
  n = PD.nParameters;
  bnds = PD.PiecewiseSurrogateModelTree(node).domainRange;
  
  % Generate a random point from a uniform distribute over bnds.
  startPt = rand(n,1).*(bnds(:,2)-bnds(:,1)) + bnds(:,1);
  
  %Determine if we need to include initial points for logX variables.
  if max(activeIdx) > n
    normStartPtLinX = 2*diag(1./diff(bnds,[],2))*(startPt-bnds(:,1)) - 1;
    normStartPtLogX = 2*diag(1./log10(bnds(:,2)./bnds(:,1)))*log10(startPt./bnds(:,1)) - 1;
    normStartPt = [normStartPtLinX; normStartPtLogX];
    xinit = normStartPt(activeIdx);
  else
    normStartPtLinX = 2*diag(1./diff(bnds,[],2))*(startPt-bnds(:,1)) - 1;
    xinit = normStartPtLinX(activeIdx);
  end

  % This xinit will be feasible for the special purpose equality
  % constraints.
  
end
