function [LB gamma] = evalCMeasObj(Dset,pt,useOutputUnc)
% evalCMeasObj. 
%
%   LB = evalCMeasObj(Dset,pt) Evaluates the cost function of the
%   consistency measure optimization at pt using the true models in the
%   ResponseModel objects but ignoring their output uncertainty (if any).
%   The result is a quasilower bound on the consistency measure of Dset
%
%   LB = evalCMeasObj(Dset,pt,useOutputUnc). If useOutputUnc==true, the
%   method evaluates the cost function of the consistency measure
%   optimization at pt using the true models and accounting for any
%   uncertainty in their outputs.  The result is a true lower bound on the
%   consistency measure of Dset. If useOutputUnc==false (the default), the
%   behavior is identical to the previous case.
%
%   Inputs:
%      Dset: a DCDataset object
%      pt: the point at which to evaluate the cost function. Can be an nx1
%      vector, or an rx1 cell array of these.

%TODO: do we want the quasilower bound to be the default?

switch nargin
  case 2
    useOutputUnc = false;
  otherwise
    error(nargchk(2,3,nargin))
end   

% Very cursory input checking.

n = Dset.nParameters;
if isnumeric(pt)
  if ~isequal(size(pt),[n 1])
    error('Inputs: dimensions of pt incompatible with the dataset')
  end
  pt = {pt};
elseif iscell(pt)
  if size(pt,2) ~= 1
     error('Inputs: invalid dimensions of cell array pt')
  end
else
  error('Inputs: invalid datatype for pt')
end

LB = zeros(size(pt));
for i1 = 1:size(pt,1)
  tempPt = pt{i1,1};
  
  [y yInt] = evalResponseModels(Dset,tempPt);

  if useOutputUnc
    yR = yInt(1:2:end-1);
    yL = yInt(2:2:end);
  else
    yR = y;
    yL = y;
  end
  
  % Given that yInt(1) <= eta(x) <= yInt(2)
  %
  % If  yInt(2) <= d + u(1-gamma)
  % and d + l(1-gamma) <= yInt(1)
  %
  % Then l(1-gamma) <= eta(x) - d <= u(1-gamma),
  % So that gamma is feasible.

  % Loop through all Pairs, and create a vector of the largest gamma
  % feasible for each pair. Then we can min over these to get the largest
  % feasible gamma for the supplied vector tempPt.

  m = Dset.nPairs;
  d = vertcat(Dset.ModelAndObservationPair.observedValue);
  uVect = vertcat(Dset.ModelAndObservationPair.observationUncertaintyPlusMinus);
  uCase = vertcat(Dset.ModelAndObservationPair.uncertaintyCase);
  gamma = zeros(m,1);
  
  % Use this point to lower bound the consistency measure
  for i2 = 1:m
    switch uCase(i2)
      case 1
        gamma(i2) = min( (d(i2) + uVect(i2,2) - yR(i2))/uVect(i2,2) , (d(i2) + uVect(i2,1) - yL(i2))/uVect(i2,1) );
      case 2
        gamma(i2) = min( (d(i2)*(1 + uVect(i2,2)) - yR(i2))/(d(i1)*uVect(i2,2)) , (d(i2)*(1 + uVect(i2,1)) - yL(i2))/(d(i2)*uVect(i2,1)) );
      case 3
        % eta(x)<=yR, yR <= d10^u(1-gamma) ==>  eta(x) <= d10^u(1-gamma), 
        % yL <= eta(x), d10^el(1-gamma) <= yL ==> d10^el(1-gamma) <= eta(x)
        gamma(i2) = min( (log10(d(i2)) - log10(yR(i2)) + uVect(i2,2))/uVect(i2,2) , (log10(d(i2)) - log10(yL(i2)) + uVect(i2,1))/uVect(i2,1) );
      case 4
        % eta(x)<=yR, yR <= d^(1+u(1-gamma)) ==>  eta(x) <= d^(1+u(1-gamma)), 
        % yL <= eta(x), d^(1+el(1-gamma)) <= yL ==> d^(1+el(1-gamma)) <= eta(x)
        gamma(i2) = min( (log10(d(i2))*(1 + uVect(i2,2)) - log10(yR(i2)))/(log10(d(i2))*uVect(i2,2)) , (log10(d(i2))*(1 + uVect(i2,1)) - log10(yL(i2)))/(log10(d(i2))*uVect(i2,1)) );
        
      otherwise
        error('Code Not Complete')
    end
  end

  %  disp('Gamma values for each constraint. The smallest is a lower bound on the Consistency Measure.')
%  disp(gamma)
  LB(i1) = min(gamma);
end

