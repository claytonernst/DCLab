function multnew = addMultZeros(mult,PD,opt,activeIdx,nodeIndices)
% function multnew = addMultZeros(mult,PD,opt,activeIdx,nodeIndices)
%
% When we perform an analysis using the surrogate models, not all the
% parameters of the dataset necessarily enter the optimization. The
% function realigns everything so that the multipliers in the outputs
% of say PolyDataset/olbCons correspond to all the dataset
% parameters. It accomplishes this by shifting things around and
% adding zeros where appropriate.


r = length(nodeIndices);
n = PD.nParameters;
m = PD.nPairs;

t = {};
const = opt.constraints;
if const(1) || const(2)
  t = [t; {'linY'}];
end
if const(3) || const(4)
  t = [t; {'logY'}];
end

if length(t) == 1
  nExpMults = 2*m;
elseif length(t) == 2
  nExpMults = 4*m;
else
  error('Usage: incorrect number of experiment constraints')
end

% if ~isequal(size(mult),[r,1])
%   error('Incorrect dimensions: mult')
% end
% if ~isequal(size(activeIdx),[r,1])
%   error('Incorrect dimensions: activeIdx')
% end
% 
% multnew = cell(r,1);
% 
% for i1 = 1:r
%   if isempty(mult{i1})
%     multnew{i1} = [];
%   elseif isstruct(mult{i1})
%     multnew{i1} = mult{i1};
%     multnew{i1}.lower = zeros(n,1);
%     multnew{i1}.lower(activeIdx{i1}) = mult{i1}.lower;
%     multnew{i1}.upper = zeros(n,1);
%     multnew{i1}.upper(activeIdx{i1}) = mult{i1}.upper;
% 
%     multnew{i1}.lower = multnew{i1}.lower.*diff(vertcat(PD.FreeParameter.range),[],2);
%     multnew{i1}.upper = multnew{i1}.upper.*diff(vertcat(PD.FreeParameter.range),[],2);
%   else
%     multnew{i1} = zeros(nExpMults+n,1);
%     multnew{i1}([1:nExpMults nExpMults+activeIdx{i1}']) = mult{i1};
% 
%     %post process multipliers. when we have scaled to -1 <= x <= 1
%     %and pose the problem as x^2 <= 1, the dual inequality looks
%     %like f(xp) >= LB - mult*2*Delta, with Delta a perturbation to
%     %the normalized constraints. 2*Delta*(UB-LB) = a perturbation
%     %to the original constraints.
%     multnew{i1}(nExpMults+1:end) = multnew{i1}(nExpMults+1:end).*diff(vertcat(PD.FreeParameter.range),[],2);
%     
%   end
% end
if isempty(mult)
    multnew = [];
elseif isstruct(mult)
    multnew = mult;
    multnew.lower = zeros(n,1);
    multnew.lower(activeIdx) = mult.lower;
    multnew.upper = zeros(n,1);
    multnew.upper(activeIdx) = mult.upper;

    multnew.lower = multnew.lower.*diff(vertcat(PD.FreeParameter.range),[],2);
    multnew.upper = multnew.upper.*diff(vertcat(PD.FreeParameter.range),[],2);
else
    multnew = zeros(nExpMults+n,1);
    multnew([1:nExpMults nExpMults+activeIdx']) = mult;

    %post process multipliers. when we have scaled to -1 <= x <= 1
    %and pose the problem as x^2 <= 1, the dual inequality looks
    %like f(xp) >= LB - mult*2*Delta, with Delta a perturbation to
    %the normalized constraints. 2*Delta*(UB-LB) = a perturbation
    %to the original constraints.
    multnew(nExpMults+1:end) = multnew(nExpMults+1:end).*diff(vertcat(PD.FreeParameter.range),[],2);

end
