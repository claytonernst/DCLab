function [y, yInt] = evalResponseModels(D,x,pairs,paramList)
%EVALRESPONSEMODELS evaluates the response models of a DCDataset at x
%
%   Note: Supplying [] for any of the optional inputs will cause the
%   corresponding default input to be used.
%
%   Y = EVALRESPONSEMODELS(D,X) evaluates the response models of the 
%   DCDataset D at the n-by-nEvals matrix of horizontally concatenated
%   parameter vectors X. The order of the parameters in X must correspond 
%   to the order of the parameters in the DCDataset. Y will be an m-by-nEvals
%   matrix, where m is the number of response models in the DCDataset.
%
%   Y = EVALRESPONSEMODLES(D,X,PAIRS) will output Y for only the subset
%   PAIRS of [1:D.nPairs]. Y will be an m-by-nEvals matrix, where m is the
%   length of PAIRS. PAIRS should be a column array and the default is
%   PAIRS = [1:D.nPairs]'. 
%
%   Y = EVALRESPONSEMODELS(D,X,PAIRS,PARAMLIST) allows you supply X in
%   order the the names appear in the column cell array of chars PARAMLIST.
%   Obviously you must include all parameters that are needed by the
%   models, but you have some flexibility with order and total number of
%   parameters in X. The default is PARAMLIST = D.parameterList.
%
%   [Y YINT] = EVALRESPONSEMODLES(...) additionally returns an
%   2*m-by-nEvals matrix of "intervals" that should (through consideration
%   of any output uncertainty in the response models) contain the true
%   value of y(x). Specifically, YINT([2*i-1 2*i],j) is a 2-by-1, whose
%   elements lower and upper bound y_i(X(:,j)), respectively.
%
%   See also ResponseModel/eval, PolyDataset/evalSurrogateModels


ni = nargin;
switch ni
  case 2
    pairs = [];
    paramList = [];
  case 3
    paramList = [];
  otherwise
    error(nargchk(2,4,ni))
end

% Check/initialize inputs.

if isempty(pairs)
  pairs = (1:D.nPairs)';
elseif size(pairs,2) ~= 1 || max(pairs) > D.nPairs
  error('Inputs: PAIRS must be a column array and its elements must reference MOPairs that are present in D')  
else
  %do nothing, input should be ok
end

if isempty(paramList)
  %assign the paramList and assume x was supplied in the proper order
  paramList = D.parameterList;
  if size(x,1)~=size(paramList,1)
    error('Inputs: X must be n-by-Nevals, where n is the number of DCDataset parameters.')
  end
else
  dsetPList = D.parameterList;
  [trash idx1 idx2] = intersect(paramList,dsetPList); %#ok
  if length(dsetPList) ~= length(idx2)
    error('Supplied paramList does not contain all the dataset parameters')
  end
  if size(x,1)~=size(paramList,1)
    error('x must be n-by-Nevals, where n is the length of the supplied paramList')
  end

  %get x in the correct order for the dataset.
  tmp(idx2,:) = x(idx1,:);
  x = tmp;
  %overwrite the supplied paramList to correspond to the modified x
  paramList = dsetPList;
end

rng = vertcat(D.FreeParameter.range);
if ~DClab.issubset([min(x,[],2) max(x,[],2)],rng,1e-10,1e-10)
  error('Usage: X is not contained in the domain  described by the DCDataset''s FreeParameters.')
end

%Initialize output
N = size(x,2);
m = length(pairs);
y = zeros(m,N);
yInt = zeros(2*m,N);

for i1 = 1:m
  Pair = D.ModelAndObservationPair(i1);
  [y(i1,:) yInt(2*i1-1:2*i1,:)] = eval(Pair.ResponseModel,x,paramList);
end
