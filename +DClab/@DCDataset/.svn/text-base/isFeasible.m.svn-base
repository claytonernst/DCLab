function bool = isFeasible(Dset,x)
%ISFEASIBLE determines if a parameter vector lies in the feasible set described by a dataset.
%
%   BOOL = ISFEASIBLE(DSET,X) returns true if the n-by-1 parameter vector X
%   lies in the feasible set of the DCDataset DSET. This involves evaluating
%   the true models at X, since surrogates are only formed when we convert
%   a DCDataset to a PolyDataset. 
%
%   See also PolyDataset/isSurrogateFeasible

%TODO we should vectorize this method.

m = Dset.nPairs;
yIntModel = zeros(m,2);
yIntData = zeros(m,2);

% This call will take care of error checking on x.
[trash tmpYInt] = evalResponseModels(Dset,x);

for i1 = 1:m
  yIntModel(i1,:) = tmpYInt(2*i1-1:2*i1)';
  Pairs = Dset.ModelAndObservationPair;
  d = Pairs(i1).observedValue;
  uVect = Pairs(i1).observationUncertaintyPlusMinus;
  uncCase = Pairs(i1).ResponseObservation.uncertaintyCase;
  
  if uncCase == 1
    yIntData(i1,:) = d+uVect;
  elseif uncCase == 2
    yIntData(i1,:) = d*(1+uVect);
  elseif uncCase == 3
    yIntData(i1,:) = d*10.^uVect;    
  elseif uncCase == 4
    yIntData(i1,:) = d.^(1+uVect);
  else
    error('Internal inconsistency, condition should never occur')
  end
end

%Check if all of these intervals overlap. 
if doesIntersect(yIntData,yIntModel,1e-10,1e-10)
  bool = true;
else
  bool = false;
end

function bool = doesIntersect(set1,set2,relTol,absTol)
%DOESINTERSECT determines if two n-D rectangles intersect
%
%   BOOL = DOESINTERSECT(SET1,SET2) returns true of the n-by-2 matrix SET1
%   describes a rectangle that intersects the rectangle described by
%   the n-by-2 matrix SET2.
%
%   BOOL = ISSUBSET(SET1,SET2,RELTOL) allows the relative slop RELTOL in
%   the intersection. It is implemented by first expanding both rectangles
%   by RELTOL, and then checking intersection. For example it SET1 = [0 1]
%   and reltol = 0.1. We will convert SET1 to [-0.1 1.1] before checking
%   intersection.
%
%   BOOL = ISSUBSET(SET1,SET2,RELTOL,ABSTOL) allows an additional
%   absolution tolerance when determining the intersection. This feature is
%   again implemented by expanding both rectangles.
%
%   The toolbox uses this function because in some situations we don't want
%   to check intersection exactly due to roundoff error, 10^log10(x) ~= x,
%   etc.


ni = nargin;
switch ni
  case 2
    relTol = 0;
    absTol = 0;
  case 3
    absTol = 0;
  otherwise
    error(nargchk(2,4,ni))
end

if ~isequal(size(set1),size(set2)) || size(set1,2) ~= 2
  error('Inputs: SET1 and SET1 must be n-by-2 matrices')
end

if any(set1(:,1) > set1(:,2))
  error('Inputs: SET1 does not describe a nonempty rectangle (lower bound greater than upper bound)');
end
if any(set2(:,1) > set2(:,2))
  error('Inputs: SET2 does not describe a nonempty rectangle (lower bound greater than upper bound)');
end
  
if relTol < 0 || absTol < 0
  error('Inputs: RELTOL and ABSTOL must be nonnegative')
end



n = size(set1,1);

%Grow the two sets.
set1 = set1 + [-relTol*diff(set1,[],2) relTol*diff(set1,[],2)] + repmat([-absTol absTol],n,1);
set2 = set2 + [-relTol*diff(set2,[],2) relTol*diff(set2,[],2)] + repmat([-absTol absTol],n,1);

dimIntersects = false(n,1);
for i1 = 1:n
  % Given [a b] and [c d]. Intersection occurs in c\in [a b], d\in [a b],
  % or [a b] \subset [c d]. set1 = [a b], set2 = [c d]
  
  if (set1(i1,1) <= set2(i1,1) && set2(i1,1) <= set1(i1,2)) || ...
     (set1(i1,1) <= set2(i1,2) && set2(i1,2) <= set1(i1,2)) || ...
     (set2(i1,1) <= set1(i1,1) && set1(i1,2) <= set2(i1,2))
     dimIntersects(i1) = true;
  end
end
  
if all(dimIntersects)
  bool = true;
else
  bool = false;
end





