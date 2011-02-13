function bool = issubset(set1,set2,relTol,absTol)
%ISSUBSET determines if one n-D rectangle is contained in another
%
%   BOOL = ISSUBSET(SET1,SET2) returns true of the nx2 matrix SET1
%   describes a rectangle that is contained in the rectangle described by
%   the nx2 matrix SET2.
%
%   BOOL = ISSUBSET(SET1,SET2,RELTOL) allows the relative slop RELTOL in
%   the containment. I.e., if SET1 = [0.9 11.1], SET2 = [1,11], and RELTOL
%   = 0.01; This function will return true since 0.9 >= 1-10*0.01 and 
%   11.1 <= 11+10*0.01
%
%   BOOL = ISSUBSET(SET1,SET2,RELTOL,ABSTOL) allows an additional
%   absolution tolerance when determining the containment.
%
%   The toolbox uses this function because in some situations we don't want
%   to check containment exactly due to roundoff error, 10^log10(x) ~= x,
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

rel = relTol*diff(set2,[],2);

lowerOK = all(set1(:,1) >= set2(:,1) - rel - absTol);
upperOK = all(set1(:,2) <= set2(:,2) + rel + absTol);

if lowerOK && upperOK
  bool = true;
else
  bool = false;
end
