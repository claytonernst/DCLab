function domstruct = createDomainStructure(paramList,range)
%CREATEDOMAINSTRUCTURE combines parameter list and matrix of their ranges into a structure
%
%   DOMSTRUCT = CREATEDOMAINSTRUCTURE(PARAMLIST,RANGE) accepts the n-by-1 cell array of unique
%   chars PARAMLIST and the n-by-2 matrix of lower and upper bounds on their respective
%   ranges. The output DOMSTRUCT is an n-by-1 struct array with fields .name and .range. 

if ~iscellstr(paramList)
  error('Inputs: PARAMLIST must be a cell array of chars')
end
[n m] = size(paramList);
if m~= 1
  error('Inputs: PARAMLIST must be a column cell array')
end
if ~isequal(size(range),[n 2])
  error('Inputs: RANGE must be a n-by-2 matrix')  
end
if any(range(:,1) > range(:,2))
   error('Inputs: RANGE does not describe a nonempty rectangle (lower bound greater than upper bound)');
end

if n == 0
  range = {};
else
  range = mat2cell(range,ones(n,1),2);
end
domstruct = struct('name',paramList,'range',range);
