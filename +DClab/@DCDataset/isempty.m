function bool = isempty(obj)
% bool = isempty(obj)
%
% A DCDataset object is empty if it contains both and empty
% ModelAndObservationPair object and an empty FreeParameter object.
%
% See also DCDataset

if isempty(obj.ModelAndObservationPair) && isempty(obj.FreeParameter)
  bool=true;
else
  bool = false;
end
