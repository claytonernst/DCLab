function [newD erasedPairs erasedParams] = merge(varargin)
%DCDATASET/MERGE combines multiple Dataset objects into a single object.
%
%   NEWD = MERGE(D1,D2,...) Combines the ModelAndObservationPairs and
%   FreeParameters of the input DCDatasets D1,D2,... to form NEWD. If
%   ModelAndObservationPairs or FreeParameters with duplicate names are
%   present in  the corresponding fields of objects in the input list, they
%   will be combined or eliminated as appropriate. See
%   ModelAndObservationPair/vertcat and FreeParameter/vertcat for a more
%   complete descriptions of this behavior. 
%
%   [NEWD ERASEDPAIRS] = MERGE(D1,D2,...) additionally returns the indices
%   of which elements of builtin('vertcat',D1.MOPair,D2.MOPair,...) were
%   eliminated (due to name duplication) to form NEWD.MOPair. This output
%   is used to properly eliminate the corresponding DCSurfaces when vertcat
%   is invoked by PolyDataset/merge.
%
%   [NEWD ERASEDPAIRS ERASEDPARAMS] = MERGE(D1,D2,...) additionally returns the indices
%   of which elements of builtin('vertcat',D1.FreePara,D2.FreePara,...) were
%   eliminated to form NEWD.FreeParameter. 
%
% See also DCDataset, ModelAndObservationPair/vertcat, FreeParameter/vertcat 

error(nargoutchk(0,3,nargout));

%remove any empty inputs from the varargin list
crap = false(1,length(varargin));
for i1 = 1:length(varargin)
    crap(i1) = builtin('isempty',varargin{i1});
end
varargin(crap) = [];

%Don't call recursively, just use the ModelAndObservationPair and
%FreeParameter vertcat.
p = length(varargin);
MOP = cell(p,1);
FP = cell(p,1);
for i1 = 1:p
  if ~isa(varargin{i1},'DClab.DCDataset')
      error('Illegal concatentation of nonidentical object types')
  end
  MOP{i1} = varargin{i1}.ModelAndObservationPair;
  FP{i1} = varargin{i1}.FreeParameter;
end
[Pairs erasedPairs] = vertcat(MOP{:});
[Params erasedParams] = vertcat(FP{:});
newD = DClab.DCDataset(Pairs,Params);

