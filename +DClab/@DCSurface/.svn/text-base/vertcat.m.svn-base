function SurfObj = vertcat(varargin)
% DCSURFACE/VERTCAT implements concatenation of DCSurface
% objects. "Empty" DCSurfaces are not included in the final output.
%
% Syntax:
%   newSurfObj = [SurfObj1; SurfObj2; ... ; SurfObjn]
%
% Inputs:
%   Each input must be a DCSurface object. 
%
% See also DCSurface

error(nargoutchk(0,1,nargout));

% varargin must contain DCSurface objects or []. Check this, then
% remove all empties from varargin.

%Use logical indexing.
emptyLogical = cellfun(@isempty,varargin);
anon = @(x) isa(x,'double');
doubleLogical = cellfun(anon,varargin);
nativeEmpty = emptyLogical & doubleLogical;

anon = @(x) isa(x,'DClab.DCSurface');
DCSLogical = cellfun(anon,varargin);
if ~all(DCSLogical | nativeEmpty)
  error('Usage: concatentation of nonidentical object types') 
end
varargin(emptyLogical) = [];

% Correct for the special case when all inputs were empty.
if isempty(varargin)
  varargin = {DClab.DCSurface};
end

%At this point varargin is a cell array with each cell containing a
%DCSurface object. Combine these with MATLAB's struct vertcat
%method.

SurfObj = builtin('vertcat',varargin{:});





