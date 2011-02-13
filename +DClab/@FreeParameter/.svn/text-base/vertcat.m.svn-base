function [newFP erased] = vertcat(varargin)
%VERTCAT Vertical concatenation of FreeParameter objects.
%
%   FP = VERTCAT(FP1,FP2,...) or FP = [FP1; FP2; ...]
%   contatenates the input objects to form the multidimension output FP.
%   The components of FP will never share a common name. See "Behavior"
%   below.
%
%   [FP ERASED] = VERTCAT(FP1,FP2,...) additionally returns a vector
%   of indices indicating which components of  
%   tmp = builtin('vertcat',varargin{:}), which will not strip out
%   free parameters with the same name, were included in FP. This output
%   may be useful for diagnosing any errors. 
%
%   Behavior:
%
%   If an "Incomplete" FreeParameter is encountered, an error with be
%   thrown, as this condition can create difficulties/ambiguities with the
%   requirement that the output of this method be an object with no shared
%   names. "Empty" FreeParameters are not included in the final output.
%
%   Each input is potentially a multidimensional vector object
%   consisting of several free parameters. Consider all individual free
%   parameters present in the input list. If multiple free parameters units
%   share an identical name, the code will check if they are identical in
%   all respects. If so, only one copy will be present in the final output.
%   If nonidentical free parameters with a common name are encountered, a
%   warning will be issued indicating that free parameters of the same name
%   have been "combined". The ranges of such duplicates will be intersected
%   to produce a free parameter. Its nominal will be set to the mean of the
%   intersected ranges, its uncertaintyType will be set to 'absolute' and
%   its uncertaintyTransformation will be set to 'none'. All other fields
%   will inherit the properties of the earliest likenamed free parameter in
%   the input list.
% 
%   See also FreeParameter

error(nargoutchk(0,2,nargout));

% varargin must contain FreeParameters objects or []. Check this, then
% remove all empties from varargin.

%Use logical indexing.
emptyLogical = cellfun(@isempty,varargin);
anon = @(x) isa(x,'double');
doubleLogical = cellfun(anon,varargin);
nativeEmpty = emptyLogical & doubleLogical;

anon = @(x) isa(x,'DClab.FreeParameter');
FPLogical = cellfun(anon,varargin);
if ~all(FPLogical | nativeEmpty)
  error('Usage: concatentation of nonidentical object types') 
end
varargin(emptyLogical) = [];

% Correct for the special case when all inputs were empty.
if isempty(varargin)
  varargin = {DClab.FreeParameter};
end

%At this point varargin is a cell array with each cell containing a
%FreeParameter object. Combine these with MATLAB's struct vertcat
%method, and then check for duplicate names.

newFP = builtin('vertcat',varargin{:});

%if there are any incompletes, just bomb.  We need to remove empties
%in all cases.  Then we need to check for duplicates.  If they occur,
%we have some work to do...

if any(vertcat(newFP.status)==0)
  error('Inputs: an incomplete object was supplied to FREEPARAMETER/VERTCAT')
end

%Check for duplicates (as indicated by shared names)
plist = {newFP.name}';
pplist = unique(plist);

if length(pplist)~= length(plist)
  %We have duplicates and need to combine them. It seems most efficient
  %to look at all the given assertions with the same name, and
  %intersect ranges. 
  forcedCombination = false;
  erased = [];
  for i1 = 1:length(pplist)
    idx = strmatch(pplist(i1),plist,'exact');
    if length(idx) > 1
      
      %Check if all components with name ppllist(i1) are identical
      allsame = true;
      for i2 = 2:length(idx)
        if ~isequal(newFP(idx(1)),newFP(idx(i2)))
          forcedCombination = true;
          allsame = false;
        end
      end

      if ~allsame
        %all assertions with a given name were not idential
        %intersect ranges, etc.
        tmp = newFP(idx);
        bnds = tmp.bounds;
        newlb = max(bnds(:,1));
        newub = min(bnds(:,2));
        if newlb >= newub
          error(['Intersecting the range of all components with name ' tmp(1).name ' resulted in an empty domain'])
        end
        newnom = mean([newlb newub]);
        newFP(idx(1)).nominal = newnom;
        newFP(idx(1)).uncertainty = [newnom-newlb newubnewnom];
        newFP(idx(1)).uncertaintyType = 'absolute';
        newFP(idx(1)).uncertaintyTransformation = 'none';

      end
        
      erased = [erased; idx(2:end)];
    end
  end
  newFP(erased) = [];

  if forcedCombination
    warning off backtrace

    %Call warning with an empty as the 2nd input to force it to display a
    %formatted string.
    warning(['Concatenating parameters assertions with identical names \n '...
             'but other differences, ranges of componenets with duplicate \n '...
             'names will be intersted. The nominal will be the mean. The \n '...
             'uncertaintyType was set to ''absolute'' and the \n'...
             'uncertaintyTransformation was set to ''none''. Other object \n'...
             'properties will be set to those of the first seen with a \n'...
             'given name. Consequently some information may be lost!'],[]); %2nd input is [] so string will be formated
    warning on backtrace
  end

end





