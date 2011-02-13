function newSurfArray = changeDomain(surfArray,oldDomain,newDomain)
%CHANGEDOMAIN alters the domain a DCSurface is defined on
%
%   NEWSURF = CHANGEDOMAIN(SURF,OLDDOMAIN,NEWDOMAIN) takes the (possibly
%   multidimensional) DCSurface object SURF and performs scaling operations
%   on its algebraic surrogate models to produce a new object NEWSURF that
%   is defined over a different domain. Both OLDDOMAIN, and NEWDOMAIN
%   should be column structure arrays with fields .name and .range.
%   OLDDOMAIN must be identical to the domain that was passed as the 2nd
%   input to the DCSurface constructor. NEWDOMAIN can add additional
%   dimensions not present in OLDDOMAIN or it can lack nonactive dimensions
%   that are present in OLDDOMAIN.
%
%   LIMITATIONS: The method will not work on "incomplete" objects, unless
%   the object is empty.
%
%   This mehtod is necessary because the algebraic form stored in the
%   object is a function of normalized (between +/-1) variables. When the
%   domain is altered, a different normalization is required so the form
%   stored in the object must be rescaled. The method also allows you to
%   add or delete dimensions to the domain.
%
%   USE WITH CARE, as this function will extrapolate if NEWDOMAIN is not a
%   subset of OLDDOMAIN. TODO: should we disallow this?
%
% See also DCSurface, PiecewiseSurrogateModelTree/changeDomain

% Description:
%
% Let p be the parameter vector; let dom1 be the original parameter
% domain; let dom2 be the new parameter domain.
%
% let x be the coded variables of the original domain, and z the coded
% variables of the new domain. Since x and y are both generated from the
% natural parameters p by an affine transformation. They are related by 
% x = my+b;
%
% dom1: p = 0.5*(x+1)*diff(oldBnds) + oldBnds(:,1);
%
% dom2: p = 0.5*(y+1)*diff(newBnds) + newBnds(:,1);
%
% Thus
% m = diff(newBnds)./diff(oldBnds);
% b = inv(diff(oldBnds))*(diff(newBnds) + 2*newBnds(:,1) - 2*oldBnds(:,2) - diff(oldBnds))


%Check the input domains.
[bool message] = DClab.isValidDomain(oldDomain,'OLDDOMAIN');
if ~isempty(message)
  error(message);
end
[bool message] = DClab.isValidDomain(newDomain,'NEWDOMAIN');
if ~isempty(message)
  error(message);
end

% Special case
if isempty(surfArray)
  newSurfArray = surfArray;
  return
end

% Check if any of the surfaces are incomplete. If so, error.
if any([surfArray.status]~=1)
  error('Inputs: CHANGEDOMAIN cannot be called with a DCSurface having any incomplete elements');
end

newSurfArray = surfArray;
for i1 = 1:length(surfArray)
  
  surf = surfArray(i1);
  newsurf = surf;
  
  newNames = {newDomain.name}';
  oldNames = {oldDomain.name}';
  activeParams = {oldDomain(surf.activeParameterIndex).name}';
  if ~isempty(setdiff(activeParams,newNames))
    error('Usage: newDomain must contain all active parameters')
  end

  % Update active Parameter list. oldNames2NewNames is a n-by-1 vector,
  % where n = length(oldDomain). oldDomain(oldNames2NewNames(i)) should be
  % newDomain(i). I.e., the ith componenent of oldNames2NewNames indicates
  % the location in newDomain of the ith parameter in oldDomain. Some
  % components of oldNames2NewNames may be NaN, since newNames doesn't need
  % to be a superset of oldNames. However, any active parameters must be
  % included in newNames.
  oldNames2NewNames = NaN(size(oldDomain));
  
  [trash newNames2sortCommonNames oldNames2sortCommonNames] = intersect(newNames,oldNames); %#ok
  oldNames2NewNames(oldNames2sortCommonNames) = newNames2sortCommonNames;

  newActiveIdxOldOrder = oldNames2NewNames(surf.activeParameterIndex);
  [newActiveIdxNewOrder old2new] = sort(newActiveIdxOldOrder);
  newsurf.activeParameterIndex = [newActiveIdxNewOrder surf.activeParameterTransformation(old2new)];

  % Make the range estimate null
  newsurf.trueModelRangeEstimate = [];

  %TODO is the range estimate ever required by our software? We should be
  %able to just set it to empty. But why?
  
  % Re-sort the polynomials to have the proper order
  poly = surf.surrogateModel; 
  if isstruct(poly)
    newpoly.num = poly.num([1; old2new+1],[1; old2new+1]);
    newpoly.den = poly.den([1; old2new+1],[1; old2new+1]);
  else
    newpoly = poly([1; old2new+1],[1; old2new+1]);
  end
  newsurf.surrogateModel = newpoly;

  % Now rescale the surfaces
  newBnds = vertcat(newDomain.range);
  oldBnds = vertcat(oldDomain.range);

  oldBndsLocal = oldBnds(surf.activeParameterIndex(old2new),:);
  newBndsLocal = newBnds(newsurf.activeParameterIndex,:);

  % Use logical indexing
  log10Dim = newsurf.activeParameterTransformation == 2;

  % Check that no log10 dimensions have a negative lower bounds
  if any(newBndsLocal(log10Dim,1)< 1e-10)
    error('log10 transformation occured on one or more dimensions whose lower limit was non positive')
  end
  
  % Transform the bounds of the log10 dimenions.
  oldBndsLocal(log10Dim,:) = log10(oldBndsLocal(log10Dim,:));
  newBndsLocal(log10Dim,:) = log10(newBndsLocal(log10Dim,:));
  
  if any(diff(newBndsLocal,[],2) == 0)
    
    %TODO fix this
    disp('finish your code, you need to eliminate an active parameter cuz its domain is constant')
    keyboard
    % you need to 
    %1) find where the "constant" value fit in the old interval. 
    %2) from the old scale, remove this parameters by setting it equal to
    %this constant value
    %3) update newBndsLocal
  end

  m = diff(newBndsLocal,[],2)./diff(oldBndsLocal,[],2);
  b = (diff(newBndsLocal,[],2) + 2*newBndsLocal(:,1) - 2*oldBndsLocal(:,1) - diff(oldBndsLocal,[],2))./diff(oldBndsLocal,[],2);

  %trans = diag([1; m]);
  %trans(2:end,1) = b;
  
  if isstruct(newsurf.surrogateModel)
    newsurf.surrogateModel.den = DClab.DCSurface.composeQuadWithAffine(newsurf.surrogateModel.den,diag(m),b);
    newsurf.surrogateModel.num = DClab.DCSurface.composeQuadWithAffine(newsurf.surrogateModel.num,diag(m),b);
  else
    newsurf.surrogateModel = DClab.DCSurface.composeQuadWithAffine(newsurf.surrogateModel,diag(m),b);
  end

  newSurfArray(i1,1) = newsurf;
end

