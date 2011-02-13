function [y yInt] = eval(Surf,domain,x)
%EVAL method evaluates a DCSurface object
%
%   Y = EVAL(SURF,DOMAIN,X) evaluates the DCSurface SURF which was fit
%   over the rectangular domain described DOMAIN at the parameter vectors
%   X. DOMAIN is an n-by-1 structure array with fields .name and .range.
%   Because a DCSurface doen't know where it came from (to save memory),
%   DOMAIN must be identical to the domain that was passed as the 2nd input
%   to the DCSurface constructor. X should be a n-by-Nevals matrix of
%   horizontally concatenated parameter vectors. The returned Y
%   is 1-by-Nevals.
%
%   [Y YINT] = EVAL(...) additionally returns a 2-by-Nevals matrix
%   YINT that is determined by adding the estimated fitting error of the
%   DCSurface to Y. YINT thus describes Nevals intervals that ought to
%   contain the corresponding values of M(x), where M is the model
%   approximated by the DCSurface. 
%
% See also DCSurface

% TODO this needs to work with vector objects

if nargin ~= 3
  error(nargchk(3,3,ni))
end

% Check the supplied domain
[bool message] = DClab.isValidDomain(domain,'DOMAIN');
if ~isempty(message)
  error(message);
end

% Check the dimensions of x
if size(x,1) ~= size(domain,1)
  error('Usage: dimensions of X and DOMAIN incompatible')
end

%Check that x is containing in the supplied domain.
domrng = vertcat(domain.range);
if ~DClab.issubset([min(x,[],2) max(x,[],2)],domrng,1e-10,1e-10)
  error('Usage: X is not contained in the supplied DOMAIN')
end

if any([Surf.status]~=1)
  error('Inputs: EVAL can only be called with "complete" (status==1) objects')
end

ns = nSurfaces(Surf);
y = zeros(ns,size(x,2));
yInt = zeros(2*ns,size(x,2));

for i1 = 1:ns

  xactive = x(Surf(i1).activeParameterIndex,:);

  % Use logical indexing
  log10Dim = Surf(i1).activeParameterTransformation == 2;
  noTransDim = Surf(i1).activeParameterTransformation == 1;

  % Recall that by design, for a parameter to be active,
  % bnds(i,2)-bnds(i,1) > 0
  actdomrng = domrng(Surf(i1).activeParameterIndex,:);
  if any( actdomrng(log10Dim,1) < 1e-10 )
    error('log10 transformation occured on one or more dimensions whose lower limit was non positive')
  end

  Neval = size(x,2);

  % Scale x to be between +/- 1. Since the range of all active parameters
  % must be nonsingleton, we don't (or at least shouldn't) have to worry
  % about dividing by zero.
  scale = diag(1./diff(actdomrng(noTransDim,:),1,2));
  xactive(noTransDim,:) = 2*scale*(xactive(noTransDim,:) - repmat(actdomrng(noTransDim,1),1,Neval)) - 1;

  scale = diag(1./log10(actdomrng(log10Dim,2)./actdomrng(log10Dim,1)));
  xactive(log10Dim,:) = 2*scale*log10(xactive(log10Dim,:)./repmat(actdomrng(log10Dim,1),1,Neval)) - 1;

  xx = [ones(1,size(xactive,2)); xactive];

  if isstruct(Surf(i1).surrogateModel)
    y(i1,:) = sum(xx.*(Surf(i1).surrogateModel.num*xx),1)./sum(xx.*(Surf(i1).surrogateModel.den*xx),1);
  else
    y(i1,:) = sum(xx.*(Surf(i1).surrogateModel*xx),1);
  end

  % peakErr(1) <= S(x) - M(x) <= peakErr(2).
  % Thus S(x) - peakErr(2) <= M(x) <= S(x) - peakErr(1)
  
  peakErr = Surf(i1).peakError;
  yInt(2*i1-1,:) = y(i1,:) - peakErr(2);
  yInt(2*i1,:) = y(i1,:) - peakErr(1);

  if strcmp(Surf(i1).responseTransformation,'log10');
    y(i1,:) = 10.^y(i1,:);
    yInt(2*i1-1:2*i1,:) = 10.^yInt(2*i1-1:2*i1,:);
  end
end