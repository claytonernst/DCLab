function [x, y] = randomEval(MM, domrng, varTrans, N)
%
%  %%% Incorrect help %%%
%
%RANDOMEVALWITHSAVE returns (x,y) pairs for a given ResponseModel
%
%   [PARAMVALUES RESPONSE] = RANDOMEVALWITHSAVE(RM,DOMAINRANGES,VARTRANS,N)
%   evalulates the ResponseModel RM at N locations in the set described by
%   the n-by-2 matrix DOMAINRANGES. The N x-values are generated by a latin
%   hypercube design which emulates a uniform distribution over
%   VARTRANS(DOMAINRANGES). VARTRANS should be a nx1 cell array containing
%   'none' or 'log10'. For example, if DOMAINRANGES = [1 10; -1 1]] and
%   VARTRANS = {'log10';'none'} the x of all returned pairs will be
%   [10^zeta(1),; zeta(2)] where zeta(1) was drawn from a uniform
%   distribution over [0 1], and zeta(2) was drawn from a uniform
%   distribution over [-1 1]. PARAMVALUES contains horizontally
%   concatenated column vectors of x and will be n-by-nevals. RESPONSE
%   contains horizontally concatenated values of y and will be 1-by-nevals.
%   If RM.saveEnabled == true. The evaluation points (x,y) will be saved in
%   a subdirectory of DClabVxxx/savedEvalutations.
%
%   [PARAMVALUES RESPONSE] = RANDOMEVALWITHSAVE(RM,DOMAINRANGES,VARTRANS,N,NCOMP)
%   distributes the task of evaluating y = M(x) at the design points among
%   NCOMP different networked computers that share a common file server
%   with the local machine. Each of these computers must have the toolbox
%   installed and be running XXX. Set NCOMP = 0 to perform all evaluations
%   on the local machine. 
%
%   [PARAMVALUES RESPONSE FILECREATED] = RANDOMEVALWITHSAVE(...)
%   additionally returns that name (an integer) of the files in which the
%   (x,y) pairs were saved.
%
%   Inputs:
%   RM: ResponseModel object
%   DOMAINRANGES: Nparam-by-2 matrix indicating the domain from which to
%      select the design points.
%   VARTRANS: Nparam-by-1 cell array with each cell containing 'none' or
%      'log10'. This input dictates how the design points are distributed.
%   N: Number of design points
%   NCOMP: a nonnegative integer indicating the number of computers the
%      task of evaluating the model at the design points should be
%      distributed over. Set this to 0 to perform all computations on the
%      local machine that called this function.

ni = nargin;
error(nargchk(4,4,ni))

%Check that the model is defined on the requested domain.
modelDomain = MM(1).domainRange;
if ~DClab.issubset(domrng,modelDomain,1e-10,1e-10)
  error('The domain of the DCMetamodel does not contain the requested domain')
end

% Create latin hypercube design. But only consider dimensions whose range
% is nonsingleton
doesVaryIdx = find(diff(domrng,[],2) > 0);
nvary = length(doesVaryIdx);
%rand('state',sum(100*clock)) %reset rand state
if N >= 3*nvary
  tmp = DClab.lhsdesign(N,nvary,'smooth','off','criterion','none');
else
  tmp = rand(N,nvary);
end

% Initialize x to be N columns of the lower limits on the ranges. We
% replace the values for dimensions that are nonsingleton with values
% generated by the above design.
x = repmat(domrng(:,1),1,N);

log10Trans = strmatch('log10',char(varTrans(doesVaryIdx)),'exact');
noTrans = strmatch('none',char(varTrans(doesVaryIdx)),'exact');

if any(domrng(doesVaryIdx(log10Trans),1) < 1e-10 )
  error('log10 transformation occured on one or more dimensions whose lower limit was non positive')
end
  
scale = diag(diff(domrng(doesVaryIdx(noTrans),:),1,2));
x(doesVaryIdx(noTrans),:) = repmat(domrng(doesVaryIdx(noTrans),1),1,N) + scale*tmp(:,noTrans)';

scale = diag(log10(domrng(doesVaryIdx(log10Trans),2)./domrng(doesVaryIdx(log10Trans),1)));
x(doesVaryIdx(log10Trans),:) = 10.^(repmat(log10(domrng(doesVaryIdx(log10Trans),1)),1,N) + scale*tmp(:,log10Trans)');

%At this point x is n-by-nevals, in the model's natural coordinates, and
%the design is complete.
y = eval(MM,x);