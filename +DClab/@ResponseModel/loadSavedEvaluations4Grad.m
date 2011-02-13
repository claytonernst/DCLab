function [x, y, filesUsed, filesQueried] = ...
    loadSavedEvaluations4Grad(RM, domrng, varTrans, Nrequested, filesToSkip)
%LOADSAVEDEVALUATIONS returns saved model evaluations from
%DClabVxxx/savedEvaluations for gradient estimates
%
%   [PARAMVALUES RESPONSE] = LOADSAVEDEVALUATIONS4GRAD(RM,DOMAINRANGES)
%   considers all saved evaluations y = M(x) of the ResponseModel RM and
%   returns those pairs (x,y) where x lies in the rectangle described by
%   the nx2 matrix DOMAINRANGES. PARAMVALUES contains horizontally
%   concatenated column vectors of x and will be n-by-nevals. RESPONSE
%   contains horizontally concatenated values of y and will be 1-by-nevals.
%
%   [PARAMVALUES RESPONSE] = LOADSAVEDEVALUATIONS4GRAD(RM,DOMAINRANGES,VARTRANS)
%   imposes the additional condition on x that it must have been generated
%   from a uniform distribution over VARTRANS(DOMAINRANGES). VARTRANS
%   should be a nx1 cell array containing 'none' or 'log10'. For example,
%   if DOMAINRANGES = [1 10; -1 1]] and VARTRANS = {'log10';'none'} the x
%   of all returned pairs will be [10^zeta(1),; zeta(2)] where zeta(1) was
%   drawn from a uniform distribution over [0 1], and zeta(2) was drawn
%   from a uniform distribution over [-1 1].
%
%   [PARAMVALUES RESPONSE] = LOADSAVEDEVALUATIONS4GRAD(RM,DOMAINRANGES,VARTRANS,N)
%   returns at most N (x,y) pairs. If more than N are available, the first
%   N found in the files are returned.
%
%   [PARAMVALUES RESPONSE] = LOADSAVEDEVALUATIONS4GRAD(RM,DOMAINRANGES,VARTRANS,N,FILESTOSKIP)
%   does not return any pairs saved in the files indicated by FILESTOSKIP.
%  
%   [PARAMVALUES RESPONSE FILESQUERIED] = LOADSAVEDEVALUATIONS4GRAD(...)
%   additionally returns a row vector listing the numbers of the files that
%   were considered during the operation. These are the files (not
%   including those listed in FILESTOSKIP) that had to be examined to
%   generate N (x,y) pairs. This output may be useful at a later stage
%   (e.g., validation) because it lists files whose contents are either not
%   useful, or were already used in some way.
%
%   [PARAMVALUES RESPONSE FILESQUERIED FILESUSED] = LOADSAVEDEVALUATIONS4GRAD(...)
%   additionally returns a row vector listing the numbers of the files that
%   the (x,y) pairs of the output were drawn from.

%initialize outputs
x = [];
y = [];
filesQueried = [];
filesUsed = [];

% If there isn't the possibility of saved evaluations, make a quick exit
if ~strcmp(RM.type,'dcModel') || RM.saveEnabled~=1
  return 
end

% Determine the directory name the evaluations would be saved in
saveDirectory = fileparts(which('savedEvaluationsDir'));
dirName = fullfile(saveDirectory,func2str(RM.model)); %TODO should we just use the name property?
for i1 = 1:length(RM.additionalInputs)
  if ischar(RM.additionalInputs{i1});
    dirName = [dirName '_' RM.additionalInputs{i1}]; %#ok
  elseif isnumeric(RM.additionalInputs{i1});
    dirName = [dirName '_' num2str(RM.additionalInputs{i1})]; %#ok
  else
    error('Internal inconsistency: condition should never occur')
  end
end

dirName = strrep(dirName,'.','p');

% If there is no directory of this name, again make a quick exit
if ~exist(dirName,'dir')
  return
end

% At this point there is a directory that make contain useful evalutations.

ni = nargin;
switch ni
  case 2
    varTrans = []; 
    Nrequested = inf;
    filesToSkip = [];
  case 3
    Nrequested = inf;
    filesToSkip = [];
  case 4
    filesToSkip = [];      
  otherwise
    error(nargchk(2,5,ni))
end

if ~isempty(varTrans)
  noTransDims = strmatch('none',char(varTrans),'exact');
  log10Dims = strmatch('log10',char(varTrans),'exact');
  if any(domrng(log10Dims,1) < 1e-1000 )
    error('log10 transformation requested on one or more dimensions whose lower limit was nonpositive')
  end
else
  %These are really the dimesions that must be subset rather than
  %exact. Pardon the abuse of variable names, but we want this
  %variable below.
  noTransDims = (1:size(domrng,1))'; 
  log10Dims = [];
end

% Load the catelog file. This is a .mat file having the same name as the
% directory. It contains a structure array named 'a'. 'a' has fields
% .domainRange and .variableTransformations.

[path name] = fileparts(dirName); %#ok
%robustLoad(fullfile(dirName,[name '.mat']),0);
s = load(fullfile(dirName,[name '.mat']));
if isfield(s,'a')
  catelog = s.a;
else
  error(['Directory ' dirName ' failed to contain ' dirName '.mat'])
end

% Determine which files we're allowed to look at.
potentialFiles = setdiff(1:length(catelog),filesToSkip);

% Determine which files have points in the gradient bunching:
potentialFiles = potentialFiles([catelog(potentialFiles).forGrad]);

% if empty, set potentialFiles to -1 so we don't enter the while loop
%if isempty(potentialFiles)
%  potentialFiles = -1;
%end

% Loop through these files until we're done.
i1 = 1;
while length(y) < Nrequested && i1<=length(potentialFiles)
  fileIdx = potentialFiles(i1);
  filesQueried = [filesQueried; fileIdx];
  
  % If varTrans is nonempty, only used files for which
  % catelog(i).variableTransformations matches varTrans
  if isempty(varTrans) || isequal(varTrans,catelog(fileIdx).variableTransformations)
    
    % Try to load data from fileIdx.mat
    try

      %See if fileIdx.mat contains any useful points. For dimenions that
      %were log10 transformed, domrng and catelog(fileIdx).domainRange must
      %correspond exactly. For dimensions that were untransformed,
      %catelog(idx).domainRange must only contain domrng. We allow so
      %relative and absolute slack in these requirements to account for
      %roundoff. Since x values are not saved in a normalized form,
      %roundoff can occur, especially if a dimension of the parameter
      %domain is singleton. 
      
      absTol = 0.1; %much bigger because we don't care for those near edge bits
      relTol = 1e-10;
      noTransOK = DClab.issubset(domrng(noTransDims,:),catelog(fileIdx).domainRange(noTransDims,:),relTol,absTol);
      log10TransOK = DClab.issubset(domrng(log10Dims,:),catelog(fileIdx).domainRange(log10Dims,:),relTol,absTol);
        
      % If a file contains useful data, load it up.

      if noTransOK && log10TransOK

        %robustLoad([fullfile(dirName,num2str(i1)) '.mat'],0);
        s = load([fullfile(dirName,num2str(fileIdx)) '.mat']);
        loadedy = s.a.responses;
        loadedx = s.a.paramValues;

        %Find the x values that are inside the current domain of interest.
        %Allow the same tolerances as before.
        %We only need to consider the noTrans dimensions. The log10Trans
        %dimensions should correspond exactly.
        nfound = length(loadedy);
        
        lowerLimit = zeros(size(domrng,1),1);
        upperLimit = zeros(size(domrng,1),1);

        lowerLimit(noTransDims) = domrng(noTransDims,1) - relTol*diff(domrng(noTransDims,:),[],2) - absTol;
        upperLimit(noTransDims) = domrng(noTransDims,2) + relTol*diff(domrng(noTransDims,:),[],2) + absTol;
        
        lowerLimit(log10Dims) = domrng(log10Dims,1) - relTol*diff(domrng(log10Dims,:),[],2) - absTol;
        upperLimit(log10Dims) = domrng(log10Dims,2) + relTol*diff(domrng(log10Dims,:),[],2) + absTol;
        
        % Use logical indexing
        %keepIdx = all(repmat(lowerLimit,1,nfound) <= loadedx,1) & ...
        %  all(repmat(upperLimit,1,nfound) >= loadedx,1);
        keepIdx = true(1,nfound); %just load 'em all

        y = [y loadedy(keepIdx)]; %#ok
        x = [x loadedx(:,keepIdx)]; %#ok
        if any(keepIdx)
          filesUsed = [filesUsed; fileIdx];
        end
      end  
    catch
      disp(['error loading ' fullfile(dirName,num2str(fileIdx)) ': skipping']);
      %keyboard
    end
  end % if isempty(varTrans) || isequal(varTrans,catelog(idx).variableTransformations)

  i1 = i1+1;
end %while

% Truncate off any extras
if length(y) > Nrequested
  x = x(:,1:Nrequested);
  y = y(1:Nrequested);
end
