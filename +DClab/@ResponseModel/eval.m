function [y yInt] = eval(RM,x,paramList)
%MODELEDRESPONCE/EVAL evaluates the model at supplied parameter vectors
%
%  Y = EVAL(RM,X) evaluates the model of the MODELEDRESPONSE object RM at
%  parameter vectors given taken from the columns of the Nparam-by-Nevals
%  matrix X. The parameter values in X must be in the order expected by the
%  model. I.e., Y(1) = M(X(1,1),...X(Nparam,1)). Y will be 1-by-Nevals.
%
%  Y = EVAL(RM,X,PARAMLIST) assumes the parameter values in each
%  columns of X are in the order indicated by the column cell array of
%  chars PARAMLIST. PARAMLIST must be a superset of the parameters needed
%  by the model. This function will then extract the needed rows of X and
%  reorder them as necessary before evaluating the model. 
%
%  [Y YINT] = EVAL(RM,X,...) additionally returns a 2-by-Nevals matrix
%  YINT that is determined by adding the appropriate components of the
%  outputUncertainty present in the object to Y. YINT thus describes Nevals
%  intervals that ought to contain the true value y(x) when M(x) exhibits
%  some bounded undermodeling. 
%
%  Inputs:
%    RM: A ResponseModel object
%    X: An Nparam-by-Nevals matrix supplying horizontally
%       concatenated parameter vectors at which to evaluate the model.
%    PARAMLIST[opt]: a Nparam-by-1 cell array of strings indicating
%       the parameter names that correspond to each row of X. This
%       is used to extract and order the parameters necessary for
%       the model.
%
%  Outputs:
%    Y: a 1-by-Nevals vector containing the model outputs at each
%      supplied input vector.
%    YINT: 2-by-Nevals matrix containing the endpoints of the
%      interval predicted to contain the true output. This is
%      obtained by added or subtracting the appropriate model
%      outputUncertainty.
%
%   See also ResponseModel

ni = nargin;
no = nargout;
if ni == 2
  paramList = RM.parameterList;
  if size(x,1) ~= length(paramList)
      error('Usage: X and RM.parameterList have incompatable dimensions')
  end
  modelx = x;
elseif ni == 3
  if size(x,1) ~= length(paramList)
      error('Usage: X and PARAMLIST have incompatable dimensions')
  end
  [sortedList idx1 idx2] = intersect(paramList,RM.parameterList); %#ok
  if ~isequal(sortedList,sort(RM.parameterList))
    error('Usage: paramList contains insufficient parameters for the ResponseModel')
  end
  modelx(idx2,:) = x(idx1,:);
else
  error(nargchk(2,3,ni))
end

%Check that the model is defined on the requested domain.
modelDomain = RM.domain;
modelDomain = vertcat(modelDomain.range);
if ~DClab.issubset([min(modelx,[],2) max(modelx,[],2)],modelDomain,1e-10,1e-10)
  error('Usage: X is not contained in the domain of the ResponseModel')
end

switch RM.type
 case 'dcModel'
  try
    y = RM.model('simulate',modelx,RM.additionalInputs{:});
  catch
    disp('Error evaluating dcModel in ResponseModel/eval')
    keyboard
  end

  %See if multiple responses is enabled. If so, then
  %the model may have returned a fat matrix rather than a 1xNevals row vector. 
  if RM.multipleResponsesEnabled
    %if enabled both featPrec and featList should be nonempty and
    %addlInputs should be at least length 1
    respSimSeg = RM.responseSimulationSegment;
    respList = RM.responseList;
    idx1 = strmatch(RM.additionalInputs{1},char(respList),'exact');
    if ~isempty(respSimSeg)
      computedResponses = respList(respSimSeg <= respSimSeg(idx1));
      idx2 = strmatch(RM.additionalInputs{1},char(computedResponses),'exact');
      y = y(idx2,:);
    else
      y = y(idx1,:);
    end
    
  else
    %do nothing, out should be a row vector.
  end

 case 'quadraticModel'
   log10Trans = strmatch('log10',char(RM.variableTransformations),'exact');
   modelx(log10Trans,:) = log10(modelx(log10Trans,:));
   xx = [ones(1,size(modelx,2)); modelx];
   y = sum(xx.*(RM.model*xx),1);
   if strcmp(RM.responseTransformation,'log10')
     y = 10.^y;
   end
 case 'linearModel'
    log10Trans = strmatch('log10',char(RM.variableTransformations),'exact');
    modelx(log10Trans,:) = log10(modelx(log10Trans,:));
    xx = [ones(1,size(modelx,2)); modelx];
    y = RM.model'*xx;
    if strcmp(RM.responseTransformation,'log10')
      y = 10.^y;
    end
  otherwise
  error('Internal Inconsistency: condition should never occur');
end

% Make sure the model doesn't violate any of the assumptions implied by the
% form of its output uncertainty.
uncCase = RM.uncertaintyCase;
[ymin idx] = min(y);
if (uncCase == 2 || uncCase == 3) && ymin <= 0
  error(['Invalid object: its output uncertainty form requires M(x) > 0; column ' num2str(idx) ' of X violated this'])
end
if uncCase == 4 && ymin <= 1
  error(['Invalid object: its output uncertainty form requires M(x) > 1; column ' num2str(idx) ' of X violated this'])
end

if no == 2
  uncVect = RM.outputUncertaintyPlusMinus'; %will be 1x2
  if uncCase == 1
    yInt = [y+uncVect(1); y+uncVect(2)];  
  elseif uncCase == 2
    yInt = [y*(1+uncVect(1)); y*(1+uncVect(2))]; 
  elseif uncCase== 3
    yInt = [y*10^uncVect(1); y*10^uncVect(2)]; 
  elseif uncCase == 4
    yInt = [y^(1+uncVect(1)); y^(1+uncVect(2))];
  else
    error('condition should never occur')
  end
end

  
  

