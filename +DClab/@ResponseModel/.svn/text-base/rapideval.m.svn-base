function y = rapideval(RM,x)
%MODELEDRESPONCE/RAPIDEVAL evaluates the model at supplied parameter vectors
%
%  This function does not check that the inputs are legal. See eval for a
%  more robust implementation. This function is used by
%  DCSurface/private/validate when we optimize the fitting error. The
%  standard eval is sufficient for most operations as the slow down from
%  input checking is minimal when everything is vectorized. validate is not
%  vectorized and additionally may try to used the model slightly outside
%  of its domain due to the tolerances/nature of the optimization.
%
%  Y = RAPIDEVAL(RM,X) evaluates the model of the MODELEDRESPONSE object MR
%  at parameter vectors taken from the columns of the Nparam-by-Nevals
%  matrix X. The parameter values in X must be in the order expected by the
%  model. I.e., Y(1) = M(X(1,1),...X(Nparam,1)). Y will be 1-by-Nevals.
%
%  Inputs:
%    RM: A ResponseModel object
%    X: An Nparam-by-Nevals matrix supplying horizontally
%       concatenated parameter vectors at which to evaluate the model.
%  Outputs:
%    Y: a 1-by-Nevals vector containing the model outputs at each
%      supplied input vector.
%
%   See also ResponseModel, ResponseModel/eval

modelx = x;
switch RM.type
 case 'dcModel'
  try
    y = RM.model('simulate',modelx,RM.additionalInputs{:});
  catch ME
    disp('Error evaluating dcModel in ResponseModel/eval')
    rethrow(ME)
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

  

