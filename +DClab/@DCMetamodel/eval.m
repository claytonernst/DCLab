function y = eval(MM,x,paramList)
%DCMETAMODEL/EVAL evaluates the model at supplied parameter vectors
%
%  Y = EVAL(MM,X) evaluates the model of the DCMetaModel object MM at
%  parameter vectors given taken from the columns of the Nparam-by-Nevals
%  matrix X. The parameter values in X must be in the order expected by the
%  model. I.e., Y(1) = M(X(1,1),...X(Nparam,1)). Y will be 1-by-Nevals.
%
%  Y = EVAL(MM,X,PARAMLIST) assumes the parameter values in each
%  columns of X are in the order indicated by the column cell array of
%  chars PARAMLIST. PARAMLIST must be a superset of the parameters needed
%  by the model. This function will then extract the needed rows of X and
%  reorder them as necessary before evaluating the model. 
%
%  Inputs:
%    MM: A DCMetamodel object
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
%
%   See also DCMetamodel

% Last modified 1/26/08

ni = nargin;
if ni == 2
  paramList = MM(1).parameterList;
  if size(x,1) ~= length(paramList)
      error('Usage: X and MM.parameterList have incompatable dimensions')
  end
  modelx = x;
elseif ni == 3
  if size(x,1) ~= length(paramList)
      error('Usage: X and PARAMLIST have incompatable dimensions')
  end
  [sortedList idx1 idx2] = intersect(paramList,MM(1).parameterList); %#ok
  if ~isequal(sortedList,sort(MM(1).parameterList))
    error('Usage: paramList contains insufficient parameters for the DCMetamodel')
  end
  modelx(idx2,:) = x(idx1,:);
else
  error(nargchk(2,3,ni))
end

%TODO add code for a refined MM.

logIdx = MM.variableTransformations == 2;
inactive = MM.variableTransformations == 0;
modelx(logIdx,:) = log10(modelx(logIdx,:));
modelx(inactive,:) = [];

switch MM.model.type
 case 'svmd'
   y = evalSvmModel(MM.model,modelx')';  
 otherwise
  error('Code not complete')
end

if strcmp(MM.responseTransformation,'log10')
  y = 10.^y;
end


