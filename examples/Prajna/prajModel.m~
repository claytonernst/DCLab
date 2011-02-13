function response = prajModel(executionType,paramValuesMatrix,time)
% model to predict x(t) for prajna's model. the third input is the time t

ni = nargin;
no = nargout;
%error checking the inputs:
error( nargchk(1,3,ni) );
error( nargoutchk(0,1,no) );
if ~ischar(executionType)	
  error(['First input to model simulation file must be of type ' ...		
         'char.'])
end

switch executionType	
 case 'getModelDomain';
   n = 2;
   tmp = cell(n,2);
   tmp(:,1) = {'p';'x0'};
   [tmp{:,2}] = deal([-inf inf]);
   response = cell2struct(tmp,{'name';'range'},2);
  
 case 'simulate'         %assign parameter values         
  N = size(paramValuesMatrix,2);
  response = zeros(2,N);

  %assign parameter values         
  for i1 = 1:N
    paramValues = paramValuesMatrix(:,i1);
    p = paramValues(1);
    x0 = paramValues(2);
    [t, y] = ode45(@model,[0 2 4],x0,[],p);
    response(:,i1) = [y(2);y(3)];
  end
  
%   ---isSaveEnabled is an optional condition, uncomment if required---
 case 'isSaveEnabled'
   response = 0; % Enter a value of 0 or 1. If 1, model evaluations invoked
            % by the 'simulate' flag will be saved and any inputs in 
            % varargin must be single row character arrays or scalars.

% ---isMultipleResponsesEnabled is an optional condition, uncomment if required---
 case 'isMultipleResponsesEnabled'
   response = 1; % Enter a value of 0 or 1. If 1, the output in the 'simulate'
            % case should be Nresp x size(paramMatrix,2). Additionally a
            % responseList must be defined in this file and the first element
            % of varargin must be a member of responseList.

% ---getResponseList is a required condition when isMultipleResponsesEnabled == 1 ---
 case 'getResponseList'
   response = {'2';'4'}; % Enter a column cell array of single line chars. The first
            % element of varargin must be an element of this cell array.
 otherwise
  error('Incorrect input type to model simulation file')
end
  
function xdot = model(t,x,p)

xdot = -p*x^3;
