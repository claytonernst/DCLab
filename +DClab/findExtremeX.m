function xext = findExtremeX(x,y,N,type)
% xext = findExtremeX(x,y,N,type)
%
% function returns the x vectors corresponding to the top (or bottom) n y
% values.
%
% x should be a matrix of horzcatted column vectors [x1 ... xN]
% y should be a row vector.
% N should be an integer
% type should be 'max' or 'min'

if length(y) < N 
  error('Insufficient points to determine requested number of extrema')
end
n = size(x,1);  
xext = zeros(n,N);
tmpy = y;

switch type
  case 'min'
    ub = max(tmpy);
    for i1 = 1:N
      [trash idx] = min(tmpy);
      xext(:,i1) = x(:,idx);
      tmpy(idx) = ub;
    end
  case 'max'
   lb = min(tmpy);
   for i1 = 1:N
      [trash idx] = max(tmpy);
      xext(:,i1) = x(:,idx);
      tmpy(idx) = lb;
    end
  otherwise
    error('Condition should never occur')
end
