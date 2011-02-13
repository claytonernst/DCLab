function [G,y] = planerotPackard(x)
%PLANEROT Givens plane rotation.
%   [G,Y] = PLANEROT(X), where X is a 2-component column vector,
%   returns a 2-by-2 orthogonal matrix G so that Y = G*X has Y(2) = 0.
%
%   Class support for input X:
%      float: double, single
%
%   See also QRINSERT, QRDELETE.

%   Copyright 1984-2004 The MathWorks, Inc. 
%   $Revision: 5.10.4.1 $  $Date: 2004/04/10 23:30:05 $

if x(2) ~= 0
   r = norm(x);
%    G = [x'; -x(2) x(1)]/r;
%    x = [r; 0];
   x(1) = x(1)/r;
   x(2) = x(2)/r;
   G = zeros(2,2);
   G(1) = x(1);
   G(2) = -x(2);
   G(3) = x(2);
   G(4) = x(1);
   y = zeros(2,1);
   y(1) = r;
else
   G = eye(2,class(x));
   y = x;
end
