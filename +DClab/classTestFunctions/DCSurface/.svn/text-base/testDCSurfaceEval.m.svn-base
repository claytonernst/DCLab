function passed = testDCSurfaceEval(dispBool)
%TESTDCSURFACEEVAL returns true if it appears OK
%
%   PASSED = TESTDCSURFACEEVAL returns true or false. Several messages will be
%   displayed to the screen.
%
%   PASSED = TESTDCSURFACEEVAL(DISPBOOL) will suppress all screen displays if
%   DISPBOOL==false.

% Last modified 9/4/07 rpf

% Copyright 2007 Ryan Feeley.   All rights reserved.
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
% 02110-1301, USA
%
% Bug reports may be sent to Ryan Feeley, rfeeley@me.berkeley.edu

if nargin == 0
  dispBool = true;
end
passed = true;

if dispBool
  disp('===Testing DCSurface eval method===')
end

% Create a ResponseModel to play with
coeffMatrix = [3 4 5]';
domain = struct('name',{'p1';'p2'},'range',{[-1 1];[-11 4]});
quadRM = ResponseModel(coeffMatrix*coeffMatrix',domain);

if dispBool
  disp('    Empty object')
end

try %#ok
  eval(DCSurface,struct('name','a','range',[1 2]),1.5);
  passed = false;
  if dispBool
    disp('   ERROR: attempt to call eval with an empty object should have failed');
  end
end

if dispBool
  disp('    Scalar object')
end

% Create the surface on a domain with an extra parameter and different
% order parameters
domain = struct('name',{'dummy';'p2';'p1'},'range',{[5 6];[-11 4];[-1 1]});
Surf = DCSurface(quadRM,domain); 

% Test eval with multiple parameter vectors
y1 = eval(quadRM,[0 0.5; -5 2]);
[y2 yInt] = eval(Surf,domain,[5.5 5.1; -5 2; 0 0.5]);

if any(abs(y1-y2)./abs(y1) > 1e-15) || any(yInt(1,:) > y1 + 1e-15) || any(yInt(2,:) + 1e-15 < y1)
  passed = false;
  if dispBool
    disp('   ERROR: eval produced the wrong output');
  end
end

try %#ok
  [y2 yInt] = eval(Surf,domain,[7; -5; 0]); %#ok
  passed = false;
  if dispBool
    disp('   ERROR: attempt to call eval with x outside domain should have failed');
  end
end

if dispBool
  disp('    Scalar object with nonoptimal transformations')
end

% Create the surface on a domain with an extra parameter and different
% order parameters
domain = struct('name',{'dummy';'p2';'p1'},'range',{[5 6];[1 4];[-1 1]});
Surf = DCSurface(quadRM,domain,'log10',{'log10';'log10';'none'}); 

y1 = eval(quadRM,[0; 3]);
[y2 yInt] = eval(Surf,domain,[5.5; 3; 0]);

if yInt(1) > y1 + 1e-15 || yInt(2) + 1e-15 < y1
  passed = false;
  if dispBool
    disp('   ERROR: eval produced the wrong output');
  end
end

if dispBool
  disp('    Vector object with nonoptimal transformations')
end

% Create the surface on a domain with an extra parameter and different
% order parameters
domain = struct('name',{'dummy';'p2';'p1'},'range',{[5 6];[1 4];[-1 1]});
Surf1 = DCSurface(quadRM,domain,'log10',{'log10';'log10';'none'}); 
Surf2 = DCSurface(quadRM,domain,'log10',{'log10';'none';'none'}); 

Surf = vertcat(Surf1,Surf2);

y1 = eval(quadRM,[0; 3]);
[y2 yInt] = eval(Surf,domain,[5.5; 3; 0]);

if yInt(1) > y1 + 1e-15 || yInt(2) + 1e-15 < y1 || yInt(3) > y1 + 1e-15 || yInt(4) + 1e-15 < y1
  passed = false;
  if dispBool
    disp('   ERROR: eval produced the wrong output');
  end
end



