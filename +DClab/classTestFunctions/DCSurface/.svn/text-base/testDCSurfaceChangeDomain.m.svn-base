function passed = testDCSurfaceChangeDomain(dispBool)
%TESTDCSURFACECHANGE returns true if it appears OK
%
%   PASSED = TESTDCSURFACECHANGE returns true or false. Several messages will be
%   displayed to the screen.
%
%   PASSED = TESTDCSURFACECHANGE(DISPBOOL) will suppress all screen displays if
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
  disp('===Testing DCSurface changeDomain method===')
end

% Create an ResponseModelto play with.
coeffMatrix = [11 4 5]';
domain = struct('name',{'p1';'p2'},'range',{[-1 1];[-11 4]});
quadRM = ResponseModel(coeffMatrix*coeffMatrix',domain);

if dispBool
  disp('    Creating DCSurface with no transformations')
end

% Create the surface on a domain with an extra parameter and different
% order parameters
domain = struct('name',{'dummy';'p2';'p1'},'range',{[5 6];[-11 4];[0.1 1]});
Surf = DCSurface(quadRM,domain); 

% Create a different domain. It includes additional parameters, eliminates
% a nonactive one, and shrinks the range of and active one.
domain2 = [struct('name','bb','range',[0 6]); domain(3); domain(2); struct('name','aa','range',[-6 6])];
domain2(3).range = [2 5];
domain2(2).range = [0.5 1];
Surf2 = changeDomain(Surf,domain,domain2);

y1 = eval(Surf,domain,[5.5; 3; 0.5]);
y2 = eval(Surf2,domain2,[1; 0.5; 3; 2]);

if abs(y1-y2)/abs(y1) > 1e-10
  passed = false;
  if dispBool
    disp('    ERROR: changing domain of surface changed the output')
  end
end

if dispBool
  disp('    Creating DCSurface with log10 output transformation')
end

Surf = DCSurface(quadRM,domain,'log10'); 
Surf2 = changeDomain(Surf,domain,domain2);

y1 = eval(Surf,domain,[5.5; 3; 0.5]);
y2 = eval(Surf2,domain2,[1; 0.5; 3; 2]);

if abs(y1-y2)/abs(y1) > 1e-10
  passed = false;
  if dispBool
    disp('    ERROR: changing domain of surface changed the output')
  end
end

if dispBool
  disp('    Creating DCSurface with log10 output transformation, a log10 input transformation')
end

Surf = DCSurface(quadRM,domain,'log10',{'log10';'none';'log10'}); 
Surf2 = changeDomain(Surf,domain,domain2);

y1 = eval(Surf,domain,[5.5; 3; 0.5]);
y2 = eval(Surf2,domain2,[1; 0.5; 3; 2]);

if abs(y1-y2)/abs(y1) > 1e-10
  passed = false;
  if dispBool
    disp('    ERROR: changing domain of surface changed the output')
  end
end

if dispBool
  disp('    Creating  a vector DCSurface containing the above three')
end

SurfA = DCSurface(quadRM,domain); 
SurfB = DCSurface(quadRM,domain,'log10'); 
SurfC = DCSurface(quadRM,domain,'log10',{'log10';'none';'log10'}); 

Surf = [SurfA; SurfB; SurfC];
Surf2 = changeDomain(Surf,domain,domain2);

y1 = eval(Surf,domain,[5.5; 3; 0.5]);
y2 = eval(Surf2,domain2,[1; 0.5; 3; 2]);

if any(abs(y1-y2)./abs(y1) > 1e-10)
  keyboard
  passed = false;
  if dispBool
    disp('    ERROR: changing domain of surface changed the output')
  end
end



