function passed = testDCSurfaceGetPrivate(dispBool)
%TESTDCSURFACEGETPRIVATE returns true if it appears OK
%
%   PASSED = TESTDCSURFACEGETPRIVATE returns true or false. Several messages will be
%   displayed to the screen.
%
%   PASSED = TESTDCSURFACEGETPRIVATE(DISPBOOL) will suppress all screen displays if
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
  disp('===Testing DCSurface getPrivate method===')
end

if dispBool
  disp('    Empty object')
end
if ~testObj(DCSurface,dispBool,1);
  passed = false;
end

if dispBool
  disp('    Scalar object')
end

coeffMatrix = [3 4 5]';
domain = struct('name',{'p1';'p2'},'range',{[-1 1];[-11 4]});
quadRM = ResponseModel(coeffMatrix*coeffMatrix',domain);
Surf1 = DCSurface(quadRM,domain); 
if ~testObj(Surf1,dispBool,1);
  passed = false;
end

if dispBool
  disp('    Vector object')
end

coeffMatrix = [3 8 5]';
domain = struct('name',{'p11';'p2'},'range',{[-1 1];[-11 4]});
quadRM = ResponseModel(coeffMatrix*coeffMatrix',domain);
Surf2 = DCSurface(quadRM,domain); 
if ~testObj([Surf1;Surf2],dispBool,1) || ~testObj([Surf1;Surf2],dispBool,2);
  passed = false;
end


function passed = testObj(obj,dispBool,element)

passed = true;
publicProps = {'type';'responseTransformation';'activeParameterIndex';...
  'activeParameterTransformation';'surrogateModel';...
  'trueModelRangeEstimate';'MOPairIdx';'surrogateFitInfo';'peakError';...
  'nSurfaces';'status'};

out = getPrivate(obj);

if ~isequal(sort(publicProps),sort(fieldnames(out)))
  passed = false;
  if dispBool
    disp('   ERROR: returned properties do not match public properties')
  end
end

% Get all public properties one by one
for i1 = 1:length(publicProps)
  try
    getPrivate(obj,publicProps{i1});
  catch
    passed = false;
    if dispBool
      disp(['   ERROR: failed to return ' publicProps{i1} ]);
    end
  end
end

% Try to get a nonexistant property
try %#ok
  getPrivate(obj,'trash')
  passed = false;
  if dispBool
    disp('   ERROR: attempt to obtain bogus property should have failed');
  end
end

% Eliminate 'nSurfaces'
idx = strmatch('nSurfaces',publicProps,'exact');
publicProps(idx) = [];

for i1 = 1:length(publicProps)
  try
    getPrivate(obj,publicProps{i1},element);
  catch
    passed = false;
    if dispBool
      disp(['   ERROR: failed to return ' publicProps{i1} ]);
    end
  end
end

try
  getPrivate(obj,[],element);
catch
  passed = false;
  if dispBool
    disp('   ERROR: failed 3 input calling syntax');
  end
end








