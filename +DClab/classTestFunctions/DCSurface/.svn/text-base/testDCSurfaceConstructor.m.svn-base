function passed = testDCSurfaceConstructor(dispBool)
%TESTDCSURFACECONSTRUCTOR returns true if it appears OK
%
%   PASSED = TESTDCSURFACECONSTRUCTOR returns true or false. Several messages will be
%   displayed to the screen.
%
%   PASSED = TESTDCSURFACECONSTRUCTOR(DISPBOOL) will suppress all screen displays if
%   DISPBOOL==false.
%
%   This file tests essentially all constructor syntax. The constructor
%   calls the isvalid method, so problems there may show up
%   here. Comment out the appropriate line of the constructor to isolate.
%
%   See also testDCSurface

% Last modified 8/28/07 rpf

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
  disp('===Testing DCSurface constructor===')
end

if dispBool
  disp('    No inputs')
end
try
  Surf = DCSurface; %#ok
catch
  passed = false;
  if dispBool
    disp('   ERROR: constructor call with no inputs failed')
  end
end

if dispBool
  disp('    One scalar input')
end
try
  Surf = DCSurface(5); %#ok
catch
  passed = false;
  if dispBool
    disp('   ERROR: constructor call with one input')
  end
end

if dispBool
  disp('    Bogus inputs')
end

try %#ok
  Surf = DCSurface(struct('a',1','b',2)); %#ok
  passed = false;
  if dispBool
    disp('   ERROR: should not have accepted structure as first input')
  end
end

coeffMatrix = [3 4 5]';
domain = struct('name',{'p1';'p2'},'range',{[-1 1];[-inf inf]});
quadRM = ResponseModel(coeffMatrix*coeffMatrix',domain);

try %#ok
  Surf = DCSurface(quadRM,domain); %#ok
  passed = false;
  if dispBool
    disp('   ERROR: should not have accepted domain with infinite range as 2nd input')
  end
end

domain = struct('name',{'p1';'p2'},'range',{[-10 1];[2 5]});
try %#ok
  Surf = DCSurface(quadRM,domain,{1}); %#ok
  passed = false;
  if dispBool
    disp('   ERROR: should not have accepted cell array as third input')
  end
end

try %#ok
  Surf = DCSurface(quadRM,domain,[],{1}); %#ok
  passed = false;
  if dispBool
    disp('   ERROR: should have complained about 4th input')
  end
end

try %#ok
  Surf = DCSurface(quadRM,domain,[],[],{1}); %#ok
  passed = false;
  if dispBool
    disp('   ERROR: should have complained about 5th input')
  end
end

try %#ok
  Surf = DCSurface(quadRM,domain,[],{'none';'log10'},7,{1}); %#ok
  passed = false;
  if dispBool
    disp('   ERROR: should have complained about 6th input')
  end
end

domain = struct('name',{'dummy';'p1';'p2'},'range',{[2 4];[-10 10];[-1 1]});
try %#ok
  Surf = DCSurface(quadRM,domain); %#ok
  passed = false;
  if dispBool
    disp('   ERROR: should have complained that domain is not contained in domain of RespModel')
  end
end

if dispBool
  disp('    Algebraic true model')
end
coeffMatrix = [3 4 5]';
domain = struct('name',{'p1';'p2'},'range',{[-10 1];[-5 2]});
quadRM = ResponseModel(coeffMatrix*coeffMatrix',domain);

%Build the surface with the order of the parameters switched.
tmpdom(1) = domain(2);
tmpdom(2,1) = domain(1);

try
  Surf = DCSurface(quadRM,tmpdom); %#ok
  t = struct(Surf);
  normQuad = t.surrogateModel;
  %Evaluate the normalized model at (1,-1). This should be the save as
  %evaluating quadRM at (-10,2).
  y1 = [1 1 -1]*normQuad*[1;1;-1];
  y2 = eval(quadRM,[-10;2]);
  if abs(y1-y2) > 1e-10
    passed = false;
    if dispBool
      disp('   ERROR: problem with scaling in DCSurface')
    end
  end
catch
  passed = false;
  if dispBool
    disp('   ERROR: constructor call with algebraic model failed')
  end
end

if dispBool
  disp('    3 inputs: transformed algebraic true model as ResponseModel')
end
addpath([fileparts(pwd) filesep 'ResponseModel'])
simpleRM = ResponseModel(@simpleDCModel);
domain = simpleRM.domain;
%Change [-inf inf] to finite.
domain(2).range = [-5 8];
%Add a dummy parameter to the begining and end and switch order
dummy1.name = 'dummy1';
dummy1.range = [-10 0];
dummy2.name = 'dummy2';
dummy2.range = [1 10];
domain = [dummy1; domain(2); domain(1); dummy2];

try
  Surf = DCSurface(simpleRM,domain,'log10'); %#ok
  % Examine the created surface.
  % There should should be two active parameters, both with none
  % transformation.
  t = struct(Surf);
  if ~isequal(t.activeParameterIndex,[2 1; 3 1])
    passed = false;
    if dispBool
      disp('   ERROR: incorrect active parameters were determined')
    end
  end
  normQuad = t.surrogateModel;
  %Evaluate the normalized model at (1,-1). Since the parameter order is
  %different, this should be the same as evaluating simpleRM at (-1,8).
  y1 = [1 1 -1]*normQuad*[1;1;-1];
  y1 = 10^y1;
  y2 = eval(simpleRM,[-1;8]);
  peakError = Surf.peakError;
  
  % Normalize since we're comparing big numbers.
  if (y1-y2)/abs(y1) < peakError(1)/abs(y1) - 1e-10 || (y1-y2)/abs(y1) > peakError(2)/abs(y1) + 1e-10 
    passed = false;
    if dispBool
      disp('   ERROR: insufficent fitting error bounds')
    end
  end  

  % There should be essentially no fitting error.
  if max(abs(peakError)) > 1e-10
    passed = false;
    if dispBool
      disp('   ERROR: fitting error should be very small')
    end
  end
  
catch
  passed = false;
  if dispBool
    disp('   ERROR: constructor call with algebraic model failed')
  end
end

if dispBool
  disp('    6 inputs: same ResponseModel, nonoptimal transformations for surface')
end

domain = simpleRM.domain;
%Change [-inf inf] to finite.
domain(2).range = [1 8];
%Add a dummy parameter to the begining and end and switch order
dummy1.name = 'dummy1';
dummy1.range = [-10 0];
dummy2.name = 'dummy2';
dummy2.range = [1 10];
domain = [dummy1; domain(2); domain(1); dummy2];

try
  Surf = DCSurface(simpleRM,domain,'none',{'none';'log10';'none';'log10'},2,DCOptions); %#ok
  % Examine the created surface.
  % There should should be two active parameters, the 1st with a log10
  % transformation and the 2nd with none
  % transformation.
  t = struct(Surf);
  if ~isequal(t.activeParameterIndex,[2 2; 3 1])
    passed = false;
    if dispBool
      disp('   ERROR: incorrect active parameters were determined')
    end
  end
  normQuad = t.surrogateModel;
  %Evaluate the normalized model at (1,-1). Since the parameter order is
  %different, this should be the same as evaluating simpleRM at (-1,8).
  y1 = [1 1 -1]*normQuad*[1;1;-1];
  y2 = eval(simpleRM,[-1;8]);
  peakError = Surf.peakError;
  
  if (y1-y2)/abs(y1) < peakError(1)/abs(y1) - 1e-10 || (y1-y2)/abs(y1) > peakError(2)/abs(y1) + 1e-10 
    passed = false;
    if dispBool
      disp('   ERROR: insufficent fitting error bounds')
    end
  end  
  
catch
  passed = false;
  if dispBool
    disp('   ERROR: constructor call failed')
  end
end

%TODO this function needs to check a lot more things. The DCSurface code is
%complicated!
%
%we'll do more in the interactive graphical function testDCSurface. This
%function will assume all the methods are "working".
