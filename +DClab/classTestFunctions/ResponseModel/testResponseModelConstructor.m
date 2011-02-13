function passed = testResponseModelConstructor(dispBool)
%TESTRESPONSEMODELCONSTRUCTOR returns true if it appears OK
%
%   PASSED = TESTRESPONSEMODELCONSTRUCTOR returns true or false. Several
%   messages will be displayed to the screen.
%
%   PASSED = TESTRESPONSEMODELCONSTRUCTOR(DISPBOOL) will suppress all
%   screen displays if DISPBOOL==false.
%
%   This file tests essentially all constructor syntax. The constructor
%   calls the isvalid method, so problems there may show up
%   here. Comment out the appropriate line of the constructor to isolate.

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
  disp('===Testing ResponseModel constructor===')
end

if dispBool
  disp('    No inputs')
end

try
  Obj = ResponseModel; %#ok
catch
  passed = false;
  if dispBool
    disp('   ERROR: constructor call with no inputs failed')
  end 
end
  
if dispBool
  disp('    bogus inputs')
end

domain = struct('name',{'p1';'p2'},'range',{[-1 1];[-inf inf]});
coeffMatrix = [3 4 5]';

try %#ok
  Obj = ResponseModel({5},domain); %#ok
  passed = false;
  if dispBool
    disp('   ERROR: should not have accepted cell array as first input')
  end
end

try %#ok
  Obj = ResponseModel('hi'); %#ok
  passed = false;
  if dispBool
    disp('   ERROR: should not have accepted char array as first input')
  end
end

try %#ok
  Obj = ResponseModel(coeffMatrix,domain,{1}); %#ok
  passed = false;
  if dispBool
    disp('   ERROR: should not have accepted cell array as third input')
  end
end

if dispBool
  disp('    2 input constructor call with linear algebraic model')
end
coeffMatrix = [3 4 5]';
domain = struct('name',{'p1';'p2'},'range',{[-1 1];[-inf inf]});

try
  Obj = ResponseModel(coeffMatrix,domain); %#ok
catch
  passed = false;
  if dispBool
    disp('   ERROR: call failed---should have been successful')
  end
end

if dispBool
  disp('    2 input constructor call with transformed linear algebraic model')
end
% Check no transformations
model.value = coeffMatrix;
model.responseTransformation = 'none'; 
model.variableTransformations = {'none';'none'};
try
  Obj = ResponseModel(model,domain); %#ok
catch
  passed = false;
  if dispBool
    disp('   ERROR: call with no transformations failed---should have been successful')
  end
end

if dispBool
  disp('    2 input constructor call with quadratic algebraic model')
end
coeffMatrix = [1 4 5; 4 5 6; 5 6 7];

try
  Obj = ResponseModel(coeffMatrix,domain); %#ok
catch
  passed = false;
  if dispBool
    disp('   ERROR: call failed---should have been successful')
  end
end

if dispBool
  disp('    2 input constructor call with transformed quadratic algebraic model')
end

% Check mixed transformations
model.value = coeffMatrix;
model.responseTransformation = 'log10'; 
model.variableTransformations = {'none';'log10'};
domain = struct('name',{'p1';'p2'},'range',{[-1 1];[1 inf]});
try
  Obj = ResponseModel(model,domain); %#ok
catch
  passed = false;
  if dispBool
    disp('   ERROR: call failed---should have been successful')
  end
end

if dispBool
  disp('    3 input constructor call with quadratic algebraic model and scalar outputUncertainty')
end
try
  Obj = ResponseModel(coeffMatrix,domain,1); %#ok
catch
  passed = false;
  if dispBool
    disp('   ERROR: call failed---should have been successful')
  end
end
  
if dispBool
  disp('    3 input constructor call with quadratic algebraic model and vector outputUncertainty')
end
try
  Obj = ResponseModel(coeffMatrix,domain,[-0.1 1]); %#ok
catch
  passed = false;
  if dispBool
    disp('   ERROR: call failed---should have been successful')
  end
end

if dispBool
  disp('    3 input constructor call with quadratic algebraic model and tranformed vector outputUncertainty')
end
outputUnc.value = [-0.1 1];
outputUnc.type = 'relative';
outputUnc.transformation = 'log10';
try
  Obj = ResponseModel(coeffMatrix,domain,outputUnc); %#ok
catch
  passed = false;
  if dispBool
    disp('   ERROR: call failed---should have been successful')
  end
end

if dispBool
  disp('    4 input constructor call with quadratic algebraic model empty guiDisplayCallback')
end
try
  Obj = ResponseModel(coeffMatrix,domain,outputUnc,[]); %#ok
catch
  passed = false;
  if dispBool
    disp('   ERROR: call failed---should have been successful')
  end
end

if dispBool
  disp('    4 input constructor call with quadratic algebraic model')
end
try
  Obj = ResponseModel(coeffMatrix,domain,outputUnc,@disp); %#ok
catch
  passed = false;
  if dispBool
    disp('   ERROR: call failed---should have been successful')
  end
end

if dispBool
  disp('    simple dcModel constructor call')
end
try
  simpleDCObj = ResponseModel(@simpleDCModel); %#ok
catch
  passed = false;
  if dispBool
    disp('   ERROR: call failed---should have been successful')
  end
end

if dispBool
  disp('    less simple dcModel constructor call')
end
try
  lsimpleDCObj = ResponseModel(@lessSimpleDCModel); %#ok
catch
  passed = false;
  if dispBool
    disp('   ERROR: call failed---should have been successful')
  end
end

if dispBool
  disp('    complex dcModel constructor call')
end
try
  cplxDCObj = ResponseModel(@complexDCModel,'resp1'); %#ok
catch
  passed = false;
  if dispBool
    disp('   ERROR: call failed---should have been successful')
  end
end
