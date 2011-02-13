function passed = testPiecewiseSurrogateModelTreeConstructor(dispBool)
%TESTPIECEWISESURROGATEMODELCONSTRUCTOR returns true if it appears OK
%
%   PASSED = TESTPIECEWISESURROGATEMODELCONSTRUCTOR returns true or false. Several
%   messages will be displayed to the screen.
%
%   PASSED = TESTPIECEWISESURROGATEMODELCONSTRUCTOR(DISPBOOL) will suppress all
%   screen displays if DISPBOOL==false.

% This file tests essentially all constructor syntax.

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
  disp('===Testing PiecewiseSurrogateModelTree constructor===')
end

if dispBool
  disp('    No inputs')
end
try
  Tree = PiecewiseSurrogateModelTree; %#ok
catch
  passed = false;
  if dispBool
    disp('   ERROR: constructor call with no inputs failed')
  end
end

if dispBool
  disp('    Two inputs, two ResponseModels')
end

%Create a few RM objects
coeffMatrix = [3 4 5]';
domain = struct('name',{'p1';'p2'},'range',{[-1 1];[-1 1]});
quadRM = ResponseModel(coeffMatrix*coeffMatrix',domain);

addpath([fileparts(pwd) filesep 'ResponseModel'])
simpleRM = ResponseModel(@simpleDCModel);

domain(1).name = 'extra1';
domain(1).range = [-1 4];
domain(2,1).name = 'p2';
domain(2).range = [-0.9 1];
domain(3,1).name = 'p1';
domain(3).range = [-0.5 0.5];
domain(4,1).name = 'extra';
domain(4).range = [-10 10];

try
  Tree = PiecewiseSurrogateModelTree({quadRM;simpleRM},domain); %#ok
catch
  passed = false;
  if dispBool
    disp('   ERROR: constructor call with two inputs failed')
  end
end


