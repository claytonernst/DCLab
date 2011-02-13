function passed = testPiecewiseSurrogateModelTreeSubdivideDomain(dispBool)
%TESTPIECEWISESURROGATEMODELSUBDIVIDEDOMAIN returns true if it appears OK
%
%   PASSED = TESTPIECEWISESURROGATEMODELSUBDIVIDEDOMAIN returns true or
%   false. Several messages will be displayed to the screen.
%
%   PASSED = TESTPIECEWISESURROGATEMODELSUBDIVIDEDOMAIN(DISPBOOL) will
%   suppress all screen displays if DISPBOOL==false.

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
  disp('===Testing PiecewiseSurrogateModelTree subdivideDomain method===')
end

if dispBool
  disp('    Splitting on nonactive dimension')
end

% Create an object to play with
%Create a few RM objects
coeffMatrix = [3 4 5]';
domain = struct('name',{'p1';'p2'},'range',{[-1 1];[-1 1]});
quadRM = ResponseModel(coeffMatrix*coeffMatrix',domain);
 
addpath([fileparts(pwd) filesep 'ResponseModel'])
simpleRM = ResponseModel(@simpleDCModel);

RMCell = {quadRM;simpleRM};

domain(1).name = 'extra1';
domain(1).range = [-1 4];
domain(2,1).name = 'p2';
domain(2).range = [-0.9 1];
domain(3,1).name = 'p1';
domain(3).range = [-0.5 0.5];
domain(4,1).name = 'extra';
domain(4).range = [-10 10];

Tree = PiecewiseSurrogateModelTree(RMCell,domain);
try
  Tree2 = subdivideDomain(Tree,RMCell,1,'extra',0); %#ok
catch
  passed = false;
  if dispBool
    disp('    ERROR: subdivideDomain call failed')
  end
end
    
if dispBool
  disp('    Splitting on active dimension')
end

try
  Tree2 = subdivideDomain(Tree,RMCell,1,'p2',0.5);
catch
  passed = false;
  if dispBool
    disp('    ERROR: subdivideDomain call failed')
  end
end

if dispBool
  disp('    Splitting further')
end

try
  Tree3 = subdivideDomain(Tree2,RMCell,3,'p1',0.15); %#ok
catch
  passed = false;
  if dispBool
    disp('    ERROR: subdivideDomain call failed')
  end
end

if dispBool
  disp('    Splitting further, with one surface inherited')
end

try
  Tree4 = subdivideDomain(Tree3,RMCell,5,'p1',0.25,[1; 0]); %#ok
catch
  passed = false;
  if dispBool
    disp('    ERROR: subdivideDomain call failed')
  end
end

if dispBool
  disp('    Trying to split nonexistent dimension, node, and location outside of domain')
end

try %#ok
  Tree5 = subdivideDomain(Tree4,RMCell,11,'p1',0.20,[1; 0]); %#ok
  passed = false;
  if dispBool
    disp('    ERROR: attempt to split nonexistant node should have failed')
  end
end

try %#ok
  Tree5 = subdivideDomain(Tree4,RMCell,7,'p1',0.20,[1; 0]); %#ok
  passed = false;
  if dispBool
    disp('    ERROR: attempt to split outside of domain of supplied node should have failed')
  end
end

try %#ok
  Tree5 = subdivideDomain(Tree4,RMCell,7,'p7',0.20,[1; 0]); %#ok
  passed = false;
  if dispBool
    disp('    ERROR: attempt to split nonexistent dimension should have failed')
  end
end

try %#ok
  Tree5 = subdivideDomain(Tree4,RMCell,7,'p7',0.20,[1 0]); %#ok
  passed = false;
  if dispBool
    disp('    ERROR: attempt to split which inproper inherit dimension should have failed')
  end
end
