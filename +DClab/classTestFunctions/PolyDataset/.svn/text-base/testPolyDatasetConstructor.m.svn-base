function passed = testPolyDatasetConstructor(dispBool)
%TESTPOLYDATASETCONSTRUCTOR returns true if it appears OK
%
%   PASSED = TESTPOLYDATASETCONSTRUCTOR returns true or false. Several
%   messages will be displayed to the screen.
%
%   PASSED = TESTPOLYDATASETCONSTRUCTOR(DISPBOOL) will suppress all
%   screen displays if DISPBOOL==false.
%
%   This file tests essentially all constructor syntax.

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
  disp('===Testing PolyDataset constructor===')
end

if dispBool
  disp('    No inputs')
end
try
  PD = PolyDataset; %#ok
catch
  passed = false;
  if dispBool
    disp('   ERROR: constructor call with no inputs failed')
  end
end

if dispBool
  disp('    One input')
end

%Create a few MOPairs objects
coeffMatrix = [3 4 5]';
domain = struct('name',{'p1';'p2'},'range',{[-1 1];[-1 1]});
quadRM = ResponseModel(coeffMatrix*coeffMatrix',domain);
RO = ResponseObservation(5,1);
Pair1 = ModelAndObservationPair(RO,quadRM,'pair1');

addpath([fileparts(pwd) filesep 'ResponseModel'])
simpleRM = ResponseModel(@simpleDCModel);

RO = ResponseObservation(11,3);
Pair2 = ModelAndObservationPair(RO,simpleRM,'pair2');

P1 = FreeParameter('p1',0,0.5);
P2 = FreeParameter('p2',0.1,0.5);
P3 = FreeParameter('extra',0,3);
Dset = DCDataset([Pair1;Pair2],[P1;P2;P3]); %#ok

try
  PD = PolyDataset(Dset); %#ok
catch
  passed = false;
  if dispBool
    disp('   ERROR: constructor call with one input failed')
  end
end

if dispBool
  disp('    Two inputs')
end
opt = DCOptions('maxBranchBoundIter',2);

try 
  PD = PolyDataset(Dset,opt); %#ok
catch
  passed = false;
  if dispBool
    disp('   ERROR: constructor call with two inputs failed')
  end
end
