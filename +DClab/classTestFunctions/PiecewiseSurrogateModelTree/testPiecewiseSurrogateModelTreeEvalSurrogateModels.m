function passed = testPiecewiseSurrogateModelTreeEvalSurrogateModels(dispBool)
%testPiecewiseSurrogateModelTreeEvalSurrogateModels returns true if it appears OK
%
%   PASSED = TESTPIECEWISESURROGATEMODELTREEEVALSURROGATEMODELS returns true or
%   false. Several messages will be displayed to the screen.
%
%   PASSED = TESTPIECEWISESURROGATEMODELTREEEVALSURROGATEMODELS(DISPBOOL) will
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
  disp('===Testing PiecewiseSurrogateModelTree evalSurrogateModels method===')
end

% Create a partitioned object to play with.
% Create an object to play with
%Create a few RM objects
coeffMatrix = [3 4 5]';
domain = struct('name',{'p1';'p2'},'range',{[-1 1];[-1 1]});
quadRM = ResponseModel(coeffMatrix*coeffMatrix',domain);
addpath([fileparts(pwd) filesep 'ResponseModel'])
simpleRM = ResponseModel(@simpleDCModel);
RMCell = {quadRM;simpleRM};

domain = struct('name',{'extra1';'p2';'p1';'extra'},'range',{[-1 4];[-0.9 1];[-0.5 0.5];[-10 10]});

if dispBool
  disp('    testing with unsubdivided tree')
end


Tree = PiecewiseSurrogateModelTree(RMCell,domain);

% Pick some random points on the domain.
Neval = 30;
domrng = vertcat(domain.range);
x = diag(diff(domrng,[],2))*rand(size(domain,1),Neval) + repmat(domrng(:,1),1,Neval);

ya = eval(RMCell{1},x,{domain.name}');
yb = eval(RMCell{2},x,{domain.name}');

try
  [y1 yInt1] = evalSurrogateModels(Tree,x);
  
  if ~doesIntersect([ya' ya'],yInt1([1 2],:)',1e-10,1e-10) || ...
         ~doesIntersect([yb' yb'],yInt1([3 4],:)',1e-10,1e-10)
    passed = false;
    if dispBool
      disp('    ERROR: evalSurrogateModels on undivided tree produced values that appear incorrect')
    end
  end
catch
  passed = false;
  if dispBool
    disp('    ERROR: evalSurrogateModels on undivided Tree failed')
  end
end

if dispBool
  disp('    testing with only 2nd RM as output')
end

try
  [y1 yInt1] = evalSurrogateModels(Tree,x,2);
  
  if size(yInt1,1) ~= 2 || ~doesIntersect([yb' yb'],yInt1([1 2],:)',1e-10,1e-10)
    passed = false;
    if dispBool
      disp('    ERROR: evalSurrogateModels produced values that appear incorrect')
    end
  end
catch
  passed = false;
  if dispBool
    disp('    ERROR: evalSurrogateModels on undivided Tree failed')
  end
end

if dispBool
  disp('    testing case where supplied x has more dimensions than tree domain')
end

domain2 = [struct('name','asdf','range',[4 6]); domain; struct('name','qwer','range',[4 6])];

% Pick some random points on the domain.
Neval = 30;
domrng = vertcat(domain2.range);
x = diag(diff(domrng,[],2))*rand(size(domain2,1),Neval) + repmat(domrng(:,1),1,Neval);

ya = eval(RMCell{1},x,{domain2.name}');
yb = eval(RMCell{2},x,{domain2.name}');

try
  [y1 yInt1] = evalSurrogateModels(Tree,x,[],{domain2.name}');
  
  if ~doesIntersect([ya' ya'],yInt1([1 2],:)',1e-10,1e-10) || ...
         ~doesIntersect([yb' yb'],yInt1([3 4],:)',1e-10,1e-10)
    passed = false;
    if dispBool
      disp('    ERROR: evalSurrogateModels produced values that appear incorrect')
    end
  end
catch
  passed = false;
  if dispBool
    disp('    ERROR: evalSurrogateModels on undivided Tree failed')
  end
end


if dispBool
  disp('    testing on subdivided tree')
end

Tree2 = subdivideDomain(Tree,RMCell,1,'extra',0); %#ok
Tree3 = subdivideDomain(Tree2,RMCell,2,'p2',0.5);
Tree4 = subdivideDomain(Tree3,RMCell,3,'p1',0.15); %#ok
Tree5 = subdivideDomain(Tree4,RMCell,5,'p1',0.25,[1; 0]); %#ok

% Pick some random points on the domain.
Neval = 30;
domrng = vertcat(domain.range);
x = diag(diff(domrng,[],2))*rand(size(domain,1),Neval) + repmat(domrng(:,1),1,Neval);

ya = eval(RMCell{1},x,{domain.name}');
yb = eval(RMCell{2},x,{domain.name}');

try
  [y1 yInt1] = evalSurrogateModels(Tree,x);
  [y2 yInt2] = evalSurrogateModels(Tree5,x);
  
  if ~doesIntersect([ya' ya'],yInt1([1 2],:)',1e-10,1e-10) || ...
         ~doesIntersect([yb' yb'],yInt1([3 4],:)',1e-10,1e-10) || ...
         ~doesIntersect([ya' ya'],yInt2([1 2],:)',1e-10,1e-10) || ...
         ~doesIntersect([yb' yb'],yInt2([3 4],:)',1e-10,1e-10)
    passed = false;
    if dispBool
      disp('    ERROR: eval on subdivided tree produced values that appear incorrect')
    end
  end
catch
  passed = false;
  if dispBool
    disp('    ERROR: evalSurrogateModels on subdivided Tree failed')
  end
end


