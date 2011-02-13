function passed = testPiecewiseSurrogateModelTreeChangeDomain(dispBool)
%testPiecewiseSurrogateModelTreeChangeDomain returns true if it appears OK
%
%   PASSED = TESTPIECEWISESURROGATEMODELTREECHANGEDOMAIN returns true or
%   false. Several messages will be displayed to the screen.
%
%   PASSED = TESTPIECEWISESURROGATEMODELTREECHANGEDOMAIN(DISPBOOL) will
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
  disp('===Testing PiecewiseSurrogateModelTree changeDomain method===')
end

% Create simple object to play with. Use quadratic response models to
% make it easier to verify that everything works. The don't have that dang
% fitting error.
tmp = rand(3);
coeffMatrix = tmp+tmp';
domain = struct('name',{'p1';'p2'},'range',{[-1 1];[-1 1]});
quadRM1 = ResponseModel(coeffMatrix*coeffMatrix',domain);

tmp = rand(3);
coeffMatrix = tmp+tmp';
domain = struct('name',{'p1';'p2'},'range',{[-10 10];[-2 2]});
quadRM2 = ResponseModel(coeffMatrix*coeffMatrix',domain);

RMCell = {quadRM1;quadRM2};
domain = struct('name',{'extra1';'p2';'p1';'extra2'},'range',{[-1 4];[-0.9 1];[-0.5 0.5];[-10 10]});
origTree = PiecewiseSurrogateModelTree(RMCell,domain);

if dispBool
  disp('    Adding a dimension to an unsubdivided tree')
end

domain = [domain; struct('name','extra3','range',[3 4])];

try
  Tree = changeDomain(origTree,domain); %#ok
  
  % Pick some random points on the domain.
  Neval = 3000;
  domrng = vertcat(domain.range);
  x = diag(diff(domrng,[],2))*rand(size(domain,1),Neval) + repmat(domrng(:,1),1,Neval);
  
  y1 = evalSurrogateModels(origTree,x,[],{domain.name}');
  y2 = evalSurrogateModels(Tree,x);
  if any(any(y1-y2 > 1e-10))
    passed = false;
    if dispBool
      disp('    ERROR: eval on modified tree produced values that appear incorrect')
    end
  end
catch
  passed = false;
  if dispBool
    disp('    ERROR: call to changeDomain failed')
  end
end

if dispBool
  disp('    Changing the order of the dimensions of an unsubdivided tree')
end

domain = domain([3 1 5 4 2]);

try
  Tree = changeDomain(origTree,domain); %#ok
  
  % Pick some random points on the domain.
  Neval = 3000;
  domrng = vertcat(domain.range);
  x = diag(diff(domrng,[],2))*rand(size(domain,1),Neval) + repmat(domrng(:,1),1,Neval);
  
  y1 = evalSurrogateModels(origTree,x,[],{domain.name}');
  y2 = evalSurrogateModels(Tree,x);
  if any(any(y1-y2 > 1e-10))
    passed = false;
    if dispBool
      disp('    ERROR: eval on modified tree produced values that appear incorrect')
    end
  end
catch
  passed = false;
  if dispBool
    disp('    ERROR: call to changeDomain failed')
  end
end

if dispBool
  disp('    Changing the domain so that a branch will be removed from a subdivided tree with no inherited surfaces')
end

Tree2 = subdivideDomain(origTree,RMCell,1,'extra1',0); %#ok
Tree3 = subdivideDomain(Tree2,RMCell,2,'p2',0.5);
Tree4 = subdivideDomain(Tree3,RMCell,3,'p1',0.15); %#ok
Tree5 = subdivideDomain(Tree4,RMCell,5,'p1',0.25); %#ok
Tree6 = subdivideDomain(Tree5,RMCell,4,'p1',-0.25); %#ok

% Originally p2 was [-0.9 1]. We cut it at 0.5 to produces nodes 4 and 5
% from 2. We change the domain so that node 4 and its children are
% eliminated.

newdomain = struct('name',{'extra1';'p2';'extra2';'p1'},'range',{[-0.9 3.9];[0.6 1];[-10 10]; [-0.5 0.5]});

try
  Tree7 = changeDomain(Tree6,newdomain); %#ok
  
  % Pick some random points on the domain.
  Neval = 3000;
  domrng = vertcat(newdomain.range);
  x = diag(diff(domrng,[],2))*rand(size(newdomain,1),Neval) + repmat(domrng(:,1),1,Neval);
  
  y1 = evalSurrogateModels(origTree,x,[],{newdomain.name}');
  y2 = evalSurrogateModels(Tree7,x);
  if any(any(y1-y2 > 1e-10))
    passed = false;
    if dispBool
      disp('    ERROR: eval on modified tree produced values that appear incorrect')
    end
  end
catch
  passed = false;
  if dispBool
    disp('    ERROR: call to changeDomain failed')
  end
end

if dispBool
  disp('    Changing the domain so that a branch will be removed from a subdivided tree with inherited surfaces')
end

% Create nodes 2 and 3 by splitting 1
Tree2 = subdivideDomain(origTree,RMCell,1,'extra1',0); 
% Create nodes 4 and 5 by splitting 2
Tree3 = subdivideDomain(Tree2,RMCell,2,'p2',0.5,[0; 1]);
% Create nodes 6 and 7 by splitting 3
Tree4 = subdivideDomain(Tree3,RMCell,3,'p1',0.15,[0; 1]); 
% Create nodes 8 and 9 by splitting 5
Tree5 = subdivideDomain(Tree4,RMCell,5,'p1',0.25,[0; 1]);
% Create nodes 10 and 11 by splitting 4
Tree6 = subdivideDomain(Tree5,RMCell,4,'p1',-0.25,[0; 1]);

% Originally p2 was [-0.9 1]. We cut it at 0.5 to produces nodes 4 and 5
% from 2. We change the domain so that node 4 and its descendants (10 and
% 11) are eliminated. node5 will become node2 pictorally (actually element
% 3 of the structure). 

newdomain = struct('name',{'extra1';'p2';'extra2';'p1'},'range',{[-0.9 3.9];[0.6 1];[-10 10]; [-0.5 0.5]});

try
  Tree7 = changeDomain(Tree6,newdomain); %#ok
  
  % Pick some random points on the domain.
  Neval = 3000;
  domrng = vertcat(newdomain.range);
  x = diag(diff(domrng,[],2))*rand(size(newdomain,1),Neval) + repmat(domrng(:,1),1,Neval);
  
  y1 = evalSurrogateModels(origTree,x,[],{newdomain.name}');
  y2 = evalSurrogateModels(Tree7,x);
  if any(any(y1-y2 > 1e-10))
    passed = false;
    if dispBool
      disp('    ERROR: eval on modified tree produced values that appear incorrect')
    end
  end
catch
  passed = false;
  if dispBool
    disp('    ERROR: call to changeDomain failed')
  end
end

% TODO since the reponse models are quadratic. This code doesn't really
% test as much as it should, since every subdivision has the same quadratic
% model (scaled for the different subdomains)


