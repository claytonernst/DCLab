function Tree6 = createPiecewiseSurrogateModelTree

addpath([fileparts(which('sourceDir')) filesep 'classTestFunctions' ...
         filesep 'ResponseModel'])

% Create an object that is as complicated as possible.
tmp = rand(3);
coeffMatrix = tmp+tmp';
modelStruct.value = coeffMatrix;
modelStruct.responseTransformation = 'log10';
modelStruct.variableTransformations = {'log10';'log10'};

domain = struct('name',{'p1';'p2'},'range',{[0.1 1];[0.1 1]});
quadRM1 = DClab.ResponseModel(modelStruct,domain);

tmp = rand(3);
coeffMatrix = tmp+tmp';
domain = struct('name',{'p2';'p3'},'range',{[-10 10];[-2 2]});
quadRM2 = DClab.ResponseModel(coeffMatrix*coeffMatrix',domain);

coeffMatrix = rand(5,1);
domain = struct('name',{'p4';'p1';'p3';'p5'},'range',{[-inf inf];[-1 1];[-2 1.4];[2 5]});
linRM = DClab.ResponseModel(coeffMatrix,domain);
simpleRM = DClab.ResponseModel(@simpleDCModel);
cplxRM = DClab.ResponseModel(@complexDCModel,'resp1');

RMCell = {quadRM1;quadRM2;simpleRM;cplxRM;linRM};

names = {'extra1';'p2';'p1';'extra2';'p3';'p5';'p4'};
ranges = {[-1 4];[0.3 1];[0.2 0.5];[-10 10];[-2 1];[3 4];[-20 20]};
domain = struct('name',names,'range',ranges);

origTree = DClab.PiecewiseSurrogateModelTree(RMCell,domain);

Tree2 = subdivideDomain(origTree,RMCell,1,'extra1',0); %#ok
Tree3 = subdivideDomain(Tree2,RMCell,2,'p2',0.5);
Tree4 = subdivideDomain(Tree3,RMCell,3,'p1',0.25); %#ok
Tree5 = subdivideDomain(Tree4,RMCell,5,'p4',0.25,[1; 1; 0; 0; 1]); %#ok
Tree6 = subdivideDomain(Tree5,RMCell,8,'p1',0.40,[1; 0; 1; 1; 1]); %#ok
