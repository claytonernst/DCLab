function PD6 = createPolyDataset

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
quadRM3 = DClab.ResponseModel(coeffMatrix*coeffMatrix',domain);
simpleRM = DClab.ResponseModel(@simpleDCModel);
cplxRM = DClab.ResponseModel(@complexDCModel,'resp1');

%RMCell = {quadRM1;quadRM2;simpleRM;cplxRM;linRM};
RMCell = {quadRM1;quadRM2;simpleRM;cplxRM;quadRM3};



names = {'extra1';'p2';'p1';'extra2';'p3';'p5';'p4'};
ranges = {[-1 4];[0.3 1];[0.2 0.5];[-10 10];[-2 1];[3 4];[-20 20]};
domain = struct('name',names,'range',ranges);

% To get some (hopefully) consistent data, evaluate each at the center
% point, and two opposite corners.
nRM = size(RMCell,1);
y = zeros(nRM,3);
rng = vertcat(ranges{:});
for i1 = 1:nRM
  y(i1,1) = eval(RMCell{i1},rng(:,1),{domain.name}');
  y(i1,2) = eval(RMCell{i1},mean(rng,2),{domain.name}');
  y(i1,3) = eval(RMCell{i1},rng(:,2),{domain.name}');
end

% Create random data and uncertainty that probably makes us consistent.
d = zeros(nRM,1);
u = zeros(nRM,1);
for i1 = 1:nRM
  mx = max(y(i1,:));
  mn = min(y(i1,:));
  d(i1) = mn + rand(1)*(mx-mn);
  u(i1) = 0.75*(mx-mn);
end

% Create the MOPairs
Pairs = DClab.ModelAndObservationPair(nRM);
for i1 = 1:nRM
  Pairs(i1) = DClab.ModelAndObservationPair(ResponseObservation(d(i1),u(i1)),RMCell{i1},num2str(i1));
end

%Create the FreeParameters.
names = {'extra1';'p2';'p1';'extra2';'p3';'p5';'p4'};
ranges = {[-1 4];[0.3 1];[0.2 0.5];[-10 10];[-2 1];[3 4];[-20 20]};
FParam = FreeParameter(length(names));
for i1 = 1:length(names)
  nom = mean(ranges{i1});
  unc = nom-ranges{i1}(1,1);
  FParam(i1) = FreeParameter(names{i1},nom,unc);
end

Dset = DCDataset(Pairs,FParam);
PD = PolyDataset(Dset);

% Do some subdividing with various inheritance states.
PD2 = subdivideDomain(PD,1,'extra1',0);
PD3 = subdivideDomain(PD2,2,'p2',0.5);
PD4 = subdivideDomain(PD3,3,'p1',0.25);
PD5 = subdivideDomain(PD4,5,'p4',0.25,[1; 1; 0; 0; 1]);
PD6 = subdivideDomain(PD5,8,'p1',0.40,[1; 0; 1; 1; 1]);



