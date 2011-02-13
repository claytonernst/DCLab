addpath([fileparts(which('sourceDir')) filesep 'classTestFunctions'])
PD = createPolyDataset;

tmp = rand(4);
coeffMatrix = tmp+tmp';
domain = struct('name',{'p2';'p5';'p3'},'range',{[0 10];[-2 5];[-4 1]});
quadRM = ResponseModel(coeffMatrix*coeffMatrix',domain);

keyboard

obj = ResponsePrediction(quadRM,PD);