addpath([fileparts(which('sourceDir')) filesep 'classTestFunctions'])
PD = createPolyDataset;

tmp = rand(4);
coeffMatrix = tmp+tmp';
domain = struct('name',{'p2';'p5';'p3'},'range',{[0 10];[-2 5];[-4 1]});
quadRM = ResponseModel(coeffMatrix*coeffMatrix',domain);

domain = DClab.createDomainStructure(PD.FreeParameter.name,PD.FreeParameter.range);

objPSMTree = PiecewiseSurrogateModelTree({quadRM},domain);
[objPSMTree trash] = mergePartitionStructures(objPSMTree,PD.PiecewiseSurrogateModelTree);

curDir = pwd;
cd([fileparts(which('sourceDir')) filesep '@ResponsePrediction' filesep 'private'])
obj = makeNqcqp(objPSMTree,4,1,PD,'approx');
cd(curDir)

lowerBnd(obj)
upperBnd(obj)
