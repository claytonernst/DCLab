function newTree = Tree4RM(oldTree,modelIdx)
%TREE4RM PiecewiseSurrogateModelTree for specific model
%
%  newTree = Tree4RM(oldTreeObj,modelIdx) take the
%  PiecewiseSurrogateModelTree oldTreeObj and creates a new
%  PiecewiseSurrogateModelTree that only corresponds to the single
%  ResponseModel indexed by modelIdx.

%initialize new object
n = nNodes(oldTree);
newTree(n) = DClab.PiecewiseSurrogateModelTree;

%fill object, but only appropriate surfaces
for i=1:n
    
    newTree(i).parameterList = oldTree(i).parameterList;
    newTree(i).domainRange = oldTree(i).domainRange;
    newTree(i).parentNode = oldTree(i).parentNode;
    newTree(i).childNodes = oldTree(i).childNodes;
    newTree(i).cutParameter = oldTree(i).cutParameter;
    newTree(i).cutLocation = oldTree(i).cutLocation;
    newTree(i).criticalRange = oldTree(i).criticalRange;
    newTree(i).trainingRange = oldTree(i).trainingRange;
    RM2LocSurf =oldTree(i).responseModel2LocationOfDCSurface(modelIdx);
    surfIdx = RM2LocSurf{1}(:,2);
    RM2LocSurf{1}(:,2) = (1:size(RM2LocSurf{1},1))';
    newTree(i).responseModel2LocationOfDCSurface = RM2LocSurf;
    newTree(i).DCSurface = oldTree(i).DCSurface(surfIdx);
    newTree(i).fittingSettings = oldTree(i).fittingSettings;
    newTree(i).displaySettings = oldTree(i).displaySettings;

end