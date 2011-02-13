function Tree = subdivideDomain(Tree,RMCell,node,paramName,loc,inherit,opt)
%SUBDIVIDEDOMAIN subdivides the domain of the surrogate models
%
%   NEWTREE = SUBDIVIDEDOMAIN(TREE,RMCELL,NODE,PARAMNAME,LOC) subdivides
%   the subdivision represented by the TREE(NODE) (which must be a leaf) of
%   the subdivided domain of the PiecewiseSurrogateModelTree TREE. The
%   domain of TREE(NODE) will be split in the dimension PARAMNAME at the
%   location LOC. This causes new DCSurfaces to be constructed over the new
%   subdivisions for each of the response models in the column cell array
%   RMCELL, effectively refining the surrogate models for these response
%   models. For NEWTREE to make sense, RMCELL must be the same collection
%   (in the same order) of ResponseModels that were originally used to
%   create TREE. The display and fitting options used to construct the new
%   DCSurfaces will be inherited from those of TREE(NODE).
%
%   NEWTREE = SUBDIVIDEDOMAIN(TREE,RMCELL,NODE,PARAMNAME,LOC,INHERIT)
%   allows you to specify if the DCSurfaces defined over the resulting
%   subdivisions should be recomputed or if they should be inherited from
%   their parents. See input description below.
%
%   NEWTREE = SUBDIVIDEDOMAIN(TREE,RMCELL,NODE,PARAMNAME,LOC,INHERIT,OPT)
%   allows you to supply a DCOptions object OPT that will dictate the
%   display and fitting options used to construct the new DCSurfaces.
%
%   Inputs
%      TREE: A PiecewiseSurrogateTree object.
%      RMCELL: a column cell array of ResponseModel objects. This must be
%         the same cell array that was originally used to create TREE. If
%         INHERIT=TRUE, RMCELL=[] is acceptible.
%      NODE: An integer that specifies which subdivision (i.e., element of
%         TREE) will be subdivided. It must be a member of
%         get(TREE,'leafNodes').
%      PARAMNAME: This is a char array specifying the dimension that the
%         partitioning plane passes though. PARAMNAME must be the name of a
%         parameter in get(TREE,'parameterList').
%      LOC: This scalar specifies the location of the partitioning plane
%         that passes though the dimension indicated by PARAMNAME.
%      INHERIT: A column vector consisting of ones and zeros. Its length
%         must be equal to the number of DCSurfaces that are defined on
%         NODE. I.e., equal to length(getSurfaces(TREE,NODE)). All
%         DCSurfaces that exist on NODE and correspond to components of
%         INHERIT that equal 1 will be inherited by the new subdivisions.
%         For components of INHERIT that equal 0, this function determines
%         which elements of RMCELL they correspond to and then creates
%         DCSurfaces over the new subdivisions for these response models.
%         This process might involve the construction of multiple
%         DCSurfaces for any given response model, and their I/O
%         transformations may differ from those of the previous
%         noninherited DCSurfaces. The sum total of these events mean that
%         a given ResponseModel may have several DCSurfaces over a new
%         subdivision, some of which were simply inherited, and some of
%         which were newly created. In the rare case where I/O
%         transformations for which a DCSurface was inherited are also used
%         by a newly created DCSurface, the inherited surface will be
%         thrown away and only the new surface will remain. If INHERIT is a
%         scalar 0, it will be interpreted as zeros(nSurf,1).
%      OPT: A DCOptions object. 
%
%   Outputs
%      NEWTREE: The resulting PiecewiseSurrogateTree object.
%
% See also PiecewiseSurrogateModelTree


%TODO would we ever want to subdivide a node that is not a leaf? Currently
%such an attempt will bomb, but if so, we could call this method
%recursively.

error(nargoutchk(0,1,nargout));
ni = nargin;
switch ni
  case 5
    inherit = [];
    opt = [];
  case 6
    opt = [];
  otherwise
    error(nargchk(5,7,ni));
end

% Input error checking and initialization of empty inputs
[message inherit opt] = checkAndInitializeInputs(Tree,RMCell,node,paramName,loc,inherit,opt);
if ~isempty(message)
  error(message)
end

% Begin actual code

N = nNodes(Tree);
dim = strmatch(paramName,char(Tree(1).parameterList),'exact');

% Update this leaf to inform it of its children (to be created)
Tree(node).childNodes = [N+1 N+2];
Tree(node).cutParameter = paramName;
Tree(node).cutLocation = loc;

%We need to find which RMs need surfaces to be recomputed for them.
%RMidx will be nSurf-by-1 and its elements indicate which RM the
%corresponding surface is for.
RMidx = Tree(node).surface2CorrespondingResponseModelIndex;

% This variable will look like
% Tree(node).responseModel2LocationOfDCSurfaces, except its cell contents
% will only include rows that reference surfaces the child leaves will
% inherit. 
responseModel2LocationOfInheritedDCSurface = cell(size(Tree(node).responseModel2LocationOfDCSurface));

% This variable will contain the indices of the response models that need
% new surfaces fit over the child leaves.
RMneedingSurfs = [];

tmp = vertcat(Tree(node).responseModel2LocationOfDCSurface{:});
for i1 = 1:length(inherit)
  if inherit(i1) == 0
    RMneedingSurfs = [RMneedingSurfs; RMidx(i1)];
  else
    % Fill out responseModel2LocationOfDCSurface referencing the inherited
    % surfaces
    responseModel2LocationOfInheritedDCSurface{RMidx(i1)} = [responseModel2LocationOfInheritedDCSurface{RMidx(i1)}; tmp(i1,:)];
  end
end
RMneedingSurfs = unique(RMneedingSurfs);

% Let's do the spliting. We create two new DCBinaryTreeNodes. One called
% lesser and one called greater.
criticalRange = Tree(1).criticalRange;
trainingRange = Tree(1).trainingRange;
% MMCell = Tree(1).DCMetamodelCell;

[displaySettings fittingSettings] = decompose(opt);

%Create the lesser child node
lesser = DClab.PiecewiseSurrogateModelTree;
domrng = Tree(node).domainRange;
domrng(dim,2) = loc;
lesser.domainRange = domrng;

domain = DClab.createDomainStructure(Tree(1).parameterList,lesser.domainRange);

lesser.parentNode = node;
[lesser.DCSurface lesser.responseModel2LocationOfDCSurface] = makeSurfs(domain,RMCell,criticalRange,trainingRange,N+1,RMneedingSurfs,responseModel2LocationOfInheritedDCSurface,opt);
lesser.fittingSettings = fittingSettings;

%Create the greater child node
greater = DClab.PiecewiseSurrogateModelTree;
domrng = Tree(node).domainRange;
domrng(dim,1) = loc;
greater.domainRange = domrng;

domain = DClab.createDomainStructure(Tree(1).parameterList,greater.domainRange);

greater.parentNode = node;
[greater.DCSurface greater.responseModel2LocationOfDCSurface] = makeSurfs(domain,RMCell,criticalRange,trainingRange,N+2,RMneedingSurfs,responseModel2LocationOfInheritedDCSurface,opt);
greater.fittingSettings = fittingSettings;

Tree = builtin('vertcat',Tree,lesser,greater);

% Update the display settings. To save memory, these are only stored in the
% first element of the Tree.
Tree(1).displaySettings = displaySettings;

%TODO go through all surfs for a given RM and eliminate any inherited ones
%that share transformations with a new one.

% ====== local functions ======

function [message inheritOUT optOUT] = checkAndInitializeInputs(Tree,RMCell,node,paramName,loc,inherit,opt)

message = '';
inheritOUT = [];
optOUT = [];

if isempty(Tree)
  message = 'Inputs: an empty TREE cannot be subdivided';
  return
end

% We check RMCell last, after inheritOUT is established

if ~isscalar(node) || ~isnumeric(node) || node ~= round(node)
  message = 'Inputs: NODE is not a scalar index';
  return
end

if ~ismember(node,leafNodes(Tree))
  message = ['Inputs: NODE: ' num2str(node) ' is not a member of TREE.leafNodes'];
  return
end

if ~ischar(paramName) || size(paramName,1) ~= 1
  message = 'Inputs: PARAMNAME must be a single row char';
  return
end

if ~ismember(paramName,Tree(1).parameterList)
  message = ['Inputs: PARAMNAME: ' paramName ' is not a member of TREE.parameterList'];
  return
end

idx = strmatch(paramName,char(Tree(1).parameterList),'exact');
if loc < Tree(node).domainRange(idx,1) || loc > Tree(node).domainRange(idx,2)
  message = 'Inputs: LOC is outside of the domain';
  return
end

nSurf = Tree(node).nSurfaces;
if isempty(inherit) || isequal(inherit,0)
  inheritOUT = zeros(nSurf,1);  
elseif ~isequal(inherit == 1 | inherit == 0,ones(nSurf,1))
  message = 'When ~= 0, INHERIT must be a column vector of 0''s and 1''s with length = get(TREE,''nSurfaces'',NODE)';
  return
else
  inheritOUT = inherit;
end

if isempty(opt)
  dispOpts = Tree(1).displaySettings;
  fitOpts = Tree(node).fittingSettings;
  opts = cell2struct([struct2cell(fitOpts);struct2cell(dispOpts)],[fieldnames(fitOpts);fieldnames(dispOpts)]);
  optOUT = gui(opts);
elseif ~isa(opt,'DClab.DCOptions')
  message = 'Inputs: OPT must be a DCOptions object';
  return
else
  optOUT = opt;
end

if isempty(RMCell)
  if ~all(inheritOUT)
    message = 'Inputs: RMCELL can only be empty if INHERIT==true';
  else
    %do nothing, when inherit is true, RMCell=[] is ok
  end
elseif ~all(cellfun('isclass',RMCell,'DClab.ResponseModel')) || size(RMCell,2)~=1
  message = 'Inputs: RMCELL must be a column cell array of ResponseModels';
else
  %do nothing
end

%   dispOpts = PD.dispOpts;
%   opts = cell2struct([struct2cell(fitOpts);struct2cell(dispOpts)],[fieldnames(fitOpts);fieldnames(dispOpts)]);

function [surf responseModel2LocationOfDCSurface] = makeSurfs(domain,RMCell,critRng,trnRng,node,RMneedingSurfs,responseModel2LocationOfInheritedDCSurface,opt)

mm = length(RMneedingSurfs);

% Determine the I/0 transformations to use for the new lesser DCSurfaces.
trans = struct('response',cell(mm,1),'variable',cell(mm,1));

% The contents of the ith element of nSurfs4RM indicate how many surfaces
% (each with different I/O transformations) will be constructed for the
% response model indicated by the ith element of RMneedingSurfs. 
nSurfs4RM = zeros(mm,1);

paramList = {domain.name}';
for i1 = 1:mm
  mIdx = RMneedingSurfs(i1);
  
%   if strcmp(RMCell{mIdx}.type,'dcModel') && ~strcmp(opt.analysisMode,'original')
%     % We need a metamodel for analysis. Construct it, and then just use
%     % the same variables transformation that are employed by the
%     % metamodel in the response surfaces.
% 
%     [trash dsetNames2sortedSharedNames modelNames2sortedSharedNames] = intersect(paramList,RMCell{mIdx}.parameterList);
%     dsetNames2ModelNames = zeros(size(trash));
%     dsetNames2ModelNames(modelNames2sortedSharedNames,1) = dsetNames2sortedSharedNames;
% 
%     % Determine the variable transformations.
%     trans(i1).response = {MMCell{mIdx}.responseTransformation};
%     trans(i1).variable = repmat({'none'},size(domain));
%  
%     
% %    disp('hacked surface transformations in /subdivideDomain')
%     %ucomment next to restore.
%     MMvarTransNum = MMCell{i1}.variableTransformations;
%     MMvarTrans = repmat({'none'},size(MMvarTransNum));
%     MMvarTrans(MMvarTransNum==2) = {'log10'};
%     trans(i1).variable(dsetNames2ModelNames,1) = MMvarTrans;
% 
%     nSurfs4RM(i1) = 1;
% 
%   else
    trans(i1) = findSurrogateIOTransformations(RMCell{mIdx},domain,opt);
    nSurfs4RM(i1) = size(trans(i1).response,2);
%   end
end

% Initialize a multidimensional DCSurface object.
if sum(nSurfs4RM)==0
  surf = DClab.DCSurface; %TODO DCSurface, and DCSurface(1) return the same thing. Is this ok? Should DCSurface(0) do something?
else
  surf = DClab.DCSurface(sum(nSurfs4RM));
end
responseModel2LocationOfDCSurface = responseModel2LocationOfInheritedDCSurface;
for i1 = 1:mm
  if i1 == 1
    lastSurfIdx = 0;
  else
    lastSurfIdx = sum(nSurfs4RM(1:i1-1));
  end
  
  mIdx = RMneedingSurfs(i1);
  
  for i2 = 1:nSurfs4RM(i1)
    surf(lastSurfIdx+i2) = DClab.DCSurface(domain,RMCell{mIdx},critRng(mIdx,:),trnRng(mIdx,:),trans(i1).response{1,i2},trans(i1).variable(:,i2),mIdx,opt);
%     surf(lastSurfIdx+i2) =
%     DClab.DCSurface(domain,RMCell{mIdx},MMCell{mIdx},critRng(mIdx,:),trnRng(mIdx,:),trans(i1).response{1,i2},trans(i1).variable(:,i2),mIdx,opt);
  end

  newSurfsIdx = lastSurfIdx + (1:nSurfs4RM(i1))';
  responseModel2LocationOfDCSurface{mIdx} = [responseModel2LocationOfDCSurface{mIdx}; [repmat(node,nSurfs4RM(i1),1) newSurfsIdx]];  
end


