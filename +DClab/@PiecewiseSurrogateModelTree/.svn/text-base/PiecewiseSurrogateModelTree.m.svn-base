classdef PiecewiseSurrogateModelTree < DClab.DCObject
    %PiecewiseSurrogateModelTree: Construct PiecewiseSurrogateModelTree object
    %
    %   %%%% HELP IS OUTDATED %%%%
    %
    %   domain2 is slightly larger. used for training metamodel
    %
    %
    %   TREE = PiecewiseSurrogateModelTree(DOMAIN,RMCELL,MMCELL) creates an
    %   object implementing a binary tree. The tree leaves represent a
    %   subdivisions of DOMAIN, and each leaf will contain DCSurface objects
    %   for the ResponseModels provided in the m-by-1 cell array RMCELL. DOMAIN
    %   is a nx1 structure array with fields .NAME and .RANGE. MMCELL is an
    %   m-by-1 cell array that contains a DCMetamodel
    %
    %   TREE = PiecewiseSurrogateModelTree(RMCELL,DOMAIN,OPT) uses settings
    %   from the supplied DCOptions object OPT to create the surrogate models,
    %   rather then the default settings produced by the DCOptions constructor.
    %
    %   TREE = PiecewiseSurrogateModelTree(RMCELL,DOMAIN,OPT,DISABLEERRCHK)
    %   omits the check of the validity of the input arguments when
    %   DISABLEERRCHK == 1.
    %
    %   See also PolyDataset, DCSurface
    %
    %   Constructor reference page in Help browser
    %      <a href="matlab:web(fullfile(fileparts(which('DCsetup')), ...
    %      'docs', 'html', 'functions', 'DCDataset.html'), ...
    %      '-helpbrowser')">DCDataset constructor</a>
    %
    %   Object reference page in Help browser
    %      <a href="matlab:web(fullfile(fileparts(which('DCsetup')), ...
    %      'docs', 'html', 'obj_dset.html'), ...
    %      '-helpbrowser')">DCDataset object</a>

    properties %public properties
        parameterList = {};
        domainRange;
        parentNode;
        childNodes;
        cutParameter = '';
        cutLocation;
%         DCMetamodelCell = {};
        criticalRange;
        trainingRange;
        DCSurface = DClab.DCSurface;
        responseModel2LocationOfDCSurface = {};
        fittingSettings;
        displaySettings;
    end
    
    properties (Dependent)
        nSurfacesPerResponseModel
        nSurfaces
        surface2CorrespondingResponseModelIndex
        descendants
    end

    methods
        function Tree = PiecewiseSurrogateModelTree(RMCell,domain1,domain2,criticalRange,trainingRange,opt,disableErrChk)
            ni = nargin;

            switch ni
                case 0
                    %use defaults
                    return
                case {1,2,3,4}
                    error('Inputs: PiecewiseSurrogateModelTree cannot be called with the supplied number of inputs')
                case 5
                    opt = DClab.DCOptions;
                    disableErrChk = false;
                case 6
                    disableErrChk = false;
                otherwise
                    error(nargchk(0,7,ni))
            end

            % Input error checking
            if ~disableErrChk
                message = inputErrorChecking(RMCell,domain1,criticalRange,trainingRange,opt);
                if ~isempty(message)
                    error(['Inputs: ' message])
                end
            end

            m = size(RMCell,1);

            % The root node is special. It contains the parameterList and the
            % DCMetamodel (if exists).

            Tree(1).parameterList = {domain1.name}';
            Tree.domainRange = vertcat(domain1.range);
            Tree.parentNode = -1; %This element is root
            Tree.childNodes = [];
            Tree.cutParameter = '';
            Tree.cutLocation = [];
%             Tree(1).DCMetamodelCell = cell(m,1);
            Tree(1).criticalRange = criticalRange;
            Tree.trainingRange = trainingRange;

            if m == 0
                surfCell = {DClab.DCSurface};
                responseModel2LocationOfDCSurface = {};
            else

                % Construct metamodels (if needed) and determine the I/0 transformations
                % used for the quadratic fits.

                trans = struct('response',cell(m,1),'variable',cell(m,1));

                % The contents of the ith element of nSurfs4RM indicate how many surfaces
                % (each with different I/O transformations) will be constructed for the
                % response model in the ith element of RMCell;
                nSurfs4RM = zeros(m,1);

                for i1 = 1:m

%                     if strcmp(RMCell{i1}.type,'dcModel') && ~strcmp(opt.analysisMode,'original')
%                         % We need a metamodel for analysis. Construct it, and then just use
%                         % the same variables transformation that are employed by the
%                         % metamodel in the response surfaces.
% 
%                         [trash dsetNames2sortedSharedNames modelNames2sortedSharedNames] = intersect(Tree(1).parameterList,RMCell{i1}.parameterList);
%                         dsetNames2ModelNames = zeros(size(trash));
%                         dsetNames2ModelNames(modelNames2sortedSharedNames,1) = dsetNames2sortedSharedNames;
% 
%                         mmDomain1 = domain1(dsetNames2ModelNames);
%                         mmDomain2 = domain2(dsetNames2ModelNames);
% 
%                         % Create the metamodel. 2nd domain is slightly larger for training.
%                         MM = DClab.DCMetamodel(RMCell{i1},mmDomain1,mmDomain2,criticalRange(i1,:),trainingRange(i1,:),opt);
% 
%                         % Determine the variable transformations.
%                         trans(i1).response = {MM.responseTransformation};
%                         trans(i1).variable = repmat({'none'},size(domain1));
% 
%                         MMvarTransNum = MM.variableTransformations;
%                         MMvarTrans = repmat({'none'},size(MMvarTransNum));
%                         MMvarTrans(MMvarTransNum==2) = {'log10'};
% 
%                         trans(i1).variable(dsetNames2ModelNames,1) = MMvarTrans;
%                         nSurfs4RM(i1) = 1;
%                         Tree(1).DCMetamodelCell{i1} = MM;
%                         Tree(1).DCMetamodelCell{i1}.userData = [criticalRange(i1,:); trainingRange(i1,:)];
%                     else
                        trans(i1) = findSurrogateIOTransformations(RMCell{i1},domain1,opt);
                        nSurfs4RM(i1) = size(trans(i1).response,2);
                        if nSurfs4RM(i1)>1
                            error('code not complete')
                        end
%                         Tree(1).DCMetamodelCell{i1} = DClab.DCDataset;
%                         Tree(1).DCMetamodelCell{i1}.userData = [criticalRange(i1,:); trainingRange(i1,:)];

%                     end

                end

                % Let NODE = PiecewiseSurrogateModelTree(end).
                % NODE.responseModel2LocationOfDCSurface is a m-by-1 cell array. Each cell
                % contains a nSurfs4RM(i)-by-2 array, call it 'findSurfs{i}', where
                % nSurfs4RM(i) is the number of distinct surfaces that are constructed for
                % the i^th element of RMCell over NODE. Let [j k] be the contents of a row
                % of findSurfs{i}. j indicates a element of the multidimensional
                % PiecewiseSurrogateModelTree that contains a surface for the i^th response
                % model over NODE. j may not equal the index of NODE, since surfaces can be
                % "inherited", rather then recomputed, when the domain is partitioned. k
                % indicates which element of PiecewiseSurrogateModelTree(j).DCSurface
                % contains this surface. I.e., PiecewiseSurrogateModelTree(j).DCSurface(k)
                % is a response surface over NODE for the i^th element of RMCell. Since
                % this constructor creates surrogates having an unpartitioned domain, the
                % elements of findSurfs are pretty simple since PiecewiseSurrogateModelTree
                % is scalar and no surfaces are inherited.

                responseModel2LocationOfDCSurface = cell(m,1);

                try

                    % Initialize a multidimensional DCSurface object.
                    surfCell = cell(sum(nSurfs4RM),1);
                    for i1 = 1:m

                        if i1 == 1
                            lastSurfIdx = 0;
                        else
                            lastSurfIdx = sum(nSurfs4RM(1:i1-1));
                        end

                        for i2 = 1:nSurfs4RM(i1)
                            surfCell{lastSurfIdx+i2} = DClab.DCSurface(domain1,RMCell{i1},criticalRange(i1,:),trainingRange(i1,:),trans(i1).response{1,i2},trans(i1).variable(:,i2),i1,opt,disableErrChk);
%                             surfCell{lastSurfIdx+i2} =
%                             DClab.DCSurface(domain1,RMCell{i1},Tree(1).DCMetamodelCell{i1},criticalRange(i1,:),trainingRange(i1,:),trans(i1).response{1,i2},trans(i1).variable(:,i2),i1,opt,disableErrChk);
                        end

                        newSurfsIdx = lastSurfIdx + (1:nSurfs4RM(i1))';
                        responseModel2LocationOfDCSurface{i1} = [repmat(1,nSurfs4RM(i1),1) newSurfsIdx];

                    end
                catch
                    disp('Problem making DCSurfaces')
                    keyboard
                end

            end

            Tree.DCSurface = vertcat(surfCell{:});
            Tree.responseModel2LocationOfDCSurface = responseModel2LocationOfDCSurface;

            [dispOpt fitOpt] = decompose(opt);
            Tree.fittingSettings = fitOpt;
            Tree(1).displaySettings = dispOpt;
        end

        %hard scalar properties
%         function out = get.parameterList(obj)
%             if obj(1).parentNode ~=-1
%                 error('Property of the root node only')
%             end
%             out = obj(1).parameterList; %does this work?
%         end
%         function out = get.DCMetamodelCell(obj)
%             if obj(1).parentNode ~=-1
%                 error('Property of the root node only')
%             end
%             out = obj(1).DCMetamodelCell; %does this work?
%         end
%         function out = get.criticalRange(obj)
%             if obj(1).parentNode ~=-1
%                 error('Property of the root node only')
%             end
%             out = obj(1).criticalRange; %does this work?
%         end
%         function out = get.displaySettings(obj)
%             if obj(1).parentNode ~=-1
%                 error('Property of the root node only')
%             end
%             out = obj(1).displaySettings; %does this work?
%         end
        
        %scalar dependent
        function out = nNodes(obj)
            if isempty(obj)
                out = 0;
            else
                out = builtin('length',obj);
            end
        end
        function out = leafNodes(obj)
            out = sort(findLeaves(obj,1));
        end
        function out = nSubdivisions(obj)
            out = length(findLeaves(obj,1));
        end
        function out = nResponseModels(obj)
            out = size(obj(1).responseModel2LocationOfDCSurface,1);
        end

        %dependent properties ...
        function out = get.nSurfacesPerResponseModel(obj)
            NS = size(vertcat(obj.responseModel2LocationOfDCSurface{:}),1);
            out = zeros(NS,1);
            for i2 = 1:NS
                out(i2) = size(obj.responseModel2LocationOfDCSurface{i2},1);
            end
        end
        function out = get.nSurfaces(obj)
            out = size(vertcat(obj.responseModel2LocationOfDCSurface{:}),1);
        end
        function out = get.surface2CorrespondingResponseModelIndex(obj)
            nSurf = size(vertcat(obj.responseModel2LocationOfDCSurface{:}),1);
            RMidx = [];
            for i2 = 1:nSurf
                RMidx = [RMidx; i2*ones(size(obj.responseModel2LocationOfDCSurface{i2},1),1)];
            end
            out = RMidx;
        end
        function out = get.descendants(Tree)
            out = sort(findDescendants(Tree,1));
        end
        
        function out = surfaces(obj,node)
            out = getSurfs(obj,node);
        end

        function bool = isempty(obj)
            error(nargoutchk(0,1,nargout))
            error(nargchk(1,1,nargin))

            if isempty(obj(1).parentNode)
                bool = true;
            else
                bool = false;
            end
        end
        function out = length(obj)
            error('LENGTH is not a valid method. Use nNodes(Tree).')
        end
        function out = size(obj,dim)
            error('SIZE is not a valid method. Use nNodes(Tree).')
        end
        function out = vertcat(varargin)
            out = builtin('vertcat',varargin{:});
        end
        function [list,sz] = displayProps(obj)
            list = [];
            sz = sprintf('%d-leaf',length(leafNodes(obj)));
        end
        function PSMTree = uminus(PSMTree)
            for i1 = 1:nNodes(PSMTree)
                PSMTree(i1).DCSurface = -PSMTree(i1).DCSurface;
            end

%             for i1 = 1:length(PSMTree(1).DCMetamodelCell)
%                 if ~isempty(PSMTree(1).DCMetamodelCell{i1})
%                     PSMTree(1).DCMetamodelCell{i1} = -PSMTree(1).DCMetamodelCell{i1};
%                 end
%             end

            PSMTree(1).criticalRange = [-PSMTree(1).criticalRange(2) -PSMTree(1).criticalRange(1)];
            PSMTree(1).trainingRange = [-PSMTree(1).trainingRange(2) -PSMTree(1).trainingRange(1)];
        end

    end %public methods

end %classdef


%%%%%%%%%%%%%%%%%LOCAL FUNCTIONS%%%%%%%%%%%%%%%%%%
function message = inputErrorChecking(RMCell,domain,criticalRange,trainingRange,opt)

message = '';

if isempty(RMCell)
    %do nothing, this is allowed
elseif ~iscell(RMCell) || ~all(cellfun('isclass',RMCell,'DClab.ResponseModel')) || size(RMCell,2)~=1
    message = 'RMCELL must be a column cell array of ResponseModels';
    return
else
    % OK
end

[bool tmp] = DClab.isValidDomain(domain,'DOMAIN');
if ~bool
    message = tmp;
    return
elseif any(any(isinf([domain.range])))
    message = 'DOMAIN must describe a bounded set';
    return
else
    %do nothing, domain should be OK
end

m = size(RMCell,1);
if m == 0
    if ~isempty(criticalRange) || ~isempty(trainingRange)
        message = 'When RMCELL is empty, criticalRange and trainingRange must be as well';
        return
    end
else
    if ~isequal(size(criticalRange),[m 2]) || ~isequal(size(trainingRange),[m 2])
        message = 'CRITICALRANGE and TRAININGRANGE must be m-by-2';
        return
    end
    if any(criticalRange(:,2) <= criticalRange(:,1))
        message = 'In each row of CRITICALRANGE, the 2nd element must be strictly greater than the 1st.';
        return
    end
    if any(trainingRange(:,2) <= trainingRange(:,1))
        message = 'In each row of TRAININGRANGE, the 2nd element must be strictly greater than the 1st.';
        return
    end
end

if ~isa(opt,'DClab.DCOptions')
    message = 'OPT must be of class DCOptions';
end
end

function Surf = getSurfs(Tree,node)
%GETSURFS returns a DCSurface object containing all surfaces defined on a node
%
%   SURF = GETSURFACES(TREE,NODE) returns a multidimensional DCSurface
%   object that contains the DCSurfaces that constitute the surrogate
%   models over the node NODE of the PiecewiseSurrogateModelTree TREE.

RM2LocOfSurf = Tree(node).responseModel2LocationOfDCSurface;
nSurfs4RM = Tree(node).nSurfacesPerResponseModel;

% Initialize a multidimensional DCSurface object.
surfCell = cell(sum(nSurfs4RM),1);

paramList = Tree(1).parameterList;
newDomain = DClab.createDomainStructure(paramList,Tree(node).domainRange);

for i1 = 1:length(nSurfs4RM)
    if i1 == 1
        lastSurfIdx = 0;
    else
        lastSurfIdx = sum(nSurfs4RM(1:i1-1));
    end

    tmp = RM2LocOfSurf{i1};
    for i2 = 1:nSurfs4RM(i1)
        j = tmp(i2,1);
        k = tmp(i2,2);
        % Get the most recent surface
        surf = Tree(j).DCSurface(k);
        if j == node;
            surfCell{lastSurfIdx+i2} = surf;
        else
            % It wasn't fit over the domain of node, so rescale it.
            oldDomain = DClab.createDomainStructure(paramList,Tree(j).domainRange);
            surfCell{lastSurfIdx+i2} = changeDomain(surf,oldDomain,newDomain);
        end
    end
end
Surf = vertcat(surfCell{:});
end

function leaves = findLeaves(Tree,idx)
% This function is recursive. To save time, the final output will not be
% 'sorted'. You can do one sort with the output of this function if
% desired.

if isempty(Tree(idx).childNodes)
    leaves = idx;
else
    leaves = [findLeaves(Tree,Tree(idx).childNodes(1)); findLeaves(Tree,Tree(idx).childNodes(2))];
end
end

function desc = findDescendants(Tree,idx)
% This function is recursive. To save time, the final output will not be
% 'sorted'. You can do one sort with the output of this function if
% desired.

if isempty(Tree(idx).childNodes)
    desc = [];
else
    childNodes = Tree(idx).childNodes';
    desc = [childNodes ; findDescendants(Tree,childNodes(1)); findDescendants(Tree,childNodes(2))];
end
end
