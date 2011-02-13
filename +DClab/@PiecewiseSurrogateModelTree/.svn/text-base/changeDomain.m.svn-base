function Tree = changeDomain(Tree,domain)
%CHANGEDOMAIN alters the domain a PiecewiseSurrogateModelTree
%
%   TREE = CHANGEDOMAIN(TREE,DOMAIN) changes the domain of the
%   PiecewiseSurrogateModelTree TREE. This will involve modifying or even
%   eliminating some of the subdivisions of the domain and performs scaling
%   operations on the algebraic surrogate models defined over any modified
%   subdivisions. DOMAIN should be column structure arrays with fields
%   .name and .range. It can change the order of the dimensions (parameter 
%   names) and/or add additional dimensions not present in the original
%   domain but cannot eliminate any existing dimensions.
%
% See also PiecewiseSurrogateModelTree, PolyDataset/changeDomain


%change bounds on root
%go through tree changing bounds where nec. an existing cut is no
%longer in the cube, remove entire subtree correpsonding to missing
%child and bump up other child to current position

ni = nargin;
no = nargout;
error(nargchk(2,2,ni));
error(nargoutchk(0,1,no));

[bool message] = DClab.isValidDomain(domain);
if bool
  error(message)
end

if any(any(isinf(vertcat(domain.range))))
  error('Inputs: DOMAIN must describe a bounded set')
end

newParamList = {domain.name}';
oldParamList = Tree(1).parameterList;

% Special Case:
if isequal(newParamList,oldParamList) && isequal(vertcat(domain.range),Tree(1).domainRange)
  % Nothing to do, make a quick exit.
  return
end

if ~isempty(setdiff(oldParamList,newParamList))
  error('The supplied {DOMAIN.name} must be a superset of TREE.parameterList')
end

[trash nPList2sort oPList2sort] = intersect(newParamList,oldParamList);

% DCMetamodelCell = Tree(1).DCMetamodelCell;
criticalRange = Tree(1).criticalRange;
trainingRange = Tree(1).trainingRange;

%recurvisely traverse the tree and fix everything
[Tree junkNodes] = changeSubdomain(Tree,1,domain,nPList2sort,oPList2sort,[]);

if length(junkNodes) ~= length(unique(junkNodes))
  error('Code has a bug')
end

% Now we need to blow away the junk nodes. This involves fixing the node
% references in the responseModel2LocationOfDCSurface to correspond to the
% components of Tree that will remain after blowing away the junk nodes and
% then blowing away the junk nodes. 

oldNodes = (1:nNodes(Tree))';
goodNodes = setdiff(oldNodes,sort(junkNodes));

% newIdx2oldIdx is a vector of old nodes, and the location of an old node
% number is its new node number. For example, newIdx2oldIdx(5) gives the
% old node number for what will become the new node 5 in the final Tree.
newIdx2oldIdx = goodNodes; 

% Note, the order that we change things is important. If we and to change
% 5-->3 and 3-->2, we must to the 2nd (earlier numerically) change first.

%TODO, an actual tree traversal would be a more efficient way to accomplish
%this updating. Fix it if its slow.

% Fix each of the good nodes one at a time.
for i1 = 1:length(goodNodes)

  if ~isempty(Tree(goodNodes(i1)).responseModel2LocationOfDCSurface)
    % Fix the node references in the .responseModel2LocationOfDCSurface
    % array.
    tmp = Tree(goodNodes(i1)).responseModel2LocationOfDCSurface;
    rowsPerCell = 0.5*cellfun('prodofsize',tmp);

    % Stack the contents into a matrix because we don't want to loop throw
    % the cells.
    tmp = vertcat(tmp{:});

    for i2 = 1:length(newIdx2oldIdx)
      % Fix the node references in the .responseModel2LocationOfDCSurface
      % array.
      idx = tmp(:,1) == newIdx2oldIdx(i2);
      tmp(idx,1) = i2;
    end

    % Recover the cell array with the first columns of the contents updated.
    tmp = mat2cell(tmp,rowsPerCell,2);
    Tree(goodNodes(i1)).responseModel2LocationOfDCSurface = tmp;
  end
    
  % Fix the parent node references. Recall the root node has parent == -1.
  % These below lines won't screw this up. If the root node is nodified,
  % changeSubdomain will have taken care of making the parent of the new
  % root -1.
  for i2 = 1:length(newIdx2oldIdx)
    if Tree(goodNodes(i1)).parentNode == newIdx2oldIdx(i2);
      Tree(goodNodes(i1)).parentNode = i2;
    end
  end
  
  % Fix the child node references.
  for i2 = 1:length(newIdx2oldIdx)
    if ~isempty(Tree(goodNodes(i1)).childNodes)
      if Tree(goodNodes(i1)).childNodes(1) == newIdx2oldIdx(i2);
        Tree(goodNodes(i1)).childNodes(1) = i2;
      end
      if Tree(goodNodes(i1)).childNodes(2) == newIdx2oldIdx(i2);
        Tree(goodNodes(i1)).childNodes(2) = i2;
      end
    end
  end
end

% Update the domain of any metamodels.
% for i1 = 1:length(DCMetamodelCell)
%   if ~isempty(DCMetamodelCell{i1})
%     DCMetamodelCell{i1} = changeDomain(DCMetamodelCell{i1},domain);
%   end
% end

% Blow away the junk nodes and update the parameter list, metamodels, etc,,
% which always resides only on the first element of Tree in order to save
% memory. 
Tree(junkNodes) = [];
Tree(1).parameterList = newParamList;
% Tree(1).DCMetamodelCell = DCMetamodelCell;
Tree(1).criticalRange = criticalRange;
Tree(1).trainingRange = trainingRange;

% ====== local function ======

function [Tree junkNodes] = changeSubdomain(Tree,currentNode,domain,nPList2sort,oPList2sort,junkNodes)
%
% Inputs:
%   TREE: PSMTree object
%   current: the node (element of TREE) that we're currently editting
%   domrng: the ranges of the dimensions of the new domain
%   nPList2sort and oPList2sort are from above. We pass these in to avoid
%   recomputing them will every recursion.
%   junkNodes: list of element of the structure array Tree that are no
%   longer needed.

n = size(domain,1);
newParamList = {domain.name}';

% Increase dimension of original node domain by giving inserted dimensions 
% bnds = [-inf inf].
nodeRangeOld = [repmat(-inf,n,1) repmat(inf,n,1)];
nodeRangeOld(nPList2sort,:) = Tree(currentNode).domainRange(oPList2sort,:);

% Intersect expanded original node subdomain with entire new domain.
newbnds = vertcat(domain.range);

lowerBnds = max(nodeRangeOld(:,1),newbnds(:,1));
upperBnds = min(nodeRangeOld(:,2),newbnds(:,2));

if any(lowerBnds > upperBnds)
  error(['The resize operation results in an empty domain for node: ' num2str(currentNode)])
end

% This is the new subdomain of currentNode.
nodeDomainNew = DClab.createDomainStructure(newParamList,[lowerBnds upperBnds]);

% Create the old subdomain of currentNode. This is used later on to change
% the domain of any surfaces defined on currentNode.
oldParamList(oPList2sort,1) = newParamList(nPList2sort);
nodeDomainOld = DClab.createDomainStructure(oldParamList,Tree(currentNode).domainRange);

% Determine if the node's cut (if exists) lies inside the new node domain.
% If not, we can remove one branch of the tree.
if ~isempty(Tree(currentNode).cutLocation)
  loc = Tree(currentNode).cutLocation;
  param = Tree(currentNode).cutParameter;
  
  temp = strmatch(param,char(newParamList),'exact');
  parent = Tree(currentNode).parentNode;
  leftChild = Tree(currentNode).childNodes(1);
  rightChild = Tree(currentNode).childNodes(2);
  
  if loc <= nodeDomainNew(temp).range(1)
    % Changing the domain resulted in the elimination of an existing cut
    % that lies to the left of the new domain. Thus we need to delete the
    % left child of currentNode and all its ancestors. Then we need to
    % replace currentNode with its (suitably modified) right child.
    
    %Find which indices of the tree structure we will wack away.
    junkNodes = [junkNodes; leftChild; Tree(leftChild).descendants];
    
    % Using the current indices (i.e., before we blow away unneeded element
    % of the structure array, replace the current node with its right child.

    % Step one, tell parent of new child

    if parent == -1
      % No parents to update.
    else
      parentsLeftChild = Tree(parent).childNodes(1);
      % Determine if the current node (which we're removing and replacing
      % with its right child) is a left or right child of its parent.
      if parentsLeftChild == currentNode
        Tree(parent).childNodes(1) = rightChild;
      else
        Tree(parent).childNodes(2) = rightChild;
      end      
    end 

    % Step two, tell right child (which will replace the current node) of
    % its new parent
    Tree(rightChild).parentNode = parent;
    
    % Step three, add the current node to the junk nodes.
    junkNodes = [junkNodes; currentNode];
        
    % Step four, add any surfaces that exist on currentNode and are
    % inherited by rightChild (note any surfaces inherited from currentNode
    % by ancestors of rightChild will also be inherited by rightChild
    % itself, so we don't need to traverse the tree to determine which
    % surfaces from currentNode are inherited). Then update indexing in
    % right child and its ancestors to reflect this.
    
    % Determine which surfaces in Tree(currentNode).DCSurface are inherited
    % by rightChild.
    tmp = vertcat(Tree(rightChild).responseModel2LocationOfDCSurface{:});
    oldIdx = find(tmp(:,1) == currentNode);
    
    if ~isempty(oldIdx)
      newIdx = size(Tree(rightChild).DCSurface,1) + (1:length(oldIdx))';
      
      %Traverse tree from rightChild on down and change [currentNode oldIdx(i)]
      %to [rightChild newIdx(i)] in any inherit references.
      Tree = updateSurfaceInheritances(Tree,rightChild,currentNode,rightChild,oldIdx,newIdx);
      
      newSurf = Tree(currentNode).DCSurface(oldIdx);
      
      % Make surfaces have the old domain of rightChild. The recursive call
      % of this function will modify the domain of all surfaces on
      % rightChild at once.
      rightDomainOld = DClab.createDomainStructure(oldParamList,Tree(rightChild).domainRange);
      newSurf = changeDomain(newSurf,nodeDomainOld,rightDomainOld);
      
      TMP = Tree(rightChild).DCSurface;
      Tree(rightChild).DCSurface = [TMP; newSurf];
    
      % TODO wierd matlab bug that this next line calls horzcat
      %Tree(rightChild).DCSurface = [Tree(rightChild).DCSurface; newSurf];
    end
    
    % Step five, recursive call to fix the domain of the right child
    [Tree junkNodes] = changeSubdomain(Tree,rightChild,domain,nPList2sort,oPList2sort,junkNodes);
    
  elseif loc >= nodeDomainNew(temp).range(2)
    % Changing the domain resulted in the elimination of an existing cut
    % that lies to the right of the new domain. Thus we need to delete the
    % right child of currentNode and all its ancestors. Then we need to
    % replace currentNode with its (suitably modified) left child.
    
    %Find which indices of the tree structure we will wack away.
    junkNodes = [junkNodes; rightChild; Tree(rightChild).descendants];
    
    % Using the current indices (i.e., before we blow away unneeded element
    % of the structure array, replace the current node with its left child.

    % Step one, tell parent of new child

    if parent == -1
      % No parents to update.
    else
      parentsLeftChild = Tree(parent).childNodes(1);
      % Determine if the current node (which we're removing and replacing
      % with its right child) is a left or right child of its parent.
      if parentsLeftChild == currentNode
        Tree(parent).childNodes(1) = leftChild;
      else
        Tree(parent).childNodes(2) = leftChild;
      end      
    end 

    % Step two, tell left child (which will replace the current node) of
    % its new parent
    Tree(leftChild).parentNode = parent;
    
    % Step three, add the current node to the junk nodes.
    junkNodes = [junkNodes; currentNode];
    
    % Step four, add any surfaces that exist on currentNode and are
    % inherited by leftChild (note any surfaces inherited from currentNode
    % by ancestors of leftChild will also be inherited by lefttChild
    % itself, so we don't need to traverse the tree to determine which
    % surfaces from currentNode are inherited). Then update indexing in
    % left child and its ancestors to reflect this.
    
    % Determine which surfaces in Tree(currentNode).DCSurface are inherited
    % by leftChild.
    tmp = vertcat(Tree(leftChild).responseModel2LocationOfDCSurface{:});
    oldIdx = find(tmp(:,1) == currentNode);
    
    if ~isempty(oldIdx)
      newIdx = size(Tree(leftChild).DCSurface,1) + (1:length(oldIdx))';
      
      %Traverse tree from leftChild on down and change [currentNode oldIdx(i)]
      %to [leftChild newIdx(i)] in any inherit references.
      Tree = updateSurfaceInheritances(Tree,leftChild,currentNode,leftChild,oldIdx,newIdx);
      
      newSurf = Tree(currentNode).DCSurface(oldIdx);
      
      % Make surfaces have the old domain of leftChild. The recursive call
      % of this function will modify the domain of all surfaces on
      % leftChild at once.
      leftDomainOld = DClab.createDomainStructure(oldParamList,Tree(leftChild).domainRange);
      newSurf = changeDomain(newSurf,nodeDomainOld,leftDomainOld);
      
      TMP = Tree(leftChild).DCSurface;
      Tree(leftChild).DCSurface = [TMP; newSurf];
    
      % TODO wierd matlab bug that this next line calls horzcat
      %Tree(leftChild).DCSurface = [Tree(leftChild).DCSurface; newSurf];
    end
    
    % Step five, recursive call to fix the domain of the right child
    [Tree junkNodes] = changeSubdomain(Tree,leftChild,domain,nPList2sort,oPList2sort,junkNodes);
    
  else
    % No children of currentNode need to be eliminated. Consequently we
    % just need to update the domainRange and surfaces of currentNode and
    % call this function recursively on both of its children.

    Tree(currentNode).DCSurface = changeDomain(Tree(currentNode).DCSurface,nodeDomainOld,nodeDomainNew);
    Tree(currentNode).domainRange = vertcat(nodeDomainNew.range);
    
    [Tree junkNodes] = changeSubdomain(Tree,leftChild,domain,nPList2sort,oPList2sort,junkNodes);
    [Tree junkNodes] = changeSubdomain(Tree,rightChild,domain,nPList2sort,oPList2sort,junkNodes);
  end
  
else
  % currrentNode has no children. Consequently we just need to update the
  % domainRange and surfaces of currentNode and exit recursion.
  
  Tree(currentNode).DCSurface = changeDomain(Tree(currentNode).DCSurface,nodeDomainOld,nodeDomainNew);
  Tree(currentNode).domainRange = vertcat(nodeDomainNew.range);  
end %if cut is not empty


function Tree = updateSurfaceInheritances(Tree,node,oldNode,newNode,oldIdx,newIdx)
% This function looks at all descendants of node. Anytime it finds
% [oldNode, element of oldIdx] in a responseModel2LocationOfDCSurface cell,
% it changes it to [newNode, corr. element of newIdx]
%
% oldIdx and newIdx must have the same dimensions

% Fix the node references in the .responseModel2LocationOfDCSurface
% array of node.
tmp = Tree(node).responseModel2LocationOfDCSurface;
rowsPerCell = 0.5*cellfun('prodofsize',tmp);

% Stack the contents into a matrix because we don't want to loop throw
% the cells.
tmp = vertcat(tmp{:});

idx1 = find(tmp(:,1) == oldNode);

% %TODO: can't we just vectorize this?
%for i1 = reshape(idx1,1,length(idx1))
%  tmp(i1,1) = newNode;
%  tmp(i1,2) = newIdx(tmp(i1,2)==oldIdx); %logical indexing
%end

tmp(idx1,1) = newNode;
tmp(idx1,2) = newIdx(tmp(idx1,2)==oldIdx); %logical indexing

% Recover the cell array with the first columns of the contents updated.
tmp = mat2cell(tmp,rowsPerCell,2);
Tree(node).responseModel2LocationOfDCSurface = tmp;

% If node is a parent, call this on its children
if ~isempty(Tree(node).childNodes)
  Tree = updateSurfaceInheritances(Tree,Tree(node).childNodes(1),oldNode,newNode,oldIdx,newIdx);
  Tree = updateSurfaceInheritances(Tree,Tree(node).childNodes(2),oldNode,newNode,oldIdx,newIdx);
end






