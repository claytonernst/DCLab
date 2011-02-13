function [PSMTree1 PSMTree2] = mergePartitionStructures(PSMTree1,PSMTree2)
%MergePartitionStructures makes the partition structure of two objects the same
%
%   [PSMTREE1 PSMTREE2] = MergePartitionStructures(PSMTREE1,PSMTREE2)
%   requires the domains of PSMTREE1 and PSMTREE2 to be the same.  It
%   applies the subdivisions of PSMTREE2 to PSMTREE1 and XXX . From a
%   storage standpoint, given two objects, one with may surfaces and/or
%   partitions, and another that is less complicated, you should make the
%   complicated object PSMTREE1. (IS THIS TRUE, OR IS THERE NO DIFFERENCE?)


% We don't want to deal with empty inputs.
if isempty(PSMTree1) || isempty(PSMTree2)
  error('Inputs: both PSMTREE1 and PSMTREE2 must be nonempty')
end

if ~isequal(PSMTree1(1).parameterList,PSMTree2(1).parameterList) || ...
       ~isequal(PSMTree1(1).domainRange,PSMTree2(1).domainRange)
  error('Inputs: the domains of PSMTREE1 and PSMTREE2 must be identical, see CHANGEDOMAIN method')
end

% To create two trees with identical partition structures we carry out the
% following steps.
%
% 1) To form the new tree2, form newTree2 from the first node of tree2.
%    This represents a unpartitioned domain and since it is the root, none
%    of its surrogates can have been inherited.
%
% 2) With inherit==true, apply each cut that partitions the domain of tree1
%    to newTree2. At this point, newTree2 and tree1 will have indentical
%    partition structures.
%
% 3) With inherit==true, produce the completed newTree1 by applying each
%    cut that partitions the domain of tree2 to each leaf of tree1.
%
% 4) This is the tricky part. As in step three, we want to apply each cut
%    that partitions the domain of tree2 to each leaf of newTree2, however,
%    we don't want to inherit the surfaces, we want to use the surfaces
%    that are present in tree2. Conceptually we need to change the domain
%    of tree2 to that of the i^th leaf of newTree2, resulting in the binary
%    tree, called say newSubtree, whose root has the same domain as the
%    i^th leaf of newTree2. Then we want to replace the i^th leaf with
%    newSubtree. This requires several substeps:
%
% 4a) Change the domain of tree2 to that of the i^th leaf of newTree2,
%     resulting in the binary tree, called say newSubtree, whose root has
%     the same domain as the i^th leaf of newTree2. Note whether the root
%     of tree2 becomes the root of newSubtree, or if changing the domain
%     lops off the root of tree2.
%
% 4b) In effect we want to replace the i^th leaf of newTree2 with the root
%     node of newSubtree. If the root of tree2 wasn't lopped off in 4a, we
%     can leave i^th leaf alone, since it currently inherits all its
%     surfaces from the root of newTree2, and these are the same surfaces
%     that are present (suitably scaled) in the root node of newSubtree. If
%     the root of tree2 was lopped off in 4a, we have to do some work. We
%     need to find which surfaces defined on the root of newSubtree were
%     not inherited from the root of tree2. These surfaces are new
%     information that is not currently in the i^th leaf, and consequently
%     these surfaces need to be placed in the i^th leaf of newTree2.
%
% 4c) Make the descendents of the root node of newSubtree new object be
%     descendents of the i^th leaf.

% Step 1:
newTree2 = PSMTree2(1);
newTree2.childNodes = []; %Remove all child info to make this a fresh unpartitioned object
newTree2.cutParameter = '';
newTree2.cutLocation = [];

% Step 2: Apply the cuts of tree1 to newTree2.

% TODO: There is surely a more efficient way to do this. However, we must apply
% the cuts of tree1 to each leaf of newTree2 in the same order that they
% were originally applied to tree1. I don't see any quick recursive
% solution to this. The bottleneck is that change domain updates the
% surfaces that are defined on it. Perhaps we could just shortcut this.
% This function won't be used much, so I really don't care.
for i1 = leafNodes(newTree2)'
  domrng = newTree2(i1).domainRange;
  dom = DClab.createDomainStructure(newTree2(1).parameterList,domrng);

  newSubtree = changeDomain(PSMTree1,dom);
  NnewNodes = nNodes(newSubtree) - 1;
  
  cutOrder = findOrderNodesWereCut(newSubtree,nNodes(newSubtree),[]);
  
  %The first column of this vector contains the node numbers in newSubtree that
  %become the second column of node numbers in newTree2. 
  idxConvert = [1 i1; 1+(1:NnewNodes)' nNodes(newTree2)+(1:NnewNodes)'];
 
  for i2 = cutOrder'
    % Divide cuttee in dimension param at location value.
    param = newSubtree(i2).cutParameter;
    value = newSubtree(i2).cutLocation;

    node2cut = idxConvert(idxConvert(:,1)==i2,2);
    newTree2 = subdivideDomain(newTree2,[],node2cut,param,value,1);
  end
end

% Step 3: Apply all cuts of tree2 to each leaf of tree1
tree1Leaves = leafNodes(PSMTree1);
for i1 = tree1Leaves'
  domrng = PSMTree1(i1).domainRange;
  dom = DClab.createDomainStructure(PSMTree1(1).parameterList,domrng);

  newSubtree = changeDomain(PSMTree2,dom);
  NnewNodes = nNodes(newSubtree) - 1;
  cutOrder = findOrderNodesWereCut(newSubtree,nNodes(newSubtree),[]);
  
  %The first column of this vector contains the node numbers in newSubtree that
  %become the second column of node numbers in PSMTree1. 
  idxConvert = [1 i1; 1+(1:NnewNodes)' nNodes(PSMTree1)+(1:NnewNodes)'];
 
  for i2 = cutOrder'
    % Divide cuttee in dimension param at location value.
    param = newSubtree(i2).cutParameter;
    value = newSubtree(i2).cutLocation;

    node2cut = idxConvert(idxConvert(:,1)==i2,2);
    PSMTree1 = subdivideDomain(PSMTree1,[],node2cut,param,value,1);
  end
end  
% PSMTree1 is now complete.

% Perform step 4 for each leaf of newTree2. We want to apply all cuts of
% tree2 to each leaf of newTree2, but we want to used the original
% noninherited surfaces of tree2 surfaces when possible.
tmp2Leaves = leafNodes(newTree2);
oldLengthTree2 = nNodes(newTree2);

for i1 = tmp2Leaves'

  domrng = newTree2(i1).domainRange;
  dom = DClab.createDomainStructure(newTree2(1).parameterList,domrng);
  
  %Step 4a: Determine which node of tree2 will become the new root after we
  %change the domain. 
  newRoot = findYoungestDescendantContaining(PSMTree2,domrng,1);
  newSubtree = changeDomain(PSMTree2,dom);
  
  % Step 4b: Use data from tree2(newRoot) to populate newTree2(i1).
  
  %Leave parameterList unchanged
  %Leave domainRange unchanged
  %Leave parentNode unchanged
  if ~isempty(PSMTree2(newRoot).childNodes)
    newTree2(i1).childNodes = oldLengthTree2 + [1 2];
  end
  newTree2(i1).cutParameter = PSMTree2(newRoot).cutParameter;
  newTree2(i1).cutLocation = PSMTree2(newRoot).cutLocation;
  
  if newRoot > 1 %Indicates root node was lopped off
  
    %If a surface was inherited from the root, we don't want to place it
    %in newTree2(i1). However, if it was created for anything decending
    %from root, we do.
    
    % Since newRoot ~= 1, we want these surfaces to be placed in
    % newTree2(i1).
    Surf = changeDomain(PSMTree2(newRoot).DCSurface,DClab.createDomainStructure(PSMTree2(1).parameterList,PSMTree2(newRoot).domainRange),dom);

    %Change references that point to newRoot (or its ancesters except for
    %root) to i1
    tmp = PSMTree2(newRoot).responseModel2LocationOfDCSurface;
    rowsPerCell = 0.5*cellfun('prodofsize',tmp);

    % Stack the contents into a matrix because we don't want to loop throw
    % the cells.
    tmp = vertcat(tmp{:});
    
    %Change references that point to newRoot to i1
    tmp((tmp(:,1) == newRoot),1) = i1;
    
    %Grab surfaces from all ancestors of newRoot except for root. 
    parent = tree2(newRoot).parent;
    while parent > 1
      currNSurf = nSurfaces(Surf);

      newSurf = changeDomain(PSMTree2(parent).DCSurface,DClab.createDomainStructure(PSMTree2(1).parameterList,PSMTree2(parent).domainRange),dom);

      % If a surface was inherited from parent, we want to place it in
      % newTree2(i1) 
      idx = find(tmp(:,1) == parent);

      tmp(idx,1) = i1;
      tmp(idx,2) = tmp(idx,2) + currNSurf;

      Surf = [Surf; newSurf];
      parent = PSMTree2(parent).parent;
    end

    % Recover the cell array with the contents updated.
    RM2LocOfSurf = mat2cell(tmp,rowsPerCell,2);

    newTree2(i1).DCSurface = Surf;
    newTree2(i1).responseModel2LocationOfDCSurface = RM2LocOfSurf;
  else
    %Do nothing. All surfaces of newTree(i1) will be inherited from
    %newTree2(1), which is what we want.
  end
  newTree2(i1).fittingSettings = PSMTree2(newRoot).fittingSettings;

  % Step 4c: Now add the children of tree2(newRoot) (i.e., newSubtree(1))
  % to newTree2;
  NnewNodes = nNodes(newSubtree) - 1;
  newTree2 = builtin('vertcat',newTree2, newSubtree(2:end));
  newLengthTree2 = builtin('length',newTree2); %newTree2 is not a valid object here. It shouldn't matter, but use builtin anyway.
    
  %The first column of this vector contains the node numbers in newSubtree that
  %become the second column of node numbers in newTree2. 
  idxConvert = [1 i1; 1+(1:NnewNodes)' nNodes(newTree2)+(1:NnewNodes)'];
  for i2 = oldLengthTree2+1:newLengthTree2
    %Fix all node references due to vertcat.
    newTree2(i2).parentNode = idxConvert(newTree2(i2).parentNode==idxConvert(:,1),2);   
    if ~isempty(newTree2(i2).childNodes)
      newTree2(i2).childNodes(1) = idxConvert(newTree2(i2).childNodes(1)==idxConvert(:,1),2);    
      newTree2(i2).childNodes(2) = idxConvert(newTree2(i2).childNodes(2)==idxConvert(:,1),2); 
    end
    
    %Change references that point to newRoot to i1
    tmp = newTree2(i2).responseModel2LocationOfDCSurface;
    rowsPerCell = 0.5*cellfun('prodofsize',tmp);

    % Stack the contents into a matrix because we don't want to loop throw
    % the cells.
    tmp = vertcat(tmp{:});
    uniq = unique(tmp(:,1));
    for i3 = uniq'
      tmp(tmp(:,1)==i3,1) = idxConvert(idxConvert(:,1)==i3,2);
    end

    % Recover the cell array with the first columns of the contents updated.
    newTree2(i2).responseModel2LocationOfDCSurface = mat2cell(tmp,rowsPerCell,2);

  end
end

PSMTree2 = newTree2;

function nodes = findOrderNodesWereCut(tree,curNode,nodes)
% This function works recursively. Initially call it with curNode =
% length(tree). The function works backwards from the end of the structure
% array to determine the order the cuts were make. The first cut was made
% at nodes(1), the 2nd at nodes(2), etc.

if tree(curNode).parentNode == -1
  %we're done, just return the input
else
  % find parent of the current node. Call recursively with the previous
  % node that does not share this parent
  par = tree(curNode).parentNode;
  if tree(curNode-1).parentNode == par
    nodes = [findOrderNodesWereCut(tree,curNode-2,nodes); par];
  else
    nodes = [findOrderNodesWereCut(tree,curNode-1,nodes); par];
  end
end

function tree1 = doCuts(tree1,tree2,cutteeIdx,cutterIdx)

%This function needs to apply the entire cut structure of tree2 (in the
%original order) to a given leaf of tree1. 

%How to determine the order of the cuts?
% Work backwards. Start at end.


%
%  The function applies the partition of the node tree2(cutterIdx) (which
%  may consist of several cuts) to the *leaf* node tree1(cutteeNode).

%do it recursively for fun
%
% check if cutterIdx of c2 has a cut.
%   no ==> return c1 unchanged
%
%  yes ==> make cut on cutteeIdx of c1
%
%        does the cut pass beneath cutteeIdx, pass above
%        cutteeIdx, or actually divide cutteeIdx
%          pass below ==> apply cut of cutterIdx.greater to cuttee
% 
%          pass above ==> apply cut of cutterIdx.lesser to cuttee
%
%          splits ==> make the cut on cuttee and apply the cut of 
%            cutterIdx.greater to greater half of cuttee and apply cut 
%            of cutterIdx.lesser to lesser half of cuttee

if isempty(tree2(cutterIdx).cutLocation)
  %do nothing
else
  % Divide cuttee in dimension param at location value.
  param = tree2(cutterIdx).cutParameter;
  value = tree2(cutterIdx).cutLocation;
  dim = strmatch(param,tree1(1).parameterList,'exact');

  % See if this cut divides cuttee.
  %
  % If it falls to the right side of cuttee (i.e., cuts empty
  % space), it does nothing, and we want to apply cuts
  % of cutter.lesser to current node.
  %
  % If it falls to the left side of cuttee (i.e., cuts empty
  % space), it does nothing, and we want to apply cuts
  % of cutter.lesser to current node.%
  %
  % If it actually divides cuttee, we want to perform the cut and then
  % apply cuts of cutter.lesser to cuttee.lesser and apply the cuts of
  % cutter.greater to cuttee.greater.
  
  temp = tree1(cutteeIdx).domainRange(dim,:);
  
  if value <= temp(1)
    % Fell to the left of cuttee, apply greater side of cutter tree to
    % current node.
    tree1 = doCuts(tree1,tree2,cutteeIdx, ...
                tree2(cutterIdx).childNodes(2));
  elseif temp(2) <= value
    % Fell to the right of cuttee, apply lesser side of cutter tree to
    % current node.
    tree1 = doCuts(tree1,tree2,cutteeIdx, ...
                tree2(cutterIdx).childNodes(1));
  else
    % The cut divides cuttee, so now lets cut each subcube
    % split is called with inherit = true so surfaces of c1 will
    % propagate down. 
    tree1 = subdivideDomain(tree1,[],cutteeIdx,param,value,1);
    
    %continue recursively
    tree1 = doCuts(tree1,tree2,tree1(cutteeIdx).childNodes(1),tree2(cutterIdx).childNodes(1));
    tree1 = doCuts(tree1,tree2,tree1(cutteeIdx).childNodes(2),tree2(cutterIdx).childNodes(2));
  end
  
end
 
function young = findYoungestDescendantContaining(tree,dom,node)
% Recursive function to determine the youngest descendant of NODE whose
% domain is a superset of the n-by-2 matrix DOM. It is assumed (but not
% verified) that the domain of NODE is a superset of DOM. This descendant
% would serve as the root node of a subtree and all cuts on that subtree
% are guarenteed to partition dom.

if isempty(tree(node).cutLocation)
  young = node;
else
  param = tree(node).cutParameter;
  value = tree(node).cutLocation;
  dim = strmatch(param,tree(1).parameterList,'exact');
  if value < dom(dim,1)
    young = findNewRoot(tree,dom,tree(node).childNodes(2));
  elseif dom(dim,2) < value
    young = findNewRoot(tree,dom,tree(node).childNodes(1));
  else
    young = node;
  end
end
  





