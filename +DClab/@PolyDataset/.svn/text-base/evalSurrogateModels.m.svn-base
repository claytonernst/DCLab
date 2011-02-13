function [y, yInt1, yInt2, node] = evalSurrogateModels(PD,x,pairs,paramList,node)
%EVALSURROGATEMODELS evaluates the surrogate models of a PolyDataset at x
%
%   Note: Supplying [] for any of the optional inputs will cause the
%   corresponding default input to be used.
%
%   Y = EVALSURROGATEMODELS(PD,X) evaluates the surrogate models of the
%   PolyDataset PD at the n-by-nEvals matrix of horizontally concatenated
%   parameter vectors X. The order of the parameters in X must correspond
%   to the order of the parameters in the PolyDataset. Y will be an m-by-1
%   cell array, where m is the number of response models in the PolyDataset.
%   Each cell will contain a p-by-nEvals matrix, where p is the number of
%   surrogate models that exist for a given response model. For example, if
%   the i'th response model has 3 surrogate models (each fit using
%   different input/output transformations) Y{i} will be 3-by-nEvals.
%
%   Y = EVALSURROGATEMODLES(PD,X,PAIRS) will output Y for only the subset
%   PAIRS of [1:PD.nPairs]. PAIRS should be a column array and the default
%   is PAIRS = [1:PD.nPairs]'.
%
%   Y = EVALSURROGATEMODELS(PD,X,PAIRS,PARAMLIST) allows you supply X in
%   order the the names appear in the column cell array of chars PARAMLIST.
%   Obviously you must include all parameters that are needed by the
%   models, but you have some flexibility with order and total number of
%   parameters in X. The default is PARAMLIST = PD.parameterList.
%
%   Y = EVALSURROGATEMODELS(PD,X,PAIRS,PARAMLIST,NODE) allows you to supple
%   the index or array of indices NODE that indicates which subdivision
%   each parameter vector (column of X) lives in. If NODE is scalar, they
%   all must live in a single subdivision. If a vector, its length must
%   equal size(X,2), i.e., one node for each parameter vector. If this
%   information is known, supply it using this input so that this
%   function does not need to compute this information.
%
%   [Y YINT1] = EVALSURROGATEMODLES(...) additionally returns an m-by-1 cell
%   of "intervals" that should (through consideration of the surrogate
%   modeling errors) enclose the output of the corresponding response model
%   M_i(X). Each cell will contain a 2*p-by-nEvals matrix, where p is the
%   number of surrogate models that exist for a given response model. For
%   example, if the i'th response model has 3 surrogate models (each fit
%   using different input/output transformations) Y{i} will be 6-by-nEvals.
%   The first two rows of this matrix supply the nEvals "intervals" formed
%   from the 1st surrogate model that should contain M_i(X). The 3rd and
%   4th rows of this matrix supply the nEvals "intervals" formed from the
%   2nd surrogate model that should contain M_i(X), etc.
%
%   [Y YINT1 YINT2] = EVALSURROGATEMODLES(...) returns an additional
%   m-by-1 cell of "intervals" that should (through consideration of the
%   surrogate modeling errors and any output uncertainty in the response
%   models) enclose the true value of corresponding response y_i(X). The
%   dimensions, etc., of YINT2 are identical to YINT1. If the response
%   models have no output uncertainty, YINT1 and YINT2 are identical.
%
%   [Y YINT1 YINT2 NODE] = EVALSURROGATEMODLES(...) additionally returns a
%   1-Nevals vector indicating the nodes at which each parameter
%   vector of X was located. This may be useful if you later evaluate the
%   surrogate model of a different ModelAndObservationPair at the same X
%   because you can pass this output as an input to this function in
%   subsequent calls, thereby avoiding the computational cost of
%   determining which subdivision of the domain of the surrogate models
%   contains contains each column of X.
%
%   See also PolyDataset, DCDataset/evalResponseModels, ResponseModel/eval


%TODO the 2nd output of getSurfaces should be the same, regardless of the
%node. In other words, all nodes of a
%PolyDataset/PiecewiseSurrogateModelTree have the same number of surrogates
%for each response model. If this is not true, the output of this function
%will need to be restructured. Is it true? I don't think it is. Should we
%just pad with zeros or last values to output matrices? We could have a
%'matrix' or 'cell' optional input.

%TODO eliminate the NODE input and output. It just complicates things to
%provide a very esoteric feature.


%error('see comments in code')
%TODO, this doesn't do what we want. There shouldn't be a cell array output
%and we should only get one point for each RM. isSurrogateFeasible relies
%on this.

ni = nargin;
switch ni
    case 2
        pairs = [];
        paramList = [];
        node = [];
    case 3
        paramList = [];
        node = [];
    case 4
        node = [];
    otherwise
        error(nargchk(2,5,ni))
end

if isempty(pairs)
    pairs = (1:PD.nPairs)';
elseif size(pairs,2) ~= 1 || max(pairs) > PD.nPairs
    error('Inputs: PAIRS must be a column array and its elements must reference MOPairs that are present in PD')
else
    %do nothing, input should be ok
end

%Initialize output
y = cell(length(pairs),1);
yInt1 = cell(length(pairs),1);
yInt2 = cell(length(pairs),1);

% nSurfacesPerPair returns an array of size(pairs). Each element contains
% the number of surfaces that exist for the corresponding pair.
mm = PD.PiecewiseSurrogateModelTree.nSurfacesPerResponseModel;

mm = mm(pairs);
N = size(x,2);
for i1 = 1:size(mm,1)
    y{i1} = zeros(mm(i1),N);
    yInt1{i1} = zeros(2*mm(i1),N);
    yInt2{i1} = zeros(2*mm(i1),N);
end

if isempty(paramList)
    %assign the paramList and assume x was supplied in the proper order
    paramList = PD.parameterList;
    if size(x,1)~=size(paramList,1)
        error('x must be n-by-Nevals, where n is the dimension of the domain of PD')
    end
else
    dsetPList = PD.parameterList;
    [trash idx1 idx2] = intersect(paramList,dsetPList); %#ok
    if length(dsetPList) ~= length(idx2)
        error('Supplied paramList does not contain all the dataset parameters')
    end
    if size(x,1)~=size(paramList,1)
        error('x must be n-by-Nevals, where n is the length of the supplied paramList')
    end

    %get x in the correct order for the dataset.
    tmp(idx2,:) = x(idx1,:);
    x = tmp;
    %overwrite the supplied paramList to correspond to the modified x
    paramList = dsetPList;

end

%If the nodes were provided, make sure they indeed contain the
%corresponding element of x.
TreeNodes = PD.PiecewiseSurrogateModelTree;
if ~isempty(node)
    if size(node,1) ~= 1
        error('Inputs: first dimension of NODE must be size 1')
    end
    if size(node,2) == 1
        xmin = min(x,[],2);
        xmax = max(x,[],2);
        rng = TreeNodes(node).domainRange;
        if ~DClab.issubset([xmin xmax],rng,1e-10,1e-10);
            error(['At least one of the supplied parameter vectors did not ' ...
                'live in the node you provided'])
        end
    elseif size(node,2) ~= N
        error(['When providing multiple nodes, you must give one for ' ...
            'each parameter vector (i.e., column of X)']);
    else

%         for i0 = 1:N
%             %Check that the nodes given are ok
%             rng = TreeNodes(node(i0).domainRange);
%             if ~DClab.issubset([x(:,i0) x(:,i0)],rng,1e-10,1e-10);
%                 error(['At least one of the nodes you provided did not match ' ...
%                     'the sudivision containing the corresponding parameter ' ...
%                     'vector'])
%             end
%         end
    end

else
    node = zeros(1,N);

    %First look for nodes that actually contain each parameter vector.
    %Having found none, allow a little bit of tolerance.
    for i1 = leafNodes(PD.PiecewiseSurrogateModelTree)
        rng = TreeNodes(i1).domainRange;
        for i0 = 1:N
            if DClab.issubset([x(:,i0) x(:,i0)],rng,0,0);
                node(i0) = i1;
            end
        end
    end
    i0 = find(node==0);
    if ~isempty(i0)
        for i1 = leafNodes(PD.PiecewiseSurrogateModelTree)
            rng = TreeNodes(i1).domainRange;
            for i2 = i0
                if DClab.issubset([x(:,i2) x(:,i2)],rng,1e-10,1e-10);
                    node(i2) = i1;
                end
            end
        end
    end

    if any(node == 0)
        error('I can not find your parameter vector in my domain')
    end
end

if size(node,2) == 1
    J = ones(1,N);
    uniqNode = node;
elseif size(node,2) ~= N
    error('You must supply as many nodes as parameter vectors')
else
    %Do an unique on nodes. We can save time by vectorizing---evaluating the
    %surrogate models of a given node at all columns of X that fall in this
    %node.
    [uniqNode I J] = unique(node); %#ok
end

for i0 = 1:length(uniqNode)

    %[surf surf2pairs] = getSurfaces(PD.PiecewiseSurrogateModelTree,uniqNode(i0),pairs);
    surf = surfaces(PD.PiecewiseSurrogateModelTree,uniqNode(i0));
    surf2pairs = PD.PiecewiseSurrogateModelTree.surface2CorrespondingResponseModelIndex;
    %surf2pairs should be the same size as surf, and contain members of pairs.

    % We did a unique on node to allow us to vectorize eval as much as
    % possible. Find which columns of x are located in uniqNode(i0).
    % Use logical indexing.
    idx1 = J == i0;

    tempx = x(:,idx1);
    try
        domain = struct('name',{PD.FreeParameter.name}','range',[]);
        for i1=1:PD.nParameters
            domain(i1).range = PD.FreeParameter(i1).range;
        end
        [ytmp bndstmp] = eval(surf,domain,tempx);
    catch
        disp('problem in PolyDataset/evalSurrogateModels')
        keyboard
    end

    %Assign output. This involves pulling off the rows of ytmp and bndstmp
    %that are relevant to a given response model, and placing these in the
    %proper location dictated by idx1. Since there are two rows in bndstmp
    %for each surface, we have this sort/cat business.
    for i1 = 1:length(pairs)
        idx2 = find(surf2pairs==pairs(i1));
        y{i1}(:,idx1) = ytmp(idx2,:);
        yInt1{i1}(:,idx1) = bndstmp(sort([idx2*2-1; idx2*2]),:);
    end
end

if nargout > 2


    % If X is a single column, (and pointwise if not)
    % yInt1(1,1) <= M(x) <= yInt1(2,1)
    %
    % The response models are such that
    % 1) M(x) + outputUnc(1) <= eta(x) <= M(x)+ outputUnc(2)
    % 2) M(x)(1 + outputUnc(1)) <= eta(x) <= M(x)(1 + outputUnc(2))
    % 3) M(x)*10^outputUnc(1) <= eta(x) <= M(x)*10^outputUnc(2)
    % 4) M(x)^(1+outputUnc(1)) <= eta(x) <= M(x)^(1+outputUnc(2))
    %
    % See pdf on uncertainties for further description

    for i1 = 1:length(pairs)
        RM = PD.ModelAndObservationPair(i1).ResponseModel;
        uncVect = RM.outputUncertaintyPlusMinus;
        uncCase = RM.uncertaintyCase;

        if uncCase == 1
            yInt2{i1} = yInt1{i1} + repmat(uncVect',[0.5 1].*size(yInt1{i1}));
        elseif uncCase == 2
            for i2 = 1:0.5*size(yInt1{i1},1)
                if uncVect(1) >= -1
                    yInt2{i1}(2*i2-1,:) = yInt1{i1}(2*i2-1,:)*(1 + uncVect(1));
                else
                    yInt2{1i}(2*i2-1,:) = yInt1{i1}(2*i2,:)*(1 + uncVect(1));
                end
                yInt2{i1}(2*i2,:) = yInt1{i1}(2*i2,:)*(1 + uncVect(2));
            end
        elseif uncCase == 3
            yInt2{i1} = yInt1{i1}.*repmat(10.^uncVect',[0.5 1].*size(yInt1{i1}));
        elseif uncCase == 4
            for i2 = 1:0.5*size(yInt1{i1},1)
                yInt2{i1}(2*i2-1,:) = max(yInt1{i1}(2*i2-1,:),0).^(1+uncVect(1));
                yInt2{i1}(2*i2,:) = yInt1{i1}(2*i2,:).^(1+uncVect(2));
            end
        else
            error('Internal inconsistency, condition should never occur')
        end
    end

end % if nargout > 2