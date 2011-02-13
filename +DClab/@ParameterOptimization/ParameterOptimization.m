classdef ParameterOptimization < DClab.DCObject
    % PARAMETEROPTIMIZATION Constructor
    %
    % obj = ParameterOptimization(Dset)
    % obj = ParameterOptimization(Dset,dcopt)
    % obj = ParameterOptimization(Dset,dcopt,norm)
    % obj = ParameterOptimization(Dset,dcopt,norm,trgWeights)
    % obj = ParameterOptimization(oldobj)
    % obj = ParameterOptimization(oldobj,dcopt)
    %
    % This function solves the optimization problem
    %
    %      min_x  ||w_i(M_i(x)-d_i)||,
    %
    % where w_i are components of the trgWeights vector, and default to
    % w_i = inv(u_i). If the uncertainties bounds are not symmetric
    % (+/-) we will shift d such that they are.
    %
    % obj = ParameterOptimization(Dset) locally optimizes the above objective function
    % with the two norm.  All options will be the defaults produced by calling
    % the DCOptions constructor with no inputs.
    %
    % obj = ParameterOptimization(Dset,dcopt) locally optimizes the above objective
    % function with the two norm using the options supplied in the DCOptions
    % object dcopt.
    %
    % obj = ParameterOptimization(Dset,dcopt,norm) locally optimizes the above objective
    % function with the given norm using the options supplied in the DCOptions
    % object dcopt. The input norm can be 'one', 'two', or 'inf'.
    %
    % obj = ParameterOptimization(Dset,dcopt,norm,trgWeights) locally optimizes the
    % above objective function with the given norm using the options supplied
    % in the DCOptions object dcopt. The input norm can be 'one', 'two', or
    % 'inf'. trgWeights is a nUnits(Dset) x 1 vector specifying the wieghts to
    % use in the objective function.
    %
    % Inputs
    %   Dset: A Dataset or PolyDataset object
    %   dcopt[optionsal]: A DCOptions object.  If not suppled, default options
    %     obtained from the DCOptions constructor are used, with the
    %     exception that maxBranchBoundIter is set to 1. This function
    %     currently does not branch and bound, so a warning is issued
    %     if maxBranchBoundIter > 1. Additionally this only uses the
    %     linXlinY transformations, so they'd better be available.
    %   norm[optional, default=two]: the norm to use in the objective function
    %     penalizing model/data mismatch. Either 'one', 'two', or 'inf'.
    %   trgWeights[default = 1/u]: A vector of weights to use in the
    %     objective function.
    %
    % Outputs
    %   obj: a ParameterOptimization object with fields .PDset, .norm, .weights,
    %     .optimInfo

    properties
        iter;
        norm = '';
        weights;
        runDate = '';
    end
    
    properties (Dependent)
        costUB
        costLB
        bestx
        wresid
        nIter
        DsetName
        Dset
        optimOpts
        surfFittingOpts
        PDset
        iterStruct
    end

    methods

        function obj = ParameterOptimization(varargin)


            ni = nargin;
            no = nargout;
            error(nargoutchk(0,1,no));

            if ni == 0
                return
            end
            if isa(varargin{1},'DClab.ParameterOptimization')
                if isempty(varargin{1})
                    warning('DClab:invalidInputs','First input to ParameterOptimization is empty, returning an empty object')
                    return
                end
                if ni == 1
                    varargin{2} = DClab.DCOptions;
                elseif ni == 2
                    if isempty(varargin{2}) && isnumeric(varargin{2})
                        varargin{2} = DClab.DCOptions;
                    end
                    assert(isa(varargin{2},'DClab.DCOptions'),'Inputs: 2nd input to ParameterOptimization of improper class.');
                else
                    error('Inputs: when 1st is a ParameterOptimization object, at most two inputs are allowed.');
                end
            elseif isa(varargin{1},'DClab.DCDataset') || isa(varargin{1},'DClab.PolyDataset')
                nU = varargin{1}.nPairs;
                if nU == 0
                    warning('DClab:invalidInputs','First input to ParameterOptimization contains no DatasetUnits, returning an empty object')
                    return
                end
                if ni==1
                    varargin{2} = DClab.DCOptions;
                    varargin{3} = 'two';
                    trgWeights = zeros(nU,1);
                    for i1 = 1:nU
                        u = varargin{1}.ModelAndObservationPair(i1).observationUncertaintyPlusMinus;
                        trgWeights(i1) = 2./diff(u);
                    end
                    varargin{4} = trgWeights;
                elseif ni==2
                    %Check the 2nd input
                    if isempty(varargin{2}) && isnumeric(varargin{2})
                        varargin{2} = DClab.DCOptions;
                    end
                    assert(isa(varargin{2},'DClab.DCOptions'),'Inputs: 2nd input to ParameterOptimization must be a DCOptions object.');
                    varargin{3} = 'two';
                    trgWeights = zeros(nU,1);
                    for i1 = 1:nU
                        u = varargin{1}.ModelAndObservationPair(i1).observationUncertainty;
                        trgWeights(i1) = 2./diff(u);
                    end
                    varargin{4} = trgWeights;
                elseif ni==3
                    %Check the 2nd input
                    if isempty(varargin{2}) && isnumeric(varargin{2})
                        varargin{2} = DClab.DCOptions;
                    end
                    assert(isa(varargin{2},'DClab.DCOptions'),'Inputs: 2nd input to ParameterOptimization must be a DCOptions object.');
                    %Check the 3rd input
                    if isempty(varargin{3}) && isnumeric(varargin{3})
                        varargin{3} = 'two';
                    end
                    assert(ischar && ismember(varargin{3},{'one','two','inf'}),'Inputs: 3rd input to ParameterOptimization must be the name of a norm.');
                    trgWeights = zeros(nU,1);
                    for i1 = 1:nU
                        u = varargin{1}.ModelAndObservationPair(i1).observationUncertainy;
                        trgWeights(i1) = 2./diff(u);
                    end
                    varargin{4} = trgWeights;
                elseif ni==4
                    %Check the 2nd input
                    if isempty(varargin{2}) && isnumeric(varargin{2})
                        varargin{2} = DClab.DCOptions;
                    end
                    assert(isa(varargin{2},'DClab.DCOptions'),'Inputs: 2nd input to ParameterOptimization must be a DCOptions object.');
                    %Check the 3rd input
                    if isempty(varargin{3}) && isnumeric(varargin{3})
                        varargin{3} = 'two';
                    end
                    assert(ischar(varargin{3}) && ismember(varargin{3},{'one','two','inf'}),'Inputs: 3rd input to ParameterOptimization must be the name of a norm.');
                    %Check the 4th input
                    if isempty(varargin{4})
                        trgWeights = zeros(nU,1);
                        for i1 = 1:nU
                            u = varargin{1}.ModelAndObservationPair(i1).observationUncertainy;
                            trgWeights(i1) = 2./diff(u);
                        end
                        varargin{4} = trgWeights;
                    else
                        assert(isnumeric(varargin{4}) && isequal(size(varargin{4}),[nU,1]),'Inputs: 4th input (trgWeights) must be a column vector of length = nUnits(Dset)')
                    end
                else
                    error(nargchk(0,4,ni));
                end
            else
                error('Inputs: 1st input to ParameterOptimization of improprer class');
            end

            %Now all inputs are legal. If the first is a ParameterOptimization object, the 2nd is
            %a DCOptions object. Otherwise the first is a DCDataset or
            %PolyDataset
            %object, the 2nd is a DCOptions object, the 3rd is the norm to use, and the
            %4th is the trgWeights to pass to the PolyDataset methods (possibly empty).

            %Initialization. Define opt, guihand, dbSaveName and iter.
            %At this point length(varargin)==2 if varargin{1} is a ParameterOptimization object,
            %and 4 otherwise.
            opt = varargin{2};
            assert(~opt.omitInnerBound,'Inputs: to use ParameterOptimization, the inner bound must be computed')
            if strcmp('linXlinY',opt.surfaceTransformation)
                opt.surfaceTransformation = 'linXlinY';
            end

            guihand = opt.guiHandle;
            if ~isempty(opt.fileName2Save)
                dbSaveName = opt.fileName2Save;
            else
                dbSaveName='';
            end
            if isa(varargin{1},'DClab.DCDataset') || isa(varargin{1},'DClab.PolyDataset')
                %We're starting from scratch.
                if isa(varargin{1},'DClab.PolyDataset')
                    PDset = varargin{1};
                else
                    if ismember(opt.display,{'iter';'all';'ALL'})
                        str = '=====Creating initial PolyDataset=====';
                        DClab.dcdispstr(str,guihand,false)
                    end
                    PDset = DClab.PolyDataset(varargin{1},opt);
                end
                norm = varargin{3};
                trgWeights = varargin{4};

                i1 = 1;
                % Initialize iter in the order we chose
                [trash junk iter.optimOpts] = decompose(opt);
                iter.PDset = PDset;
                iter.worstLBleaf = [];
                iter.dimSplit = [];
                iter.splitLoc = [];
                iter.leafInfo.leafIndices = PDset.leafNodes;


                % iter.norm = [];
                % iter.bestx = [];
                % iter.wresid = [];
                % iter.weights = [];
                % iter.cost = [];
                % iter.method = [];
                % iter.PDset = [];
                % iter.opt = opt;

                if ~isempty(dbSaveName)
                    save(dbSaveName,'iter')
                end

                if ~opt.omitOuterBound
                    if ismember(opt.display,{'iter';'all';'ALL'})
                        str = ['  Computing lower bound on the ' norm '-norm objective'];
                        DClab.dcdispstr(str,opt.guiHandle,false)
                    end
                    lb = olbFitX(iter(i1).PDset,opt,norm,trgWeights);
                    iter(i1,1).leafInfo.lower.bnd = lb.star.bnd;
                    [LB LBIdx] = min(iter(i1).leafInfo.lower.bnd);
                    iter(i1).worstLBleaf = iter(i1).leafInfo.leafIndices(LBIdx);
                else
                    % Assign default outputs to uncomputed quantities
                    LB = nan;
                end

                %The inner bound is always computed
                if ismember(opt.display,{'iter';'all';'ALL'})
                    str = ['  Optimizing model parameters with a ' norm '-norm objective'];
                    DClab.dcdispstr(str,opt.guiHandle,false)
                end
                ub = ilbFitX(iter(i1).PDset,opt,norm,trgWeights);

                iter(i1,1).leafInfo.upper.bnd = ub.pound.bnd;
                iter(i1,1).leafInfo.upper.bestx = ub.pound.xfeas;
                iter(i1,1).leafInfo.upper.wresid = ub.pound.resid;
                iter(i1,1).leafInfo.upper.method = ub.pound.method;
                UB = min(iter(i1,1).leafInfo.upper.bnd);

                if ~isempty(dbSaveName)
                    save(dbSaveName,'iter')
                end

                if opt.omitInnerBound || opt.omitOuterBound
                    if opt.maxBranchBoundIter > 1
                        if ismember(opt.display,{'final';'notify';'iter';'all';'ALL'})
                            str = 'Exiting with 1 iteration: both inner and outer bounds required to branch and bound';
                            DClab.dcdispstr(str,guihand,false)
                        end
                    end
                    if ismember(opt.display,{'final';'iter';'all';'ALL'})
                        str = ['UB: ' num2str(UB)];
                        DClab.dcdispstr(str,guihand,false)
                        str = ['LB: ' num2str(LB)];
                        DClab.dcdispstr(str,guihand,false)
                    end

                    if UB < 100*opt.tolFun && ismember(opt.display,{'notify';'iter';'all';'ALL'})
                        str = '  Warning: Tolerances on objective may be too loose for convergence, consider decreasing tolFun';
                        DClab.dcdispstr(str,opt.guiHandle,false)
                    end

                    obj.iter = iter;
                    obj.norm = norm;
                    obj.weights = trgWeights;
                    obj.runDate = date;
                    return
                end

            else
                %We recieved a ParameterOptimization object. Make sure we have something to do,
                %then proceed into the while loop.
                iter = varargin{1}.iterStruct;
                norm = varargin{1}.norm;
                trgWeights = varargin{1}.weights;
                i1 = length(iter);

                if opt.omitInnerBound || opt.omitOuterBound
                    str = 'Exiting without performing an additional iteration: both inner and outer bounds required to branch and bound';
                    DClab.dcdispstr(str,guihand,false)
                    obj = varargin{1};
                    return
                end
                if opt.maxBranchBoundIter <= length(iter)
                    if ismember(opt.display,{'final';'notify';'iter';'all';'ALL'})
                        str = 'Exiting without performing an additional iteration: maxBranchBoundIter already met in the previous iteration';
                        DClab.dcdispstr(str,guihand,false)
                    end
                    obj = varargin{1};
                    return
                end

                LB = varargin{1}.costLB;
                UB = varargin{1}.costUB;
                if isnan(LB) || isnan(UB)
                    str = 'Exiting without performing an additional iteration: both inner and outer bounds were not computed in the previous iteration';
                    DClab.dcdispstr(str,guihand,false)
                    obj = varargin{1};
                    return
                end
                if UB-LB <= opt.bbTermTol
                    str = 'Exiting without performing an additional iteration: bbTermTol already met in the previous iteration';
                    DClab.dcdispstr(str,guihand,false)
                    obj = varargin{1};
                    return
                end

                %If we made it this far, we're ready to perform another iteration.
                if ismember(opt.display,{'iter';'all';'ALL'})
                    str = 'First displayed bounds are for the old ParameterOptimization object';
                    DClab.dcdispstr(str,guihand,false)
                end

            end

            % Conditional to enter the branch and bound iteration
            while UB-LB > opt.branchBoundTermTol && i1 < opt.maxBranchBoundIter

                %display for fun
                if ismember(opt.display,{'iter';'all';'ALL'})
                    str = ['UB: ' num2str(UB)];
                    DClab.dcdispstr(str,guihand,false)
                    str = ['LB: ' num2str(LB)];
                    DClab.dcdispstr(str,guihand,false)
                end

                % Increament the iteration count
                i1 = i1+1;
                [dispOpts fitOpts iter(i1,1).optimOpts] = decompose(opt);

                % Decide where to divide prev worstLBleaf
                [dimSplit, splitLoc, inherit] = chooseDivision(iter(i1-1).PDset,iter(i1-1).worstLBleaf);
                iter(i1).dimSplit = dimSplit;
                iter(i1).splitLoc = splitLoc;

                if any(strmatch(opt.display,{'iter';'all';'ALL'}))
                    str = ['===Iteration ' num2str(i1) ': spliting dimension ' dimSplit ' of leaf ' num2str(iter(i1-1).worstLBleaf) ...
                        ' at location ' num2str(splitLoc) '==='];
                    DClab.dcdispstr(str,guihand,false)
                end

                % Subdivide the domain of PDset where indicated.
                % First provide the current display options.
                oldPDset = iter(i1-1).PDset;
                oldPDset.displaySettings=dispOpts;
                iter(i1).PDset = split(oldPDset,iter(i1-1).worstLBleaf,inherit,dimSplit,splitLoc,fitOpts);
                iter(i1).leafInfo.leafIndices = iter(i1).PDset.leafNodes;

                % If requested, save iter structure
                if ~isempty(dbSaveName)
                    save(dbSaveName,'iter')
                end

                % Find which leaves were added with that division
                newLeaves = setdiff(iter(i1).leafInfo.leafIndices,iter(i1-1).leafInfo.leafIndices);

                %we only need do optimization on the newly added leaves, then we
                %will zero out the bnd, mults, logX, from the parent of these
                %new leaves and add these quantities from the to new
                %leaves. that way we alway have available bnd, mults, and logX
                %for each leaf.
                try
                    if ismember(opt.display,{'iter';'all';'ALL'})
                        str = ['  Computing lower bound on the ' norm '-norm objective'];
                        DClab.dcdispstr(str,guihand,false)
                    end
                    lb = olbFitX(iter(i1).PDset,opt,norm,trgWeights,newLeaves);

                    if ismember(opt.display,{'iter';'all';'ALL'})
                        str = ['  Optimizing model parameters with a ' norm '-norm objective'];
                        DClab.dcdispstr(str,opt.guiHandle,false)
                    end
                    ub = ilbFitX(iter(i1).PDset,opt,norm,trgWeights,newLeaves);

                catch
                    str = 'failed call to olbFitX or ilbFitX in ParameterOptimization';
                    DClab.dcdispstr(str,guihand,false)
                    DClab.dcdispstr(lasterror,guihand,true);
                    keyboard
                end

                %do the zeroing out an updating alluded to before. three elements
                %of the iter stuct need to be updated/modified: iter.leafInfo.upper,
                %iter.leafInfo.lower1, and iter.leafInfo.lower2.
                parentIdx = find(iter(i1-1).leafInfo.leafIndices == iter(i1-1).worstLBleaf);

                % Copy previous values
                iter(i1,1).leafInfo = iter(i1-1).leafInfo;
                % Correct leaf Indices
                iter(i1,1).leafInfo.leafIndices = iter(i1).PDset.leafNodes;

                %Eliminate info relevant to the leaf that was split and append info from the two daughter leaves
                iter(i1).leafInfo.lower.bnd(parentIdx) = [];
                iter(i1).leafInfo.lower.bnd = [iter(i1).leafInfo.lower.bnd; lb.star.bnd];

                iter(i1).leafInfo.upper.bnd(parentIdx) = [];
                iter(i1).leafInfo.upper.bestx(parentIdx) = [];
                iter(i1).leafInfo.upper.wresid(parentIdx) = [];
                iter(i1).leafInfo.upper.method(parentIdx) = [];

                iter(i1).leafInfo.upper.bnd = [iter(i1).leafInfo.upper.bnd; ub.pound.bnd];
                iter(i1).leafInfo.upper.bestx = [iter(i1).leafInfo.upper.bestx; ub.pound.xfeas];
                iter(i1).leafInfo.upper.wresid = [iter(i1).leafInfo.upper.wresid; ub.pound.resid];
                iter(i1).leafInfo.upper.method = [iter(i1).leafInfo.upper.method; ub.pound.method];

                [LB LBIdx] = min(iter(i1).leafInfo.lower.bnd);
                iter(i1).worstLBleaf = iter(i1).leafInfo.leafIndices(LBIdx);
                UB = min(iter(i1).leafInfo.upper.bnd);

                if ~isempty(dbSaveName)
                    save(dbSaveName,'iter')
                end

            end

            if UB < 100*opt.tolFun && ismember(opt.display,{'notify';'iter';'all';'ALL'})
                str = '  Warning: Tolerances on objective may be too loose for convergence, consider decreasing tolFun';
                DClab.dcdispstr(str,opt.guiHandle,false)
            end

            if UB-LB > opt.tolFun && any(strmatch(opt.display,{'final';'notify';'iter';'all';'ALL'}))
                str = 'Exiting: maximum number of branch and bound iterations reached';
                DClab.dcdispstr(str,guihand,false)
            end

            if any(strmatch(opt.display,{'final';'iter';'all';'ALL'}))
                str = ['UB: ' num2str(UB)];
                DClab.dcdispstr(str,guihand,false)
                str = ['LB: ' num2str(LB)];
                DClab.dcdispstr(str,guihand,false)
                str = '=====Optimization complete=====';
                DClab.dcdispstr(str,guihand,false)
            end

            %assign output
            obj.iter = iter;
            obj.norm = norm;
            obj.weights = trgWeights;
            obj.runDate = date;
        end
        
        function out = get.costUB(obj)
            out = min(obj.iter(end).leafInfo.upper.bnd);
        end
        function out = get.costLB(obj)
            opt = obj.iter(end).optimOpts;
            if opt.omitOuterBound
                out = [];
            else
                out = min(obj.iter(end).leafInfo.lower.bnd);
            end
        end
        function out = get.bestx(obj)
            [UB UBidx] = min(obj.iter(end).leafInfo.upper.bnd);
            out = obj.iter(end).leafInfo.upper.bestx{UBidx};
        end
        function out = get.wresid(obj)
            [UB UBidx] = min(obj.iter(end).leafInfo.upper.bnd);
            out = obj.iter(end).leafInfo.upper.wresid{UBidx};
        end
        function out = get.nIter(obj)
            out = length(obj.iter);
        end
        function out = get.DsetName(obj)
            out = obj.iter(end).PDset.name;
        end
        function out = get.Dset(obj)
            %Use the first PDset in iter just because it uses less memory
            Units = obj.iter(1).PDset.ModelAndObservationPairs;
            ParamAssns = obj.iter(1).PDset.FreeParameter;
            name = obj.iter(1).PDset.name;
            out = DCDataset(Units,ParamAssns);
            out.name=name;
        end
        function out = get.optimOpts(obj)
            opt = obj.iter(end).optimOpts;
            out = opt;
        end
        function out = get.surfFittingOpts(obj)
            out = obj.iter(end).PDset.fittingSettings;
        end
        function out = get.PDset(obj)
            out = obj.iter(end).PDset;
        end
        function out = get.iterStruct(obj)
            out = obj.iter;
        end
        
        function bool = isempty(obj)
            % isempty Check whether a ParameterOptimization object is empty
            %
            % bool = isempty(obj) returns true of the object's date property is empty.
            %
            % See also ParameterOptimization

            if isempty(obj.runDate)
                bool = true;
            else
                bool = false;
            end
        end
        function [list,sz] = displayProps(obj)
            list = [];
            sz = '';
        end
    end %public methods
    
end %classdef
