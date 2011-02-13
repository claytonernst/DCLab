classdef ConsistencyTest < DClab.DCObject
    %ConsistencyTest computes the consistency measure of a dataset
    %
    %      ConsistencyTest starting from a DCDataset or PolyDataset:
    %
    %   TESTOBJ = ConsistencyTest(DSET) bounds the value of the consistency
    %   measure of the DCDataset or PolyDataset object DSET.
    %
    %   TESTOBJ = ConsistencyTest(DSET,DCOPT) allows you to supply options
    %   to the algorithm via the DCOptions object DCOPT.
    %
    %      ConsistencyTest starting from an existing ConsistencyTest
    %      object:
    %
    %   TESTOBJ = ConsistencyTest(OLDTESTOBJ) performs additional iterations of
    %   the algorithm, warm-starting with the ConsistencyTest object OLDTESTOBJ.
    %
    %   TESTOBJ = ConsistencyTest(OLDTESTOBJ,DCOPT) uses options from the
    %   DCOptions object DCOPT for the additional iterations of the algorithm.
    %
    % Inputs:
    %   DSET: A DCDataset object or a PolyDataset object previously formed
    %      from a DCDataset object. If DSET is a PolyDataset and DCOPT is not
    %      supplied, and further subdivision of the domain of PDSET will use
    %      the default fitting options.
    %   OLDTESTOBJ: An existing ConsistencyTest object.
    %   DCOPT[optional]: A DCOptions object.  If not suppled or empty, the
    %     default options obtained from the DCOptions constructor are used. If
    %     maxBranchBoundIter = 1, no subdividing will occur. If either
    %     omitInnerBound or omitOuterBound is true, you will not be able to
    %     perform additional branch-and-bound iterations at a future time by
    %     calling this constructor with the resulting TESTOBJ.
    %
    % Outputs:
    %   TESTOBJ: a ConsistencyTest object with a fields .iter and .runDate. .iter is a
    %       structure array with fields .optimOpts, .PDset, .worstUBleaf, .dimSplit,
    %       .splitLoc, and .leafInfo. These fields of the i'th
    %       element of .iter contain the following data:
    %     optimOpts: structure containing the optimization options for the i'th
    %       iteration
    %     PDSet: Partitioned PolyDataset used for the i'th iteration. Saving
    %       this object every iteration is probably unnessary since in principle
    %       the object used in earlier iterations should be recoverable from the
    %       object used in the last iteration.
    %     worstUBleaf: indicates which leaf of the current PDset gives the
    %       worst upper bound. This leaf will be divided on the next
    %       iteration. This field will be empty if opt.omitOuterBound = true.
    %     dimSplit: indicates which dimension (i.e., parameter name) of the
    %       the iter(i1-1).worstUBleaf was divided by the current iteration.
    %       iter(1).dimSplit will be empty.
    %     splitLoc: indicates where dimSplit was divided. iter(1).splitLoc will
    %       be emmpty
    %     leafInfo: big structure with fields .leafIndices, .upper.bnd,
    %       .upper.mult, .lower.bnd, .lower.mult, .lower.xfeas. that
    %       contain the results for the optimizations over each leaf listed in
    %       leafInfo.leafIndices
    
    
    properties
        iter;
        runDate;
    end
    
    properties (Dependent)
        UB;
        LB;
        LBx;
        upperBndMults;
        upperBndSens;
        nIter;
        DatasetName;
        Dataset;
        optimOpts;
        surfFittingOpts;
        PolyDataset;
        iterStruct;
    end
    
    
    methods
        
        function testObj = ConsistencyTest(varargin)
            
            % Calling syntax we will support
            % ConsistencyTest
            % ConsistencyTest(DSET)    (DSET is type DCDataset or PolyDataset)
            % ConsistencyTest(DSET,OPT)
            % ConsistencyTest(DSET,[])
            % ConsistencyTest(OLDOBJ)
            % ConsistencyTest(OLDOBJ,OPT)
            % ConsistencyTest(OLDOBJ,[])
            
            % We will warn if OLDOBJ or DSET are empty.
            
            % Our strategy is to massage the input list so that it looks like
            % ConsistencyTest(DSET,OPT) or
            % ConsistencyTest(OLDOBJ,OPT).
            
            %TODO, to save memory, we should probably delete PDSet
            %from all but the last element of iter.
            
            %TODO: all errors and warnings should be displayed with dcdispstr ??
            
            % Input error checking. If the second input is empty or not given, replace it with
            % DCOptions.
            ni = nargin;
            no = nargout;
            error(nargchk(0,2,ni));
            error(nargoutchk(0,1,no));
            
            % Special cases
            if ni == 0
                return
            end
            if isempty(varargin{1})
                warning('DClab:invalidInputs','First input to CONSISTENCYTEST is empty, returning an empty object')
                return
            end
            
            % Massage the inputs
            assert(isa(varargin{1},'DClab.ConsistencyTest') || isa(varargin{1},'DClab.DCDataset'),'Inputs: 1st input to CONSISTENCYTEST of improper class');
            if isa(varargin{1},'DClab.ConsistencyTest')
                warm = true;
                assert(ni==1 || ni==2,'Inputs: when 1st is a ConsistencyTest object, at most two inputs are allowed.');
                if ni == 1
                    varargin{2} = DClab.DCOptions; %fill in default
                elseif ni == 2
                    if isempty(varargin{2}) && isnumeric(varargin{2})
                        varargin{2} = DClab.DCOptions;
                    end
                    assert(isa(varargin{2},'DClab.DCOptions'),'Inputs: 2nd input to CONSISTENCYTEST of improper class, expecting type DCOptions.');
                end
            elseif isa(varargin{1},'DClab.DCDataset')
                % Since PolyDataset inherits from DCDataset, we only need to check isa DCDataset
                warm = false;
                assert(ni==1 || ni==2,'Inputs: when 1st is a ConsistencyTest object, at most two inputs are allowed.');
                if ni == 1
                    varargin{2} = DClab.DCOptions; %fill in default
                elseif ni == 2
                    if isempty(varargin{2}) && isnumeric(varargin{2})
                        varargin{2} = DClab.DCOptions;
                    end
                    assert(isa(varargin{2},'DClab.DCOptions'),'Inputs: 2nd input to CONSISTENCYTEST of improper class, expecting type DCOptions.');
                end
            end
            
            %Initialization. Define opt, guihand, dbSaveName, and iter.
            %At this point length(varargin)==2
            
            % Define some useful constants:
            opt = varargin{2};
            guihand = opt.guiHandle;
            
            if ~isempty(opt.fileName2Save)
                dbSaveName = opt.fileName2Save;
            else
                dbSaveName='';
            end
            
            % This variable will become true if we failed to meet the bbTermTol or if we
            % received an old TEST object and no additional iterations were performed.  When
            % true, opt.display='notify' will cause the exit message to be displayed.
            
            weirdexit = false;
            
            if ~warm
                %We're starting from scratch. Note, since all erquested bounds are
                %computed during initialization, we consider the initialization step to
                %be the first iteration.
                
                stime = clock;
                iterStruct = BBinitialization(varargin{:},guihand,dbSaveName);
                
                if ~isempty(dbSaveName)
                    save(dbSaveName,'iter')
                end
                
                % Determine where to subdivide on the next iteration
                if opt.omitOuterBound
                    UB = NaN;
                    LB = max(iterStruct.lower.bnd);
                elseif opt.omitInnerBound
                    UB = max(iterStruct.upper.bnd);
                    LB = NaN;
                else
                    [UB UBIdx] = max(iterStruct.upper.bnd);
                    LB = max(iterStruct.lower.bnd);
                    iterStruct.worstUBleaf = iterStruct.leafNodes(UBIdx);
                end
                
                % Should we proceed into the while loop? If not, determine the exit
                % condition.
                if UB-LB <= opt.branchBoundTermTol
                    term = true;
                    exitmsg = 'Exiting with 1 iteration: branchBoundTermTol met';
                elseif opt.maxBranchBoundIter <= 1
                    term = true;
                    weirdexit = true;
                    exitmsg = 'Exiting: maximum number of branch and bound iterations reached';
                elseif opt.omitInnerBound || opt.omitOuterBound
                    term = true;
                    weirdexit = true;
                    exitmsg = 'Exiting with 1 iteration: both inner and outer bounds required to branch and bound';
                else
                    term = false;
                    i1 = 1;
                end
                iterStruct.runtime = etime(clock,stime);
                
            else
                % We've been warm started
                iterStruct = varargin{1}.iter;
                UB = varargin{1}.UB;
                LB = varargin{1}.LB;
                
                % Should we proceed into the while loop? If not, determine the exit
                % condition.
                if opt.omitInnerBound || opt.omitOuterBound
                    term = true;
                    weirdexit = true;
                    exitmsg = 'Returning original object: both inner and outer bounds are required to branch and bound';
                elseif opt.maxBranchBoundIter <= length(iterStruct)
                    term = true;
                    weirdexit = true;
                    exitmsg = 'Returning original object: maxBranchBoundIter already met in the previous iteration';
                elseif isnan(UB) || isnan(LB)
                    term = true;
                    weirdexit = true;
                    exitmsg = 'Returning original object: both inner and outer bounds were not computed in the previous iteration';
                elseif UB-LB <= opt.branchBoundTermTol
                    term = true;
                    weirdexit = true;
                    exitmsg = 'Returning original object: branchBoundTermTol already met in the previous iteration';
                else
                    %If we made it this far, we're ready to perform another iteration.
                    if ismember(opt.display,{'iter';'all';'ALL'})
                        str = 'First displayed bounds are for the old ConsistencyTest object';
                        DClab.dcdispstr(str,guihand,false)
                    end
                    term = false;
                    i1 = length(iterStruct);
                end
            end
            
            % Conditional to enter the branch and bound iteration
            while ~term
                
                stime = clock;
                % Increament the iteration count
                i1 = i1+1;
                [dispOpts fitOpts iterStruct(i1,1).optimOpts] = decompose(opt);
                
                % Decide where to divide prev worstUBleaf
                
                [trash Uidx] = max(iterStruct(i1-1).upper.bnd);
                
                %                 if strcmp(opt.analysisMode,'original')
                % If we're in original mode, assume the problem is with duality.
                
                pt = [];
                mults = [];
                Q = iterStruct(i1-1).upper.Q{Uidx};
                leaf = iterStruct(i1-1).leafNodes(Uidx);
                computeGap = false;
                sigma_t = LB; %dummy value. This won't be used.
                
                %display for fun
                if ismember(opt.display,{'iter';'all';'ALL'})
                    str = ['UB: ' num2str(UB)];
                    DClab.dcdispstr(str,guihand,false)
                    str = ['LB: ' num2str(LB)];
                    DClab.dcdispstr(str,guihand,false)
                end
                
                %                 else
                %                     % If we're in one of the metabased modes, determine if the present
                %                     % issue is duality or approximation.
                %
                %                     [sigma_t sidx] = max(iter(i1-1).upper.qbnd);
                %                     CDmm = max(iter(i1-1).upper.sbnd);
                %
                %                     %display for fun
                %                     if ismember(opt.display,{'iter';'all';'ALL'})
                %                         str = ['UB: ' num2str(UB)];
                %                         DClab.dcdispstr(str,guihand,false)
                %                         str = ['t:  ' num2str(sigma_t)];
                %                         DClab.dcdispstr(str,guihand,false)
                %                         str = ['s:  ' num2str(CDmm)];
                %                         DClab.dcdispstr(str,guihand,false)
                %                         str = ['LB: ' num2str(LB)];
                %                         DClab.dcdispstr(str,guihand,false)
                %                     end
                %
                %                     if sigma_t-s <= opt.branchBoundQGapTol
                %                         % problem is duality
                %                         pt = [];
                %                         mults = [];
                %                         Q = iter(i1-1).upper.Q{Uidx};
                %                         leaf = iter(i1-1).leafNodes(Uidx);
                %                         computeGap = false;
                %                         if ismember(opt.display,{'all';'ALL'})
                %                             str = 'Branching will seek to reduce the d-gap';
                %                             DClab.dcdispstr(str,guihand,false)
                %                         end
                %                     else
                %                         % problem is fitting error
                %                         pt = iter(i1-1).upper.surrX{sidx};
                %                         mults = max(iter(i1-1).upper.pairMult{sidx},[],2);
                %                         Q = [];
                %                         leaf = iter(i1-1).leafNodes(sidx);
                %                         computeGap = true;
                %                         if ismember(opt.display,{'all';'ALL'})
                %                             str = 'Branching will seek to reduce the q-gap';
                %                             DClab.dcdispstr(str,guihand,false)
                %                         end
                %                     end
                %                 end
                
                [dimSplit, splitLoc, inherit] = DClab.ConsistencyTest.findDimension2Divide(iterStruct(i1-1).PDset,leaf,pt,mults,Q,UB,sigma_t);
                
                iterStruct(i1).dimSplit = dimSplit;
                iterStruct(i1).splitLoc = splitLoc;
                
                if any(strmatch(opt.display,{'iter';'all';'ALL'}))
                    str = ['===Iteration ' num2str(i1) ': spliting dimension ' dimSplit ' of leaf ' num2str(leaf) ...
                        ' at location ' num2str(splitLoc) '==='];
                    DClab.dcdispstr(str,guihand,false)
                end
                
                % Subdivide the domain of PDset where indicated.
                % First provide the current display options.
                oldPDset = iterStruct(i1-1).PDset;
                oldPDset.displaySettings = dispOpts;
                iterStruct(i1).PDset = subdivideDomain(oldPDset,leaf,dimSplit,splitLoc,inherit,opt);
                iterStruct(i1).leafNodes = iterStruct(i1).PDset.leafNodes;
                
                %TODO we're erasing the old pdset to save memory
                iterStruct(i1-1).PDset = [];
                
                % If requested, save iter structure
                if ~isempty(dbSaveName)
                    save(dbSaveName,'iter')
                end
                
                % Find which leaves were added with that division
                newLeaves = setdiff(iterStruct(i1).leafNodes,iterStruct(i1-1).leafNodes);
                
                %we only need do optimization on the newly added leaves, then we
                %will zero out the bnd, mults, logX, from the parent of these
                %new leaves and add these quantities from the to new
                %leaves. that way we alway have available bnd, mults, and logX
                %for each leaf.
                try
                    if ismember(opt.display,{'iter';'all';'ALL'})
                        str = '  Computing upper bound';
                        DClab.dcdispstr(str,guihand,false)
                    end
                    
                    ub = DClab.ConsistencyTest.maximizeUB(iterStruct(i1).PDset,opt,newLeaves,computeGap);
                    
%                     if ~computeGap
%                         % qbnd and sbnd were not updated. fill in with old values.
%                         val1 = max(iterStruct(i1-1).upper.qbnd);
%                         ub.q.bnd = [val1;val1];
%                         
%                         [val2 idx] = max(iterStruct(i1-1).upper.sbnd);
%                         ub.mm.bnd = [val2; val2];
%                         ub.mm.xfeas = [iterStruct(i1-1).upper.surrX(idx); iterStruct(i1-1).upper.surrX(idx)];
%                     end
                    
                    if ismember(opt.display,{'iter';'all';'ALL'})
                        str = '  Computing lower bound';
                        DClab.dcdispstr(str,guihand,false)
                    end
                    
                    %if computeGap % we're still trying to reduce q-gap
                    lb = DClab.ConsistencyTest.maximizeLB(iterStruct(i1).PDset,opt,newLeaves);
                    %else
                    %TODO: Don't bother recomputing???
                    %Actually, we probably should. Our lower bounds can't be that good.
                    [tmp idx] = max(iterStruct(i1-1).lower.bnd);
                    lb.pound.bnd = tmp*[1;1];
                    lb.pound.xfeas = [iterStruct(i1-1).lower.xfeas(idx);iterStruct(i1-1).lower.xfeas(idx)];
                    %end
                    
                catch ME
                    str = 'failed call to maximizeUB or maximizeLB in CONSISTENCYTEST';
                    DClab.dcdispstr(str,guihand,false)
                    DClab.dcdispstr(ME,guihand,true);
                end
                
                %do the zeroing out an updating alluded to before. three elements
                %of the iter stuct need to be updated/modified: iter.leafInfo.upper,
                %iter.leafInfo.lower1, and iter.leafInfo.lower2.
                parentIdx = find(iterStruct(i1-1).leafNodes == leaf);
                
                % Copy previous values
                iterStruct(i1,1).upper = iterStruct(i1-1).upper;
                iterStruct(i1,1).lower = iterStruct(i1-1).lower;
                
                % Update the leaf Indices
                iterStruct(i1,1).leafNodes = iterStruct(i1).PDset.leafNodes;
                
                %Eliminate info relevant to the leaf that was split and append info from the two daughter leaves
                
                iterStruct(i1).upper.bnd(parentIdx) = [];
                %iterStruct(i1).upper.qbnd(parentIdx) = [];
                %iterStruct(i1).upper.sbnd(parentIdx) = [];
                %iterStruct(i1).upper.surrX(parentIdx) = [];
                iterStruct(i1).upper.pairMult(parentIdx) = [];
                iterStruct(i1).upper.paramMult(parentIdx) = [];
                iterStruct(i1).upper.Q(parentIdx) = [];
                iterStruct(i1).upper.api(parentIdx) = [];
                
                iterStruct(i1).upper.bnd = [iterStruct(i1).upper.bnd; ub.star.bnd];
                %iterStruct(i1).upper.qbnd = [iterStruct(i1).upper.qbnd; ub.q.bnd];
                %iterStruct(i1).upper.sbnd = [iterStruct(i1).upper.sbnd; ub.mm.bnd];
                %iterStruct(i1).upper.surrX = [iterStruct(i1).upper.surrX; ub.mm.xfeas];
                iterStruct(i1).upper.pairMult = [iterStruct(i1).upper.pairMult; ub.star.pairMult];
                iterStruct(i1).upper.paramMult = [iterStruct(i1).upper.paramMult; ub.star.paramMult];
                iterStruct(i1).upper.Q = [iterStruct(i1).upper.Q; ub.star.Q];
                iterStruct(i1).upper.api = [iterStruct(i1).upper.api; ub.star.api];
                
                iterStruct(i1).lower.bnd(parentIdx) = [];
                iterStruct(i1).lower.xfeas(parentIdx) = [];
                iterStruct(i1).lower.bnd = [iterStruct(i1).lower.bnd; lb.pound.bnd];
                iterStruct(i1).lower.xfeas = [iterStruct(i1).lower.xfeas; lb.pound.xfeas];
                
                [UB UBIdx] = max(iterStruct(i1).upper.bnd);
                iterStruct(i1).worstUBleaf = iterStruct(i1).leafNodes(UBIdx);
                LB = max(iterStruct(i1).lower.bnd);
                
                if ~isempty(dbSaveName)
                    save(dbSaveName,'iter')
                end
                
                if UB - LB <= opt.branchBoundTermTol
                    term = true;
                    exitmsg = ['Exiting after ' num2str(i1) ' iterations: branchBoundTermTol met'];
                elseif i1 >= opt.maxBranchBoundIter
                    term = true;
                    weirdexit = true;
                    exitmsg = 'Exiting: maximum number of branch and bound iterations reached';
                else
                    % do nothing, we need to keep iterating
                end
                iterStruct(i1).runtime = etime(clock,stime);
                
            end
            
            if strcmp(opt.display,'notify') && weirdexit
                DClab.dcdispstr(exitmsg,guihand,false)
            elseif any(strmatch(opt.display,{'final';'iter';'all';'ALL'}))
                DClab.dcdispstr(exitmsg,guihand,false)
                str = ['UB: ' num2str(UB)];
                DClab.dcdispstr(str,guihand,false)
                str = ['LB: ' num2str(LB)];
                DClab.dcdispstr(str,guihand,false)
                if UB < 0
                    str = '=====The DCDataset is Inconsistent=====';
                elseif LB > 0
                    str = '=====The DCDataset is Consistent=====';
                else
                    str = '=====The consistency test was inconclusive=====';
                end
                DClab.dcdispstr(str,guihand,false)
            else
                %display nothing
            end
            
            %assign output
            testObj.iter = iterStruct;
            testObj.runDate = date;
            
        end
        
        function out = get.UB(obj)
            opt = obj.iter(end).optimOpts;
            if ~opt.omitOuterBound
                out = max(obj.iter(end).upper.bnd);
            else
                out = nan;
            end
        end
        function out = get.LB(obj)
            opt = obj.iter(end).optimOpts;
            if ~opt.omitInnerBound
                out = max(obj.iter(end).lower.bnd);
            else
                out = nan;
            end
        end
        function out = get.LBx(obj)
            opt = obj.iter(end).optimOpts;
            if ~opt.omitInnerBound
                [LB LBidx] = max(obj.iter(end).lower.bnd);
                out = obj.iter(end).lower.xfeas{LBidx};
            else
                out = [];
            end
        end
        function out = get.upperBndMults(obj)
            opt = obj.iter(end).optimOpts;
            if opt.omitOuterBound
                out = struct([]);
            else
                [UB UBidx] = max(obj.iter(end).upper.bnd);
                upperBndMults.expl = obj.iter(end).upper.pairMult{UBidx}(:,1);
                upperBndMults.expu = obj.iter(end).upper.pairMult{UBidx}(:,2);
                tmpl = obj.iter(end).upper.paramMult{UBidx}(:,1);
                tmpu = obj.iter(end).upper.paramMult{UBidx}(:,2);
                
                % Combine parameters for linX and logX.
                PDset = obj.iter(end).PDset;
                n = PDset.nParameters;
                uncCase = vertcat(PDset.FreeParameter.uncertaintyCase);
                range = vertcat(PDset.FreeParameter.range);
                if length(tmpl) > n
                    both = true;
                else
                    both = false;
                end
                
                paraml = zeros(n,1);
                paramu = zeros(n,1);
                
                %TODO: check that this math makes some sort of sense.
                for i1 = 1:n
                    if uncCase(i1)==1 || uncCase(i1)==2
                        % Return multipliers for lin X.
                        if both
                            paraml(i1) = tmpl(i1) + tmpl(i1+n)/log(10)/range(i1,1);
                            paramu(i1) = tmpu(i1) + tmpu(i1+n)/log(10)/range(i1,2);
                        else
                            paraml(i1) = tmpl(i1);
                            paramu(i1) = tmpu(i1);
                        end
                        
                    else
                        if both
                            paraml(i1) = tmpl(i1)*log(10)*range(i1,1)+tmpl(i1+n);
                            paramu(i1) = tmpu(i1)*log(10)*range(i1,2)+tmpu(i1+n);
                        else
                            paraml(i1) = tmpl(i1)*log(10)*range(i1,1);
                            paramu(i1) = tmpu(i1)*log(10)*range(i1,2);
                        end
                    end
                end
                
                upperBndMults.paraml = paraml;
                upperBndMults.paramu = paramu;
                out = upperBndMults;
                
            end
        end
        function out = get.upperBndSens(obj)
            mults = obj.upperBndMults;
            UB = obj.UB;
            out.expl = -mults.expl*(1-UB);
            out.expu = mults.expu*(1-UB);
            
            [UB UBidx] = max(obj.iter(end).upper.bnd);
            worstNode = obj.iter(end).leafNodes(UBidx);
            worstRange = obj.iter(end).PDset.PiecewiseSurrogateModelTree(worstNode).domainRange;
            totalRange = vertcat(obj.iter(end).PDset.FreeParameter.range);
            api = obj.iter(end).upper.api{UBidx};
            n = obj.iter(end).PDset.nParameters;
            if size(totalRange,1)==n
                Q = zeros(n+1);
            else
                Q = zeros(2*n+1);
            end
            Q([1;api],[1;api]) = obj.iter(end).upper.Q{UBidx};
            
            out.paraml = mults.paraml;
            out.paramu = mults.paramu;
            
            n = size(totalRange,1);
            if length(mults.paraml)==n
                for i1 = 1:size(totalRange,1)
                    if abs(worstRange(i1,1)-totalRange(i1,1)) < 10*eps
                        out.paraml(i1) = mults.paraml(i1).*(Q(1+i1,1)-1.1);
                    else
                        out.paraml(i1) = 0;
                    end
                    if abs(worstRange(i1,2)-totalRange(i1,2)) < 10*eps
                        out.paramu(i1) = mults.paramu(i1)*(Q(1+i1,1)+1.1);
                    else
                        out.paramu(i1) = 0;
                    end
                end
            else
                for i1 = 1:size(totalRange,1)
                    if abs(worstRange(i1,1)-totalRange(i1,1)) < 10*eps
                        out.paraml(i1) = mults.paraml(i1).*(0.5*(Q(1+i1,1)+Q(1+i1+n,1))-1.1);
                    else
                        out.paraml(i1) = 0;
                    end
                    if abs(worstRange(i1,2)-totalRange(i1,2)) < 10*eps
                        out.paramu(i1) = mults.paramu(i1)*(0.5*(Q(1+i1,1)+Q(1+i1+n,1))+1.1);
                    else
                        out.paramu(i1) = 0;
                    end
                end
            end
            %now convert back to original units
            out.paramu = max(out.paramu.*2./diff(totalRange,[],2),0);
            out.paraml = min(out.paraml.*2./diff(totalRange,[],2),0);
        end
        function out = get.nIter(obj)
            out = length(obj.iter);
        end
        function out = get.DatasetName(obj)
            out = obj.iter(end).PDset.name;
        end;
        function out = get.Dataset(obj)
            %Use the first PDset in iter just because it uses less memory
            Pairs = obj.iter(end).PDset.ModelAndObservationPair;
            FreeParams = obj.iter(end).PDset.FreeParameter;
            name = obj.iter(end).PDset.name;
            out = DClab.DCDataset(Pairs,FreeParams);
            out.name=name;
        end
        function out = get.optimOpts(obj)
            out = obj.iter(end).optimOpts;
        end
        function out = get.surfFittingOpts(obj)
            out = obj.iter(end).PDset.fittingSettings;
        end
        function out = get.PolyDataset(obj)
            out = obj.iter(end).PDset;
        end
        function out = get.iterStruct(obj)
            out = obj.iter;
        end
        
        function bool = isempty(obj)
            % bool = isempty(obj)
            %
            % A ConsistTest object is considered empty if obj(1).PDset is.
            %
            % See also ConsistTest
            
            if isempty(obj.runDate)
                bool = true;
            else
                bool = false;
            end
        end
        function [list,sz] = displayProps(obj)
            len = length(obj.iter);
            if len==0
                list = [];
                sz = 'empty';
            else
                list{1} = sprintf('Consistency measure bounds: [%0.4g, %0.4g]',obj.LB,obj.UB);
                if obj.LB>0
                    list{2} = 'Dataset is consistent';
                elseif obj.UB<0
                    list{2} = 'Dataset is inconsistent';
                else
                    list{2} = 'Consistency test was inconclusive:';
                    list{3} = '  More iterations may resolve sign of consistency measure.';
                end
                sz = sprintf('%d iteration',len);
            end
        end
    end %public methods
    
    methods (Access=private, Static)
        %[dim, loc, inherit] = findDimension2Divide(PDset,leaf,pt,mults,Q,UB,qbnd);
        %[obj,activeIdx] = makeCmoptim(PD,nodeIdx,approxType);
        %[s1 s2] = maximizeLB(PD,opt,nodes,computeGap,startPt);
        %[fval xfeas] = maximizeLB_meta(MMCell,d,u,uncCase,FP,domRng);
        %[s1 s2] = maximizeUB(PD,opt,nodes,computeGap,startPt);
        %[fval xfeas] = maximizeUB_meta(MMCell,d,u,uncCase,FP,domRng);
        
        function [obj,activeIdx] = makeCmoptim(PD,nodeIdx,approxType)
            %MAKECMOPTIM sets up an optimization problem
            %
            %   OBJ = MAKECMOPTIM(PD,NODEIDX,APPROXTYPE) creates an optimization
            %   problem object for the consistency test problem. The objective is the
            %   slack variable gamma.  PD is a PolyDataset that describes the
            %   constraints.
            %
            %   [OBJ,ACTIVEIDX] = MAKECMOPTIM(...) additionally returns a vector
            %   indicating which parameters in the PolyDataset became optimization
            %   variables. 1:n indicate the linX variables; n+1:n*2 indicate the logX
            %   variables. ACTIVEIDX will be a subset of 1:2*n.
            %
            %   Inputs:
            %      PD: a PolyDataset object that describes the constraints.
            %      NODEIDX: is a scalar that specifies the domain subdomain of PD for the
            %         problem.
            %      APPROXTYPE: either 'inner', 'outer', 'crude'. This indicates how
            %         the fitting errors of the surrogate models in PD should be
            %         handled. If 'inner', used them to pose the optimization over an
            %         inner approximation of the feasible set. If 'outer', used them to
            %         pose the optimization over an outer approximation of the feasible
            %         set. If 'crude', ignore them.
            %   Outputs:
            %      OBJ: an nqcqp object using surrogate models that are quadratic in
            %         either the parameters and/or log10 of the parameters.
            
            % code assumes activeIdx is the same for linX and logX constraints
            
            %error(nargchk(3,3,nargin,'struct'));
            
            % Determine which parameters are active. There are n parameters in the
            % dataset. Pretend there are 2n, with the first n untransformed and the
            % last n with a log10 transformation. The indices we produce will be
            % indices into this ficticious 2n-by-1 parameter list.
            np = PD.nParameters;
            conSurfs = surfaces(PD.PiecewiseSurrogateModelTree,nodeIdx);
            %surfs2RM = surface2CorrespondingResponseModel(PD.PiecewiseSurrogateModelTree(nodeIdx));
            nsPerRM = PD.PiecewiseSurrogateModelTree(nodeIdx).nSurfacesPerResponseModel;
            ns = length(conSurfs);
            
            surfActIdx = cell(ns,1);
            for i1 = 1:ns
                surfActIdx{i1} = conSurfs(i1).activeParameterIndex;
                trans = conSurfs(i1).activeParameterTransformation;
                surfActIdx{i1}(trans==2) = surfActIdx{i1}(trans==2) + np;
            end
            
            % The call to unique will sort activeIdx
            activeIdx = unique(vertcat(surfActIdx{:}));
            n = length(activeIdx);
            
            nQuad = sum(nsPerRM);
            Zquad = cell(nQuad,1);
            d = zeros(nQuad,1);
            u = zeros(nQuad,2);
            ou = zeros(nQuad,2);
            se = zeros(nQuad,2);
            uncCase = zeros(nQuad,1);
            respTrans = cell(nQuad,1);
            
            % Loop through constraint response models and assemble constraints
            for i1 = 1:PD.nPairs
                tmp_uncCase = PD.ModelAndObservationPair(i1).uncertaintyCase;
                tmp_d = PD.ModelAndObservationPair(i1).observedValue;
                tmp_u = PD.ModelAndObservationPair(i1).observationUncertaintyPlusMinus;
                tmp_ou = PD.ModelAndObservationPair(i1).ResponseModel.outputUncertaintyPlusMinus;
                
                % We want to "enforce" the constraints
                % L <= M(x) <= U
                %
                % Using se1 <= S(x) - M(x) <= se2 or se1 <= S(x) - log10(M(x)) <= se2
                
                if i1 == 1
                    lastSurfIdx = 0;
                else
                    lastSurfIdx = sum(nsPerRM(1:i1-1));
                end
                
                for i2 = 1:nsPerRM(i1)
                    [sortedSurfParam surfParam2SortedSurfParam] = sort(surfActIdx{lastSurfIdx+i2});
                    [trash sortedSurfParam2LocInActive] = intersect(activeIdx,sortedSurfParam);
                    surfParam2LocInActive = zeros(size(sortedSurfParam2LocInActive));
                    surfParam2LocInActive(surfParam2SortedSurfParam) = sortedSurfParam2LocInActive;
                    
                    tmp_se = conSurfs(lastSurfIdx+i2).peakError;
                    poly = conSurfs(lastSurfIdx+i2).surrogateModel;
                    nact = length(surfActIdx{lastSurfIdx+i2});
                    tmp_respTrans = conSurfs(lastSurfIdx+i2).responseTransformation;
                    
                    switch approxType
                        case 'outer'
                            % do nothing, fitErr is correct
                        case 'crude'
                            tmp_se = [0 0];
                        case 'inner'
                            tmp_se = [tmp_se(2) tmp_se(1)];
                    end
                    
                    tmp_Quad = spalloc(n+1,n+1,(nact+1)^2);
                    tmp_Quad([1;surfParam2LocInActive+1],[1;surfParam2LocInActive+1]) = poly;
                    
                    Zquad{(lastSurfIdx+i2)} = tmp_Quad;
                    d(lastSurfIdx+i2) = tmp_d;
                    u(lastSurfIdx+i2,:) = tmp_u;
                    ou(lastSurfIdx+i2,:) = tmp_ou;
                    se(lastSurfIdx+i2,:) = tmp_se;
                    uncCase(lastSurfIdx+i2) = tmp_uncCase;
                    respTrans{lastSurfIdx+i2} = tmp_respTrans;
                end
                
            end % for i1 = 1PD.nPairs
            
            % Find which variables employ both transformations
            linIdx = activeIdx <= np;
            
            [trash xActWithEqConstLin xActWithEqConstLog] = intersect(activeIdx(linIdx),activeIdx(~linIdx)-np);
            xActWithEqConstLog = xActWithEqConstLog+sum(linIdx);
            
            % For i\in xActWithEqConstLin, activeIdx(i) should be the linear version of
            % a parameter that also has a log version. Similarly for j\in
            % xActWithEqConstLog, activeIdx(j) should be the log version of a parameter
            % that also has a linear version.
            
            range = PD.PiecewiseSurrogateModelTree(nodeIdx).domainRange;
            nElim = length(xActWithEqConstLog);
            linXlogX = zeros(nElim,5);
            
            for i1 = 1:nElim
                ii = xActWithEqConstLin(i1);
                jj = xActWithEqConstLog(i1);
                
                a = range(activeIdx(ii),1);
                b = range(activeIdx(ii),2);
                
                c1 = b/a;
                c2 = a;
                c3 = b-a;
                
                linXlogX(i1,:) = [ii jj c1 c2 c3];
            end
            
            %constraint bounds on x's. since everything is normalized, this
            %is trivial.
            LB = -ones(n,1);
            UB = ones(n,1);
            
            obj = DClab.cmoptim(Zquad,d,u,LB,UB,ou,uncCase,se,respTrans,linXlogX);
        end
        
        function [dim, loc, inherit] = findDimension2Divide(PDset,leaf,pt,mults,Q,UB,qbnd)
            % function [dim, value, inherit] = findDimension2Divide(PDset,leaf,pt,mults,Q,UB,qbnd)
            %
            % Functionality: This method determines where to subdivide the hypercube.
            % It contains two algorithms.  The first is used to shrink the gap between
            % the local solution to a quadratic problem and the corresponding
            % outerbound generated by the s-procedure.  The second algorithm generates
            % the subdivision location based on the fitting error of the quadratic
            % surrogates.  If "pt" is empty, the first algorithm is used.
            %
            % Inputs
            %  PDset: a PolyDataset object with the same partition structure as PSMTree
            %  leaf: the index to the leaf of the tree that we are currently working on
            %  pt: a x vector used by the 2nd algorithm
            %  mults: multiplier vector
            %  Q: Q matrix from the rank relaxation
            %  UB: upper bound from the rank relaxation
            %  qbnd: lower bound from the local search
            %
            % Outputs:
            %   dim: a string naming the dimension to be divided
            %   loc: the location of the cutting plane in that dimension (in the
            %   original parameter range units)
            %   inherit: currently defaults to an array of zeros, indicating
            %     that new surfaces should be computed for each model upon subdivision.
            %
            % ALGORITHM 1 (REDUCE D_GAP)
            %   The s-procedure is exact if the Q matrix has rank 1.  The code
            %   implements a heuristic which subdivide the cube by analyzing the
            %   eigenvector corresponding to the 2nd largest eigenvalue of Q.  In the
            %   worst case, the algorithm is brute force, but in practice is will
            %   significantly reduce the gap in relatively few iterations.
            %
            % ALGORITHM 2 (REDUCE ?_GAP)
            %   The algorithm assumes that the local solution pt is actually the global
            %   solution.  It then computes how much just the fitting error would cost us in
            %   terms of the upper bound.
            %
            %
            % See also PolyDataset/subdivideDomain, PiecewiseSurrogateModelTree/subdivideDomain
            
            
            pList = PDset.parameterList;
            
            if isempty(pt)
                %assume we're trying to reduce d-gap
                
                % Create an optimization problem object, but don't solve it.  We tweak
                % variable bounds in the object to evaluate the potential improvement for
                % splitting different dimensions.
                [optProb activeIdx] = DClab.ConsistencyTest.makeCmoptim(PDset,leaf,'outer');
                
                [V E] = eig(Q);
                tmp = V*E; %multiplying by E serves no purpose with the current code.
                
                % In an effort to compute the dimension to divide quickly, we perform
                % from one to
                % three passes.  Each time the threshholds are adjusted so that more and
                % more dimensions are considered.
                
                %==initial pass==
                
                % Look at the 2nd through last components of the 2nd largest eigenvector of
                % Q.  Ignore the first component because of [1 x].
                irank1 = abs(tmp(2:end,end-1));
                
                % Remove the constant offset and consider any dimensions with 75% or more
                % of the largest component.  These are "worthy" of consideration.
                irank1 = irank1-min(irank1);
                worthy = find(irank1 > 0.75*max(irank1));
                
                % The previous heuristic will select one or more dimensions as "worthy".
                % Performance is typically improved by considering more than just the
                % most worthy.  The while loop ensures at least three trial subdivisions
                % are made (provided the dimenion of x is at least 3).
                cnt = 1;
                while length(worthy) < 3 && length(sum(irank1>=0.01*max(irank1))) >= 3
                    worthy = find(irank1 > 0.75^(1+cnt)*max(irank1));
                    cnt = cnt+1;
                end
                
                % Now split at each worthy dimenions and call sedumi to see how much we
                % reduced the gap.  The brute force at this stage (pre actually
                % subdividing the cube in the overall scheme) is justified because we
                % must make intelligent cuts in the overall scheme or the number of
                % subcubes explodes.  For every "bad" cut, we must make two "good" cuts
                % to recover from the bad choice.
                
                % Reduce the sedumi tolerance from our default to speed computation. We
                % don't particularly need accuracy at this point.
                opt = DClab.DCOptions('sedumiParEps',1e-7);
                
                % Grab the original cube bounds.  Then split each worthy dimension and
                % evaluate the new upper bound both left and right subcube.
                default = optProb.variableBnds;
                
                tmpBnd = zeros(length(worthy),2);
                for i2 = 1:length(worthy) 
                    bndL = default;
                    bndL(worthy(i2),2) = 0; % set upper bound of the worthy dimension to 0 to generate a left subcube
                    optProbL = optProb;
                    optProbL.variableBnds = bndL; %evaluate the upper bound on the subcube
                    bndR = default;
                    bndR(worthy(i2),1) = 0;
                    optProbR = optProb;
                    optProbR.variableBnds = bndR;
                    tmpBnd(i2,1) = -lowerBnd(optProbL,opt); %lowerBnd performs a minimization.  flip the sign.
                    tmpBnd(i2,2) = -lowerBnd(optProbR,opt);
                end
                % Dump the results to the screen.
                %max(tmpBnd,[],2) %take a max because we need to improve the bound on both the left and right subcubes.
                %min(max(tmpBnd,[],2)) %determine which trial subdivision performed best on both subcubes.
                
                % Determine which of the worthy dimensions performed the best.
                [betterUB idx] = min(max(tmpBnd,[],2));
                
                % Look for a 5% reduction in the gap.  If we fail to acheive this,
                % enlarge the pool of worthy dimensions.
                if (UB - betterUB) < 0.05*(UB-qbnd)
                    %If we failed to acheive the desired improvement,
                    % make another pass with more dimensions considered.
                    
                    %disp('==========making pass 2===========')
                    
                    % Cache the previous results to avoid recomputing for these dimesions.
                    oldWorthy = worthy;
                    oldTmpBnd = tmpBnd;
                    
                    % This time we consider any dimensions whose components of the "2nd
                    % largest eigenvector of Q" are at least 1/3 as large as the largest.
                    worthy = find(irank1 > 0.33*max(irank1));
                    
                    % The previous heuristic will select one or more dimensions as "worthy".
                    % Performance is typically improved by considering more than just the
                    % most worthy.  The while loop ensures at least four trial subdivisions
                    % are made (provided the dimenion of x is at least 4).
                    cnt = 1;
                    while length(worthy) < 4 && length(sum(irank1>=0.01*max(irank1))) >= 4
                        worthy = find(irank1 > 0.33^(1+cnt)*max(irank1));
                        cnt = cnt+1;
                    end
                    
                    tmpBnd = zeros(length(worthy),2);
                    
                    % Determine where to insert the cached results from pass 1.
                    [asdf idxx] = intersect(worthy,oldWorthy);
                    tmpBnd(idxx,:) = oldTmpBnd;
                    
                    % Indices into worthy that represent the additional dimensions
                    % considered in pass 2.
                    newWorthy = setdiff(1:length(worthy),idxx);
                    
                    for i3 = reshape(newWorthy,1,length(newWorthy)) 
                        bndL = default;
                        bndL(worthy(i3),2) = 0;
                        optProbL = optProb;
                        optProbL.variableBnds = bndL;
                        bndR = default;
                        bndR(worthy(i3),1) = 0;
                        optProbR = optProb;
                        optProbR.variableBnds = bndR;
                        tmpBnd(i3,1) = -lowerBnd(optProbL,opt);
                        tmpBnd(i3,2) = -lowerBnd(optProbR,opt);
                    end
                    % Dump the results to the screen.
                    %max(tmpBnd,[],2) %take a max because we need to improve the bound on both the left and right subcubes.
                    %min(max(tmpBnd,[],2)) %determine which trial subdivision performed best on both subcubes.
                    
                    % Determine which of the worthy dimensions performed the best.
                    [betterUB idx] = min(max(tmpBnd,[],2));
                    
                end
                
                % Look for a 5% reduction in the gap.  If we fail to acheive this,
                % enlarge the pool of worthy dimensions to include all.
                
                % TODO we may
                % want to perform another intermediate pass before we jump all the way to
                % including all.
                
                if (UB - betterUB) < 0.05*(UB-qbnd)
                    %make another pass
                    %disp('==========making pass 3===========')
                    
                    % Cache the previous results to avoid recomputing for these dimesions.
                    oldWorthy = worthy;
                    oldTmpBnd = tmpBnd;
                    
                    % Include all dimensions as worthy.
                    worthy = (1:length(irank1))';
                    
                    tmpBnd = zeros(length(worthy),2);
                    
                    % Determine where to insert the cached results from pass 2 and insert them.
                    [asdf idxx] = intersect(worthy,oldWorthy);
                    tmpBnd(idxx,:) = oldTmpBnd;
                    
                    newWorthy = setdiff(1:length(worthy),idxx);
                    
                    for i2 = reshape(newWorthy,1,length(newWorthy)) 
                        bndL = default;
                        bndL(worthy(i2),2) = 0;
                        optProbL = optProb;
                        optProbL.variableBnds = bndL;
                        bndR = default;
                        bndR(worthy(i2),1) = 0;
                        optProbR = optProb;
                        optProbR.variableBnds = bndR;
                        tmpBnd(i2,1) = -lowerBnd(optProbL,opt);
                        tmpBnd(i2,2) = -lowerBnd(optProbR,opt);
                    end
                    
                    %max(tmpBnd,[],2)
                    %min(max(tmpBnd,[],2))
                    
                    [trash idx] = min(max(tmpBnd,[],2));
                    
                end
                % After we have performed between 1 to 3 passes to determine the best
                % dimension to subdivide, we try to subdivide at a few places in addition
                % to right in the middle.
                
                % If the bounds are +/- 1 (they are) we try to subdivide at -0.33, 0, and
                % 0.33.  We have the results from splitting at 0 from the previous
                % computations.
                
                % try at 2 more places.
                newBnd=zeros(3,2);
                
                newBnd(2,:) = tmpBnd(idx,:); %Insert prev result for split at 0.
                
                bndL = default;
                bndL(worthy(idx),2) = -0.33;
                optProbL = optProb;
                optProbL.variableBnds = bndL;
                bndR = default;
                bndR(worthy(idx),1) = -0.33;
                optProbR = optProb;
                optProbR.variableBnds = bndR;
                newBnd(1,1) = -lowerBnd(optProbL);
                newBnd(1,2) = -lowerBnd(optProbR);
                
                bndL = default;
                bndL(worthy(idx),2) = +0.33;
                optProbL = optProb;
                optProbL.variableBnds = bndL;
                bndR = default;
                bndR(worthy(idx),1) = +0.33;
                optProbR = optProb;
                optProbR.variableBnds = bndR;
                newBnd(3,1) = -lowerBnd(optProbL);
                newBnd(3,2) = -lowerBnd(optProbR);
                
                [best loc] = min(max(newBnd,[],2));
                
                bestLoc = loc-2;
                % loc == 1 indicates divide at x=0.33,
                % loc == 0 indicates divide at x=0,
                % loc == -1 indicates divide at x=-0.33.
                
                bestVar = activeIdx(worthy(idx));
                
                % Determine the subdivision location in the actual parameter ranges (as
                % opposed to +/- 1)
                
                % If the length of variables in the optimization problem is twice the
                % number of parameters, we are including both lin(x) and log10(x) in the
                % optimization variables.
                if bestVar > length(pList)
                    % We must have decided to subdivide on the log10 transformation
                    % version of the parameter.
                    bestVar = bestVar-length(pList);
                    
                    rng = log10(PDset.PiecewiseSurrogateModelTree(leaf).domainRange(bestVar,:));
                    if bestLoc == -1
                        loc = rng(1)+0.33*(rng(2)-rng(1));
                    elseif bestLoc == 0
                        loc = mean(rng);
                    elseif bestLoc == 1
                        loc = rng(1)+0.66*(rng(2)-rng(1));
                    else
                        disp('Im confused: condition should never occur')
                        keyboard
                    end
                    loc = 10^loc;
                    
                else
                    % We must have decided to subdivide on the linear transformation
                    % version of the parameter.
                    rng = PDset.PiecewiseSurrogateModelTree(leaf).domainRange(bestVar,:);
                    if bestLoc == -1
                        loc = rng(1)+0.33*(rng(2)-rng(1));
                    elseif bestLoc == 0
                        loc = mean(rng);
                    elseif bestLoc == 1
                        loc = rng(1)+0.66*(rng(2)-rng(1));
                    else
                        disp('Im confused: condition should never occur')
                        keyboard
                    end
                end
                
                % Determine the parameter to subdivide.
                dim = pList{bestVar};
                inherit = zeros(PDset.PiecewiseSurrogateModelTree(leaf).nSurfaces,1);
                
                
            else
                
                % ALGORITHM 2
                
                % eval slack of metamodels.
                %  y = evalResponseModels(PDset,pt);
                MetamodelCell = PDset.PiecewiseSurrogateModelTree(1).DCMetamodelCell;
                m = length(MetamodelCell);
                y = zeros(m,1);
                yL = zeros(m,1);
                yR = zeros(m,1);
                % yL <= M(x) <= yR (M(x) is true model)
                for i1 = 1:m 
                    y(i1) = MetamodelCell{i1}.eval(pt);
                    bS = MetamodelCell{i1}.fitInfo.peakErrorFromOptimization;
                    if strcmp(MetamodelCell{i1}.responseTransformation,'log10')
                        yL(i1) = 10^(log10(y(i1))-bS(2));
                        yR(i1) = 10^(log10(y(i1))-bS(1));
                    else
                        yL(i1) = y(i1)-bS(2);
                        yR(i1) = y(i1)-bS(1);
                    end
                end
                
                d = vertcat(PDset.ModelAndObservationPair.observedValue);
                uVect = vertcat(PDset.ModelAndObservationPair.observationUncertaintyPlusMinus);
                uCase = vertcat(PDset.ModelAndObservationPair.uncertaintyCase);
                gamma = zeros(m,1);
                
                % Assume pt is the global solution, and compute the upper bound based on
                % just the fitting error
                for i2 = 1:m 
                    switch uCase(i2)
                        case 1
                            gamma(i2) = min( (d(i2) + uVect(i2,2) - yL(i2))/uVect(i2,2) , (d(i2) + uVect(i2,1) - yR(i2))/uVect(i2,1) );
                        case 2
                            gamma(i2) = min( (d(i2)*(1 + uVect(i2,2)) - yL(i2))/(d(i1)*uVect(i2,2)) , (d(i2)*(1 + uVect(i2,1)) - yR(i2))/(d(i2)*uVect(i2,1)) );
                        case 3
                            % eta(x)<=yR, yR <= d10^u(1-gamma) ==>  eta(x) <= d10^u(1-gamma),
                            % yL <= eta(x), d10^el(1-gamma) <= yL ==> d10^el(1-gamma) <= eta(x)
                            gamma(i2) = min( (log10(d(i2)) - log10(yL(i2)) + uVect(i2,2))/uVect(i2,2) , (log10(d(i2)) - log10(yR(i2)) + uVect(i2,1))/uVect(i2,1) );
                        case 4
                            % eta(x)<=yR, yR <= d^(1+u(1-gamma)) ==>  eta(x) <= d^(1+u(1-gamma)),
                            % yL <= eta(x), d^(1+el(1-gamma)) <= yL ==> d^(1+el(1-gamma)) <= eta(x)
                            gamma(i2) = min( (log10(d(i2))*(1 + uVect(i2,2)) - log10(yL(i2)))/(log10(d(i2))*uVect(i2,2)) , (log10(d(i2))*(1 + uVect(i2,1)) - log10(yR(i2)))/(log10(d(i2))*uVect(i2,1)) );
                            
                        otherwise
                            error('Code Not Complete')
                    end
                end
                gammaOpt = min(gamma);
                
                % find the constraint slack at the optimal gamma
                slack = zeros(m,2);
                for i2 = 1:m 
                    switch uCase(i2)
                        case 1
                            % s1 + el(1-gamma)<= y-d <= -s2 + u(1-gamma)
                            slack(i2,1) = yL(i2)-d(i2)-uVect(i2,1)*(1-gammaOpt);
                            slack(i2,2) = -yR(i2)+d(i2)+uVect(i2,2)*(1-gammaOpt);
                        otherwise
                            error('Code Not Complete')
                    end
                end
                
                slack = min(slack,[],2);
                
                %UtilOne = (max(gamma)-gamma)/(max(gamma)-min(gamma));
                
                % from the upper bound over the worst subdomain, fitting errors are
                % important based on the multipliers?
                
                SurfArray = surfaces(PDset.PiecewiseSurrogateModelTree,leaf);
                MetamodelCell = PDset.PiecewiseSurrogateModelTree(1).DCMetamodelCell;
                
                surfs = struct(SurfArray);
                metaErr = zeros(m,2);
                surfErr = zeros(m,2);
                RMCell = cell(m,1);
                for i2 = 1:m 
                    surfErr(i2,:) = surfs(i2).surrogateFitInfo.peakErrorFromOptimization;
                    if ~isempty(MetamodelCell{i2})
                        metaErr(i2,:) = MetamodelCell{i2}.fitInfo.peakErrorFromOptimization;
                    end
                    RMCell{i2} = PDset.ModelAndObservationPair(i2).ResponseModel;
                end
                
                %what is the fitting error due to the quadratics.
                quadErr = surfErr-metaErr;
                
                % We consider to types of "Utility" is splitting at a particular
                % dimension. One is based on constrant slack and the other uses the
                % multipliers.  Kappa is a weighting factor for combining these
                % utilities.
                
                %fitting error is costly if ...
                UtilOne = max((diff(quadErr,[],2)-slack),0)./diff(uVect,[],2);
                
                UtilOne = UtilOne/max(UtilOne);
                
                UtilTwo = mults.*diff(quadErr,[],2);
                UtilTwo = UtilTwo/max(UtilTwo);
                
                kappa = 0.5;
                
                Util = UtilOne+kappa*UtilTwo;
                Util = Util./max(Util);
                
                %keyboard
                
                domain = DClab.createDomainStructure(PDset.parameterList,PDset.PiecewiseSurrogateModelTree(leaf).domainRange);
                trainingRange = vertcat(PDset.ModelAndObservationPair.trainingRange);
                [dim,loc] = findSubdivideLocation(SurfArray,MetamodelCell,RMCell,trainingRange,domain,Util);
                
                %inherit = ~todivide;
                inherit = Util < 0.2;
                
            end
            
            
            
            function [dim,loc, todivide] = findSubdivideLocation(SurfArray,MetamodelCell,RMCell,trainingRange,domain,impact)
                
                % we assume expErr and MetamodelCell, RMCell, trainingRange are in the same order.
                
                %At some pt we need to estimate expected err improvement in dcsurface??
                
                %find current avg err for each surf
                
                %nSurfs = nSurfaces(SurfArray);
                
                %ToDo: use methods
                SurfArray = struct(SurfArray);
                surf2meta=vertcat(SurfArray.MOPairIdx);
                
                %impact = fitting error/experimental error
                
                % impact = zeros(nSurfs,1);
                % for i1 = 1:nSurfs
                %   impact(i1) = diff(SurfArray(i1).surrogateFitInfo.peakErrorFromOptimization)/expErr(surf2meta(i1));
                % end
                
                toUse=find(impact>= 3*max(impact)/4);
                nToUse=length(toUse);
                
                %TODO make sure impact comes in the same size as SurfArray.
                
                todivide = false(size(impact));
                todivide(toUse) = true;
                
                linAct = cell(nToUse,1);
                logAct = cell(nToUse,1);
                for j1 = 1:nToUse
                    linAct{toUse(i1)} = SurfArray(toUse(j1)).activeParameterIndex(SurfArray(toUse(j1)).activeParameterTransformation==1);
                    logAct{toUse(i1)} = SurfArray(toUse(j1)).activeParameterIndex(SurfArray(toUse(j1)).activeParameterTransformation==2);
                end
                allLinAct = unique(vertcat(linAct{:}));
                allLogAct = unique(vertcat(logAct{:}));
                allAct=unique([allLinAct; allLogAct]);
                nact=length(allAct);
                
                PList = {domain.name}';
                bnds = vertcat(domain.range);
                %actPList={domain(allAct).name}';
                %actDomain=domain(allAct);
                
                % algorithm parameters:
                %   N = 100*(n+1)*(n+2)/2; % number of points over the current domain that are used to determine the splitting location.
                Ncuts = 8; % number of trial cuts in each dimensions used to determine the optimal cut.
                N = 150*(nact+1)*(nact+2)/2;
                
                n = length(domain);
                randdes=DClab.lhsdesign(N,n,'criterion','none')';
                xdes=zeros(n,N);
                
                cuts = zeros(n,Ncuts);
                for j1=1:length(allLinAct)
                    xdes(allLinAct(j1),:)=randdes(allLinAct(j1),:)*diff(bnds(allLinAct(j1),:))+bnds(allLinAct(j1),1);
                    tmp = linspace(bnds(allLinAct(j1),1),bnds(allLinAct(j1),2),Ncuts+2);
                    tmp(1) = [];
                    tmp(end) = [];
                    cuts(allLinAct(j1),:) = tmp;
                end
                for j1=1:length(allLogAct)
                    if xdes(allLogAct(j1),1)==0 && xdes(allLogAct(j1),2)==0
                        %impossible for it to also have been lin
                        xdes(allLogAct(j1),:)=10.^( randdes(allLogAct(j1),:)*log10(bnds(allLogAct(j1),2)/bnds(allLogAct(j1),1))+log10(bnds(allLogAct(j1),1)) );
                        tmp = logspace(log10(bnds(allLogAct(j1),1)),log10(bnds(allLogAct(j1),2)),Ncuts+2);
                        tmp(1) = [];
                        tmp(end) = [];
                        cuts(allLogAct(j1),:) = tmp;
                    else
                        %overwrite with half distributed in log spaceimpossible for it to also have been lin
                        rprm=randperm(N);
                        rprm=rprm(1:round(N/2));
                        xdes(allLogAct(j1),rprm)=10.^( randdes(allLogAct(j1),rprm)*log10(bnds(allLogAct(j1),2)/bnds(allLogAct(j1),1))+log10(bnds(allLogAct(j1),1)) );
                        
                        % We obtain N points from a linspace, and add in the 2nd to the middle from
                        % a logspace, giving around the original N total.
                        N_ = ceil(0.66*Ncuts);
                        
                        tmpln = linspace(bnds(allLogAct(j1),1),bnds(allLogAct(j1),2),N_+2);
                        tmplg = logspace(log10(bnds(allLogAct(j1),1)),log10(bnds(allLogAct(j1),2)),N_+2);
                        tmp = sort([tmpln tmplg(2:round(N_/2))]);
                        
                        tmp(1) = [];
                        tmp(end) = [];
                        if length(tmp)~=Ncuts;
                            disp('NcutsProb')
                            keyboard
                        end
                        cuts(allLogAct(i1),:) = tmp;
                    end
                end
                
                
                
                xT=cell(nToUse,1);
                yT=cell(nToUse,1);
                
                %get training
                trnRng = zeros(nToUse,2);
                
                for j1=1:nToUse 
                    MM=MetamodelCell{surf2meta(toUse(j1))};
                    trnRng(j1,:) = trainingRange(surf2meta(toUse(j1)),:);
                    
                    if isempty(MM)
                        yT{j1} = RMCell{surfmeta(toUse(j1))}.eval(xdes,PList);
                        if strcmp(SurfArray(toUse),'log10')
                            yT{j1} = log10(yT{j1});
                            trnRng(j1,:) = log10(trnRng(j1,:));
                        end
                    else
                        yT{j1} = MM.eval(xdes,PList);
                        if strcmp(MM.responseTransformation,'log10')
                            yT{j1} = log10(yT{j1});
                            trnRng(j1,:) = log10(trnRng(j1,:));
                        end
                    end
                    tmpX = xdes;
                    la = linAct{toUse(j1)};
                    ga = logAct{toUse(j1)};
                    tmpX(la,:) = 2*diag(1./diff(bnds(la,:),[],2))*(xdes(la,:)-repmat(bnds(la,1),1,N)) - 1;
                    tmpX(ga,:) = 2*diag(1./diff(log10(bnds(ga,:)),[],2))*(log10(xdes(ga,:))-repmat(log10(bnds(ga,1)),1,N)) - 1;
                    junk = setdiff((1:n)',[la;ga]);
                    tmpX(junk,:) = [];
                    xT{i1} = tmpX;
                    
                end
                
                % Determine the performance that is possible if we
                
                % Split each dimension in half and fit the residual.
                %MX = DClab.buildQuadX(Xtest');
                %coeffA = MX\yhat';
                %possibleresid = (MX*coeffA)' - yhat;
                %stdposs = std(possibleresid);
                %swingposs = max(possibleresid)-min(possibleresid);
                
                
                
                %  % Only normalize once. The problem shouldn't get too ill-conditioned if
                %  % we don't renormalize over each trial partion.
                %  bnds = tree(worstleaf).bnds;
                %  normx = 2*(x-repmat(bnds(:,1),1,N))'*diag(1./diff(bnds,[],2))-1;
                
                % Cycle thru the dimensions.
                
                stdnew = zeros(nact,Ncuts);
                extentnew = zeros(nact,Ncuts);
                for j1 = 1:nact
                    for j2 = 1:Ncuts
                        
                        % Partition the x data into the two trial subdomains.
                        idxL = xdes(allAct(j1),:) <= cuts(allAct(j1),j2);
                        idxR = xdes(j1,:) > cuts(allAct(j1),j2);
                        
                        for j3 = 1:nToUse 
                            
                            % Fit with quadratic.
                            
                            yL = yT{j3}(idxL)';
                            xL = xT{j3}(:,idxL)';
                            junkIdx = yL <= trnRng(j3,1) | yL >= trnRng(j3,2);
                            yL = yL(~junkIdx);
                            xL = xL(~junkIdx,:);
                            
                            if length(yL) > (size(xL,2)+1)*(size(xL,2)+2)/2;
                                [coeffL errL] = DClab.twoNormRegression(DClab.buildQuadX(xL),yL);
                            else
                                errL = 0;
                            end
                            
                            yR = yT{j3}(idxR)';
                            xR = xT{j3}(:,idxR)';
                            junkIdx = yR <= trnRng(j3,1) | yR >= trnRng(j3,2);
                            yR = yR(~junkIdx);
                            xR = xR(~junkIdx,:);
                            
                            if length(yR) > (size(xR,2)+1)*(size(xR,2)+2)/2;
                                [coeffR errR] = DClab.twoNormRegression(DClab.buildQuadX(xR),yR);
                            else
                                errR = 0;
                            end
                            
                            newerr = [errL; errR];
                            stdnew(j1,j2) = stdnew(j1,j2) + ( impact(toUse(j3))*std(newerr) )^2; %minimize the two norm
                            extentnew(j1,j2) = extentnew(j1,j2) + ( impact(toUse(j3))*(max(newerr)-min(newerr)) )^2;
                        end
                        %      newerr = [errL; errR];
                        %      stdnew(i1,i2) = std(newerr);
                        %      extentnew(i1,i2) = max(newerr)-min(newerr);
                    end
                end
                
                [loweststd whichCut] = min(stdnew,[],2);
                [crap whichDim] = min(loweststd);
                
                loc = cuts(allAct(whichDim),whichCut(whichDim));
                dim = PList{allAct(whichDim)};
                
            end
        end
        
        function [s1 s2] = maximizeLB(PD,opt,nodes,computeGap,startPt)
            %MAXIMIZELB computes an lower bound on the consistency measure of PD
            %
            %   S1 = MAXIMIZELB(PD) determines an upper bound on the consistency measure
            %   of the dataset PD.
            %
            %   S1 = MAXIMIZELB(PD,OPT) allows you to supply options for the
            %   optimization with the DClab.DCOptions object OPT.
            %
            %   S1 = MAXIMIZELB(PD,OPT,NODES) allows you specify which leaf
            %   nodes to optimize over. The default is NODES = PDSET.leafNodes.
            %
            %   S1 = MAXIMIZELB(PD,OPT,NODES,COMPUTEGAP) will perform a local
            %   search for the maximum so you may estimate how well the SDP relaxation
            %   of the quadratic problem did.
            %
            %   S1 = MAXIMIZELB(PD,OPT,NODES,COMPUTEGAP,STARTPNT) allows you to
            %   supply an initial seed for the local search. STARTPNT must be an rx1
            %   cell array, where r = length(NODES).
            %
            %   [S1 S2] = MAXIMIZELB(...) returns a second output that contains much
            %   more information from the optimization(s)
            %
            %   Inputs:
            %      PD: A PolyDataset object describing the constraints.
            %      OPT[optional]: A DCOptions object.  Used fields are .tolFun,
            %         .tolCon, .sedumiParEps, and .display
            %      NODES[optional]: rx1 vector of indicies indicating the subdivisions
            %         of the domain over which to perform the optimization. Must be a
            %         subset of PSMTREE.leafNodes.
            %      COMPUTEGAP[optional]: A boolean. If COMPUTEGAP==1, any nqcqps used to
            %         lower bound the minimum of Mo will also be subjected to a local
            %         search. Except in the special case where the true models are
            %         algebraic, it is not clear this is useful. See below.
            %      STARTPNT[optional]: A rx1 cell array of nx1 vectors of starting
            %         points for the local searchs (applicable if computeGap=true). The
            %         order should correspond to the order of the parameters in
            %         PSMTREE.
            
            %   OLD HELP BELOW, possibly inrelevant help below
            %
            % Inputs:
            %   MnotPD: A PolyDataset object for the predictor model. Currently,
            %     this must be a length one object and its containing model assertion
            %     must have length(featureList)<=1;
            %   PD: A PolyDataset object describing the contraints.
            %   opt[optional]: A DCOptions object.  Used fields are .tolFun, .tolCon,
            %     .sedumiParEps, .constraints, and .display
            %   nodeIndices[optional]: rx1 vector of indicies into the surrogateTree of the
            %     merged PolyDataset mPD = merge(MnotPD,PD).
            %   computeGap[optional]: A boolean. If computeGap = 1, any nqcqps used to
            %     lower bound the minimum of Mo will also be subjected to a local
            %     search. Except in the special case where the true models are
            %     algebraic, it is not clear this is useful. See below.
            %   startPt[optional]: A rx1 cell array of nx1 vectors of starting points
            %     for the local searchs (applicable if computeGap=true). The order
            %     should correspond to the order of the parameters in the merged
            %     PolyDataset mPD = merge(MnotPD,PD).
            %
            % Outputs:
            %   s1: structure summarizing the results, further detail below
            %   s2: structure containing all results, further detail below
            %
            % Description
            %
            % The purpose of this function is to obtain an upper bound on the solution
            % of the mathematical program
            %
            % min Mo(x)   subj. to: |M_e(x)-d_e| \leq u_e*(1-\gamma)           (eq1)
            %                        x in the nodeIndices subdiv of H
            %
            % where Mo is the ModelAssertion of the PolyDataset MnotPD, M_e are the
            % true models present in the ModelAssertions of the PolyDataset PD and d_e,
            % u_e are the data and uncertainty present in the corresponding
            % ExperimentAssertions. (actually it is not quite that simple since the
            % uncertainty u_e can be asymmetric and the models Mo,M_e may have output
            % uncertainty).
            %
            % Depending on the surrogate models that are available, up to 8
            % approximations to (eq1) are solved, (only up to four if no fitting error)
            % each which takes the form of a nonconvex quadratically constrained
            % quadratic program (nqcqp). If computeGap = true, the S-procedure is used
            % to determine a global lower bound on the optimum for each nqcqp. If the
            % gap between this lower bound and the solution found by local search is
            % zero, the local search found the global minimum of the nqcqp. However,
            % since each nqcqp is just an approximation to eq1 (except in the special
            % case where M_e is quadratic and has no output uncertainty), it is not
            % clear that the gap being zero has any practical meaning. Consequently,
            % computeGap = false is the default behavior.
            %
            % Each method to approximate (eq1) determines a 'good' x at which to
            % evaluate the objective function, producing an upper bound on the optimum.
            % Finding such good x's is difficult. The first 4 methods consider
            % optimizations over a quadratically described inner approximation of the
            % feasible set realized by subtracting off the fitting error. Fit errors
            % are often large, causing this inner approximation to be empty. In this
            % case, inf used as an upper bound on the optimum. These four methods
            % differ only in the transformations used to construct the surrogate
            % models. The nqcqp's solved are tersely presented next.
            %
            % method1. here we have surrogate models S_e that approximate M_e such that
            %       fiterr_e(1)+M_e(x) <= S_e(x) <= fiterr_e(2)+M_e(x)
            %   and possibly
            %       fiterr'_e(1)+log10(M_e(x)) <= S'_e(x) <= fiterr'_e(2)+log10(M_e(x)) (prime denotes model in logY)
            %
            %   we solve the nqcqp
            %   min S_o(x)  subj. to:
            %      d_e-u_e + fiterr_e(2) <= S_e(x) <= d_e + u_e +  fiterr_e(1)
            %      log10(d_e-u_e) + fiterr'_e(2) <= S'_e(x) <= log10(d_e+u_e) + fiterr'_e(1)  (if linXlogY surrogates exist)
            %      x in the nodeIndices subdiv of H
            %
            % method2. here we have surrogate models S_e that approximate M_e such that
            %       fiterr_e(1)+M_e(x) <= S_e(log10(x)) <= fiterr_e(2)+M_e(x)
            %   and possibly
            %       fiterr'_e(1)+log10(M_e(x)) <= S'_e(log10(x)) <= fiterr'_e(2)+log10(M_e(x)) (prime denotes model in logY)
            %
            %   we solve the nqcqp
            %   min S_o(x')  subj. to:
            %      d_e-u_e + fiterr_e(2) <= S_e(x') <= d_e + u_e +  fiterr_e(1)
            %      log10(d_e-u_e) + fiterr'_e(2) <= S'_e(x') <= log10(d_e+u_e) + fiterr'_e(1)  (if linXlogY surrogates exist)
            %      where x' = log10(x) and x in the nodeIndices subdiv of H
            %
            % method3. here we have surrogate models S_e that approximate M_e such that
            %       fiterr'_e(1)+log10(M_e(x)) <= S'_e(x) <= fiterr'_e(2)+log10(M_e(x))  (prime denotes model in logY)
            %   and possibly
            %       fiterr_e(1)+M_e(x) <= S_e(x) <= fiterr_e(2)+M_e(x)
            %
            %   we solve the nqcqp
            %   min S'_o(x)  subj. to:
            %      log10(d_e-u_e) + fiterr'_e(2) <= S'_e(x) <= log10(d_e+u_e) + fiterr'_e(1)
            %      d_e-u_e + fiterr_e(2) <= S_e(x) <= d_e + u_e +  fiterr_e(1) (if linXlinY surrogates exist)
            %      x in the nodeIndices subdiv of H
            %
            % method4. here we have surrogate models S_e that approximate M_e such that
            %       fiterr'_e(1)+log10(M_e(x)) <= S'_e(log10(x)) <= fiterr'_e(2)+log10(M_e(x)) (prime denotes model in logY)
            %   and possibly
            %       fiterr_e(1)+M_e(x) <= S_e(log10(x)) <= fiterr_e(2)+M_e(x)
            %
            %   we solve the nqcqp
            %   min S'_o(x')  subj. to:
            %      log10(d_e-u_e) + fiterr'_e(2) <= S'_e(x') <= log10(d_e+u_e) + fiterr'_e(1)
            %      d_e-u_e + fiterr_e(2) <= S_e(x') <= d_e + u_e +  fiterr_e(1) (if linXlinY surrogates exist)
            %      where x' = log10(x) and x in the nodeIndices subdiv of H
            %
            % The next 4 methods again consider the above four optimizations with the
            % exception that the fitting error is ignored.  This results in an
            % optimization over a quadratically described crude approximation of the
            % feasible set. As before, the four methods are identical with the
            % exception of the surrogate model transformations. These methods suffer
            % from the problem that the optimizer determined over this crude
            % approximation may not be feasible for the actual problem. The approach
            % taken is now described. We begin by finding a point x that is in the
            % 'center' of the crude approximation be minimizing a constraint slack
            % variable. This minimum of this slack variable (call it delta) will be
            % negative if the crude approximation is nonempty, and hopefully the x
            % found is in the true feasible set. Then perform we minimize the quadratic
            % objective subject to
            %
            % -0.95*delta -d_e +u_e(1) <= S_e(x) <= d_e + u_e(2) + 0.95*delta
            %
            % -0.8*delta -d_e +u_e(1) <= S_e(x) <= d_e + u_e(2) + 0.8*delta
            %
            % -0.6*delta -d_e +u_e(1) <= S_e(x) <= d_e + u_e(2) + 0.6*delta
            %
            % -0.4*delta -d_e +u_e(1) <= S_e(x) <= d_e + u_e(2) + 0.4*delta
            %
            % -0.2*delta -d_e +u_e(1) <= S_e(x) <= d_e + u_e(2) + 0.2*delta
            %
            % TODO for this to work as desired, the u_e's must be about the same
            % magnitude. Someday we should place additional info in the nqcqp's to make
            % it easy to rescale so this occurs.
            %
            % This creates a short sequence of x's each somewhat improving the value of
            % the objective function. Likely only some of these points will be feasible
            % for the constraints of (eq1). The upper bounds on the optimum of eq1 from
            % these last four methods are obtained by evaluating the true objective
            % function at each FEASIBLE x in this short sequence.
            %
            % Note, if the objective in in linXlinY, constraints in both linXlinY and
            % linXlogY are used if available. Other cases are similar. Essentially we
            % used as many available constraints as possible.
            %
            % Description of the structure array s2
            %
            % s2.m1.star.bnd contains an rx1 vector of lower bounds from the
            %                s-proc (if computeGap = true), where r=length(nodeIndices)
            % s2.m1.star.mult contains an rx1 cell array of vectors containing
            %                 the corresponding multipliers from the s-proc
            %                 optimization. It is not clear that these have any
            %                 meaning. For consistency with the object PD,
            %                 postprocessing occurs so that if a parameter of PD was
            %                 not a decision variable in the optimization (due to not
            %                 being an active parameter) a "multipler" is still
            %                 returned corresponding to that parameter. It is just set
            %                 to zero.
            % s2.m1.pound.bnd contains an rxk matrix of upper bounds to eq2. Each
            %                 component is determined by evaluating the objective in
            %                 eq1 at the local search solutions of the nqcqps. If the
            %                 local search solution is infeasible for eq1, the bound is
            %                 inf. Here k = 1+opt.nRestart.
            % s2.m1.pound.xfeas contains an rxk cell array of parameter vectors
            %                 corresponding to the local minimum found by local
            %                 search. Post-processing occurs so that these values are
            %                 nontransformed, i.e., they are in whatever units the
            %                 original parameter assertions were provided in. If a
            %                 parameter of PD was not a decision variable in the
            %                 optimization (due to not being an active parameter), the
            %                 "nominal value" of that parameter is inserted in the
            %                 proper location in s2.m1.pound.xfeas.
            % s2.m1.pound.mult is similar to ...star.mult, it is just the
            %                 multipliers from the local search. It is not
            %                 clear that these are useful in the nonquadratic
            %                 case, because they have nothing to do with the
            %                 optimization we want to solve (eq1).
            %
            % s2.m2 through s2.m8 contain the above info from the corresponding methods.
            %
            % Description of the structure array s1
            %
            % s1.pound.bnd    contains an rx1 vector of the best overall (from all
            %                 the restarts and all the methods) upper bound on the
            %                 solution to eq1 for each of the r subdivisions.
            % s1.pound.method contains an rx1 vector of elements in {1,...,8} reporting
            %                 which method produced the best answer
            % s1.pound.mult   contains an rx1 cell of the multipliers
            %                 corresponding to the best overall optimizations. As in
            %                 the above, it is not clear this are meaningful in the
            %                 nonquadratic case
            % s1.pound.xfeas  contains the parameter vector that produced the
            %                 best bound in s1.pound.bound
            % s1.pound.gap if computeGap=true, this rx1 vector contains the
            %                 difference between s1.pound.bnd and lowerbound
            %                 solutions of the ncqcps that produced the bestx's in
            %                 local search. In the nonquad case this seems
            %                 meaningless...
            %
            % See also PolyDataset, Prediction, PolyDataset/olbPred
            
            
            %TODO fix help
            
            %TODO, since we evaluate the true models, the multipliers and gap are
            %essentially meaningless. Should we even return them?
            
            
            %TODO: update the help.
            
            ni = nargin;
            switch ni
                case 1
                    opt = [];
                    nodes = [];
                    computeGap = 0;
                    startPt = [];
                case 2
                    nodes = [];
                    computeGap = 0;
                    startPt = [];
                case 3
                    computeGap = 0;
                    startPt = [];
                case 4
                    startPt = [];
                otherwise
                    error(nargchk(1,5,ni))
            end
            
            %===input checking/empty input initialization
            %TODO, since this is a private method, we can probably do away with most of
            %the input checks once it has been tested.
            
            leafNodes = PD.leafNodes;
            if isempty(nodes)
                nodes = leafNodes;
                r = length(nodes);
            else
                r = length(nodes);
                if ~isequal(size(nodes),[r 1]) || ~isnumeric(nodes)
                    error('Inputs: NODES must be a column vector')
                end
                if ~isempty(setdiff(nodes,leafNodes))
                    error('Inputs: NODES must be a subset of PDSet.leafNodes.')
                end
            end
            
            if isempty(opt)
                opt = DClab.DCOptions;
            end
            if isempty(computeGap)
                computeGap = false;
            end
            if ~isempty(startPt) && ~iscell(startPt) && ~isequal(size(startPt),[r 1])
                error('STARTPNT must be an rx1 cell array, where r is the number of subdivisions being optimized over.')
            end
            %===end input checking/empty input initialization
            
            % Unless there is no fitting error, we will solve using two methods. This
            % first solves on optimization over the inner approximation of the feasible
            % set (which may be empty). The second uses a crude approximation of the
            % feasible set to find points that are hopefully feasible and hopefully
            % give small costs.
            
            % Useful constants
            nNodes = length(nodes);
            n = PD.nParameters;
            nTries = opt.nRestart+1;
            
            %Initialize some variables.
            s2.m1.pound.bnd = cell(nNodes,1);
            s2.m1.pound.xfeas = cell(nNodes,1);
            s2.m1.pound.pairMult = cell(nNodes,1);
            s2.m1.pound.paramMult = cell(nNodes,1);
            
            s2.m2.pound.bnd = cell(nNodes,1);
            s2.m2.pound.xfeas = cell(nNodes,1);
            s2.m2.pound.pairMult = cell(nNodes,1);
            s2.m2.pound.paramMult = cell(nNodes,1);
            
            s2.m3.pound.bnd = cell(nNodes,1);
            s2.m3.pound.xfeas = cell(nNodes,1);
            s2.m3.pound.pairMult = cell(nNodes,1);
            s2.m3.pound.paramMult = cell(nNodes,1);
            
            s1.pound.bnd = zeros(nNodes,1);
            s1.pound.xfeas = cell(nNodes,1);
            s1.pound.pairMult = cell(nNodes,1);
            s1.pound.paramMult = cell(nNodes,1);
            
            if computeGap
                s2.m1.star.bnd = zeros(nNodes,1);
                s2.m1.star.pairMult = cell(nNodes,1);
                s2.m1.star.paramMult = cell(nNodes,1);
                
                s1.star.bnd = zeros(nNodes,1);
                s1.star.pairMult = cell(nNodes,1);
                s1.star.paramMult = cell(nNodes,1);
                
                s1.gap = zeros(nNodes,1);
            end
            
            %analMode = opt.analysisMode;
            
            for i1 = 1:nNodes
                
                %if strcmp(analMode,'original') || strcmp(analMode,'metamodelBasedA')
                
                
                [optProb.m1 activeIdx] = DClab.ConsistencyTest.makeCmoptim(PD,nodes(i1),'inner');
                optProb.m2 = DClab.ConsistencyTest.makeCmoptim(PD,nodes(i1),'crude');
                
                if computeGap
                    [bnd,mult,Q] = lowerBnd(optProb.m1,opt);
                    s2.m1.star.bnd(i1) = -bnd;
                    
                    %Determine if we need to include multiplier for logX variables.
                    if ~isempty(mult)
                        if max(activeIdx) > n
                            pmult = zeros(2*n,2);
                        else
                            pmult = zeros(n,2);
                        end
                        pmult(activeIdx,:) = [mult.lower mult.upper];
                        emult = [mult.quadineq(1:2:end-1) mult.quadineq(2:2:end)];
                        s2.m1.star.pairMult{i1} = emult;
                        s2.m1.star.paramMult{i1} = pmult;
                    end
                end
                
                %TODO: the logY multipliers will be screwed up...need to postprocess some
                %day
                if ~exist('Q','var') %see if we can use the s-proc results for a seed point
                    Q = [];
                end
                
                for i3 = 1:nTries
                    %if ~isempty(startPt{i1}), findStartPt will massage startPt{i1} so
                    %it has the correct dimensions for the optimization. The second
                    %time through the start point will be based on Q. All subsequent
                    %times it will be randomly generated.
                    %
                    %if isempty(startPt{i1}), The first time through a start point will
                    %be based on Q; all subsequent times it will be randomly generated.
                    
                    if isempty(startPt) || isempty(startPt{i1})
                        if i3 == 1
                            xinit = DClab.findStartPt(PD,nodes(i1),activeIdx,[],Q);
                        else
                            xinit = DClab.findStartPt(PD,nodes(i1),activeIdx,[],[]);
                        end
                    else
                        if i3 == 1
                            xinit = DClab.findStartPt(PD,nodes(i1),activeIdx,startPt{i1},Q);
                        elseif i3 == 2
                            xinit = DClab.findStartPt(PD,nodes(i1),activeIdx,[],Q);
                        else
                            xinit = DClab.findStartPt(PD,nodes(i1),activeIdx,[],[]);
                        end
                    end
                    
                    %refresh gui
                    if strcmpi(opt.display,'all')
                        str = [blanks(4) 'method 1: fmincon local search try ' num2str(i3) ' of ' num2str(nTries)];
                        DClab.dcdispstr(str,opt.guiHandle,false)
                    else
                        drawnow
                    end
                    
                    [fval,xfeas,mult] = upperBnd(optProb.m1,xinit,opt);
                    
                    xfeas =  DClab.transformAndAddNominals(xfeas,PD,activeIdx,nodes(i1));
                    
                    % TODO: if xfeas is nonempty, it will be in the true feasible set
                    % provided our error estimates are ok. Should we check to make sure?
                    % Calling isFeasible(PD,xfeas{i3}) would do the trick.
                    
                    if isempty(xfeas)
                        s2.m1.pound.bnd{i1}(1,i3) = -inf; %set to inf if no feasible point found
                    else
                        s2.m1.pound.bnd{i1}(1,i3) = evalCMeasObj(PD,xfeas);
                    end
                    
                    s2.m1.pound.xfeas{i1}{1,i3} =  xfeas;
                    
                    %Determine if we need to include multiplier for logX variables.
                    if max(activeIdx) > n
                        pmult = zeros(2*n,2);
                    else
                        pmult = zeros(n,2);
                    end
                    pmult(activeIdx,:) = [mult.lower mult.upper];
                    emult = [mult.quadineq(1:2:end-1) mult.quadineq(2:2:end)];
                    s2.m1.pound.pairMult{i1}{1,i3} = emult;
                    s2.m1.pound.paramMult{i1}{1,i3} = pmult;
                end %for i3 = 1:nTries
                
                %TODO is isequal a sufficiently fast way to check if there is any
                %surrogate modeling error?
                if ~isequal(optProb.m1,optProb.m2)
                    if computeGap
                        [bnd,mult,Q] = lowerBnd(optProb.m2,opt);
                        s2.m2.star.bnd(i1) = -bnd;
                        
                        %Determine if we need to include multiplier for logX variables.
                        if ~isempty(mult)
                            if max(activeIdx) > n
                                pmult = zeros(2*n,2);
                            else
                                pmult = zeros(n,2);
                            end
                            pmult(activeIdx,:) = [mult.lower mult.upper];
                            emult = [mult.quadineq(1:2:end-1) mult.quadineq(2:2:end)];
                            s2.m2.star.pairMult{i1} = emult;
                            s2.m2.star.paramMult{i1} = pmult;
                        end
                    end
                    
                    %TODO: the logY multipliers will be screwed up...need to postprocess some
                    %day
                    if ~exist('Q','var') %see if we can use the s-proc results for a seed point
                        Q = [];
                    end
                    
                    for i3 = 1:nTries
                        %if ~isempty(startPt{i1}), findStartPt will massage startPt{i1} so
                        %it has the correct dimensions for the optimization. The second
                        %time through the start point will be based on Q. All subsequent
                        %times it will be randomly generated.
                        %
                        %if isempty(startPt{i1}), The first time through a start point will
                        %be based on Q; all subsequent times it will be randomly generated.
                        
                        if isempty(startPt) || isempty(startPt{i1})
                            if i3 == 1
                                xinit = DClab.findStartPt(PD,nodes(i1),activeIdx,[],Q);
                            else
                                xinit = DClab.findStartPt(PD,nodes(i1),activeIdx,[],[]);
                            end
                        else
                            if i3 == 1
                                xinit = DClab.findStartPt(PD,nodes(i1),activeIdx,startPt{i1},Q);
                            elseif i3 == 2
                                xinit = DClab.findStartPt(PD,nodes(i1),activeIdx,[],Q);
                            else
                                xinit = DClab.findStartPt(PD,nodes(i1),activeIdx,[],[]);
                            end
                        end
                        
                        %refresh gui
                        if strcmpi(opt.display,'all')
                            str = [blanks(4) 'method 2: fmincon local search try ' num2str(i3) ' of ' num2str(nTries)];
                            DClab.dcdispstr(str,opt.guiHandle,false)
                        else
                            drawnow
                        end
                        
                        [fval,xfeas,mult] = upperBnd(optProb.m2,xinit,opt);
                        xfeas =  DClab.transformAndAddNominals(xfeas,PD,activeIdx,nodes(i1));
                        
                        % TODO: if xfeas is nonempty, it will be in the true feasible set
                        % provided our error estimates are ok. Should we check to make sure?
                        % Calling isFeasible(PD,xfeas{i3}) would do the trick.
                        
                        if isempty(xfeas)
                            s2.m2.pound.bnd{i1}(1,i3) = -inf; %set to inf if no feasible point found
                        else
                            %DO STUFF TO GET LOWER BOUND
                            s2.m2.pound.bnd{i1}(1,i3) = evalCMeasObj(PD,xfeas);
                        end
                        
                        s2.m2.pound.xfeas{i1}{1,i3} =  xfeas;
                        
                        %Determine if we need to include multiplier for logX variables.
                        if max(activeIdx) > n
                            pmult = zeros(2*n,2);
                        else
                            pmult = zeros(n,2);
                        end
                        pmult(activeIdx,:) = [mult.lower mult.upper];
                        emult = [mult.quadineq(1:2:end-1) mult.quadineq(2:2:end)];
                        s2.m2.pound.pairMult{i1}{1,i3} = emult;
                        s2.m2.pound.paramMult{i1}{1,i3} = pmult;
                    end %for i3 = 1:nTries
                    
                end %~isequal(optProb.m1,optProb.m2)
                
                %end
                
                %   if strcmp(analMode,'metamodelBasedA') || strcmp(analMode,'metamodelBasedB')
                %     % If there are metamodels. Lets try to optimize with them as well.
                %     MMCell = PD.PiecewiseSurrogateModelTree(1).DCMetamodelCell;
                %     if ~isempty(MMCell{1})
                %
                %       d = vertcat(PD.ModelAndObservationPair.observedValue);
                %       u = vertcat(PD.ModelAndObservationPair.observationUncertaintyPlusMinus);
                %       uncCase = vertcat(PD.ModelAndObservationPair.uncertaintyCase);
                %
                %       domRng = PD.PiecewiseSurrogateModelTree(nodes(i1)).domainRange;
                %       FP = PD.FreeParameter;
                %
                %       for i3 = 1:nTries
                %         [fval xfeas] = DClab.ConsistencyTest.maximizeLB_meta(MMCell,d,u,uncCase,FP,domRng);
                %
                %         if isempty(xfeas)
                %           s2.m3.pound.bnd{i1}(1,i3) = -inf; %set to inf if no feasible point found
                %         else
                %           %DO STUFF TO GET LOWER BOUND
                %           s2.m3.pound.bnd{i1}(1,i3) = evalCMeasObj(PD,xfeas);
                %         end
                %         s2.m3.pound.xfeas{i1}{1,i3} =  xfeas;
                %         s2.m3.pound.pairMult{i1}{1,i3} = [];
                %         s2.m3.pound.paramMult{i1}{1,i3} = [];
                %       end
                %     end
                %   end
                
                if ~isempty(s2.m1.pound.bnd{i1})
                    m1max = max(s2.m1.pound.bnd{i1});
                else
                    m1max = -inf;
                end
                if ~isempty(s2.m2.pound.bnd{i1})
                    m2max = max(s2.m2.pound.bnd{i1});
                else
                    m2max = -inf;
                end
                if ~isempty(s2.m3.pound.bnd{i1})
                    m3max = max(s2.m3.pound.bnd{i1});
                else
                    m3max = -inf;
                end
                
                if m1max >= m2max && m1max >= m3max
                    bestmethod = 'm1';
                elseif m2max >= m1max && m2max >= m3max
                    bestmethod = 'm2';
                else
                    bestmethod = 'm3';
                end
                if strcmpi(opt.display,'all')
                    str = [blanks(4) 'method ' bestmethod(2) ' gave the best lower bound in maximizeLB'];
                    DClab.dcdispstr(str,opt.guiHandle,false)
                end
                
                [rowmax best2] = max(s2.(bestmethod).pound.bnd{i1},[],2);
                [trash best1] = max(rowmax);
                
                s1.pound.bnd(i1) = s2.(bestmethod).pound.bnd{i1}(best1,best2);
                s1.pound.xfeas{i1} = s2.(bestmethod).pound.xfeas{i1}{best1,best2};
                s1.pound.pairMult{i1} = s2.(bestmethod).pound.pairMult{i1}{best1,best2};
                s1.pound.paramMult{i1} = s2.(bestmethod).pound.paramMult{i1}{best1,best2};
                
                if computeGap
                    s1.star.bnd(i1) = s2.(bestmethod).star.bnd(i1);
                    s1.star.pairMult{i1} = s2.(bestmethod).star.pairMult{i1};
                    s1.star.paramMult{i1} = s2.(bestmethod).star.paramMult{i1};
                    s1.gap(i1,1) = s1.pound.bnd(i1)-s1.star.bnd(i1);
                end
                
            end
        end
        function [s1 s2] = maximizeUB(PD,opt,nodes,computeGap,startPt)
            %MAXIMIZEUB computes an upper bound on the consistency measure of PD
            %
            %   S1 = MAXIMIZEUB(PD) determines an upper bound on the consistency measure
            %   of the dataset PD.
            %
            %   S1 = MAXIMIZEUB(PD,OPT) allows you to supply options for the
            %   optimization with the DCOptions object OPT.
            %
            %   S1 = MAXIMIZEUB(PD,OPT,NODES) allows you specify which leaf
            %   nodes to optimize over. The default is NODES = PDSET.leafNodes.
            %
            %   S1 = MAXIMIZEUB(PD,OPT,NODES,COMPUTEGAP) will perform a local
            %   search for the maximum so you may estimate how well the SDP relaxation
            %   of the quadratic problem did.
            %
            %   S1 = MAXIMIZEUB(PD,OPT,NODES,COMPUTEGAP,STARTPNT) allows you to
            %   supply an initial seed for the local search. STARTPNT must be an rx1
            %   cell array, where r = length(NODES).
            %
            %   [S1 S2] = MAXIMIZEUB(...) returns a second output that contains much
            %   more information from the optimization(s)
            %
            %   Inputs:
            %      PD: A PolyDataset object describing the constraints.
            %      OPT[optional]: A DCOptions object.  Used fields are .tolFun,
            %         .tolCon, .sedumiParEps, and .display
            %      NODES[optional]: rx1 vector of indicies indicating the subdivisions
            %         of the domain over which to perform the optimization. Must be a
            %         subset of PSMTREE.leafNodes.
            %      COMPUTEGAP[optional]: A boolean. If COMPUTEGAP==1, any nqcqps used to
            %         lower bound the minimum of Mo will also be subjected to a local
            %         search. Except in the special case where the true models are
            %         algebraic, it is not clear this is useful. See below.
            %      STARTPNT[optional]: A rx1 cell array of nx1 vectors of starting
            %         points for the local searchs (applicable if computeGap=true). The
            %         order should correspond to the order of the parameters in
            %         PSMTREE.
            %
            % Outputs:
            %   S1: structure summarizing the results, further detail below
            %   S2: structure containing all results, further detail below
            %
            % Description
            %
            % The purpose of this function is to obtain an upper bound on the solution
            % of the mathematical program
            %
            % max \gamma subj. to: |M_e(x)-d_e| \leq u_e*(1-\gamma) (eq1)
            %                        x in the nodeIndices subdiv of H
            %
            % where M_e are the true models present in the ModelAssertions of the
            % PolyDataset PD and d_e, u_e are the data and uncertainty present in the
            % corresponding ExperimentAssertions. (actually it is not quite that
            % simple since the uncertainty u_e can be asymmetric and the model M_e may
            % have output uncertainty).
            %
            % An approximation to the optimization problem (eq2) is considered that
            % takes the form of a nonconvex quadratically constrained quadratic program
            % (nqcqp). The s-procedure is used to lower bound the optimum of each
            % nqcqp. If computeGap = true, a local search is used to determine an upper
            % bound on the optimum for each nqcqp. If the gap between this upper bound
            % and the solution found by local search is zero, the local search found
            % the global minimum of the nqcqp. However, since each nqcqp is just an
            % approximation to eq2 (except in the special case where M_e is quadratic
            % and has no output uncertainty), it is not clear that the gap being zero
            % has any practical meaning. Consequently, computeGap = false is the
            % default behavior.
            %
            % If the surrogate models approximate y = M(x) and/or y = M(log10(x)), then
            % a nqcqp is constructed and solved in both the x and logx decision
            % variables. If the surrogate models approximate log10(y) = M(x) and/or
            % log10(y) = M(log10(x)), an iterative approach is taken in which a
            % sequence of nqcqp's are solved. See nqcqpe.
            %
            % See also maximizeLB, ConsistencyTest
            
            
            %TODO, where should we post-process the multipliers? Should this function
            %return multipliers for both linX and logX? What about when there are
            %multiple surrogates (with different I/0 transformations) for each
            %response model. Where do we pare this down?
            %
            %Since this function acts on PolyDatasets which have surrogate models, in
            %both linX and logX variables, the multipliers from this function should be
            %for everything.
            %
            %Lets return a ns-by-2 vector for the surrogates, with multipliers for the
            %upper and lower bound. Also a 2*n-by-2 vector for parameters, with
            %multipliers for the upper and lower bound.
            
            %TODO: it is a little weird the we want startPt to be n-by-1, but the
            %multipliers and any returned points from the local search are 2n-by-1.
            
            % Since a polydataset knows about all its surfaces, I think its ok to
            % return multipliers for the constraints imposed by each surfaces. However,
            % since it just has n parameters, maybe it would be ok to only return
            % multipliers for the parameter bounds. This will require massaging.
            
            %TODO: update the help.
            
            ni = nargin;
            switch ni
                case 1
                    opt = [];
                    nodes = [];
                    computeGap = 0;
                    startPt = [];
                case 2
                    nodes = [];
                    computeGap = 0;
                    startPt = [];
                case 3
                    computeGap = 0;
                    startPt = [];
                case 4
                    startPt = [];
                otherwise
                    error(nargchk(1,5,ni))
            end
            
            %===input checking/empty input initialization
            %TODO, since this is a private method, we can probably do away with most of
            %the input checks once it has been tested.
            
            leafNodes = PD.leafNodes;
            if isempty(nodes)
                nodes = leafNodes;
                r = length(nodes);
            else
                r = length(nodes);
                if ~isequal(size(nodes),[r 1]) || ~isnumeric(nodes)
                    error('Inputs: NODES must be a column vector')
                end
                if ~isempty(setdiff(nodes,leafNodes))
                    error('Inputs: NODES must be a subset of PDSet.leafNodes.')
                end
            end
            
            if isempty(opt)
                opt = DClab.DCOptions;
            end
            if isempty(computeGap)
                computeGap = false;
            end
            if ~isempty(startPt) && ~iscell(startPt) && ~isequal(size(startPt),[r 1])
                error('STARTPNT must be an rx1 cell array, where r is the number of subdivisions being optimized over.')
            end
            %===end input checking/empty input initialization
            
            %s2.star.bnd is an r-by-1 cell array of column vectors. The
            %length(star.bnd{i}) is equal to the number of response surfaces that were
            %created for the objective function over the i^th element of NODE.
            %s2.star.pairMult, s2.star.paramMult, etc are r-by-1 cell arrays of cells.
            
            %s1.star.bnd is an r-by-1 vector
            %s1.star.pairMult is an r-by-1 cell array of ns-by-2 matrices
            %s1.star.paramMult is an r-by-1 cell array of 2n-by-2 matrices
            
            nNodes = length(nodes);
            n = PD.nParameters;
            nTries = opt.nRestart+1;
            
            s2.star.bnd = zeros(nNodes,1);
            s2.star.pairMult = cell(nNodes,1);
            s2.star.paramMult = cell(nNodes,1);
            s2.star.Q = cell(nNodes,1);
            s2.star.api = cell(nNodes,1);
            
            s1.star.bnd = zeros(nNodes,1);
            s1.star.pairMult = cell(nNodes,1);
            s1.star.paramMult = cell(nNodes,1);
            s1.star.Q = cell(nNodes,1);
            s1.star.api = cell(nNodes,1);
            
            if computeGap
                s2.q.bnd = cell(nNodes,1);
                s2.q.xfeas = cell(nNodes,1);
                s2.q.pairMult = cell(nNodes,1);
                s2.q.paramMult = cell(nNodes,1);
                
                s2.mm.bnd = cell(nNodes,1);
                s2.mm.xfeas = cell(nNodes,1);
                
                s1.q.bnd = zeros(nNodes,1);
                s1.q.xfeas = cell(nNodes,1);
                s1.q.pairMult = cell(nNodes,1);
                s1.q.paramMult = cell(nNodes,1);
                s1.mm.bnd = zeros(nNodes,1);
                s1.mm.xfeas = cell(nNodes,1);
                
                s1.qgap = zeros(nNodes,1);
                s1.mmgap = zeros(nNodes,1);
            end
            
            % analMode = opt.analysisMode;
            
            for i1 = 1:nNodes
                
                [optProb activeIdx] = DClab.ConsistencyTest.makeCmoptim(PD,nodes(i1),'outer');
                [bnd,mult,Q] = lowerBnd(optProb,opt);
                bnd = -bnd;
                
                %%%%works fair
                %   [trash idx] = max(mult.upper+mult.lower);
                %   var = activeIdx(idx);
                %   if var > n
                %     var = var-n;
                %   end
                %%%end works fair
                
                %   ttt = struct(optProb);
                %   tmp = zeros(size(Q));
                %   for i2 = 1:2:length(ttt.ZquadFinal)/2
                %     tmp = tmp+ttt.ZquadFinal{i2};
                %   end
                %   keyboard
                
                %   % determine which variable looks bad to Q
                %   dQ = diag(Q);
                %   dQ = dQ(2:end);
                %   [trash idx] = max(dQ);
                %   var = activeIdx(idx);
                %   if var > n
                %     var = var-n;
                %   end
                
                %   [U S V] = svd(Q([2:end],[2:end]));
                %   nn = length(U);
                %   tmp = zeros(nn,1);
                %   for i2 = 1:nn
                %     tmp = tmp + U(:,i2)*sqrt(S(i2,i2));
                %   end
                %   [trash idx] = max(abs(tmp));
                %   var = activeIdx(idx);
                %   if var > n
                %      var = var-n;
                %   end
                
                %  [U S V] = svd(Q);
                %  S = diag(S);
                %  idx = S > 1e-2;
                %  tmp = zeros(sum(idx),1);
                %  for i2 = 1:sum(idx)
                %    tmp(i2) = U(1,i2)^2*S(i2);
                %  end
                
                %  keyboard
                
                % find which is nearest to 1.
                %  [keep idx] = min(abs(tmp-1));
                
                %  nn = length(Q);
                %  tmp = zeros(nn,1);
                %  for i2 = setdiff(1:nn,idx)
                %    tmp = tmp + U(:,i2)*sqrt(S(i2));
                %  end
                %  % 1st element of tmp is hopefully pretty small.
                %
                %  [trash idx] = max(abs(tmp(2:end)));
                %  var = activeIdx(idx);
                %  if var > n
                %    var = var-n;
                %  end
                
                %  tmp = diag(Q);
                %  idx = max(tmp(2:end));
                %  [trash idx] = max(abs(tmp(2:end)));
                %  var = activeIdx(idx);
                %  if var > n
                %    var = var-n;
                %  end
                
                s2.star.bnd(i1) = bnd;
                
                %TODO The multipliers are to the constraints L+se(1) <= S(x) <=U+se(2).
                %Somewhere outside of this function additional post processing will
                %need to occur to realize the desired "sensitivies".
                
                %TODO, we only want to return one set of multipliers per parameter, so
                %if there was a linX and logX version of a parameter in the
                %optimization, we need to fix this here.
                
                %Determine if we need to include multiplier for logX variables.
                if ~isempty(mult)
                    if max(activeIdx) > n
                        pmult = zeros(2*n,2);
                    else
                        pmult = zeros(n,2);
                    end
                    pmult(activeIdx,:) = [mult.lower mult.upper];
                    emult = [mult.quadineq(1:2:end-1) mult.quadineq(2:2:end)];
                    s2.star.pairMult{i1} = emult;
                    s2.star.paramMult{i1} = pmult;
                end
                
                s2.star.Q{i1} = Q;
                s2.star.api{i1} = activeIdx;
                
                if computeGap
                    s2.q.bnd{i1} = zeros(1,nTries);
                    s2.q.xfeas{i1} = cell(1,nTries);
                    s2.q.pairMult{i1} = cell(1,nTries);
                    s2.q.paramMult{i1} = cell(1,nTries);
                    
                    for i3 = 1:nTries
                        %if ~isempty(startPt{i1}), findStartPt will massage startPt{i1} so
                        %it has the correct dimensions for the optimization. The second
                        %time through the start point will be based on Q. All subsequent
                        %times it will be randomly generated.
                        %
                        %if isempty(startPt{i1}), The first time through a start point will
                        %be based on Q; all subsequent times it will be randomly generated.
                        
                        if isempty(startPt) || isempty(startPt{i1})
                            if i3 == 1
                                xinit = DClab.findStartPt(PD,nodes(i1),activeIdx,[],Q);
                            else
                                xinit = DClab.findStartPt(PD,nodes(i1),activeIdx,[],[]);
                            end
                        else
                            if i3 == 1
                                xinit = DClab.findStartPt(PD,nodes(i1),activeIdx,startPt{i1},Q);
                            elseif i3 == 2
                                xinit = DClab.findStartPt(PD,nodes(i1),activeIdx,[],Q);
                            else
                                xinit = DClab.findStartPt(PD,nodes(i1),activeIdx,[],[]);
                            end
                        end
                        
                        %refresh gui
                        if strcmpi(opt.display,'all')
                            str = [blanks(4) 'fmincon local search try ' num2str(i3) ' of ' num2str(nTries)];
                            DClab.dcdispstr(str,opt.guiHandle,false)
                        else
                            drawnow
                        end
                        
                        [fval,xfeas,mult] = upperBnd(optProb,xinit,opt,true);
                        
                        s2.q.bnd{i1}(1,i3) = -fval;
                        s2.q.xfeas{i1}{1,i3} =  DClab.transformAndAddNominals(xfeas,PD,activeIdx,nodes(i1));
                        
                        if ~isempty(mult)
                            %Determine if we need to include multiplier for logX variables.
                            if max(activeIdx) > n
                                pmult = zeros(2*n,2);
                            else
                                pmult = zeros(n,2);
                            end
                            pmult(activeIdx,:) = [mult.lower mult.upper];
                            emult = [mult.quadineq(1:2:end-1) mult.quadineq(2:2:end)];
                            s2.q.pairMult{i1}{1,i3} = emult;
                            s2.q.paramMult{i1}{1,i3} = pmult;
                        end
                        
                    end
                    
                    %     if strcmp(analMode,'metamodelBasedA') || strcmp(analMode,'metamodelBasedB')
                    %       s2.mm.bnd{i1} = zeros(1,nTries);
                    %       s2.mm.xfeas{i1} = cell(1,nTries);
                    %
                    %       MMCell = PD.PiecewiseSurrogateModelTree(1).DCMetamodelCell;
                    %       if ~isempty(MMCell{1})
                    %
                    %         d = vertcat(PD.ModelAndObservationPair.observedValue);
                    %         u = vertcat(PD.ModelAndObservationPair.observationUncertaintyPlusMinus);
                    %         uncCase = vertcat(PD.ModelAndObservationPair.uncertaintyCase);
                    %
                    %         domRng = PD.PiecewiseSurrogateModelTree(nodes(i1)).domainRange;
                    %         FP = PD.FreeParameter;
                    %
                    %         for i3 = 1:nTries
                    %           [fval xfeas] = DClab.ConsistencyTest.maximizeUB_meta(MMCell,d,u,uncCase,FP,domRng);
                    %
                    %           if isempty(xfeas)
                    %             s2.mm.bnd{i1}(1,i3) = -inf; %set to inf if no feasible point found
                    %           else
                    %             s2.mm.bnd{i1}(1,i3) = fval;
                    %           end
                    %           s2.mm.xfeas{i1}{1,i3} =  xfeas;
                    %         end
                    %       end
                    %     end
                end %if compute gap
                
                s1.star = s2.star;
                if computeGap
                    [trash best2] = max(s2.q.bnd{i1}(1,:));
                    s1.q.bnd(i1) = s2.q.bnd{i1}(1,best2);
                    s1.q.xfeas{i1} = s2.q.xfeas{i1}{1,best2};
                    s1.q.pairMult{i1} = s2.q.pairMult{i1}{1,best2};
                    s1.q.paramMult{i1} = s2.q.paramMult{i1}{1,best2};
                    s1.qgap(i1,1) = s1.q.bnd(i1)-s1.star.bnd(i1);
                    
                    %     if strcmp(analMode,'metamodelBasedA') || strcmp(analMode,'metamodelBasedB')
                    %       [trash best2] = max(s2.mm.bnd{i1}(1,:));
                    %       s1.mm.bnd(i1) = s2.mm.bnd{i1}(1,best2);
                    %       s1.mm.xfeas{i1} = s2.mm.xfeas{i1}{1,best2};
                    %       s1.mmgap(i1,1) = s1.mm.bnd(i1)-s1.q.bnd(i1);
                    %     else
                    s1.mm.bnd(i1) = s1.q.bnd(i1);
                    s1.mm.xfeas{i1} = s1.q.xfeas{i1};
                    s1.mmgap(i1,1) = 0;
                    %     end
                end
            end
        end
        
        
    end % private methods
    
end %classdef



%LOCAL functions

function iter = BBinitialization(Dset,opt,guihand,dbSaveName)

%Initialize the fields of iter in the desired order.
fn = {'optimOpts','PDset','worstUBleaf','dimSplit',...
    'splitLoc','leafNodes','upper','lower','runtime'};
iter = cell2struct(cell(9,1),fn);

iter.upper = struct('bnd',[],'qbnd',[],'sbnd',[],'surrX',[],'pairMult',[],'paramMult',[],'Q',[]);
iter.lower = struct('bnd',[],'xfeas',[]);

%Define initial PDset
if ~isa(Dset,'DClab.PolyDataset')
    if ismember(opt.display,{'iter';'all';'ALL'})
        str = '=====Creating initial PolyDataset from DCDataset=====';
        DClab.dcdispstr(str,guihand,false)
    end
    PDset = DClab.PolyDataset(Dset,opt);
else
    PDset = Dset;
end

%At this point PDset should be defined. Begin to populate iter.

[trash junk iter.optimOpts] = decompose(opt);
iter.PDset = PDset;

% Tell leafInfo structure the leaves of the bounds we're about to compute
iter.leafNodes = PDset.leafNodes;

if ~isempty(dbSaveName)
    save(dbSaveName,'iter')
end

if ~opt.omitOuterBound
    
    % Compute outer bounds
    try
        if ismember(opt.display,{'iter';'all';'ALL'})
            str = '  Computing upper bound';
            DClab.dcdispstr(str,guihand,false)
        end
        
        ub = DClab.ConsistencyTest.maximizeUB(PDset,opt,iter.leafNodes,false);
    catch ME
        str = 'call to mximizeUB failed during upper bound computation of ConsistencyTest.';
        DClab.dcdispstr(str,guihand,false)
        DClab.dcdispstr(ME,guihand,true)
    end
    
    %  % Correct outer upper bound for minimizing a negative function (either
    %  % -M(x) or -log10(M(x)) = log10(1/M(x))) instead of maximizing
    %  for i2 = 1:length(oub.star.bnd)
    %    if strcmp(oub.objTrans{i2},'log10')
    %      oub.star.bnd(i2) = 1/oub.star.bnd(i2);
    %    else
    %      oub.star.bnd(i2) = -oub.star.bnd(i2);
    %    end
    %  end
    
    iter.upper.bnd = ub.star.bnd;
    %iter.upper.qbnd = ub.q.bnd;
    %iter.upper.sbnd = ub.mm.bnd;
    %iter.upper.surrX = ub.mm.xfeas;
    iter.upper.pairMult = ub.star.pairMult;
    iter.upper.paramMult = ub.star.paramMult;
    iter.upper.Q = ub.star.Q;
    iter.upper.api = ub.star.api;
end

if ~opt.omitInnerBound
    % Compute lower bound
    try
        if ismember(opt.display,{'iter';'all';'ALL'})
            str = '  Computing lower bound';
            DClab.dcdispstr(str,guihand,false)
        end
        lb = DClab.ConsistencyTest.maximizeLB(PDset,opt,iter.leafNodes);
    catch ME
        str = 'call to maximizeLB failed in lower bound computation of ConsistencyTest';
        DClab.dcdispstr(str,guihand,false)
        DClab.dcdispstr(ME,guihand,true)
    end
    % Correct the inner upper bound for minimizing a negative: consider the
    % uncertain model eta(x): M(x)+ou(1)<=eta(x)<=M(x)+ou(2). minimizeUB
    % returns M(xopt)+ou(2), as this is an upper bound on min eta(x). However
    % for the upper inner bound, we want M(xopt)+ou(1).
    
    iter.lower.bnd = lb.pound.bnd;
    iter.lower.xfeas = lb.pound.xfeas;
end
end

