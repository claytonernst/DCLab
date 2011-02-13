classdef PosteriorObservation < DClab.DCObject
    %TODO: Update help once good name picked
    
    %RESPONSEPREDICTION predicts the range of a ResponseModel subject to a DCDataset
    %
    %      ResponsePrediction starting from a ResponseModel constrained by a
    %      Dataset or PolyDataset:
    %
    %   PRED = ResponsePrediction(RM,DSET) bounds the extremal values of the
    %   model contained in the ResponseModel object RM over the feasible
    %   set
    %   described by the DCDataset or PolyDataset object DSET.
    %
    %   PRED = ResponsePrediction(RM,DSET,DCOPT) allows you to supply options
    %   to the algorithm via the DCOptions object DCOPT.
    %
    %      ResponsePrediction starting from a ResponseModel and a
    %      PiecewiseSurrogateModelTree for that ResponseModel:
    %
    %   PRED = ResponsePrediction(RM,SMTREE,DSET) bounds the extremal values
    %   of the model contained in the ResponseModel object RM over the feasible
    %   set described by the DCDataset or PolyDataset object DSET. Computational
    %   time is reduced by supplying an existing PiecewiseSurrogateModelTree
    %   SMTREE for RM. SMTREE could perhaps have been obtained from a
    %   previous ResponsePrediction object for RM. Any further subdividing of
    %   SMTREE will be done using the default fitting options.
    %
    %   PRED = ResponsePrediction(RM,SMTREE,DSET,DCOPT) allows you to supply
    %   options to the algorithm via the DCOptions object DCOPT.
    %
    %      ResponsePrediction starting from an existing ResponsePrediction:
    %
    %   PRED = ResponsePrediction(OLDPRED) performs additional iterations of
    %   the algorithm, warm-starting with the ResponsePrediction object
    %   OLDPRED.
    %
    %   PRED = ResponsePrediction(OLDPRED,DCOPT) uses options from the
    %   DCOptions object DCOPT for the additional iterations of the algorithm.
    %
    %   Inputs:
    %      RM: A ResponseModel Object
    %      DSET: A DCDataset object or a PolyDataset object previously formed
    %         from a DCDataset object. If DSET is a PolyDataset and DCOPT is not
    %         supplied, and further subdivision of the domain of PDSET will use
    %         the default fitting options.
    %      SMTREE: A PiecewiseSurrogateModelTree previously formed for RM. If
    %         DCOPT not supplied, any further subdivision of PSMTREE will use
    %         the default fitting options.
    %      OLDPRED: An existing ResponsePrediction object.
    %      DCCOPT: A DCOptions object.  If not suppled, default options
    %         obtained from the DCOptions constructor are used. If
    %         maxBranchBoundIter = 1, no subdividing will occur. If either
    %         omitInnerBound or omitOuterBound is true, you will not be able to
    %         perform additional branch-and-bound iterations at a future time by
    %         calling this constructor with the resulting output PRED.
    %
    %   Outputs:
    %      PRED: a Prediction object with fields .iter and .runDate. .iter is a
    %         structure array with fields .optimOpts, .MnotPD, .PDset,
    %         .worstOBLeaf, .dimSplit, .splitLoc, and .leafInfo. These fields
    %         of the i'th element of .iter contain the following data:
    %         optimOpts: structure containing the optimization options for
    %            the i'th iteration
    %         RM:
    %         PSMTree4RM: PiecewiseSurrogateModelTree containing the
    %            surrogate model for the RM subject to prediction used for the
    %            i'th iteration.
    %         PDset: Partitioned PolyDataset used for the i'th iteration.
    %            Saving this object every iteration is probably unnessary since
    %            in principle the object used in earlier iterations should be
    %            recoverable from the object used in the last iteration.
    %         worstOBLeaf: indicates which leaf of the current PDset gives the
    %            worst outer bound. This leaf will be divided on the next
    %            iteration. This field will be empty if
    %            opt.omitOuterBound==true.
    %         dimSplit: indicates which dimension (i.e., parameter name) of the
    %            the iter(i1-1).worstOBLeaf was divided by the current
    %            iteration. iter(1).dimSplit will be empty.
    %         splitLoc: indicates where dimSplit was divided. iter(1).splitLoc
    %            will be emmpty
    %         .leafNodes: node numbers of each leaf that correspond to the info
    %            in the arrays in .upper.output.bnd, etc.
    %         .upper.outer.bnd: nLeaves-by-1 array of upper output bnds over
    %            the corresponding leaves listed in leafNodes.
    %         .upper.outer.mult:
    %         .upper.inner.bnd:
    %         .upper.inner.mult:
    %         .upper.inner.xfeas:
    %         the same for .lower
    
   
    properties
        runDate;
        iter;
    end
    
    properties (Dependent)
        postOuterBnds
        postInnerBnds
        postBndsx
        postObsRanges
        outerBndMults
        outerBndSens
        nIter
        ResponseModel
        DatasetName
        Dataset
        PolyDataset
        PSMTree4RM
        optimOpts
        surfFittingOpts
    end
    
    methods
        function obj = PosteriorObservation(varargin)
            
            % Calling syntax we will support
            % PosteriorObservation
            % PosteriorObservation(DSET)    (DSET is type DCDataset or PolyDataset)
            % PosteriorObservation(DSET,OPT)
            % PosteriorObservation(POSTOBS)
            % PosteriorObservation(POSTOBS,OPT)
            % PosteriorObservation(CONSISTTEST)
            % PosteriorObservation(CONSISTTEST,OPT)
            
            % We will warn if POSTOBS, or DSET are empty.
            
            %TODO, to save memory, we should probably delete PDSet, RM, and PSMTree4RM
            %from all but the last element of iter.
            
            %TODO: all errors and warnings should be displayed with DClab.dcdispstr ??
            
            ni = nargin;
            no = nargout;
            error(nargchk(0,4,ni));
            error(nargoutchk(0,1,no));
            
            % Special cases (one more is around line 200 below):
            if ni == 0
                return
            end
            if isempty(varargin{1})
                warning('DClab:invalidInputs','First input to RESPONSEPREDICTION is empty, returning an empty object')
                return
            end
            
            % Massage the inputs
            if isa(varargin{1},'DClab.PosteriorObservation')
                warm = true;
            elseif isa(varargin{1},'DClab.DCDataset') || isa(varargin{1},'DClab.ConsistencyTest')
                warm = false;
            else
                error('Inputs: 1st input to POSTERIOROBSERVATION of improper class');
            end
            if ni == 1
                varargin{2} = DClab.DCOptions; %fill in default
            elseif ni == 2
                assert(isa(varargin{2},'DClab.DCOptions'),'Inputs: 2nd input to POSTERIOROBSERVATION of improper class, expecting DCOptions.');
            else
                error('Inputs: when 1st is a PosteriorObservation object, at most two inputs are allowed.');
            end
            
            % At this point, varargin is of length 2 or length 4. When of length 2, it
            % contains a ResponsePrediction and DCOptions object. When of length 4, it
            % contains a ResponseModel,PiecewiseSurrogateModelTree,DCDataset/PolyDataset,
            % and DCOptions object.
            
            
            % The last special case
            if ~warm && isempty(varargin{1})
                warning('DClab:invalidInputs','DSET was empty, returning an empty object')
                return
            end
            
            % Define some useful constants:
            opt = varargin{2};
            guihand = opt.guiHandle;
            if ~isempty(opt.fileName2Save)
                dbSaveName = opt.fileName2Save;
            else
                dbSaveName='';
            end
            
            if opt.omitInnerBound && opt.omitOuterBound
                str = 'Either outer bounds or inner bounds must be computed.';
                DClab.dcdispstr(str,guihand,true)
            end
            
            % This variable will become true if we failed to meet the bbTermTol or if we
            % received a PRED object and no additional iterations were performed.  When
            % true, opt.display='notify' will cause the exit message to be displayed.
            weirdexit = false;
            
            if ~warm
                %We're starting from scratch. Note, since all requested bounds are
                %computed during initialization, we consider the initialization step to
                %be the first iteration.
                
                tic
                iterStruct = BBinitialization(varargin{:},guihand,dbSaveName);
                m=iterStruct(end).PDset.nPairs;
                
                if ~isempty(dbSaveName)
                    save(dbSaveName,'iter')
                end
                
                UBo = zeros(m,1);
                LBo = zeros(m,1);
                LBi = zeros(m,1);
                UBi = zeros(m,1);
                for modelIdx = 1:m
                    % Determine where to subdivide on the next iteration
                    if opt.omitOuterBound && opt.maxBranchBoundIter<=1
                        UBo(modelIdx) = NaN;
                        LBo(modelIdx) = NaN;
                        UBi(modelIdx) = max(iterStruct.upper(modelIdx).inner.bnd);
                        LBi(modelIdx) = min(iterStruct.lower(modelIdx).inner.bnd);
                    elseif opt.omitInnerBound && opt.maxBranchBoundIter<=1
                        UBo(modelIdx) = max(iterStruct.upper(modelIdx).outer.bnd);
                        LBo(modelIdx) = min(iterStruct.lower(modelIdx).outer.bnd);
                        UBi(modelIdx) = NaN;
                        LBi(modelIdx) = NaN;
                    else
                        [UBo(modelIdx) UBoIdx] = max(iterStruct.upper(modelIdx).outer.bnd);
                        [LBo(modelIdx) LBoIdx] = min(iterStruct.lower(modelIdx).outer.bnd);
                        UBi(modelIdx) = max(iterStruct.upper(modelIdx).inner.bnd);
                        LBi(modelIdx) = min(iterStruct.lower(modelIdx).inner.bnd);
                        %if UBo(modelIdx)-UBi(modelIdx) > LBi(modelIdx) - LBo(modelIdx)
                        %    iter.worstOBleaf(modelIdx) = iter.leafNodes(UBoIdx);
                        %else
                        %    iter.worstOBleaf(modelIdx) = iter.leafNodes(LBoIdx);
                        %end
                    end
                    if UBo(modelIdx)<UBi(modelIdx)
                        UBo(modelIdx) = UBi(modelIdx);
                    end
                    if LBo>LBi
                        LBo(modelIdx) = LBi(modelIdx);
                    end
                end
                
                % Should we proceed into the while loop? If not, determine the exit
                % condition.
                if all(UBo-UBi <= opt.branchBoundTermTol) && all(LBi-LBo <= opt.branchBoundTermTol)
                    term = true;
                    exitmsg = 'Exiting with 1 iteration: branchBoundTermTol met';
                elseif opt.maxBranchBoundIter <= 1
                    term = true;
                    weirdexit = true;
                    exitmsg = 'Exiting: maximum number of branch and bound iterations reached';
                %elseif opt.omitInnerBound || opt.omitOuterBound
                %    term = true;
                %    weirdexit = true;
                %    exitmsg = 'Exiting with 1 iteration: both inner and outer bounds required to branch and bound';
                else
                    term = false;
                    i1 = 1;
                end
                
                iterStruct.runtime = toc;
                
            else
                % We've been warm started
                iterStruct = varargin{1}.iter;
                m = iterStruct(end).PDset.nPairs;
                ob = varargin{1}.postOuterBnds;
                ib = varargin{1}.postInnerBnds;
                UBo = ob(:,2);
                UBi = ib(:,2);
                LBi = ib(:,1);
                LBo = ob(:,1);
                for modelIdx=1:m
                    if UBo(modelIdx)<UBi(modelIdx)
                        UBo(modelIdx) = UBi(modelIdx);
                    end
                    if LBo(modelIdx)>LBi(modelIdx)
                        LBo(modelIdx) = LBi(modelIdx);
                    end
                end
                
                %make sure each type of bound is computed if continuing to
                %iterate
                if opt.maxBranchBoundIter > length(iterStruct)
                    if any(isnan(LBi))
                        if ismember(opt.display,{'iter';'all';'ALL'}) && opt.omitInnerBound
                            str = 'Performing a one-time inner bound calculation for iteration convergence';
                            DClab.dcdispstr(str,guihand,false)
                        end
                        for modelIdx = find(isnan(LBi))'
                            PSMTree4RM = Tree4RM(iterStruct(end).PDset.PiecewiseSurrogateModelTree,modelIdx);
                            RM = iterStruct(end).PDset.ModelAndObservationPair(modelIdx).ResponseModel;
                            try
                                ilb = DClab.PosteriorObservation.minimizeUB(PSMTree4RM,RM,iterStruct(end).PDset,'min',opt,iterStruct(end).leafNodes);
                                iub = DClab.PosteriorObservation.minimizeUB(PSMTree4RM,RM,iterStruct(end).PDset,'max',opt,iterStruct(end).leafNodes);
                            catch ME
                                str = 'Error during warm start one-time inner bound calculation.';
                                DClab.dcdispstr(str,guihand,false)
                                DClab.dcdispstr(ME,guihand,true)
                            end
                            for i2 = 1:length(iub.pound.bnd)
                                iub.pound.bnd(i2) = -iub.pound.bnd(i2);
                            end
                            iterStruct(end).upper(modelIdx).inner.bnd = iub.pound.bnd;
                            iterStruct(end).upper(modelIdx).inner.xfeas = iub.pound.xfeas;
                            
                            iterStruct(end).lower(modelIdx).inner.bnd = ilb.pound.bnd;
                            iterStruct(end).lower(modelIdx).inner.xfeas = ilb.pound.xfeas;
                            UBi(modelIdx) = max(iterStruct(end).upper(modelIdx).inner.bnd);
                            LBi(modelIdx) = min(iterStruct(end).lower(modelIdx).inner.bnd);
                        end
                    end
                    if any(isnan(LBo))
                        if ismember(opt.display,{'iter';'all';'ALL'}) && opt.omitOuterBound
                            str = 'Performing a one-time outer bound calculation for iteration convergence';
                            DClab.dcdispstr(str,guihand,false)
                        end
                        for modelIdx = find(isnan(LBo))'
                            PSMTree4RM = Tree4RM(iterStruct(end).PDset.PiecewiseSurrogateModelTree,modelIdx);
                            try
                                olb = DClab.PosteriorObservation.minimizeLB(PSMTree4RM,iterStruct(end).PDset,'min',opt,iterStruct(end).leafNodes,false);
                                oub = DClab.PosteriorObservation.minimizeLB(PSMTree4RM,iterStruct(end).PDset,'max',opt,iterStruct(end).leafNodes,false);
                            catch ME
                                str='Error warm starting one-time outer bound';
                                DClab.dcdispstr(str,guihand,false)
                                DClab.dcdispstr(ME,guihand,true)
                            end
                            
                            for i2 = 1:length(oub.star.bnd)
                                if strcmp(oub.objTrans{i2},'log10')
                                    oub.star.bnd(i2) = 1/oub.star.bnd(i2);
                                else
                                    oub.star.bnd(i2) = -oub.star.bnd(i2);
                                end
                            end
                            iterStruct(end).upper(modelIdx).outer.bnd = oub.star.bnd;
                            iterStruct(end).upper(modelIdx).outer.pairMult = oub.star.pairMult;
                            iterStruct(end).upper(modelIdx).outer.paramMult = oub.star.paramMult;
                            iterStruct(end).upper(modelIdx).outer.Q = oub.star.Q;
                            iterStruct(end).upper(modelIdx).outer.api = oub.star.api;
                            
                            iterStruct(end).lower(modelIdx).outer.bnd = olb.star.bnd;
                            iterStruct(end).lower(modelIdx).outer.pairMult = olb.star.pairMult;
                            iterStruct(end).lower(modelIdx).outer.paramMult = olb.star.paramMult;
                            iterStruct(end).lower(modelIdx).outer.Q = olb.star.Q;
                            iterStruct(end).lower(modelIdx).outer.api = olb.star.api;
                            UBo(modelIdx) = max(iterStruct(end).upper(modelIdx).outer.bnd);
                            LBo(modelIdx) = min(iterStruct(end).lower(modelIdx).outer.bnd);
                        end
                    end
                end
                
                % Should we proceed into the while loop? If not, determine the exit
                % condition.
                %if opt.omitInnerBound || opt.omitOuterBound
                %    term = true;
                %    weirdexit = true;
                %    exitmsg = 'Returning original object: both inner and outer bounds are required to branch and bound';
                %else
                if opt.maxBranchBoundIter <= length(iterStruct)
                    term = true;
                    weirdexit = true;
                    exitmsg = 'Returning original object: maxBranchBoundIter already met in the previous iteration';
                elseif any(isnan(UBo)) || any(isnan(UBi)) || any(isnan(LBi)) || any(isnan(LBo))
                    term = true;
                    weirdexit = true;
                    exitmsg = 'Returning original object: both inner and outer bounds were not computed in the previous iteration';
                elseif max(UBo-UBi) <= opt.branchBoundTermTol && max(LBi-LBo) <= opt.branchBoundTermTol
                    term = true;
                    weirdexit = true;
                    exitmsg = 'Returning original object: branchBoundTermTol already met in the previous iteration';
                else
                    %If we made it this far, we're ready to perform another iteration.
                    if ismember(opt.display,{'iter';'all';'ALL'})
                        str = 'First displayed bounds are for the old ResponsePrediction object';
                        DClab.dcdispstr(str,guihand,false)
                    end
                    term = false;
                    i1 = length(iterStruct);
                end
            end
            
            % Conditional to enter the branch and bound iteration
            while ~term
                
                tic
                % Increament the iteration count
                i1 = i1+1;
                [dispOpts fitOpts iterStruct(i1,1).optimOpts] = decompose(opt);
                %iter(i1).RM = iter(i1-1).RM;
                
                bndU = zeros(m,1);
                Uidx = zeros(m,1);
                bndL = zeros(m,1);
                Lidx = zeros(m,1);
                for modelIdx=1:m
                    [bndU(modelIdx) Uidx(modelIdx)] = max(iterStruct(i1-1).upper(modelIdx).outer.bnd);
                    [bndL(modelIdx) Lidx(modelIdx)] = min(iterStruct(i1-1).lower(modelIdx).outer.bnd);
                end
                
                %upper or lower bound is worse?
                if max(bndU-UBi) > max(LBi-bndL)
                    isMax=true;
                    [M whichMod] = max(bndU-UBi);
                    UB = UBo(whichMod);
                    whichLeaf = Uidx(whichMod);
                else
                    isMax=false;
                    [M whichMod] = max(LBi-bndL);
                    UB = LBo(whichMod);
                    whichLeaf = Lidx(whichMod);
                end
                
                %                 if strcmp(opt.analysisMode,'original')
                %                     % If we're in original mode, assume the problem is with duality.
                if isMax
                    Q = iterStruct(i1-1).upper(whichMod).outer.Q{Uidx(whichMod)};
                    sigma_t = UBi; % shouldn't be used in this case;
                else
                    Q = iterStruct(i1-1).lower(whichMod).outer.Q{Lidx(whichMod)};
                    sigma_t = LBi; % shouldn't be used in this case;
                end
                leaf = iterStruct(i1-1).leafNodes(whichLeaf);
                computeGap = false;
                
                %display for fun
                if ismember(opt.display,{'iter';'all';'ALL'})
                    for modelIdx = 1:m
                        str = sprintf('Outer: [%5.4f %5.4f]  Inner: [%5.4f %5.4f]',LBo(modelIdx),UBo(modelIdx),LBi(modelIdx),UBi(modelIdx));
                        DClab.dcdispstr(str,guihand,false)
                    end
                end
                
                %                 else
                %                     % If we're in one of the metabased modes, determine if the present
                %                     % issue is duality or approximation.
                %                     try
                %                         CDmm = max(iter(i1-1).upper.outer.sbnd);
                %                     catch
                %                         error('Response Prediction can''t handle meta models right now')
                %                     end
                %
                %                     %display for fun
                %                     if ismember(opt.display,{'iter';'all';'ALL'})
                %                         str = sprintf('Outer: [%5.4f %5.4f]  Inner: [%5.4f %5.4f]',LBo(modelIdx),UBo(modelIdx),LBi(modelIdx),UBi(modelIdx));
                %                         DClab.dcdispstr(str,guihand,false)
                %                     end
                %
                %                     if sigma_t-CDmm <= opt.branchBoundQGapTol
                %                         % problem is duality
                %                         if isMax
                %                             Q = iter(i1-1).upper.outer.Q{Uidx};
                %                             sigma_t = UBi; % shouldn't be used in this case;
                %                         else
                %                             Q = iter(i1-1).lower.outer.Q{Lidx};
                %                             sigma_t = LBi; % shouldn't be used in this case;
                %                         end
                %                         leaf = iter(i1-1).leafNodes(Uidx);
                %                         computeGap = false;
                %                         if ismember(opt.display,{'all';'ALL'})
                %                             str = 'Branching will seek to reduce the d-gap';
                %                             DClab.dcdispstr(str,guihand,false)
                %                         end
                %                     else
                %                         % problem is fitting error
                %                         if isMax
                %                             [sigma_t sidx] = max(iter(i1-1).upper.outer.qbnd);
                %                         else
                %                             [sigma_t sidx] = min(iter(i1-1).lower.outer.qbnd);
                %                         end
                %                         Q = [];
                %                         leaf = iter(i1-1).leafNodes(sidx);
                %                         computeGap = true;
                %                         if ismember(opt.display,{'all';'ALL'})
                %                             str = 'Branching will seek to reduce the q-gap';
                %                             DClab.dcdispstr(str,guihand,false)
                %                         end
                %                     end
                %                 end
                
                % Decide where to divide worstOBleaf
                PSMTree4RM = Tree4RM(iterStruct(i1-1).PDset.PiecewiseSurrogateModelTree,whichMod);
                if opt.omitOuterBound
                    if isMax
                        [dimSplit,dimSplit] = max(abs(iterStruct(i1-1).upper(whichMod).inner.xfeas{leaf}));
                    else
                        [dimSplit,dimSplit] = max(abs(iterStruct(i1-1).lower(whichMod).inner.xfeas{leaf}));
                    end
                    splitLoc = mean(iterStruct(i1-1).PDset.PiecewiseSurrogateModelTree(leaf).domainRange(dimSplit,:));
                    dimSplit = iterStruct(i1-1).PDset.FreeParameter(dimSplit).name;
                    %inherit1 = zeros(PSMTree4RM(leaf).nSurfaces,1);
                    inherit2 = zeros(iterStruct(i1-1).PDset.PiecewiseSurrogateModelTree(leaf).nSurfaces,1);
                else
                    [dimSplit, splitLoc, inherit1 inherit2] = DClab.PosteriorObservation.findDimension2Divide(PSMTree4RM,leaf,1,iterStruct(i1-1).PDset,UB,sigma_t,Q,isMax);
                end
                %[dim, loc, inherit] =
                %findDimension2Divide_pred(PSMTree,leaf,surfIdx,PDset,UB,qbnd,Q,isMax)
                iterStruct(i1).dimSplit = dimSplit;
                iterStruct(i1).splitLoc = splitLoc;
                iterStruct(i1-1).worstOBLeaf = leaf;
                
                if any(strmatch(opt.display,{'iter';'all';'ALL'}))
                    str = ['===Iteration ' num2str(i1) ': spliting dimension ' dimSplit ' of leaf ' num2str(leaf) ' at location ' num2str(splitLoc) '==='];
                    DClab.dcdispstr(str,guihand,false)
                end
                
                % Subdivide the domain of PDset and PSMTree where indicated.
                
                % TODO: Do we need to provide the current display options here? See the
                % analogous section of ConsistencyTest. Perhaps not. I'm not sure that
                % there is any reason for the displaySettings fields of a PolyDataset.
                
                %iter(i1,1).PSMTree4RM = subdivideDomain(iter(i1-1).PSMTree4RM,{iter(i1-1).RM},leaf,dimSplit,splitLoc,inherit1,opt);
                iterStruct(i1,1).PDset = subdivideDomain(iterStruct(i1-1).PDset,leaf,dimSplit,splitLoc,inherit2,opt);
                
                iterStruct(i1,1).leafNodes = leafNodes(iterStruct(i1).PDset.PiecewiseSurrogateModelTree);
                
                %The PolyDataset in the new information has all this
                %info...
                iterStruct(i1-1,1).PDset = [];
                
                % If requested, save iter structure
                if ~isempty(dbSaveName)
                    save(dbSaveName,'iter')
                end
                
                % Find which leaves were added by the subdivision
                newLeaves = setdiff(iterStruct(i1).leafNodes,iterStruct(i1-1).leafNodes);
                
                %We only need to do optimization on the newly added leaves. We fix copy
                %the previous info into the current idx of iter. Then we
                %zero only the bnd, mults, etc from the parent of these
                %just added leaves and replace them with the info from the two new
                %leaves. that way we alway have bnd, mults, and xfeas available for each
                %leaf with out having to dig through the structure.
                
                %TODO: does this just complicate things unnecessarily?
                
                %Copy previous values
                iterStruct(i1).upper = iterStruct(i1-1).upper;
                iterStruct(i1).lower = iterStruct(i1-1).lower;
                
                for modelIdx=1:m
                    PSMTree4RM = Tree4RM(iterStruct(i1).PDset.PiecewiseSurrogateModelTree,modelIdx);
                    RM = iterStruct(i1).PDset.ModelAndObservationPair(modelIdx).ResponseModel;
                    if ismember(opt.display,{'iter';'all';'ALL'})
                        str = [sprintf('% 2d',modelIdx) ': ' RM.name];
                        DClab.dcdispstr(str,guihand,false)
                    end
                    if ~opt.omitOuterBound
                        try
                            if ismember(opt.display,{'iter';'all';'ALL'})
                                str = '      Computing outer lower bound';
                                DClab.dcdispstr(str,guihand,false)
                            end
                            olb = DClab.PosteriorObservation.minimizeLB(PSMTree4RM,iterStruct(i1).PDset,'min',opt,newLeaves,computeGap);
                            if ismember(opt.display,{'iter';'all';'ALL'})
                                str = '      Computing outer upper bound';
                                DClab.dcdispstr(str,guihand,false)
                            end
                            oub = DClab.PosteriorObservation.minimizeLB(PSMTree4RM,iterStruct(i1).PDset,'max',opt,newLeaves,computeGap);
                        catch ME
                            str = 'Call in ResponsePrediction to minimizeLB failed';
                            DClab.dcdispstr(str,guihand,false)
                            DClab.dcdispstr(ME,guihand,true)
                        end
                        % Correct outer upper bound for minimizing a negative function (either
                        % -M(x) or -log10(M(x)) = log10(1/M(x))) instead of maximizing
                        for i2 = 1:length(oub.star.bnd)
                            if strcmp(oub.objTrans{i2},'log10')
                                oub.star.bnd(i2) = 1/oub.star.bnd(i2);
                            else
                                oub.star.bnd(i2) = -oub.star.bnd(i2);
                            end
                        end
                    end
                    
                    if ~opt.omitInnerBound
                        try
                            if ismember(opt.display,{'iter';'all';'ALL'})
                                str = '      Computing inner lower bound';
                                DClab.dcdispstr(str,guihand,false)
                            end
                            ilb = DClab.PosteriorObservation.minimizeUB(PSMTree4RM,RM,iterStruct(i1).PDset,'min',opt,newLeaves,computeGap);
                            if ismember(opt.display,{'iter';'all';'ALL'})
                                str = '      Computing inner upper bound';
                                DClab.dcdispstr(str,guihand,false)
                            end
                            iub = DClab.PosteriorObservation.minimizeUB(PSMTree4RM,RM,iterStruct(i1).PDset,'max',opt,newLeaves,computeGap);
                        catch ME
                            str = 'Call in ResponsePrediction to minimizeUB failed';
                            DClab.dcdispstr(str,guihand,false)
                            DClab.dcdispstr(ME,guihand,true)
                        end
                        % Correct the upper inner bound for minimizing a negative: This is less
                        % complicated than the outer bound case because the bound came from
                        % evaluating the true model.
                        for i2 = 1:length(iub.pound.bnd)
                            iub.pound.bnd(i2) = -iub.pound.bnd(i2);
                        end
                    end
                    
                    %do the zeroing out and updating mentioned to before.
                    parentIdx = find(iterStruct(i1-1).leafNodes == iterStruct(i1-1).worstOBLeaf);
                    
                    
                    %Update the leaf Indices
                    iterStruct(i1,1).leafNodes = leafNodes(iterStruct(i1).PDset.PiecewiseSurrogateModelTree);
                    %Eliminate info relevant to the parent leaf that was split and append
                    %info from the two daughter leaves
                    if ~opt.omitOuterBound
                        iterStruct(i1).upper(modelIdx).outer.bnd(parentIdx) = [];
                        iterStruct(i1).upper(modelIdx).outer.pairMult(parentIdx) = [];
                        iterStruct(i1).upper(modelIdx).outer.paramMult(parentIdx) = [];
                        iterStruct(i1).upper(modelIdx).outer.Q(parentIdx) = [];
                        iterStruct(i1).upper(modelIdx).outer.api(parentIdx) = [];
                        iterStruct(i1).upper(modelIdx).outer.bnd = [iterStruct(i1).upper(modelIdx).outer.bnd; oub.star.bnd];
                        iterStruct(i1).upper(modelIdx).outer.pairMult = [iterStruct(i1).upper(modelIdx).outer.pairMult; oub.star.pairMult];
                        iterStruct(i1).upper(modelIdx).outer.paramMult = [iterStruct(i1).upper(modelIdx).outer.paramMult; oub.star.paramMult];
                        iterStruct(i1).upper(modelIdx).outer.Q = [iterStruct(i1).upper(modelIdx).outer.Q; oub.star.Q];
                        iterStruct(i1).upper(modelIdx).outer.api = [iterStruct(i1).upper(modelIdx).outer.api; oub.star.api];
                        
                        iterStruct(i1,1).lower(modelIdx).outer.bnd(parentIdx) = [];
                        iterStruct(i1,1).lower(modelIdx).outer.pairMult(parentIdx) = [];
                        iterStruct(i1,1).lower(modelIdx).outer.paramMult(parentIdx) = [];
                        iterStruct(i1,1).lower(modelIdx).outer.Q(parentIdx) = [];
                        iterStruct(i1,1).lower(modelIdx).outer.api(parentIdx) = [];
                        iterStruct(i1).lower(modelIdx).outer.bnd = [iterStruct(i1).lower(modelIdx).outer.bnd; olb.star.bnd];
                        iterStruct(i1).lower(modelIdx).outer.pairMult = [iterStruct(i1).lower(modelIdx).outer.pairMult; olb.star.pairMult];
                        iterStruct(i1).lower(modelIdx).outer.paramMult = [iterStruct(i1).lower(modelIdx).outer.paramMult; olb.star.paramMult];
                        iterStruct(i1).lower(modelIdx).outer.Q = [iterStruct(i1).lower(modelIdx).outer.Q; olb.star.Q];
                        iterStruct(i1).lower(modelIdx).outer.api = [iterStruct(i1).lower(modelIdx).outer.api; olb.star.api];
                    else
                        iterStruct(i1).upper(modelIdx).outer.bnd = [iterStruct(i1).upper(modelIdx).outer.bnd; iterStruct(i1-1).upper(modelIdx).outer.bnd(end)];
                        iterStruct(i1).upper(modelIdx).outer.pairMult = [iterStruct(i1).upper(modelIdx).outer.pairMult; iterStruct(i1-1).upper(modelIdx).outer.pairMult(end)];
                        iterStruct(i1).upper(modelIdx).outer.paramMult = [iterStruct(i1).upper(modelIdx).outer.paramMult; iterStruct(i1-1).upper(modelIdx).outer.paramMult(end)];
                        iterStruct(i1).upper(modelIdx).outer.Q = [iterStruct(i1).upper(modelIdx).outer.Q; iterStruct(i1-1).upper(modelIdx).outer.Q(end)];
                        iterStruct(i1).upper(modelIdx).outer.api = [iterStruct(i1).upper(modelIdx).outer.api; iterStruct(i1-1).upper(modelIdx).outer.api(end)];
                        
                        iterStruct(i1).lower(modelIdx).outer.bnd = [iterStruct(i1).lower(modelIdx).outer.bnd; iterStruct(i1-1).lower(modelIdx).outer.bnd(end)];
                        iterStruct(i1).lower(modelIdx).outer.pairMult = [iterStruct(i1).lower(modelIdx).outer.pairMult; iterStruct(i1-1).lower(modelIdx).outer.pairMult(end)];
                        iterStruct(i1).lower(modelIdx).outer.paramMult = [iterStruct(i1).lower(modelIdx).outer.paramMult; iterStruct(i1-1).lower(modelIdx).outer.paramMult(end)];
                        iterStruct(i1).lower(modelIdx).outer.Q = [iterStruct(i1).lower(modelIdx).outer.Q; iterStruct(i1-1).lower(modelIdx).outer.Q(end)];
                        iterStruct(i1).lower(modelIdx).outer.api = [iterStruct(i1).lower(modelIdx).outer.api; iterStruct(i1-1).lower(modelIdx).outer.api(end)];
                    end
                    
                    if ~opt.omitInnerBound
                        iterStruct(i1).upper(modelIdx).inner.bnd(parentIdx) = [];
                        iterStruct(i1).upper(modelIdx).inner.xfeas(parentIdx) = [];
                        iterStruct(i1).upper(modelIdx).inner.bnd = [iterStruct(i1).upper(modelIdx).inner.bnd; iub.pound.bnd];
                        iterStruct(i1).upper(modelIdx).inner.xfeas = [iterStruct(i1).upper(modelIdx).inner.xfeas; iub.pound.xfeas];
                        
                        iterStruct(i1).lower(modelIdx).inner.bnd(parentIdx) = [];
                        iterStruct(i1).lower(modelIdx).inner.xfeas(parentIdx) = [];
                        iterStruct(i1).lower(modelIdx).inner.bnd = [iterStruct(i1).lower(modelIdx).inner.bnd; ilb.pound.bnd];
                        iterStruct(i1).lower(modelIdx).inner.xfeas = [iterStruct(i1).lower(modelIdx).inner.xfeas; ilb.pound.xfeas];
                    else
                        iterStruct(i1).upper(modelIdx).inner.bnd = [iterStruct(i1).upper(modelIdx).inner.bnd; iterStruct(i1-1).upper(modelIdx).inner.bnd(end)];
                        iterStruct(i1).upper(modelIdx).inner.xfeas = [iterStruct(i1).upper(modelIdx).inner.xfeas; iterStruct(i1-1).upper(modelIdx).inner.xfeas(end)];
                        
                        iterStruct(i1).lower(modelIdx).inner.bnd = [iterStruct(i1).lower(modelIdx).inner.bnd; iterStruct(i1-1).lower(modelIdx).inner.bnd];
                        iterStruct(i1).lower(modelIdx).inner.xfeas = [iterStruct(i1).lower(modelIdx).inner.xfeas; iterStruct(i1-1).lower(modelIdx).inner.xfeas(end)];
                    end
                    
                    %if computeGap
                    %    iterStruct(i1).upper(modelIdx).outer.qbnd(parentIdx) = [];
                    %    iterStruct(i1,1).lower(modelIdx).outer.qbnd(parentIdx) = [];
                    %    iterStruct(i1).upper(modelIdx).outer.qbnd = [iterStruct(i1).upper(modelIdx).outer.qbnd; oub.pound.bnd];
                    %    iterStruct(i1).lower(modelIdx).outer.qbnd = [iterStruct(i1).lower(modelIdx).outer.qbnd; olb.pound.bnd];
                    %end
                    
                    UBi(modelIdx) = max(iterStruct(i1).upper(modelIdx).inner.bnd);
                    LBi(modelIdx) = min(iterStruct(i1).lower(modelIdx).inner.bnd);
                    
                    [UBo(modelIdx) UBoIdx(modelIdx)] = max(iterStruct(i1).upper(modelIdx).outer.bnd);
                    [LBo(modelIdx) LBoIdx(modelIdx)] = min(iterStruct(i1).lower(modelIdx).outer.bnd);
                    if UBo(modelIdx)<UBi(modelIdx)
                        UBo(modelIdx) = UBi(modelIdx);
                    end
                    if LBo(modelIdx)>LBi(modelIdx)
                        LBo(modelIdx) = LBi(modelIdx);
                    end
                    
                    %if UBo(modelIdx)-UBi(modelIdx) > LBi(modelIdx) - LBo(modelIdx)
                    %    iter(i1).worstOBleaf(modelIdx) = iter(i1).leafNodes(UBoIdx(modelIdx));
                    %else
                    %    iter(i1).worstOBleaf(modelIdx) = iter(i1).leafNodes(LBoIdx(modelIdx));
                    %end
                end
                
                if ~isempty(dbSaveName)
                    save(dbSaveName,'iter')
                end
                
                if max(UBo-UBi) <= opt.branchBoundTermTol && max(LBi-LBo) <= opt.branchBoundTermTol
                    term = true;
                    exitmsg = ['Exiting after ' num2str(i1) ' iterations: branchBoundTermTol met'];
                elseif i1 >= opt.maxBranchBoundIter
                    term = true;
                    weirdexit = true;
                    exitmsg = 'Exiting: maximum number of branch and bound iterations reached';
                else
                    % do nothing, we need to keep iterating
                end
                
                iterStruct(i1).runtime = toc;
                
            end
            
            if strcmp(opt.display,'notify') && weirdexit
                DClab.dcdispstr(exitmsg,guihand,false)
            elseif any(strmatch(opt.display,{'final';'iter';'all';'ALL'}))
                DClab.dcdispstr(exitmsg,guihand,false)
                str = '=====Final Bounds=====';
                DClab.dcdispstr(str,guihand,false)
                for modelIdx = 1:m
                        str = sprintf('Outer: [%5.4f %5.4f]  Inner: [%5.4f %5.4f]',LBo(modelIdx),UBo(modelIdx),LBi(modelIdx),UBi(modelIdx));
                        DClab.dcdispstr(str,guihand,false)
                end
            else
                %display nothing
            end
            
            %assign output
            obj.iter = iterStruct;
            obj.runDate = date;
        end
        
        function bool = isempty(obj)
            % bool = isempty(obj)
            %
            % Prediction object is considered empty if obj(1).MnotPD is.
            %
            % See also Prediction
            
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
                m=obj.iter(end).PDset.nPairs;
                list = cell(m,1);
                out = obj.postOuterBnds;
                in = obj.postInnerBnds;
                for modelIdx = 1:m
                    list{modelIdx} = sprintf('Outer: [%5.4f %5.4f]  Inner: [%5.4f %5.4f]',out(modelIdx,:),in(modelIdx,:));
                end
                sz = sprintf('%d iteration',len);
            end
        end
        
        function out = get.postOuterBnds(obj)
            m = obj.Dataset.nPairs;
            if isempty(obj.iter(end).upper(1).outer)
                out = nan(m,2);
            else
                out = zeros(m,2);
                for i=1:m
                    out(i,1) = min(obj.iter(end).lower(i).outer.bnd);
                    out(i,2) = max(obj.iter(end).upper(i).outer.bnd);
                end
            end
        end
        function out = get.postInnerBnds(obj)
            m = obj.Dataset.nPairs;
            if isempty(obj.iter(end).upper(1).inner)
                out = nan(m,2);
            else
                out = zeros(m,2);
                for i=1:m
                    out(i,1) = min(obj.iter(end).lower(i).inner.bnd);
                    out(i,2) = max(obj.iter(end).upper(i).inner.bnd);
                end
            end
        end
        function out = get.postBndsx(obj)
            opt = obj.iter(end).optimOpts;
            if opt.omitInnerBound && opt.maxBranchBoundIter>1
                idx = 1;
            else
                idx = length(obj.iter);
            end
            m = obj.Dataset.nPairs;
            out = cell(m,2);
            for i=1:m
                [LB,LBi] = min(obj.iter(idx).lower(i).inner.bnd);
                out(i,1) = obj.iter(idx).lower(i).inner.xfeas{LBi};
                [UB,UBi] = max(obj.iter(idx).upper(i).inner.bnd);
                out(i,2) = obj.iter(idx).upper(i).inner.xfeas{UBi};
            end
        end
        function out = get.outerBndMults(obj)
            opt = obj.iter(end).optimOpts;
            if opt.omitOuterBound
                out = struct([]);
            else
                m = obj.iter(end).PDset.nPairs;
                fn = [{'expl','expu','paraml','paramu'}; repmat({cell(m,1)},1,4)];
                out.upper = struct(fn{:});
                out.lower = struct(fn{:});
                for i = 1:m
                    [UB UBoidx] = max(obj.iter(end).upper(i).outer.bnd);
                    [LB LBoidx] = min(obj.iter(end).lower(i).outer.bnd);
                    
                    out.upper.expl = obj.iter(end).upper(i).outer.pairMult{UBoidx}(:,1);
                    out.upper.expu = obj.iter(end).upper(i).outer.pairMult{UBoidx}(:,2);
                    out.lower.expl = obj.iter(end).lower(i).outer.pairMult{LBoidx}(:,1);
                    out.lower.expu = obj.iter(end).lower(i).outer.pairMult{LBoidx}(:,2);
                    
                    
                    upper.tmpl = obj.iter(end).upper(i).outer.paramMult{UBoidx}(:,1);
                    upper.tmpu = obj.iter(end).upper(i).outer.paramMult{UBoidx}(:,2);
                    lower.tmpl = obj.iter(end).lower(i).outer.paramMult{LBoidx}(:,1);
                    lower.tmpu = obj.iter(end).lower(i).outer.paramMult{LBoidx}(:,2);
                    
                    % Combine parameters for linX and logX.
                    PDset = obj.iter(end).PDset;
                    n = PDset.nParameters;
                    uncCase = vertcat(PDset.FreeParameter.uncertaintyCase);
                    range = vertcat(PDset.FreeParameter.range);
                    if length(lower.tmpl) > n
                        both = true;
                    else
                        both = false;
                    end
                    
                    upper.paraml = zeros(n,1);
                    upper.paramu = zeros(n,1);
                    lower.paraml = zeros(n,1);
                    lower.paramu = zeros(n,1);
                    
                    %TODO: check that this math makes some sort of sense.
                    for i1 = 1:n
                        if uncCase(i1)==1 || uncCase(i1)==2
                            % Return multipliers for lin X.
                            if both
                                upper.paraml(i1) = upper.tmpl(i1) + upper.tmpl(i1+n)/log(10)/range(i1,1);
                                upper.paramu(i1) = upper.tmpu(i1) + upper.tmpu(i1+n)/log(10)/range(i1,2);
                                lower.paraml(i1) = lower.tmpl(i1) + lower.tmpl(i1+n)/log(10)/range(i1,1);
                                lower.paramu(i1) = lower.tmpu(i1) + lower.tmpu(i1+n)/log(10)/range(i1,2);
                            else
                                upper.paraml(i1) = upper.tmpl(i1);
                                upper.paramu(i1) = upper.tmpu(i1);
                                lower.paraml(i1) = lower.tmpl(i1);
                                lower.paramu(i1) = lower.tmpu(i1);
                            end
                            
                        else
                            if both
                                upper.paraml(i1) = upper.tmpl(i1)*log(10)*range(i1,1)+upper.tmpl(i1+n);
                                upper.paramu(i1) = upper.tmpu(i1)*log(10)*range(i1,2)+upper.tmpu(i1+n);
                                lower.paraml(i1) = lower.tmpl(i1)*log(10)*range(i1,1)+lower.tmpl(i1+n);
                                lower.paramu(i1) = lower.tmpu(i1)*log(10)*range(i1,2)+lower.tmpu(i1+n);
                            else
                                upper.paraml(i1) = upper.tmpl(i1)*log(10)*range(i1,1);
                                upper.paramu(i1) = upper.tmpu(i1)*log(10)*range(i1,2);
                                lower.paraml(i1) = lower.tmpl(i1)*log(10)*range(i1,1);
                                lower.paramu(i1) = lower.tmpu(i1)*log(10)*range(i1,2);
                            end
                        end
                    end
                    
                    out.upper(i).paraml = upper.paraml;
                    out.upper(i).paramu = upper.paramu;
                    out.lower(i).paraml = lower.paraml;
                    out.lower(i).paramu = lower.paramu;
                end
                
            end
        end
        function out = get.outerBndSens(obj)
            mults = obj.outerBndMults;
            if isempty(mults)
                out = [];
            else
                m = obj.iter(end).PDset.nPairs;
                fn = [{'expl','expu','paraml','paramu'}; repmat({cell(m,1)},1,4)];
                out.upper = struct(fn{:});
                out.lower = struct(fn{:});
                for i=1:m
                    %get transform variables for param sensitivity
                    [UBoidx UBoidx] = max(obj.iter(end).upper(i).outer.bnd);
                    [LBoidx LBoidx] = min(obj.iter(end).lower(i).outer.bnd);
                    Qu = obj.iter(end).upper.outer.Q{UBoidx};
                    Ql = obj.iter(end).lower.outer.Q{LBoidx};
                    worstNodeU = obj.iter(end).leafNodes(UBoidx);
                    worstNodeL = obj.iter(end).leafNodes(LBoidx);
                    worstRangeU = obj.iter(end).PDset.PiecewiseSurrogateModelTree(worstNodeU).domainRange;
                    worstRangeL = obj.iter(end).PDset.PiecewiseSurrogateModelTree(worstNodeL).domainRange;
                    rng = vertcat(obj.iter(end).PDset.FreeParameter.range);
                    
                    %experiment sensitivity
                    out.lower(i).expl = mults.lower.expl;
                    out.lower(i).expu = -mults.lower.expu;
                    out.upper(i).expl = -mults.upper.expl;
                    out.upper(i).expu = mults.upper.expu;
                    
                    %parameter sensitivity
                    out.lower(i).paraml = -mults.lower.paraml;
                    out.lower(i).paramu = -mults.lower.paramu;
                    out.upper(i).paraml = mults.upper.paraml;
                    out.upper(i).paramu = mults.upper.paramu;
                    
                    n = size(totalRange,1);
                    if size(rng)>n
                        for j=1:size(rng,1)
                            if abs(worstRangeL(j,1)-rng(j,1)) < 10*eps
                                out.lower(i).paraml(j) = -mults.lower(i).paraml(j).*(Ql(1+j,1)-1.1);
                            else
                                out.lower(i).paraml(j) = 0;
                            end
                            if abs(worstRangeL(j,2)-rng(j,2)) < 10*eps
                                out.lower(i).paramu(j) = -mults.lower(i).paramu(j).*(Ql(1+j,1)+1.1);
                            else
                                out.lower(i).paramu(j) = 0;
                            end
                            if abs(worstRangeU(j,1)-rng(j,1)) < 10*eps
                                out.upper(i).paraml(j) = mults.upper(i).paraml(j).*(Qu(1+j,1)-1.1);
                            else
                                out.upper(i).paraml(j) = 0;
                            end
                            if abs(worstRangeU(j,2)-rng(j,2)) < 10*eps
                                out.upper(i).paramu(j) = mults.upper(i).paramu(j).*(Qu(1+j,1)+1.1);
                            else
                                out.upper(i).paramu(j) = 0;
                            end
                        end
                    else
                        for j=1:size(rng,1)
                            if abs(worstRangeL(j,1)-rng(j,1)) < 10*eps
                                out.lower(i).paraml(j) = -mults.lower(i).paraml(j).*(0.5*(Ql(1+j,1)+Ql(1+j+n,1))-1.1);
                            else
                                out.lower(i).paraml(j) = 0;
                            end
                            if abs(worstRangeL(j,2)-rng(j,2)) < 10*eps
                                out.lower(i).paramu(j) = -mults.lower(i).paramu(j).*(0.5*(Ql(1+j,1)+Ql(1+j+n,1))+1.1);
                            else
                                out.lower(i).paramu(j) = 0;
                            end
                            if abs(worstRangeU(j,1)-rng(j,1)) < 10*eps
                                out.upper(i).paraml(j) = mults.upper(i).paraml(j).*(0.5*(Qu(1+j,1)+Qu(1+j+n,1))-1.1);
                            else
                                out.upper(i).paraml(j) = 0;
                            end
                            if abs(worstRangeU(j,2)-rng(j,2)) < 10*eps
                                out.upper(i).paramu(j) = mults.upper(i).paramu(j).*(0.5*(Qu(1+j,1)+Qu(1+j+n,1))+1.1);
                            else
                                out.upper(i).paramu(j) = 0;
                            end
                        end
                    end
                    %now convert back to original units
                    out.upper(i).paramu = out.upper(i).paramu.*2./diff(rng,[],2);
                    out.upper(i).paraml = out.upper(i).paraml.*2./diff(rng,[],2);
                    out.lower(i).paramu = out.lower(i).paramu.*2./diff(rng,[],2);
                    out.lower(i).paraml = out.lower(i).paraml.*2./diff(rng,[],2);
                end
            end
        end
        function out = get.nIter(obj)
            out = length(obj.iter);
        end
        function out = get.ResponseModel(obj)
            out = obj.iter(end).RM;
        end
        function out = get.DatasetName(obj)
            out = obj.iter(end).PDset.name;
        end
        function out = get.Dataset(obj)
            %Use the first PDset in iter just because it uses less memory
            Pairs = obj.iter(end).PDset.ModelAndObservationPair;
            FreeParams = obj.iter(end).PDset.FreeParameter;
            name = obj.iter(end).PDset.name;
            out = DClab.DCDataset(Pairs,FreeParams);
            out.name=name;
        end
        function out = get.PolyDataset(obj)
            out = obj.iter(end).PDset;
        end
        function out = get.PSMTree4RM(obj)
            out = obj.iter(end).PSMTree4RM;
        end
        function out = get.optimOpts(obj)
            opt = obj.iter(end).optimOpts;
            out = opt;
        end
        function out = get.surfFittingOpts(obj)
            out = obj.iter(end).PDset.fittingSettings;
        end
        
    end %public methods
    
    
    methods (Access=private, Static)
        
        function [dim, loc, inherit1, inherit2] = findDimension2Divide(PSMTree,leaf,surfIdx,PDset,UB,qbnd,Q,isMax)
            % function [dim, value, inherit] = findDimension2Divide(PDset,leaf,UB,qbnd,Q)
            %
            % Functionality: This method determines where to subdivide the hypercube.
            % to shrink the gap between
            % the local solution to a quadratic problem and the corresponding
            % outerbound generated by the s-procedure.
            %
            % Inputs
            %      PSMTREE: a PiecewiseSurrogateModelTree that provides the objective
            %         function.
            %      LEAF: is a scalar that specifies the domain subdivision (of both
            %         PSMTREE and PD) for the problem. A DCSurface fit over this
            %         subdivision (suitably shifted to account for fitting error)
            %         becomes the objective function.
            %      SURFIDX: a scalar that specifies which DCSurface for the objective fit
            %         over the NODEIDX subdivision should be used. SURFIDX is typically 1.
            %         It may not be if multiple DCSurfaces (each with different I/O
            %         transformations) were fit over the NODEIDX subdivision.
            %  PDset: a PolyDataset object with the same partition structure as PSMTree
            %  UB: outer bound from the rank relaxation
            %  qbnd: inner bound from the local search
            %  Q: Q matrix from the rank relaxation
            %  isMax: booleen indicating whether we are maximizing (working on the
            %    upper bound)
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
            %
            % See also PolyDataset/subdivideDomain, PiecewiseSurrogateModelTree/subdivideDomain
            
            
            pList = PDset.parameterList;
            
            %assume we're trying to reduce d-gap
            
            % Create an optimization problem object, but don't solve it.  We tweak
            % variable bounds in the object to evaluate the potential improvement for
            % splitting different dimensions.
            [optProb activeIdx] = DClab.PosteriorObservation.makeNqcqp(PSMTree,leaf,surfIdx,PDset,'outer');
            
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
                optProbL.variableBnds = bndL;
                bndR = default;
                bndR(worthy(i2),1) = 0;
                optProbR = optProb;
                optProbR.variableBnds = bndR;
                
                if isMax
                    %must be working on the upper bound
                    tmpBnd(i2,1) = -lowerBnd(optProbL,opt); %lowerBnd performs a minimization.  flip the sign.
                    tmpBnd(i2,2) = -lowerBnd(optProbR,opt);
                else
                    tmpBnd(i2,1) = lowerBnd(optProbL,opt); %lowerBnd performs a minimization.  flip the sign.
                    tmpBnd(i2,2) = lowerBnd(optProbR,opt);
                end
            end
            
            % Dump the results to the screen.
            if isMax
                %max(tmpBnd,[],2) %take a max because we need to improve the bound on both the left and right subcubes.
                %min(max(tmpBnd,[],2)) %determine which trial subdivision performed best on both subcubes.
                
                % Determine which of the worthy dimensions performed the best.
                [betterUB idx] = min(max(tmpBnd,[],2));
            else
                %min(tmpBnd,[],2)
                %max(min(tmpBnd,[],2))
                [betterUB idx] = max(min(tmpBnd,[],2));
            end
            
            % Look for a 5% reduction in the gap.  If we fail to acheive this,
            % enlarge the pool of worthy dimensions.
            if abs(UB - betterUB) < 0.05*abs(UB-qbnd)
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
                
                for i2 = reshape(newWorthy,1,length(newWorthy))
                    bndL = default;
                    bndL(worthy(i2),2) = 0;
                    optProbL = optProb;
                    optProbL.variableBnds = bndL;
                    bndR = default;
                    bndR(worthy(i2),1) = 0;
                    optProbR = optProb;
                    optProbR.variableBnds = bndR;
                    if isMax
                        tmpBnd(i2,1) = -lowerBnd(optProbL,opt);
                        tmpBnd(i2,2) = -lowerBnd(optProbR,opt);
                    else
                        tmpBnd(i2,1) = lowerBnd(optProbL,opt);
                        tmpBnd(i2,2) = lowerBnd(optProbR,opt);
                    end
                end
                
                % Dump the results to the screen.
                if isMax
                    %max(tmpBnd,[],2) %take a max because we need to improve the bound on both the left and right subcubes.
                    %min(max(tmpBnd,[],2)) %determine which trial subdivision performed best on both subcubes.
                    
                    % Determine which of the worthy dimensions performed the best.
                    [betterUB idx] = min(max(tmpBnd,[],2));
                else
                    %min(tmpBnd,[],2)
                    %max(min(tmpBnd,[],2))
                    [betterUB idx] = max(min(tmpBnd,[],2));
                end
                
            end
            
            % Look for a 5% reduction in the gap.  If we fail to acheive this,
            % enlarge the pool of worthy dimensions to include all.
            
            % TODO we may
            % want to perform another intermediate pass before we jump all the way to
            % including all.
            
            if abs(UB - betterUB) < 0.05*abs(UB-qbnd)
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
                    
                    if UB > qbnd
                        %must be working on the upper bound
                        tmpBnd(i2,1) = -lowerBnd(optProbL,opt); %lowerBnd performs a minimization.  flip the sign.
                        tmpBnd(i2,2) = -lowerBnd(optProbR,opt);
                    else
                        tmpBnd(i2,1) = lowerBnd(optProbL,opt); %lowerBnd performs a minimization.  flip the sign.
                        tmpBnd(i2,2) = lowerBnd(optProbR,opt);
                    end
                end
                
                % Dump the results to the screen.
                if isMax
                    %max(tmpBnd,[],2) %take a max because we need to improve the bound on both the left and right subcubes.
                    %min(max(tmpBnd,[],2)) %determine which trial subdivision performed best on both subcubes.
                    
                    % Determine which of the worthy dimensions performed the best.
                    [betterUB idx] = min(max(tmpBnd,[],2));
                else
                    %min(tmpBnd,[],2)
                    %max(min(tmpBnd,[],2))
                    [betterUB idx] = max(min(tmpBnd,[],2));
                end
                
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
            
            % Don't recompute any fits.  We're trying to reduce the gap in the
            % S-procedure at this time
            inherit1 = zeros(PSMTree(leaf).nSurfaces,1);
            inherit2 = zeros(PDset.PiecewiseSurrogateModelTree(leaf).nSurfaces,1);
        end
        
        function [obj,activeIdx] = makeNqcqp(PSMTree,nodeIdx,surfIdx,PD,approxType)
            %MAKENQCQP sets up an optimization problem
            %
            %   OBJ = MAKENQCQP(PSMTREE,NODEIDX,SURFIDX,PD,APPROXTYPE) creates an
            %   optimization problem object. The objective is taken from the SURFIDX
            %   surface defined on the subdivision NODEIDX of PSMTREE. PD is a
            %   PolyDataset that describes the constraints.
            %
            %   [OBJ,ACTIVEIDX] = MAKENQCQP(...) additionally returns a vector
            %   indicating which parameters in the PolyDataset became optimization
            %   variables. 1:n indicate the linX variables; n+1:n*2 indicate the logX
            %   variables. ACTIVEIDX will be a subset of 1:2*n.
            %
            %   Inputs:
            %      PSMTREE: a PiecewiseSurrogateModelTree that provides the objective
            %         function.
            %      NODEIDX: is a scalar that specifies the domain subdivision (of both
            %         PSMTREE and PD) for the problem. A DCSurface fit over this
            %         subdivision (suitably shifted to account for fitting error)
            %         becomes the objective function.
            %      SURFIDX: a scalar that specifies which DCSurface for the objective fit
            %         over the NODEIDX subdivision should be used. SURFIDX is typically 1.
            %         It may not be if multiple DCSurfaces (each with different I/O
            %         transformations) were fit over the NODEIDX subdivision.
            %      PD: a PolyDataset object that describes the constraints. The
            %         partition structure of its domain must be identical to that of
            %         PSMTREE.
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
            
            %error(nargchk(5,5,nargin,'struct'));
            
            objSurfs = surfaces(PSMTree,nodeIdx);
            objSurf = objSurfs(surfIdx);
            
            % Determine which parameters are active. There are n parameters in the
            % dataset. Pretend there are 2n, with the first n untransformed and the
            % last n with a log10 transformation. The indices we produce will be
            % indices into this ficticious 2n-by-1 parameter list.
            
            np = PD.nParameters;
            objActiveIdx = objSurf.activeParameterIndex;
            trans = objSurf.activeParameterTransformation;
            objActiveIdx(trans==2) = objActiveIdx(trans==2) + np;
            
            conSurfs = surfaces(PD.PiecewiseSurrogateModelTree,nodeIdx);
            
            if isempty(conSurfs)
                activeIdx = sort(objActiveIdx);
                n = length(activeIdx);
                const = [];
            else
                
                nsPerRM = PD.PiecewiseSurrogateModelTree(nodeIdx).nSurfacesPerResponseModel;
                ns = length(conSurfs);
                
                surfActIdx = cell(ns,1);
                for i1 = 1:ns
                    surfActIdx{i1} = conSurfs(i1).activeParameterIndex;
                    trans = conSurfs(i1).activeParameterTransformation;
                    surfActIdx{i1}(trans==2) = surfActIdx{i1}(trans==2) + np;
                end
                
                % The call to unique will sort activeIdx
                activeIdx = unique([objActiveIdx; vertcat(surfActIdx{:})]);
                n = length(activeIdx);
                
                nConst = 2*sum(nsPerRM);
                const = cell(nConst,1);
                
                L = vertcat(PD.ModelAndObservationPair.constraintLowerLimit);
                U = vertcat(PD.ModelAndObservationPair.constraintUpperLimit);
                
                % Loop through constraint response models and assemble constraints
                for i1 = 1:PD.nPairs
                    
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
                        
                        se = conSurfs(lastSurfIdx+i2).peakError;
                        poly = conSurfs(lastSurfIdx+i2).surrogateModel;
                        nact = length(surfActIdx{lastSurfIdx+i2});
                        respTrans = conSurfs(lastSurfIdx+i2).responseTransformation;
                        
                        switch approxType
                            case 'outer'
                                % do nothing, fitErr is correct
                            case 'crude'
                                se = [0 0];
                            case 'inner'
                                se = [se(2) se(1)];
                        end
                        %allocate for (n+1)x(n=1) matrices
                        tempL = spalloc(n+1,n+1,(nact+1)^2); %constraint from S(x) <=
                        tempR = tempL; %constraint from <= S(x)
                        
                        tempL([1;surfParam2LocInActive+1],[1;surfParam2LocInActive+1]) = -poly;
                        tempR([1;surfParam2LocInActive+1],[1;surfParam2LocInActive+1]) = poly;
                        
                        if strcmp(respTrans,'log10')
                            %The constraints are
                            %  right constraint: S(x) <= log10(M(x)) + se_2 <= log10(U) + se_2
                            %  left constraint: log10(L) + se_1 <= log10(M(x)) + se_1 <= S(x)
                            if L(i1) < 0
                                error('Surrogate in log y should not have been allowed')
                            elseif L(i1) == 0
                                % S(x) is unbounded below. Use blkdiag(-1,zeros(n)) as the
                                % constraint. TODO: does the choice of the (1,1) entry affect
                                % anything?
                                tempL = spalloc(n+1,n+1,1);
                                tempL(1,1) = -1;
                            else
                                tempL(1,1) = tempL(1,1) + log10(L(i1)) + se(1);
                            end
                            if isinf(U(i1))
                                %S(x) is unbounded above. Use blkdiag(-1,zeros(n)) as the
                                % constraint. TODO: does the choice of the (1,1) entry affect
                                % anything?
                                tempR = spalloc(n+1,n+1,1);
                                tempR(1,1) = -1;
                            else
                                tempR(1,1) = tempR(1,1) - log10(U(i1)) - se(2);
                            end
                            
                        else
                            if isinf(L(i1))
                                % S(x) is unbounded below. Use blkdiag(-1,zeros(n)) as the
                                % place-holder constraint. TODO: does the choice of the (1,1) entry affect
                                % anything?
                                tempL = spalloc(n+1,n+1,1);
                                tempL(1,1) = -1;
                            else
                                tempL(1,1) = tempL(1,1) + L(i1) + se(1);
                            end
                            
                            if isinf(U(i1))
                                %S(x) is unbounded above. Use blkdiag(-1,zeros(n)) as the
                                % constraint. TODO: does the choose of the (1,1) entry affect
                                % anything?
                                tempR = spalloc(n+1,n+1,1);
                                tempR(1,1) = -1;
                            else
                                tempR(1,1) = tempR(1,1) - U(i1) - se(2);
                            end
                        end
                        
                        const{(lastSurfIdx+i2)*2-1} = tempL;
                        const{(lastSurfIdx+i2)*2} = tempR;
                    end
                    
                end % for i1 = 1:PD.nPairs
            end
            
            % Find which variables employ both transformations
            linIdx = activeIdx <= np;
            [trash xActWithEqConstLin xActWithEqConstLog] = intersect(activeIdx(linIdx),activeIdx(~linIdx)-np);
            xActWithEqConstLog = xActWithEqConstLog+sum(linIdx);
            
            % For i\in xActWithEqConstLin, activeIdx(i) should be the linear version of
            % a parameter that also has a log version. Similarly for j\in
            % xActWithEqConstLog, activeIdx(j) should be the log version of a parameter
            % that also has a linear version.
            
            range = PSMTree(nodeIdx).domainRange;
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
            
            %constrain bounds on x's. since everything is normalized, this
            %is trivial.
            optimBnds = [-ones(n,1) ones(n,1)];
            
            %build objective function
            fitErr = objSurf.peakError;
            poly = objSurf.surrogateModel;
            
            switch approxType
                case 'outer'
                    % do nothing
                case 'crude'
                    fitErr = [0 0];
                case 'inner'
                    fitErr = [fitErr(2) fitErr(1)];
            end
            
            mm = length(objActiveIdx);
            [sortedSurfParam surfParam2SortedSurfParam] = sort(objActiveIdx);
            [trash sortedSurfParam2LocInActive] = intersect(activeIdx,sortedSurfParam);
            surfParam2LocInActive = zeros(size(sortedSurfParam2LocInActive));
            surfParam2LocInActive(surfParam2SortedSurfParam) = sortedSurfParam2LocInActive;
            
            %allocate for (n+1)x(n=1) matrices
            Z0 = spalloc(n+1,n+1,(mm+1)^2);
            Z0([1;surfParam2LocInActive+1],[1;surfParam2LocInActive+1]) = poly;
            
            % S(x)-M(x) \leq fitErr(2) ==> S(x) - fitErr(2) <= M(x)
            Z0(1) = Z0(1) - fitErr(2);
            
            obj = DClab.nqcqp(Z0,const,optimBnds(:,1),optimBnds(:,2),linXlogX,0);
        end
        
        function [s1 s2] = minimizeLB(PSMTree,PD,extrema,opt,nodes,computeGap,startPt)
            %MINIMIZELB computes a lower bound on the minimim of the PSMTree function
            %
            %   S1 = MINIMIZELB(PSMTREE,PD) determines a lower bound on the minimim of
            %   the function described by the PiecewiseSurrogateModelTree PSMTREE
            %   subject to the constraints implied by the PolyDataset PD. PSDTREE and
            %   PD.PiecewiseSurrogateModelTree must have identical subdivided domains.
            %   Use PiecewiseSurrogateModelTree/mergePartitionStructures to acheive
            %   this. The output S1 is a structure array containing information
            %   described below.
            %
            %   S1 = MINIMIZELB(PSMTREE,PD,EXTREMA) allows you to use this code to
            %   effectively determine an upper bound on the maximum of the PSMTree
            %   function. If EXTREMA='min' (default), it functions as before. If
            %   EXTREMA='max', a lower bound on the minimum of -RM(x) will be returned.
            %
            %   S1 = MINIMIZELB(PSMTREE,PD,EXTREMA,OPT) allows you to supply options
            %   for the optimization with the DCOptions object OPT.
            %
            %   S1 = MINIMIZELB(PSMTREE,PD,EXTREMA,OPT,NODES) allows you specify which
            %   leaf nodes to optimize over. The default is NODES =
            %   PSMTREE.leafNodes.
            %
            %   S1 = MINIMIZELB(PSMTREE,PD,EXTREMA,OPT,NODES,COMPUTEGAP) will perform a
            %   local search for the minimim so you may estimate how well the SDP
            %   relaxation of the quadratic problem did.
            %
            %   S1 = MINIMIZELB(PSMTREE,PD,EXTREMA,OPT,NODES,COMPUTEGAP,STARTPNT)
            %   allows you to supply an initial seed for the local search. STARTPNT
            %   must be an rx1 cell array, where r = length(NODES).
            %
            %   [S1 S2] = MINIMIZELB(...) returns a second output that contains much
            %   more information from the optimization(s)
            %
            %   Inputs:
            %      PSMTREE: A PiecewiseSurrogateModelTree that specifies the objective
            %         function.
            %      PD: A PolyDataset object describing the constraints. Its subdivided
            %         domain must be identical to that of PSMTREE.
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
            % The purpose of this function is to obtain an lower bound on the solution
            % of the mathematical program
            %
            % min Mo(x)   subj. to: |M_e(x)-d_e| \leq u_e                (eq1)
            %                        x in the NODES(i) subdiv of H
            %
            % where Mo is the ResponseModel of the PolyDataset MnotPD, M_e are the
            % true models present in the ResponseModels of the PolyDataset PD and d_e,
            % u_e are the data and uncertainty present in the corresponding
            % ResponseObservations. (actually it is not quite that simple since the
            % uncertainty u_e can be asymmetric and the models Mo, M_e may have output
            % uncertainty).
            %
            % An approximation to the optimization problem (eq1) is considered that
            % takes the form of a nonconvex quadratically constrained quadratic program
            % (nqcqp). The s-procedure is used to lower bound the optimum of each
            % nqcqp. If computeGap = true, a local search is used to determine an upper
            % bound on the optimum for each nqcqp. If the gap between this upper bound
            % and the solution found by local search is zero, the local search found
            % the global minimum of the nqcqp. However, since each nqcqp is just an
            % approximation to eq2 (except in the special case where Mo & M_e's are
            % quadratic and have no output uncertainty), it is not clear that the gap
            % being zero has any practical meaning. Consequently, computeGap = false is
            % the default behavior.
            %
            % Depending on the surrogate models that are available, up to four nqcqps
            % may be constructed using linXlinY, logXlinY, linXlogY, and/or logXlogY
            % transformations.
            %
            % See also PolyDataset, ResposnePrediction, PolyDataset/ilbPred
            
            %TODO, where should we post-process the multipliers? Should this function
            %return multipliers for both linX and logX? What about when there are
            %multiplier surrogates (with different I/0 transformations) for each
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
                case 2
                    extrema = 'min';
                    opt = [];
                    nodes = [];
                    computeGap = 0;
                    startPt = [];
                case 3
                    opt = [];
                    nodes = [];
                    computeGap = 0;
                    startPt = [];
                case 4
                    nodes = [];
                    computeGap = 0;
                    startPt = [];
                case 5
                    computeGap = 0;
                    startPt = [];
                case 6
                    startPt = [];
                otherwise
                    error(nargchk(2,7,ni))
            end
            
            %===input checking/empty input initialization
            %TODO, since this is a private method, we can probably do away with most of
            %the input checks once it has been tested.
            
            if nResponseModels(PSMTree) ~= 1
                error('Inputs: 1st input must be a PiecewiseSurrogateModelTree formed from a single ResponseModel.')
            end
            
            if strcmp(extrema,'max')
                PSMTree = -PSMTree;
            end
            
            %Verify that PSMTree and PD have the same subdivied domains.
            leaf1 = leafNodes(PSMTree);
            leaf2 = PD.leafNodes;
            
            if ~isequal(leaf1,leaf2)
                error('Inputs: PSMTREE and PD do not have domains with identical subdivisions.')
            elseif ~isequal(PSMTree(1).parameterList,PD.parameterList)
                error('Inputs: PSMTREE and PD do not have domains with identical subdivisions.')
            else
                PDPSMTree = PD.PiecewiseSurrogateModelTree;
                for i1 = 1:length(leaf1)
                    if ~isequal(PSMTree(leaf1(i1)).domainRange,PDPSMTree(leaf1(i1)).domainRange)
                        error('Inputs: PSMTREE and PD do not have domains with identical subdivisions.')
                    end
                end
            end
            
            if isempty(nodes)
                nodes = leaf1;
                r = length(nodes);
            else
                r = length(nodes);
                if ~isequal(size(nodes),[r 1]) || ~isnumeric(nodes)
                    error('Inputs: NODES must be a column vector')
                end
                if ~isempty(setdiff(nodes,leaf1))
                    error('Inputs: NODES must be a subset of PSMTREE.leafNodes.')
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
            
            % Can we have multiple transformations on the objective? I don't see why
            % not.
            %
            % So we should loop through the nodes, loop through the different
            % objective function tranformations, and create an optimization for each.
            
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
            
            %s2.objTrans{i}{j} is the jth objective transformation for the ith node
            %Similar for .star.pairMult, .star.paramMult
            s2.objTrans = cell(nNodes,1);
            
            %s2.star.bnd{i1}(j) is the lower bound from the jth objective
            %transformation for the ith node.
            s2.star.bnd = cell(nNodes,1);
            s2.star.pairMult = cell(nNodes,1);
            s2.star.paramMult = cell(nNodes,1);
            s2.star.Q = cell(nNodes,1);
            s2.star.api = cell(nNodes,1);
            
            %s1.objTrans{i} is the objective transformation for the ith node that gave
            %the best, i.e., largest, lower bound among the different objective
            %transformations.
            s1.objTrans = cell(nNodes,1);
            s1.star.bnd = zeros(nNodes,1);
            s1.star.pairMult = cell(nNodes,1);
            s1.star.paramMult = cell(nNodes,1);
            s1.star.Q = cell(nNodes,1);
            s1.star.api = cell(nNodes,1);
            
            if computeGap
                s2.pound.bnd = cell(nNodes,1);
                s2.pound.xfeas = cell(nNodes,1);
                s2.pound.pairMult = cell(nNodes,1);
                s2.pound.paramMult = cell(nNodes,1);
                
                s1.pound.bnd = zeros(nNodes,1);
                s1.pound.xfeas = cell(nNodes,1);
                s1.pound.pairMult = cell(nNodes,1);
                s1.pound.paramMult = cell(nNodes,1);
                
                s1.gap = zeros(nNodes,1);
            end
            
            for i1 = 1:nNodes
                objSurfs = surfaces(PSMTree,i1);
                %nObjSurfs may differ depending on the node.
                nObjSurfs = length(objSurfs);
                s2.star.bnd{i1} = zeros(nObjSurfs,1);
                s2.star.pairMult{i1} = cell(nObjSurfs,1);
                s2.star.paramMult{i1} = cell(nObjSurfs,1);
                s2.objTrans{i1} = cell(nObjSurfs,1);
                
                for i2 = 1:nObjSurfs
                    [optProb activeIdx] = DClab.PosteriorObservation.makeNqcqp(PSMTree,nodes(i1),i2,PD,'outer');
                    s2.objTrans{i1}{i2} = objSurfs(i2).responseTransformation;
                    
                    %Set the is feasible property at some point if we know it is.
                    %I.e., if it was feasible for the first surf, it should be feasible for
                    %the subsequent surfaces.
                    
                    drawnow
                    [bnd,mult,Q] = lowerBnd(optProb,opt);
                    
                    if strcmp(s2.objTrans{i1}{i2},'log10')
                        s2.star.bnd{i1}(i2) = 10.^bnd;
                    else
                        s2.star.bnd{i1}(i2) = bnd;
                    end
                    
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
                        s2.star.pairMult{i1}{i2} = emult;
                        s2.star.paramMult{i1}{i2} = pmult;
                    end
                    
                    s2.star.Q{i1} = Q;
                    s2.star.api{i1} = activeIdx;
                    
                    if computeGap
                        s2.pound.bnd{i1} = zeros(nObjSurfs,nTries);
                        s2.pound.xfeas{i1} = cell(nObjSurfs,nTries);
                        s2.pound.pairMult{i1} = cell(nObjSurfs,nTries);
                        s2.pound.paramMult{i1} = cell(nObjSurfs,nTries);
                        
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
                            
                            if strcmp(strcmp(s2.objTrans{i1}{i2},'log10'),'log10')
                                s2.pound.bnd{i1}(i2,i3) = 10.^-fval;
                            else
                                s2.pound.bnd{i1}(i2,i3) = -fval;
                            end
                            
                            s2.pound.xfeas{i1}{i2,i3} =  DClab.transformAndAddNominals(xfeas,PD,activeIdx,nodes(i1));
                            
                            if ~isempty(mult)
                                %Determine if we need to include multiplier for logX variables.
                                if max(activeIdx) > n
                                    pmult = zeros(2*n,2);
                                else
                                    pmult = zeros(n,2);
                                end
                                pmult(activeIdx,:) = [mult.lower mult.upper];
                                emult = [mult.quadineq(1:2:end-1) mult.quadineq(2:2:end)];
                                s2.pound.pairMult{i1}{i2,i3} = emult;
                                s2.pound.paramMult{i1}{i2,i3} = pmult;
                            end
                        end
                    end
                    
                end
                
                [trash best] = min(s2.star.bnd{i1});
                s1.objTrans{i1} = s2.objTrans{i1}{best};
                
                s1.star.bnd(i1) = s2.star.bnd{i1}(best);
                s1.star.pairMult{i1} = s2.star.pairMult{i1}{best};
                s1.star.paramMult{i1} = s2.star.paramMult{i1}{best};
                s1.star.Q{i1} = s2.star.Q{i1};
                s1.star.api{i1} = s2.star.api{i1};
                
                if computeGap
                    [trash best2] = min(s2.pound.bnd{i1}(best,:));
                    s1.pound.bnd(i1) = s2.pound.bnd{i1}(best,best2);
                    s1.pound.xfeas{i1} = s2.pound.xfeas{i1}{best,best2};
                    s1.pound.pairMult{i1} = s2.pound.pairMult{i1}{best,best2};
                    s1.pound.paramMult{i1} = s2.pound.paramMult{i1}{best,best2};
                    s1.gap(i1,1) = s1.pound.bnd(i1)-s1.star.bnd(i1);
                end
            end
            
            % t = {};
            % if opt.constraints(1)
            %   t = [t; {'linXlinY'}];
            % end
            % if opt.constraints(2)
            %   t = [t; {'logXlinY'}];
            % end
            % if opt.constraints(3)
            %   t = [t; {'linXlogY'}];
            % end
            % if opt.constraints(4)
            %   t = [t; {'logXlogY'}];
            % end
            %
            % % Construct the nqcqp's
            % [nqcqp.linXlinY, nqcqp.logXlinY, nqcqp.linXlogY, nqcqp.logXlogY, activeIdx] = pred_nqcqp(MnotPD,PD,opt,nodeIndices,'outer');
            %
            % %solve the optimizations
            % for i1 = 1:length(t)
            %   [bnd,mult,Q] = lb(nqcqp.(t{i1}),opt);
            %   if strcmp(t{i1},'linXlogY') || strcmp(t{i1},'logXlogY')
            %     s2.(t{i1}).star.bnd = 10.^bnd;
            %   else
            %     s2.(t{i1}).star.bnd = bnd;
            %   end
            %   %TODO: the logY multipliers will be screwed up...need to postprocess some
            %   %day
            %   s2.(t{i1}).star.mult = addMultZeros(mult,PD,opt,activeIdx,nodeIndices);
            %
            %   if computeGap
            %     s2.pound.paramList = MnotPD.parameterList;
            %     for i2 = 1:opt.nRestart+1
            %       xinit = findStartPtPred(MnotPD,activeIdx,nodeIndices,startPt,t{i1},Q,i2);
            %       [fval,xfeas,mult] = ub(nqcqp.(t{i1}),xinit,opt);
            %
            %       if strcmp(t{i1},'linXlogY') || strcmp(t{i1},'logXlogY')
            %         s2.(t{i1}).pound.bnd(:,i2) = 10.^fval;
            %       else
            %         s2.(t{i1}).pound.bnd(:,i2) = fval;
            %       end
            %       %TODO note: the logY multipliers will be screwed up
            %       s2.(t{i1}).pound.xfeas(:,i2) =  DClab.transformAndAddNominals(xfeas,PD,activeIdx,nodeIndices,t{i1});
            %       s2.(t{i1}).pound.mult(:,i2) = addMultZeros(mult,PD,opt,activeIdx,nodeIndices);
            %     end
            %   end
            % end
            %
            % % Determine which method gives better results (big is better). The
            % % contents of the structure s2 are, for example, s2.linXlinY.star.bnd is a
            % % r-by-1 array, s2.linXlinY.star.mult is a cell array of the same
            % % dimensions.
            % tmpBnd = repmat(-inf,r,4);
            % if ismember('linXlinY',t)
            %   tmpBnd(:,1) = s2.linXlinY.star.bnd; %get results from linXlinY
            % end
            % if ismember('logXlinY',t)
            %   tmpBnd(:,2) = s2.logXlinY.star.bnd; %get results from logXlinY
            % end
            % if ismember('linXlogY',t)
            %   tmpBnd(:,3) = s2.linXlogY.star.bnd;
            % end
            % if ismember('logXlogY',t)
            %   tmpBnd(:,4) = s2.logXlogY.star.bnd;
            % end
            % [s1.star.bnd s1.star.trans] = max(tmpBnd,[],2); %#ok
            %
            % % Summarize the relevant results in the structure s1
            % for i1 = 1:r
            %   if s1.star.trans(i1) == 1
            %     best = 'linXlinY';
            %   elseif s1.star.trans(i1) == 2
            %     best = 'logXlinY';
            %   elseif s1.star.trans(i1) == 3
            %     best = 'linXlogY';
            %   else
            %     best = 'logXlogY';
            %   end
            %
            %   % display which was best if there were choices
            %   if length(t) > 1 && strcmpi(opt.display,'all')
            %     str = [blanks(4) best ' transformation gave the best outer bound in olbPred'];
            %     DClab.dcdispstr(str,opt.guiHandle,false)
            %   end
            %
            %   s1.star.mult(i1,1) = s2.(best).star.mult(i1);
            %   if computeGap,
            %     s1.gap(i1,1) = min(s2.(best).pound(i1,:)) - s1.star.bnd(i1);
            %   end
            % end
        end
        
        function [s1 s2] = minimizeUB(PSMTree,RM,PD,extrema,opt,nodes,computeGap,startPt)
            %MINIMIZEUB computes an upper bound on the minimim of the PSMTree function
            %
            %   S1 = MINIMIZEUB(PSMTREE,RM,PD) determines an upper bound on the minimim
            %   of the function described by the ResponseModel RM subject to the
            %   constraints implied by the PolyDataset PD. PSDTREE and
            %   PD.PiecewiseSurrogateModelTree must have identical subdivided domains.
            %   Use PiecewiseSurrogateModelTree/mergePartitionStructures to acheive
            %   this. The output S1 is a structure array containing information
            %   described below.
            %
            %   S1 = MINIMIZEUB(PSMTREE,RM,PD,EXTREMA) allows you to use this code to
            %   effectively determine a lower bound on the maximum of the PSMTree
            %   function. If EXTREMA='min' (default), it functions as before. If
            %   EXTREMA='max', a upper bound on the minimum of -RM(x) will be returned.
            %
            %   S1 = MINIMIZEUB(PSMTREE,RM,PD,EXTREMA,OPT) allows you to supply options
            %   for the optimization with the DCOptions object OPT.
            %
            %   S1 = MINIMIZEUB(PSMTREE,RM,PD,EXTREMA,OPT,NODES) allows you specify
            %   which leaf nodes to optimize over. The default is NODES =
            %   PSMTREE.leafNodes.
            %
            %   S1 = MINIMIZEUB(PSMTREE,RM,PD,EXTREMA,OPT,NODES,COMPUTEGAP) will use
            %   SDP to compute a lower bound on XXX
            %
            %   S1 = MINIMIZEUB(PSMTREE,RM,PD,EXTREMA,OPT,NODES,COMPUTEGAP,STARTPNT)
            %   allows you to supply an initial seed for the local search. STARTPNT
            %   must be an rx1 cell array, where r = length(NODES).
            %
            %   [S1 S2] = MINIMIZEUB(...) returns a second output that contains much
            %   more information (all retries) from the optimization(s)
            %
            %   ALL MULTIPLIERS IN THE OUTPUT DON'T MEAN MUCH. THE FINAL BOUNDS ARE
            %   OBTAINED FROM EVALUATING THE TRUE MODELS.
            %
            %   OLD HELP BELOW
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
            
            ni = nargin;
            switch ni
                case 3
                    extrema = 'min';
                    opt = [];
                    nodes = [];
                    computeGap = 0;
                    startPt = [];
                case 4
                    opt = [];
                    nodes = [];
                    computeGap = 0;
                    startPt = [];
                case 5
                    nodes = [];
                    computeGap = 0;
                    startPt = [];
                case 6
                    computeGap = 0;
                    startPt = [];
                case 7
                    startPt = [];
                otherwise
                    error(nargchk(2,7,ni))
            end
            
            %===input checking/empty input initialization
            %TODO, since this is a private method, we can probably do away with most of
            %the input checks once it has been tested.
            
            if nResponseModels(PSMTree) ~= 1
                error('Inputs: 1st input must be a PiecewiseSurrogateModelTree formed from a single ResponseModel.')
            end
            
            if strcmp(extrema,'max')
                PSMTree = -PSMTree;
            end
            
            %Verify that PSMTree and PD have the same subdivied domains.
            leaf1 = leafNodes(PSMTree);
            leaf2 = PD.leafNodes;
            if ~isequal(leaf1,leaf2)
                error('Inputs: PSMTREE and PD do not have domains with identical subdivisions.')
            elseif ~isequal(PSMTree(1).parameterList,PD.parameterList)
                error('Inputs: PSMTREE and PD do not have domains with identical subdivisions.')
            else
                PDPSMTree = PD.PiecewiseSurrogateModelTree;
                for i1 = 1:length(leaf1)
                    if ~isequal(PSMTree(leaf1(i1)).domainRange,PDPSMTree(leaf1(i1)).domainRange)
                        error('Inputs: PSMTREE and PD do not have domains with identical subdivisions.')
                    end
                end
            end
            
            if isempty(nodes)
                nodes = leaf1;
                r = length(nodes);
            else
                r = length(nodes);
                if ~isequal(size(nodes),[r 1]) || ~isnumeric(nodes)
                    error('Inputs: NODES must be a column vector')
                end
                if ~isempty(setdiff(nodes,leaf1))
                    error('Inputs: NODES must be a subset of PSMTREE.leafNodes.')
                end
            end
            
            if isempty(opt)
                opt = DClab.DCOptions;
            end
            if isempty(computeGap)
                computeGap = false;
            end
            if ~isempty(startPt) && ~iscell(startPt) && ~isequal(size(startPt),[r 1])
                error('STARTPNT must be an rx1 cell array, where r is the number of subdomains being optimized over.')
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
            paramList = PD.parameterList;
            
            %Initialize some variables.
            %s2.objTrans{i}{j} is the jth objective transformation for the ith node
            %Similar for .star.pairMult, .star.paramMult
            s2.objTrans = cell(nNodes,1);
            
            %s2.m1.pound.bnd{i}(j) is the upper bound from method 1, from the jth
            %objective transformation for the ith node.
            s2.m1.pound.bnd = cell(nNodes,1);
            s2.m1.pound.xfeas = cell(nNodes,1);
            s2.m1.pound.pairMult = cell(nNodes,1);
            s2.m1.pound.paramMult = cell(nNodes,1);
            
            s2.m2.pound.bnd = cell(nNodes,1);
            s2.m2.pound.xfeas = cell(nNodes,1);
            s2.m2.pound.pairMult = cell(nNodes,1);
            s2.m2.pound.paramMult = cell(nNodes,1);
            
            s1.objTrans = cell(nNodes,1);
            s1.pound.bnd = zeros(nNodes,1);
            s1.pound.xfeas = cell(nNodes,1);
            s1.pound.pairMult = cell(nNodes,1);
            s1.pound.paramMult = cell(nNodes,1);
            
            if computeGap
                s2.m1.star.bnd = cell(nNodes,1);
                s2.m1.star.pairMult = cell(nNodes,1);
                s2.m1.star.paramMult = cell(nNodes,1);
                
                s1.star.bnd = zeros(nNodes,1);
                s1.star.pairMult = cell(nNodes,1);
                s1.star.paramMult = cell(nNodes,1);
                
                s1.gap = zeros(nNodes,1);
            end
            
            for i1 = 1:nNodes
                
                objSurfs = surfaces(PSMTree,i1);
                
                %nObjSurfs may differ depending on the node.
                nObjSurfs = length(objSurfs);
                s2.objTrans{i1} = cell(nObjSurfs,1);
                
                s2.m1.pound.bnd{i1} = zeros(nObjSurfs,nTries);
                s2.m1.pound.pairMult{i1} = cell(nObjSurfs,nTries);
                s2.m1.pound.paramMult{i1} = cell(nObjSurfs,nTries);
                
                for i2 = 1:nObjSurfs
                    
                    s2.objTrans{i1}{i2} = objSurfs(i2).responseTransformation;
                    
                    [optProb.m1 activeIdx] = DClab.PosteriorObservation.makeNqcqp(PSMTree,nodes(i1),i2,PD,'inner');
                    optProb.m2 = DClab.PosteriorObservation.makeNqcqp(PSMTree,nodes(i1),i2,PD,'crude');
                    
                    if computeGap
                        
                        s2.m1.star.bnd{i1} = zeros(nObjSurfs,1);
                        s2.m1.star.pairMult{i1} = cell(nObjSurfs,1);
                        s2.m1.star.paramMult{i1} = cell(nObjSurfs,1);
                        
                        drawnow
                        [bnd,mult,Q] = lowerBnd(optProb.m1,opt);
                        if strcmp(s2.objTrans{i1}{i2},'log10')
                            s2.m1.star.bnd{i1}(i2) = 10.^bnd;
                        else
                            s2.m1.star.bnd{i1}(i2) = bnd;
                        end
                        
                        %Determine if we need to include multiplier for logX variables.
                        if max(activeIdx) > n
                            pmult = zeros(2*n,2);
                        else
                            pmult = zeros(n,2);
                        end
                        pmult(activeIdx,:) = [mult.lower mult.upper];
                        emult = [mult.quadineq(1:2:end-1) mult.quadineq(2:2:end)];
                        s2.m1.star.pairMult{i1}{i2} = emult;
                        s2.m1.star.paramMult{i1}{i2} = pmult;
                        
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
                        
                        if isempty(xfeas)
                            s2.m1.pound.bnd{i1}(i2,i3) = inf; %set to inf if no feasible point found
                        else
                            [y yInt] = RM.eval(xfeas,paramList);
                            if strcmp(extrema,'min')
                                s2.m1.pound.bnd{i1}(i2,i3) = yInt(2);
                            else
                                s2.m1.pound.bnd{i1}(i2,i3) = -yInt(1);
                            end
                        end
                        
                        %This whole thing is kinda odd. How do we handle the output
                        %uncertainty for inner bounds? Unless we ignore it or subtract it,
                        %using the output uncertainty in the prediction doesn't make much
                        %sense. More work needs to be done in makeNqcqp!!
                        
                        %note: mults for logY are not scaled
                        
                        s2.m1.pound.xfeas{i1}{i2,i3} =  xfeas;
                        
                        %Determine if we need to include multiplier for logX variables.
                        if max(activeIdx) > n
                            pmult = zeros(2*n,2);
                        else
                            pmult = zeros(n,2);
                        end
                        pmult(activeIdx,:) = [mult.lower mult.upper];
                        emult = [mult.quadineq(1:2:end-1) mult.quadineq(2:2:end)];
                        s2.m1.pound.pairMult{i1}{i2,i3} = emult;
                        s2.m1.pound.paramMult{i1}{i2,i3} = pmult;
                    end %for i3 = 1:nTries
                    
                    %TODO is isequal a sufficiently fast why to check if there is any
                    %surrogate modeling error?
                    if ~isequal(optProb.m1,optProb.m2)
                        % Find feasible points of objects
                        state = warning('off','MATLAB:class:cannotUpdateClass');
                        optProb.m2 = struct(findFeasPnt(optProb.m2,opt));
                        warning(state);
                        
                        %There is not optimization, hence no gap, for method 2.
                        
                        xfeas = optProb.m2.xfeas; %should be a 1-by-5 cell array
                        N = size(xfeas,2);
                        
                        for i3 = 1:N
                            xfeas{i3} = DClab.transformAndAddNominals(xfeas{i3},PD,activeIdx,nodes(i1));
                            if isempty(xfeas{i3})
                                feasBool = false;
                            else
                                feasBool = isFeasible(PD,xfeas{i3});
                            end
                            if ~feasBool
                                s2.m2.pound.bnd{i1}(i2,i3) = inf;
                            else
                                [y yInt] = RM.eval(xfeas{i3},paramList);
                                if strcmp(extrema,'min')
                                    s2.m2.pound.bnd{i1}(i2,i3) = yInt(2);
                                else
                                    s2.m2.pound.bnd{i1}(i2,i3) = -yInt(1);
                                end
                            end
                        end
                        s2.m2.pound.xfeas{i1}(i2,:) =  xfeas;
                        s2.m2.pound.pairMult{i1}(i2,:) = cell(size(xfeas)); %no mults for method 2
                        s2.m2.pound.paramMult{i1}(i2,:) = cell(size(xfeas));
                    end
                end %for i2 = 1:nObjSurfs
                
                m1min = min(min(s2.m1.pound.bnd{i1}));
                if ~isempty(s2.m2.pound.bnd{i1})
                    m2min = min(min(s2.m2.pound.bnd{i1}));
                    if m1min <= m2min
                        bestmethod = 'm1';
                    else
                        bestmethod = 'm2';
                    end
                    if strcmpi(opt.display,'all')
                        str = [blanks(4) 'method ' bestmethod(2) ' gave the best inner bound in minimizeUB'];
                        DClab.dcdispstr(str,opt.guiHandle,false)
                    end
                else
                    bestmethod = 'm1';
                end
                
                [rowmin best2] = min(s2.(bestmethod).pound.bnd{i1},[],2);
                [trash best1] = min(rowmin);
                
                s1.objTrans{i1} = s2.objTrans{i1}{best1};
                s1.pound.bnd(i1) = s2.(bestmethod).pound.bnd{i1}(best1,best2);
                s1.pound.xfeas{i1} = s2.(bestmethod).pound.xfeas{i1}{best1,best2};
                s1.pound.pairMult{i1} = s2.(bestmethod).pound.pairMult{i1}{best1,best2};
                s1.pound.paramMult{i1} = s2.(bestmethod).pound.paramMult{i1}{best1,best2};
                
                if computeGap && strcmp(bestmethod,'m1')
                    s1.star.bnd(i1) = s2.(bestmethod).star.bnd{i1}(best1);
                    s1.star.pairMult{i1} = s2.(bestmethod).star.pairMult{i1}{best1};
                    s1.star.paramMult{i1} = s2.(bestmethod).star.paramMult{i1}{best1};
                    s1.gap(i1,1) = s1.pound.bnd(i1)-s1.star.bnd(i1);
                end
            end
            
        end
    end
    
end %classdef


function iter = BBinitialization(Dset,opt,guihand,dbSaveName)

%Initialize the fields of iter in the desired order.
fn = {'optimOpts','PDset','worstOBLeaf','dimSplit',...
    'splitLoc','leafNodes','upper','lower','runtime'};
iter = cell2struct(cell(9,1),fn);



%Define initial PDset
if isa(Dset,'DClab.ConsistencyTest')
    PDset = Dset.iter(end).PDset;
elseif ~isa(Dset,'DClab.PolyDataset')
    if ismember(opt.display,{'iter';'all';'ALL'})
        str = '=====Creating initial PolyDataset from DCDataset=====';
        DClab.dcdispstr(str,guihand,false)
    end
    PDset = DClab.PolyDataset(Dset,opt);
else
    PDset = Dset;
end

m = PDset.nPairs;

fn = {'outer',cell(m,1),'inner',cell(m,1)};
iter.upper = struct(fn{:});
iter.lower = struct(fn{:});

% Make sure the partition structures are identical
% dsetTree = PDset.PiecewiseSurrogateModelTree;
% [dsetTree PSMTree4RM] = mergePartitionStructures(dsetTree,PSMTree4RM);
% PDset.PiecewiseSurrogateModelTree = dsetTree;

%At this point, RM, PSMTree, and PDset should be defined. Begin
%to populate iter.

[trash junk iter.optimOpts] = decompose(opt);
%iter.RM = RM;
%iter.PSMTree4RM = PSMTree4RM;
iter.PDset = PDset;

% Tell leafInfo structure the leaves of the bounds we're about to compute
iter.leafNodes = leafNodes(PDset.PiecewiseSurrogateModelTree);

if ~isempty(dbSaveName)
    save(dbSaveName,'iter')
end

for modIdx = 1:m
    PSMTree4RM = Tree4RM(PDset.PiecewiseSurrogateModelTree,modIdx);
    RM = PDset.ModelAndObservationPair(modIdx).ResponseModel;
    if ismember(opt.display,{'iter';'all';'ALL'})
        str = [sprintf('% 2d',modIdx) ': ' RM.name];
        DClab.dcdispstr(str,guihand,false)
    end
    if ~opt.omitOuterBound || opt.maxBranchBoundIter>1
        
        % Compute outer bounds
        try
            %currently never computeGap...
            %if opt.maxBranchBoundIter>1
            %    computeGap = true;
            %else
            computeGap = false;
            %end
            if ismember(opt.display,{'iter';'all';'ALL'}) && opt.omitOuterBound
                str = '      Computing one-time outer lower bound';
                DClab.dcdispstr(str,guihand,false)
            elseif ismember(opt.display,{'iter';'all';'ALL'})
                str = '      Computing outer lower bound';
                DClab.dcdispstr(str,guihand,false)
            end
            olb = DClab.PosteriorObservation.minimizeLB(PSMTree4RM,PDset,'min',opt,iter.leafNodes,computeGap);
            if ismember(opt.display,{'iter';'all';'ALL'}) && opt.omitOuterBound
                str = '      Computing one-time outer upper bound';
                DClab.dcdispstr(str,guihand,false)
            elseif ismember(opt.display,{'iter';'all';'ALL'})
                str = '      Computing outer upper bound';
                DClab.dcdispstr(str,guihand,false)
            end
            oub = DClab.PosteriorObservation.minimizeLB(PSMTree4RM,PDset,'max',opt,iter.leafNodes,computeGap);
        catch ME
            str = 'call to minimizeLB failed during outer bound computation of ResponsePrediction.';
            DClab.dcdispstr(str,guihand,false)
            DClab.dcdispstr(ME,guihand,true)
        end
        
        % Correct outer upper bound for minimizing a negative function (either
        % -M(x) or -log10(M(x)) = log10(1/M(x))) instead of maximizing
        for i2 = 1:length(oub.star.bnd)
            if strcmp(oub.objTrans{i2},'log10')
                oub.star.bnd(i2) = 1/oub.star.bnd(i2);
            else
                oub.star.bnd(i2) = -oub.star.bnd(i2);
            end
        end
        
        iter.upper(modIdx).outer.bnd = oub.star.bnd;
        iter.upper(modIdx).outer.pairMult = oub.star.pairMult;
        iter.upper(modIdx).outer.paramMult = oub.star.paramMult;
        iter.upper(modIdx).outer.Q = oub.star.Q;
        iter.upper(modIdx).outer.api = oub.star.api;
        
        iter.lower(modIdx).outer.bnd = olb.star.bnd;
        iter.lower(modIdx).outer.pairMult = olb.star.pairMult;
        iter.lower(modIdx).outer.paramMult = olb.star.paramMult;
        iter.lower(modIdx).outer.Q = olb.star.Q;
        iter.lower(modIdx).outer.api = olb.star.api;
        
        %if computeGap
        %    iter.upper(modIdx).outer.qbnd = oub.pound.bnd;
        %    iter.lower(modIdx).outer.qbnd = olb.pound.bnd;
        %else
            iter.upper(modIdx).outer.qbnd = [];
            iter.lower(modIdx).outer.qbnd = [];
        %end
        
    end
    
    if ~opt.omitInnerBound || opt.maxBranchBoundIter>1
        % Compute inner bounds
        try
            if ismember(opt.display,{'iter';'all';'ALL'}) && opt.omitInnerBound
                str = '      Computing one-time inner lower bound';
                DClab.dcdispstr(str,guihand,false)
            elseif ismember(opt.display,{'iter';'all';'ALL'})
                str = '      Computing inner lower bound';
                DClab.dcdispstr(str,guihand,false)
            end
            ilb = DClab.PosteriorObservation.minimizeUB(PSMTree4RM,RM,PDset,'min',opt,iter.leafNodes);
            if ismember(opt.display,{'iter';'all';'ALL'}) && opt.omitInnerBound
                str = '      Computing one-time inner upper bound';
                DClab.dcdispstr(str,guihand,false)
            elseif ismember(opt.display,{'iter';'all';'ALL'})
                str = '      Computing inner upper bound';
                DClab.dcdispstr(str,guihand,false)
            end
            iub = DClab.PosteriorObservation.minimizeUB(PSMTree4RM,RM,PDset,'max',opt,iter.leafNodes);
        catch ME
            str = 'call to minimizeUB failed in inner bound computation of ResponsePrediction';
            DClab.dcdispstr(str,guihand,false)
            DClab.dcdispstr(ME,guihand,true) %display error
        end
        % Correct the upper inner bound for minimizing a negative: This is less
        % complicated than the outer bound case because the bound came from
        % evaluating the true model.
        for i2 = 1:length(iub.pound.bnd)
            iub.pound.bnd(i2) = -iub.pound.bnd(i2);
        end
        
        iter.upper(modIdx).inner.bnd = iub.pound.bnd;
        iter.upper(modIdx).inner.xfeas = iub.pound.xfeas;
        
        iter.lower(modIdx).inner.bnd = ilb.pound.bnd;
        iter.lower(modIdx).inner.xfeas = ilb.pound.xfeas;
    end
end
end



