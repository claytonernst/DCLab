function [s1 s2] = olbFitX(PD,opt,norm,trgWeights,nodeIndices,computeGap,startPt)
% [s1 s2] = olbFitX(PD,opt,norm,trgWeights,nodeIndices,computeGap,startPt)
%
% olbFitX is used to lower bound the norm cost
% on the residual vector.
%
% Inputs:
%   PD: A PolyDataset object.
%   opt: A DCOptions object or structure with fields .display,
%     .sedumiParsEps, .tolFun, .tolCon, .nRestart, .constraints, and .guihandle
%   norm: either 'one','two',or 'inf'.
%   trgWeights: A mx1 vector of weights w_e in the above
%     optimization.
%   nodeIndices[optional]: rx1 vector of indicies into the surrogateTree
%     of PD. This specifies which r subdivisions of H we will approximate the
%     above optimization over.
%   computeGap[optional]: A boolean. If computeGap = 1, the s-procedure
%     will be used to lower bound the solution to all nqcqp's. See below.
%   startPt[optional]: A rx1 cell array of nx1 vectors of starting
%     points for the optimizations. The order should correspond to the
%     order of the parameters in PD.
%
% Outputs:
%   s1: structure summarizing the results, further detail below
%   s2: structure containing all results, further detail below
%
%
% See also PolyDataset, ParameterOptimization, PolyDataset/ilbFitX

ni = nargin;

switch ni
    case 4
        nodeIndices = [];
        computeGap = 0;
        startPt = [];
    case 5
        computeGap = 0;
        startPt = [];
    case 6
        startPt = [];
    otherwise
        error(nargchk(4,7,ni))
end

%input checking/empty input initialization
if isempty(nodeIndices)
    nodeIndices = PD.leafNodes;
end
if isempty(computeGap)
    computeGap = false;
end
r = length(nodeIndices);
if ~isempty(startPt) && ~iscell(startPt) && ~isequal(size(startPt),[r 1])
    error('startPt must be a cell array')
end
%end input checking/empty input initialization

%we can chose to solve the problem by two methods:
% m1: use outer approximation in linX
% m2: use outer approximation in logX

%currently we only can handle linXlinY transformations, and we require their
%existence
if ~opt.constraints(1)
    error('Inputs, currently can not handle your transformations')
end
t = {'m1'};

switch norm
    case 'one'
        [optProb.m1, optProb.m2, activeIdx] = optim1norm_nqcqp(PD,opt,trgWeights,'outer',nodeIndices);
    case 'two'
        [optProb.m1, optProb.m2, activeIdx] = optim2norm_nqcqp(PD,opt,trgWeights,'outer',nodeIndices);
    case 'inf'
        [optProb.m1, optProb.m2, activeIdx] = optimInfnorm_nqcqp(PD,opt,trgWeights,'outer',nodeIndices);
    otherwise
        error('Internal inconsistency, check code')
end

for i1 = 1:length(t)

    bnd = zeros(r,1);
    mult = cell(r,1);
    Q = cell(r,1);
    for i2=1:r
        [bnd(i2),mult{i2},Q{i2}] = lowerBnd(optProb.(t{i1}){i2},opt);
    end
%     [bnd,mi] = min(bnd);
%     mult = mult{mi};
%     Q = Q{mi};
    switch norm
        case {'one','inf'}
            s2.(t{i1}).star.bnd = bnd;
        case 'two'
            bnd = max(bnd,0);
            s2.(t{i1}).star.bnd = sqrt(bnd);
        otherwise
            error('Internal inconsistency, condition should never occur')
    end

    if computeGap
        s2.paramList = PD.parameterList;
        for i2 = 1:opt.nRestart+1
            switch norm
                case {'one','two'}
                    xinit = findStartPtOptim(PD,activeIdx,nodeIndices,startPt,t{i1},Q,i2);
                case 'inf'
                    xinit = findStartPtCons(PD,activeIdx,nodeIndices,startPt,t{i1},Q,i2); %like consistency, just one slack variable
                otherwise
                    error('Internal inconsistency, condition should never occur')
            end

            %refresh gui
            if strcmpi(opt.display,'all')
                str = [blanks(4) 'fmincon local search try ' num2str(i1) ' of ' num2str(opt.nRestart+1)];
                DClab.dcdispstr(str,opt.guiHandle,false)
            else
                drawnow
            end
            r = length(xinit);
            fval = zeros(r,1);
            xfeas = cell(r,1);
            mult = cell(r,1);
            for i3=1:r
                fval(i3) = qpub_FMINCON(optProb.(t{i1}){i3},xinit{i3},opt);
            end

            if strcmp(norm,'two')
                fval = min(fval,0);
                fval = sqrt(fval);
            end
            %Remove last m xfeas and multiplier that correspond to
            %gammas. We don't want to see these in our final output.
            %       switch norm
            %         case {'one','two'}
            %           m = nUnits(PD);
            %           for i3 = 1:r
            %             xfeas{i3}(end-m+1:end) = [];
            %             mult{i3}.lower(end-m+1:end) = [];
            %             mult{i3}.upper(end-m+1:end) = [];
            %           end
            %         case 'inf'
            %           for i3 = 1:r
            %             xfeas{i3}(end) = [];
            %             mult{i3}.lower(end) = [];
            %             mult{i3}.upper(end) = [];
            %           end
            %         otherwise
            %           error('Internal inconsistency, condition should never occur')
            %       end

            s2.(t{i1}).pound.bnd(:,i2) = fval;
            %Convert xfeas to original units, and if any DCDataset variables
            %were not optimization variables, insert their nominal values
            %into the returned output.
            %      s2.(t{i1}).pound.xfeas(:,i2) =  DClab.transformAndAddNominals(xfeas,PD,activeIdx,nodeIndices,t{i1});
            %Add multipliers of 0 to any dataset parameter that where not
            %optimization variables.
            %      s2.(t{i1}).pound.mult(:,i2) = addMultZeros(mult,PD,opt,activeIdx,nodeIndices);
        end
    end
end

% Determine which method gives better results (since we're lower bounding,
% big is better). The contents of the structure s2 are, for example,
% s2.m1.pound.bnd is a r-by-nRestarts+1 array, s2.m1.pound.xfeas is a cell
% array of the same dimensions, and so is s2.m1.pound.mult.
tmpBnd = repmat(-inf,r,2);
if ismember('m1',t)
    tmpBnd(:,1) = s2.m1.star.bnd; %get results from linX
end
if ismember('m2',t)
    tmpBnd(:,2) = s2.m2.star.bnd; %get results from logX
end
[s1.star.bnd s1.star.method] = max(tmpBnd,[],2); %#ok

% Summarize the relevant results in the structure s1
for i1 = 1:r
    if s1.star.method(i1) == 1
        method = 'm1';
        trans = 'linX';
    else
        method = 'm2';
        trans = 'logX';
    end

    % display which was best if there were choices
    if length(t) > 1 && strcmpi(opt.display,'all')
        str = [blanks(4) trans ' transformation gave the best lower bound in oneNormFitXlb'];
        DClab.dcdispstr(str,opt.guiHandle,false)
    end

    %  s1.star.mult(i1,1) = s2.(best).star.mult(i1);
    if computeGap,
        s1.gap(i1,1) = min(s2.(method).pound(i1,:)) - s1.star.bnd(i1);
    end
end
