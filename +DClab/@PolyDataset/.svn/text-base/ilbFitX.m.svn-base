function [s1 s2] = ilbFitX(PD,opt,norm,trgWeights,nodeIndices,computeGap,startPt)
% [s1 s2] = ilbFitX(PD,opt,norm,trgWeights,nodeIndices,computeGap,startPt)
%
% ilbFitX is used to fit model parameters to data using a specified norm cost
% on the residual vector. It performs up to four different approximations
% to the optimization
%
% min \sum_e w_e*|M_e(x)-d_e|   subj to: x in some subdivision of H (eq1)
%
% Inputs:
%   PD: A PolyDataset object.
%   opt: A DCOptions object or structure with fields .display,
%     .sedumiParsEps, .tolFun, .tolCon, .nRestart, .constraints, and .guihandle
%   norm: either 'one','two',or 'inf'.
%   trgWeights: A mx1 vector of weights w_e in the above optimization.
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
% Description for one norm objective
%
% The optimization in (eq1) is equivalent to
%
% min \sum \gamma_e   subj. to: w_e*|M_e(x)-d_e| \leq \gamma_e (eq2)
%                               x in the nodeIndices subdiv of H
%
% where M_e are the true models present in the ModelAssertions of the
% PolyDataset PD. If output uncertainty is present in M_e,
% mean([M_e(x)+outputUnc(1) M_e(x)+outputUnc(2)]) replaces M_e in the above
% expression. See the 2nd output of ModelAssertion/eval or the
% ModelAssertion constructor for more information. If the uncertainty of
% the e'th ExperimentAssertion is asymmetric, d_e is the above expression is
% not the data value in the assertion, but rather is centered, i.e.,
% d_e = EA.data - mean(EA.uncVect). Of course, if the uncertainty is
% symmetric (typical), d_e = EA.data.
%
% w_e are weights, which are supplied by the input trgWeights.
%
% Each approximation to the optimization problem (eq2) takes the form of a
% nonconvex quadratically constrained quadratic program (nqcqp). If
% computeGap = true, the S-procedure is used to determine a global lower
% bound on the optimum for each nqcqp. If the gap between this lower bound
% and the solution found by local search is zero, the local search found
% the global minimum of the nqcqp. However, since each nqcqp is just an
% approximation to eq1 (except in the special case where M_e is quadratic
% and has no output uncertainty), it is not
% clear that the gap being zero has any practical meaning. Consequently,
% computeGap = false is the default behavior.
%
% Each of the four methods used to find suboptimal solutions to eq1 seek to
% find good vectors x at which to evaluate the objective function in eq1
% (the one norm of the weighted residual vector). Their behavior is
% described below. Since each of these is a local search, they are each
% restarted opt.nRestarts times, where opt is a DCOptions object contained
% in the PD.as follows.
%
% method1: here surrogate models S_e quadratic in x are given that
%   approximate M_e such that
%       fiterr_e(1)+M_e(x) <= S_e(x) <= fiterr_e(2)+M_e(x)
%
%   we then solve the optimization
%   min \sum \gamma_e   subj. to:
%    -\gamma_e + w_e*fiterr_e(2) <= w_e*(S_e(x)-d_e) <= w_e*fiterr_e(1) + \gamma_e   (m1)
%    x in the nodeIndices subdiv of H
%
%   notice the optimal gamma's for this problem are alway larger than the optimal gammas for eq2.
%   we then evaluate the objective function in eq1 at the best x from (m1)
%
% method2: here surrogate models S_e quadratic in LOG10(x) are given that
%   approximate M_e such that
%       fiterr_e(1)+M_e(x) <= S_e(log10(x)) <= fiterr_e(2)+M_e(x)
%
%   we then solve the optimization
%   min \sum \gamma_e   subj. to:
%    -\gamma_e + w_e*fiterr_e(2) <= w_e*(S_e(x')-d_e) <= w_e*fiterr_e(1) + \gamma_e   (m2)
%    where x' = log10(x) and x in the nodeIndices subdiv of H
%
%   notice the optimal gamma's for this problem are alway larger than the optimal gammas for eq2.
%   we then evaluate the objective function in eq1 at the best x from (m2)
%
% method3: here we have the same surrogate models S_e as in method1, but
%   ignore the fitting error in formulating the optimization. I.e., we solve
%   min \sum \gamma_e   subj. to:
%    -\gamma_e <= w_e*(S_e(x)-d_e) <= \gamma_e   (m3)
%    x in the nodeIndices subdiv of H
%
%   notice the optimal gamma's may be bigger or smaller than the optimal gammas for eq2.
%   we then evaluate the objective function in eq1 at the best x from (m3)
%
% method4: here we have the same surrogate models S_e as in method2, but
%   again ignore the fitting error in formulating the optimization. I.e., we solve
%   min \sum \gamma_e   subj. to:
%    -\gamma_e <= w_e*(S_e(x')-d_e) <= \gamma_e   (m3)
%    where x' = log10(x) and x in the nodeIndices subdiv of H
%
%   notice the optimal gamma's may be bigger or smaller than the optimal gammas for eq2.
%   we then evaluate the objective function in eq1 at the best x from (m4)
%
% It is unclear which of these methods is superior (if in fact one is). The
% output s2 of this function contains the solutions to all these problems (including restarts).
% The output s1 contains only the best solution (in terms of acheiving the
% lowest cost in eq1).
%
% Description for two norm objective:
%
% The optimization in (eq1) is equivalent to
%
% min \sqrt\sum \gamma_e^2   subj. to: w_e*|M_e(x)-d_e| \leq \gamma_e (eq2)
%                               x in the nodeIndices subdiv of H
%
% where M_e are the true models present in the ModelAssertions of the
% PolyDataset PD. If output uncertainty is present in M_e,
% mean([M_e(x)+outputUnc(1) M_e(x)+outputUnc(2)]) replaces M_e in the above
% expression. See the 2nd output of ModelAssertion/eval or the
% ModelAssertion constructor for more information. If the uncertainty of
% the e'th ExperimentAssertion is asymmetric, d_e is the above expression is
% not the data value in the assertion, but rather is centered, i.e.,
% d_e = EA.data - mean(EA.uncVect). Of course, if the uncertainty is
% symmetric (typical), d_e = EA.data.
%
% w_e are weights, which may optionally be supplied by the input
% trgWeights. The default is w_e = inv(diff(EA.uncVect)/2), where EA is the
% e'th ExperimentAssertion in the PolyDataset PD.
%
% Each approximation to the optimization problem (eq2) takes the form of a
% nonconvex quadratically constrained quadratic program (nqcqp). If
% computeGap = true, the S-procedure is used to determine a global lower
% bound on the optimum for each nqcqp. If the gap between this lower bound
% and the solution found by local search is zero, the local search found
% the global minimum of the nqcqp. However, since each nqcqp is just an
% approximation to eq1 (except in the special case where M_e is quadratic
% and has no output uncertainty), it is not
% clear that the gap being zero has any practical meaning. Consequently,
% computeGap = false is the default behavior.
%
% Each of the four methods used to find suboptimal solutions to eq1 seek to
% find good vectors x at which to evaluate the objective function in eq1
% (the two norm of the weighted residual vector). Their behavior is
% described below. Since each of these is a local search, they are each
% restarted opt.nRestarts times, where opt is a DCOptions object contained
% in the PD.as follows.
%
% method1: here surrogate models S_e quadratic in x are given that
%   approximate M_e such that
%       fiterr_e(1)+M_e(x) <= S_e(x) <= fiterr_e(2)+M_e(x)
%
%   we then solve the optimization
%   min \sqrt\sum \gamma_e^2   subj. to:
%    -\gamma_e + w_e*fiterr_e(2) <= w_e*(S_e(x)-d_e) <= w_e*fiterr_e(1) + \gamma_e   (m1)
%    x in the nodeIndices subdiv of H
%
%   notice the optimal gamma's for this problem are alway larger than the optimal gammas for eq2.
%   we then evaluate the objective function in eq1 at the best x from (m1)
%
% method2: here surrogate models S_e quadratic in LOG10(x) are given that
%   approximate M_e such that
%       fiterr_e(1)+M_e(x) <= S_e(log10(x)) <= fiterr_e(2)+M_e(x)
%
%   we then solve the optimization
%   min \sqrt\sum \gamma_e^2   subj. to:
%    -\gamma_e + w_e*fiterr_e(2) <= w_e*(S_e(x')-d_e) <= w_e*fiterr_e(1) + \gamma_e   (m2)
%    where x' = log10(x) and x in the nodeIndices subdiv of H
%
%   notice the optimal gamma's for this problem are alway larger than the optimal gammas for eq2.
%   we then evaluate the objective function in eq1 at the best x from (m2)
%
% method3: here we have the same surrogate models S_e as in method1, but
%   ignore the fitting error in formulating the optimization. I.e., we solve
%   min \sqrt\sum \gamma_e^2   subj. to:
%    -\gamma_e <= w_e*(S_e(x)-d_e) <= \gamma_e   (m3)
%    x in the nodeIndices subdiv of H
%
%   notice the optimal gamma's may be bigger or smaller than the optimal gammas for eq2.
%   we then evaluate the objective function in eq1 at the best x from (m3)
%
% method4: here we have the same surrogate models S_e as in method2, but
%   again ignore the fitting error in formulating the optimization. I.e., we solve
%   min \sqrt\sum \gamma_e^2   subj. to:
%    -\gamma_e <= w_e*(S_e(x')-d_e) <= \gamma_e   (m3)
%    where x' = log10(x) and x in the nodeIndices subdiv of H
%
%   notice the optimal gamma's may be bigger or smaller than the optimal gammas for eq2.
%   we then evaluate the objective function in eq1 at the best x from (m4)
%
% It is unclear which of these methods is superior (if in fact one is). The
% output s2 of this function contains the solutions to all these problems (including restarts).
% The output s1 contains only the best solution (in terms of acheiving the
% lowest cost in eq1).
%
% Description for inf-norm objective:
%
% The optimization in (eq1) is equivalent to
%
% min \gamma   subj. to: w_e*|M_e(x)-d_e| \leq \gamma (eq2)
%                        x in the nodeIndices subdiv of H
%
% where M_e are the true models present in the ModelAssertions of the
% PolyDataset PD. If output uncertainty is present in M_e,
% mean([M_e(x)+outputUnc(1) M_e(x)+outputUnc(2)]) replaces M_e in the above
% expression. See the 2nd output of ModelAssertion/eval or the
% ModelAssertion constructor for more information. If the uncertainty of
% the e'th ExperimentAssertion is asymmetric, d_e is the above expression is
% not the data value in the assertion, but rather is centered, i.e.,
% d_e = EA.data - mean(EA.uncVect). Of course, if the uncertainty is
% symmetric (typical), d_e = EA.data.
%
% w_e are weights, which may optionally be supplied by the input
% trgWeights. The default is w_e = inv(diff(EA.uncVect)/2), where EA is the
% e'th ExperimentAssertion in the PolyDataset PD.
%
% Each approximation to the optimization problem (eq2) takes the form of a
% nonconvex quadratically constrained quadratic program (nqcqp). If
% computeGap = true, the S-procedure is used to determine a global lower
% bound on the optimum for each nqcqp. If the gap between this lower bound
% and the solution found by local search is zero, the local search found
% the global minimum of the nqcqp. However, since each nqcqp is just an
% approximation to eq1 (except in the special case where M_e is quadratic
% and has no output uncertainty), it is not
% clear that the gap being zero has any practical meaning. Consequently,
% computeGap = false is the default behavior.
%
% Each of the four methods used to find suboptimal solutions to eq1 seek to
% find good vectors x at which to evaluate the objective function in eq1
% (the inf norm of the weighted residual vector). Their behavior is
% described below. Since each of these is a local search, they are each
% restarted opt.nRestarts times, where opt is a DCOptions object contained
% in the PD.as follows.
%
% method1: here surrogate models S_e quadratic in x are given that
%   approximate M_e such that
%       fiterr_e(1)+M_e(x) <= S_e(x) <= fiterr_e(2)+M_e(x)
%
%   we then solve the optimization
%   min \gamma   subj. to:
%    -\gamma + w_e*fiterr_e(2) <= w_e*(S_e(x)-d_e) <= w_e*fiterr_e(1) + \gamma   (m1)
%    x in the nodeIndices subdiv of H
%
%   notice the optimal gamma's for this problem are alway larger than the optimal gammas for eq2.
%   we then evaluate the objective function in eq1 at the best x from (m1)
%
% method2: here surrogate models S_e quadratic in LOG10(x) are given that
%   approximate M_e such that
%       fiterr_e(1)+M_e(x) <= S_e(log10(x)) <= fiterr_e(2)+M_e(x)
%
%   we then solve the optimization
%   min \gamma   subj. to:
%    -\gamma + w_e*fiterr_e(2) <= w_e*(S_e(x')-d_e) <= w_e*fiterr_e(1) + \gamma   (m2)
%    where x' = log10(x) and x in the nodeIndices subdiv of H
%
%   notice the optimal gamma's for this problem are alway larger than the optimal gammas for eq2.
%   we then evaluate the objective function in eq1 at the best x from (m2)
%
% method3: here we have the same surrogate models S_e as in method1, but
%   ignore the fitting error in formulating the optimization. I.e., we solve
%   min \gamma   subj. to:
%    -\gamma <= w_e*(S_e(x)-d_e) <= \gamma   (m3)
%    x in the nodeIndices subdiv of H
%
%   notice the optimal gamma's may be bigger or smaller than the optimal gammas for eq2.
%   we then evaluate the objective function in eq1 at the best x from (m3)
%
% method4: here we have the same surrogate models S_e as in method2, but
%   again ignore the fitting error in formulating the optimization. I.e., we solve
%   min \sum \gamma   subj. to:
%    -\gamma <= w_e*(S_e(x')-d_e) <= \gamma   (m3)
%    where x' = log10(x) and x in the nodeIndices subdiv of H
%
%   notice the optimal gamma's may be bigger or smaller than the optimal gammas for eq2.
%   we then evaluate the objective function in eq1 at the best x from (m4)
%
% It is unclear which of these methods is superior (if in fact one is). The
% output s2 of this function contains the solutions to all these problems (including restarts).
% The output s1 contains only the best solution (in terms of acheiving the
% lowest cost in eq1).
%
% Description of the structure array s2
%
% s2.m1.star.bnd contains an rx1 vector of lower bounds from the
%                s-proc (if computeGap = true), where r=length(nodeIndices)
% s2.m1.star.mult contains an rx1 cell array of vectors containing
%                 the corresponding multipliers from the s-proc
%                 optimization. For consistency with the object PD,
%                 if a parameter of PD was not a variable in the
%                 optimization (due to not being an active
%                 parameter) a "multipler" is still returned
%                 corresponding to that parameter. It is just set
%                 to zero. It is not clear that these have any meaning.
% s2.m1.pound.bnd contains an rxk matrix of upper bounds to eq1. Each
%                 component is determined by evaluating the objective
%                 in eq1 at the local search solutions of the
%                 nqcqps. Here k = 1+opt.nRestart.
% s2.m1.pound.xfeas contains an rxk cell array of parameter vectors
%                 corresponding to the local minimum found by local
%                 search. These values are nontransformed, i.e.,
%                 they are in whatever units the original parameter
%                 assertions were provided in. If a parameter of PD
%                 was not a variable in the optimization (due to
%                 not being an active parameter), the "nominal
%                 value" of that parameter is inserted in the
%                 proper location in s2.m1.pound.xfeas.
% s2.m1.pound.mult is similar to ...star.mult, it is just the
%                 multipliers from the local search. It is not
%                 clear that these are useful in the nonquadratic
%                 case, because they have nothing to do with the
%                 optimization we want to solve (eq1).
%
% s2.m2, s2.m3, and s2.m4 contain the above info from method2,
% method3, and method4 respectively.
%
% Description of the structure array s1
%
% s1.pound.bnd    contains an rx1 vector of the best overall (from all
%                 the restarts and all the methods) upper bound on the
%                 solution to eq1 for each of the r subdivisions.
% s1.pound.method contains an rx1 vector of 1,2,3,or 4 reporting
%                 which method produced the best answer
% s1.pound.mult   contains an rx1 cell of the multipliers
%                 corresponding to the best overall optimizations.
% s1.pound.xfeas  contains the parameter vector that produced the
%                 best bound in s1.pound.bound
% s1.pound.resid  contains an rx1 cell array of residual vectors
%                 M(x)-d evaluated at x=s1.pound.xfeas
% s1.pound.gap if computeGap=true, this rx1 vector contains the
%                 difference between s1.pound.bnd and lowerbound
%                 solutions of the ncqcps that produced the bestx's in
%                 local search. In the nonquad case this seems
%                 meaningless...
%
% See also PolyDataset, ParameterOptimization, PolyDataset/olbFitX

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

%we can chose to solve the problem by four methods:
% m1: use inner approximation in linX
% m2: use inner approximation in logX
% m3: use crude approximation in linX
% m4: use crude approximation in logX

%currently we only can handle linX transformations, and we require their
%existence
t = {'m1'};
if ~isPerfectFit(PD)
    t = [t;{'m3'}]; %m3 and m4 are redundant if there is no fitting error
end

switch norm
    case 'one'
        if any(ismember({'m1';'m2'},t))
            [optProb.m1, optProb.m2, activeIdx] = optim1norm_nqcqp(PD,opt,trgWeights,'inner',nodeIndices);
        end
        if any(ismember({'m3';'m4'},t))
            [optProb.m3, optProb.m4, activeIdx] = optim1norm_nqcqp(PD,opt,trgWeights,'crude',nodeIndices);
        end
    case 'two'
        if any(ismember({'m1';'m2'},t))
            [optProb.m1, optProb.m2, activeIdx] = optim2norm_nqcqp(PD,opt,trgWeights,'inner',nodeIndices);
        end
        if any(ismember({'m3';'m4'},t))
            [optProb.m3, optProb.m4, activeIdx] = optim2norm_nqcqp(PD,opt,trgWeights,'crude',nodeIndices);
        end
    case 'inf'
        if any(ismember({'m1';'m2'},t))
            [optProb.m1, optProb.m2, activeIdx] = optimInfnorm_nqcqp(PD,opt,trgWeights,'inner',nodeIndices);
        end
        if any(ismember({'m3';'m4'},t))
            [optProb.m3, optProb.m4, activeIdx] = optimInfnorm_nqcqp(PD,opt,trgWeights,'crude',nodeIndices);
        end
    otherwise
        error('Internal inconsistency, check code')
end

s2.paramList = PD.parameterList;
for i1 = 1:length(t)
    if computeGap
        bnd = zeros(r,1);
        mult = cell(r,1);
        Q = cell(r,1);
        for i2=1:r
            [bnd(i2),mult{i2},Q{i2}] = lowerBound(optProb.(t{i1}){i2},opt);
        end
        [bnd,mi] = min(bnd);
        mult = mult{mi};
        Q = Q{mi};
        switch norm
            case {'one','inf'}
                s2.(t{i1}).star.bnd = bnd;
            case 'two'
                bnd = max(bnd,0);
                s2.(t{i1}).star.bnd = sqrt(bnd);
            otherwise
                error('Internal inconsistency, condition should never occur')
        end
    end
    if ~exist('Q','var') %see if we can use the s-proc results for a seed point
        Q = [];
    end

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
            str = [blanks(4) 'fmincon local search try ' num2str(i2) ' of ' num2str(opt.nRestart+1)];
            DClab.dcdispstr(str,opt.guiHandle,false)
        else
            drawnow
        end

        r = length(xinit);
        fval = zeros(r,1);
        xfeas = cell(r,1);
        mult = cell(r,1);
        for i3=1:r
            [fval(i3),xfeas{i3},mult{i3}] = qpub_FMINCON(optProb.(t{i1}){i3},xinit{i3},opt);
        end

        %Remove last m xfeas and multiplier that correspond to
        %gammas. We don't want to see these in our final output.
        switch norm
            case {'one','two'}
                m = PD.nPairs;
                for i3 = 1:r
                    xfeas{i3}(end-m+1:end) = [];
                    mult{i3}.lower(end-m+1:end) = [];
                    mult{i3}.upper(end-m+1:end) = [];
                end
            case 'inf'
                for i3 = 1:r
                    xfeas{i3}(end) = [];
                    mult{i3}.lower(end) = [];
                    mult{i3}.upper(end) = [];
                end
            otherwise
                error('Internal inconsistency, condition should never occur')
        end

        s2.(t{i1}).pound.bnd(:,i2) = inf(r,1); %overwrite in a few lines

        %find best among different leaves
%         [fval,bestidx] = min(fval);
%         xfeas = xfeas{bestidx};
%         mult = mult{bestidx};

        %Convert xfeas to original units, and if any Dataset variables
        %were not optimization variables, insert their nominal values
        %into the returned output.
        s2.(t{i1}).pound.xfeas(:,i2) = cell(r,1);
        s2.(t{i1}).pound.mult(:,i2) = cell(r,1);
        for i3=1:r
            s2.(t{i1}).pound.xfeas{i3,i2} =  DClab.transformAndAddNominals(xfeas{i3},PD,activeIdx{i3},nodeIndices);
            %Add multipliers of 0 to any dataset parameter that where not
            %optimization variables.
            s2.(t{i1}).pound.mult{i3,i2} = DClab.addMultZeros(mult{i3},PD,opt,activeIdx{i3},nodeIndices);
        end
    end

    %searching on the surrogates found a (hopefully) good point.
    %now eval the true models to actually get an upper bound on the
    %solution to eq1 above.
    [s2.(t{i1}).pound.bnd s2.(t{i1}).pound.resid] = evalCost(PD,s2.(t{i1}).pound.xfeas,norm,trgWeights);
end

% Determine which method gives better results (small is better). The
% contents of the structure s2 are, for example, s2.m1.pound.bnd is a
% r-by-nRestarts+1 array, s2.m1.pound.xfeas is a cell array of the same
% dimensions, and so is s2.m1.pound.mult.
tmpBnd = repmat(inf,r,4);
tmpIdx = zeros(r,4);
if ismember('m1',t)
    [tmpBnd(:,1) tmpIdx(:,1)] = min(s2.m1.pound.bnd,[],2); %get best of the method1 retrys for each partion
end
if ismember('m2',t)
    [tmpBnd(:,2) tmpIdx(:,2)] = min(s2.m2.pound.bnd,[],2); %get best of the method2 retrys for each partion
end
if ismember('m3',t)
    [tmpBnd(:,3) tmpIdx(:,3)] = min(s2.m3.pound.bnd,[],2);
end
if ismember('m4',t)
    [tmpBnd(:,4) tmpIdx(:,4)] = min(s2.m4.pound.bnd,[],2);
end
[s1.pound.bnd s1.pound.method] = min(tmpBnd,[],2); %#ok

% Summarize the relevant results in the structure s1
for i1 = 1:r
    best = ['m' num2str(s1.pound.method(i1))];
    idx = tmpIdx(i1,s1.pound.method(i1));

    s1.pound.mult(i1,1) = s2.(best).pound.mult(i1,idx);
    s1.pound.xfeas(i1,1) = s2.(best).pound.xfeas(i1,idx);
    s1.pound.resid(i1,1) = s2.(best).pound.resid(i1,idx);
    if computeGap
        s1.gap(i1,1) = s1.pound.bnd(i1) - s2.(best).star.bnd(i1);
    end
end

function [ibnd resid] = evalCost(PD,pt,norm,trgWeights)
%function to evaluate the objective function in eq1 above at a
%given x value.
%
%
% PD is a PolyDataset object, pt is a rxk cell array, trgWeights is
% an mx1 double array.
%
% ibnd is a matrix, resid is a cell array.
if isnumeric(pt)
    pt = {pt};
end

ibnd = zeros(size(pt));
resid = cell(size(pt));
for i1 = 1:size(pt,1) %over divisions
    for i2 = 1:size(pt,2) %over restarts
        resid{i1,i2} = zeros(PD.nPairs,1);
        tempPt = pt{i1,i2};
        %now use this point to lower bound the consistency measure
        for i3 = 1:PD.nPairs
            [trash yInt] = eval(PD.ModelAndObservationPair(i3).ResponseModel,tempPt,PD.parameterList); %#ok
            u = PD.ModelAndObservationPair(i3).observationUncertaintyPlusMinus;
            d = PD.ModelAndObservationPair(i3).observedValue;

            %optimization is min sum gammai subj to wi*|Mi(x)-di| <= gammai
            d = d-mean(u); %center d
            y = mean(yInt);
            resid{i1,i2}(i3) = trgWeights(i3)*abs(y-d);
        end
        switch norm
            case 'one'
                ibnd(i1,i2) = sum(resid{i1,i2});
            case 'two'
                ibnd(i1,i2) = sqrt(sum(resid{i1,i2}.^2));
            case 'inf'
                ibnd(i1,i2) = max(resid{i1,i2});
            otherwise
                error('Internal inconsistency, condition should never occur')
        end
    end
end
