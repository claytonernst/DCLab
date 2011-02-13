classdef DCMetamodel < DClab.DCObject
    % model = DCMetamodel(RM,domain,criticalRange,trainingRange,opt)
    %
    % Domain is expected to be in the same variable order as the domain of RM.
    % range of domain2 is slightly larger. used for training

    properties
        parameterList = {};
        domainRange;
        responseTransformation;
        variableTransformations;
        parentNode;
        childNodes;
        model;
        fitInfo;
    end

    methods
        function model = DCMetamodel(RM,domain1,domain2,criticalRange,trainingRange,opt)


            %
            % Fields:
            % parameterList (just on root node to save mem)
            % domainRange
            % responseTransformation
            % variableTransformations (0 for nonactive, 1 for none, 2 for log10
            % parentNode
            % childNodes
            % model
            % modelError

            % TODO add error checking

            %Check that the model is defined on the requested domain.

            %modelDomain = RM.domain;
            %modelDomain = vertcat(modelDomain.range);
            %if ~issubset(domrng,modelDomain,1e-10,1e-10)
            %  error('The domain of the ResponseModel does not contain the requested domain')
            %end

            ni = nargin;

            switch ni
                case 0
                    return
                case {1,2}
                    error('Incorrect number of input arguments')
                case 3
                    criticalRange = [-inf inf];
                    trainingRange = [-inf inf];
                    opt = DClab.DCOptions;
                case 4
                    trainingRange = [-inf inf];
                    opt = DClab.DCOptions;
                case 5
                    opt = DClab.DCOptions;
                otherwise
                    error( nargchk(0,6,ni) );
            end

            n = length(domain1);

            [respTrans varTransNum varTransChar] = determineTransformations(RM,domain1);

            domrng = vertcat(domain1.range);

            model.parameterList = {domain1.name}';
            model.domainRange = domrng;
            model.responseTransformation = respTrans;
            model.variableTransformations = varTransNum; % use integer representation to save mem.
            model.parentNode = -1;
            model.childNodes = [];

            fitInfo.fittingTime = [];
            fitInfo.optimErrTime = [];
            fitInfo.nSamplePoints4Fit = [];
            fitInfo.nSamplePoints4Validation = [];
            fitInfo.errorType = 'absolute';
            fitInfo.averageErrorOnSamplePoints4Fit = [];
            fitInfo.averageErrorOnSamplePoints4Validation = [];
            fitInfo.peakErrorOnSamplePoints4Fit = [];
            fitInfo.peakErrorOnSamplePoints4Validation = [];
            fitInfo.peakErrorFromOptimization = [];
            fitInfo.optimFailed = [];
            fitInfo.xGivingErrMin = [];
            fitInfo.xGivingErrMax = [];

            type = 'svmd';
            na = sum(varTransNum~=0);
            N = opt.nPntsPerCoeff4Metamodel*(na+1)*(na+2)/2;
            t1 = cputime;

            [model.model xT yT filesUsed] = constructMetamodel(RM,domain2,respTrans,varTransChar,varTransNum,N,type,trainingRange,opt);
            fitInfo.fittingTime = cputime-t1;

            yHatT = evalPrivate(model.model,xT,respTrans,varTransNum);
            if strcmp(respTrans,'log10')
                errT = log10(yHatT)-log10(yT);
            else
                errT = yHatT-yT;
            end
            fitInfo.nSamplePoints4Fit = N;
            fitInfo.averageErrorOnSamplePoints4Fit = mean(abs(errT));
            fitInfo.peakErrorOnSamplePoints4Fit = [min(errT) max(errT)];
            %TODO add field for error on critical sample points

            %toss out training points that lie outside of domain1
            keep = true(size(yT));
            for i1 = 1:n
                keep = keep & domrng(i1,1) <= xT(i1,:) & xT(i1,:) <= domrng(i1,2);
            end
            xT = xT(:,keep);
            yT = yT(keep);
            errT = errT(keep);

            % Create validation points. These should be distributed in all variables.
            % Distribute nonactive variables according to a uniform transformation.

            Nv = 75*(na+1)*(na+2)/2;

            [xV yV] = loadSavedEvaluations(RM,domrng,varTransChar,Nv,filesUsed);

            Nsamples = length(yV);

            if Nsamples < Nv
                [newXV, newYV] = randomEvalWithSave(RM,domrng,varTransChar,Nv-Nsamples,opt.nComputer);
                xV = [xV newXV];
                yV = [yV newYV];
            end

            % Reject data that didn't lie in the critical range.
            if ~any(isinf(criticalRange))
                junkDataT = isnan(yT) | yT > criticalRange(2) | yT < criticalRange(1);
                errT = errT(~junkDataT);
                xT = xT(:,~junkDataT);

                junkDataV = isnan(yV) | yV > criticalRange(2) | yV < criticalRange(1);
                yV = yV(~junkDataV);
                xV = xV(:,~junkDataV);
            elseif isinf(criticalRange(1))
                junkDataT = isnan(yT) | yT > criticalRange(2);
                errT = errT(~junkDataT);
                xT = xT(:,~junkDataT);

                junkDataV = isnan(yV) | yV > criticalRange(2);
                yV = yV(~junkDataV);
                xV = xV(:,~junkDataV);
            elseif isinf(criticalRange(2))
                junkDataT = isnan(yT) | yT < criticalRange(1);
                errT = errT(~junkDataT);
                xT = xT(:,~junkDataT);

                junkDataV = isnan(yV) | yV < criticalRange(1);
                yV = yV(~junkDataV);
                xV = xV(:,~junkDataV);
            else
                junkDataT = isnan(yT);
                errT = erT(~junkDataT);
                xT = xT(:,~junkDataT);
                junkDataV = isnan(yV);
                yV = yV(~junkDataV);
                xV = xV(:,~junkDataV);
            end

            if length(yV) < 0.05*Nv
                disp('More than 95 %% of the validation data was rejected.')
                keyboard
            end

            % Determine the error on the validation points.
            yHatV = evalPrivate(model.model,xV,respTrans,varTransNum);
            if strcmp(respTrans,'log10')
                errV = log10(yHatV)-log10(yV);
            else
                errV = yHatV-yV;
            end

            fitInfo.nSamplePoints4Validation = Nv;
            fitInfo.averageErrorOnSamplePoints4Validation = mean(abs(errV));
            fitInfo.peakErrorOnSamplePoints4Validation = [min(errV) max(errV)];

            % Optimizes the fitting error.

            xAll = [xT xV];
            errAll = [errT errV];

            %if any(isnan(yV)) && ( isempty(criticalRange) || isempty(trainingRange) )
            %  error('Response Model give NaN outputs, but either criticalRange or trainingRange was undefined.')
            %end

            % Reject data that didn't lie in the critical range.
            %if ~isempty(criticalRange)
            %  junkData = isnan(yV) | yV > criticalRange(2) | yV < criticalRange(1);
            %  yV = yV(~junkData);
            %  xV = xV(:,~junkData);
            %  if isempty(yV)
            %    disp('No validation data in the critical range')
            %    keyboard
            %  end
            %end
            %
            % yHatV = evalPrivate(model.model,xV,respTrans,varTransNum);
            % errV = yHatV-yV;
            %
            % fitInfo.nSamplePoints4Validation = Nv;
            % fitInfo.averageErrorOnSamplePoints4Validation = mean(abs(errV));
            % fitInfo.peakErrorOnSamplePoints4Validation = [min(errV) max(errV)];
            %
            % % Optimize fitting error:
            %
            % % Reject data that didn't lie in the critical range.
            % if ~isempty(criticalRange)
            %   junkData = isnan(yT) | yT > criticalRange(2) | yT < criticalRange(1);
            %   errT = errT(~junkData);
            %   xT = xT(:,~junkData);
            % end
            %
            % xAll = [xT xV];
            % errAll = [errT errV];

            % TODO: Number of local searches for metamodel. This should be an option and can
            % be zero. Currently 0 will bomb
            n2Opt = opt.nLocalValidationSearches4Metamodel;

            if n2Opt < 1
                fitInfo.optimErrTime = 0;
                fitInfo.optimFailed = [false false];
                fitInfo.peakErrorFromOptimization = [0 0];
                fitInfo.xGivingErrMin = [];
                fitInfo.xGivingErrMax = [];
            else

                % Create a handle to the ResponseModel in coded variables.
                xCoded2xSTR = '';
                for i1 = 1:n
                    if varTransNum(i1) == 2
                        xCoded2xSTR = [xCoded2xSTR '10.^( 0.5*log10(domrng(' num2str(i1) ',2)/domrng(' num2str(i1) ',1))*(x(' num2str(i1) ',:) + 1) + log10(domrng(' num2str(i1) ',1)) );']; %#ok
                    else
                        xCoded2xSTR = [xCoded2xSTR '0.5*(diff(domrng(' num2str(i1) ',:)))*(x(' num2str(i1) ',:) + 1) + domrng(' num2str(i1) ',1);']; %#ok
                    end
                end
                xCoded2xSTR(end) = ']';
                xCoded2xSTR = ['[' xCoded2xSTR];

                eval(['M_A = @(x) rapideval(RM,' xCoded2xSTR ');']);

                % Handle to a modified model that is everywhere defined that can be used
                % with fmincon to assess the fitting error without difficulty. If the
                % ResponseModel evaluates to NaN, that value will be replaced with
                % criticalRange +/- 50;

                % TODO: we need to verify that the user's criticalRange(1) is posivite if
                % we decide upon a log10 transformation for Y.
                if ~isinf(criticalRange(2))
                    M_A_clipped = @(x) clipy(x,M_A,criticalRange(2)+50);
                elseif ~isinf(criticalRange(1))
                    % Manipulate to avoid problems with log10 transform.
                    if strcmp(respTrans,'log10')
                        M_A_clipped = @(x) clipy(x,M_A,0.5*criticalRange(1));
                    else
                        M_A_clipped = @(x) clipy(x,M_A,criticalRange(1)-50);
                    end
                else
                    M_A_clipped = M_A; %assume the model is everywhere defined if no crit range given.
                    if any(junkDataT) || any(junkDataV)
                        disp('NaN model outputs detected, but the critical range is infinite. Error optimization may encounter difficulties.')
                    end
                end

                % Create handle to metamodel in coded variables.
                eval(['meta = @(x) evalPrivate(model.model,' xCoded2xSTR ',respTrans,varTransNum);']);

                fmopt = optimset('fmincon');
                fmopt = optimset(fmopt,'display','off','RelLineSrchBnd',0.05,'RelLineSrchBndDuration',20,'DiffMinChange',1e-4,'maxfunevals',600,'largescale','off');
                fmopt2 = optimset(fmopt,'RelLineSrchBnd',0.005); %We reduce the linesearch step size if we end up infeasible. Since we start feasible, this should help...

                eFlagH = zeros(1,n2Opt);
                eFlagL = zeros(1,n2Opt);
                ERRH = zeros(1,n2Opt);
                ERRL = zeros(1,n2Opt);
                XH = zeros(n,n2Opt);
                XL = zeros(n,n2Opt);

                % Initial seeds:
                xHseed = x2Coded(DClab.findExtremeX(xAll,errAll,n2Opt,'max'),domrng,varTransNum);
                xLseed = x2Coded(DClab.findExtremeX(xAll,errAll,n2Opt,'min'),domrng,varTransNum);

                % Objective functions
                if strcmp(respTrans,'log10')
                    objfunH = @(x) log10(M_A_clipped(x))-log10(meta(x));
                    objfunL = @(x) log10(meta(x))-log10(M_A_clipped(x));
                else
                    objfunH = @(x) M_A_clipped(x)-meta(x);
                    objfunL = @(x) meta(x)-M_A_clipped(x);
                end

                % Constrant function
                if all(isinf(criticalRange))
                    confun = [];
                elseif isinf(criticalRange(1))
                    confun = @(x) criticalRangeUConstFun(x,M_A_clipped,criticalRange(2));
                elseif isinf(criticalRange(2))
                    confun = @(x) criticalRangeLConstFun(x,M_A_clipped,criticalRange(1));
                else
                    confun = @(x) criticalRangeConstFun(x,M_A_clipped,criticalRange);
                end

                t1 = cputime;
                for i1 = 1:n2Opt
                    [XH(:,i1) ERRH(i1) eFlagH(i1)] = fmincon(objfunH,xHseed(:,i1),[],[],[],[],-ones(n,1),ones(n,1),confun,fmopt);
                    [XL(:,i1) ERRL(i1) eFlagL(i1)] = fmincon(objfunL,xLseed(:,i1),[],[],[],[],-ones(n,1),ones(n,1),confun,fmopt);

                    % If the optimization exited with an infeasible solution, try again
                    % with a smaller line search stepsize.
                    if eFlagH(i1) == -2
                        [XH(:,i1) ERRH(i1) eFlagH(i1)] = fmincon(objfunH,xHseed(:,i1),[],[],[],[],-ones(n,1),ones(n,1),confun,fmopt2);
                    end
                    if eFlagL(i1) == -2
                        [XL(:,i1) ERRL(i1) eFlagL(i1)] = fmincon(objfunL,xLseed(:,i1),[],[],[],[],-ones(n,1),ones(n,1),confun,fmopt2);
                    end

                    % If maxfunevals occurred but the final value is (essentially) feasible, keep it.

                    if eFlagH(i1) == 0
                        if isempty(confun)
                            eFlagH(i1) = 1;
                        else
                            if  max(confun(XH(:,i1))) <= 1e-4
                                eFlagH(i1) = 1;
                            end
                        end
                    end
                    if eFlagL(i1) == 0
                        if isempty(confun)
                            eFlagL(i1) = 1;
                        else
                            if  max(confun(XL(:,i1))) <= 1e-4
                                eFlagL(i1) = 1;
                            end
                        end
                    end
                end
                fitInfo.optimErrTime = cputime-t1;

                ERRH = -ERRH;
                ERRH = ERRH(eFlagH>0);
                ERRL = ERRL(eFlagL>0);
                XH = XH(:,eFlagH>0);
                XL = XL(:,eFlagL>0);

                if isempty(ERRL)
                    fitInfo.optimFailed(1,1) = true;
                    fitInfo.peakErrorFromOptimization(1,1) = min(errAll);
                    fitInfo.xGivingErrMin = xLseed(:,1);
                else
                    fitInfo.optimFailed(1,1) = false;
                    [fitInfo.peakErrorFromOptimization(1,1) idx] = min(ERRL);
                    fitInfo.xGivingErrMin = XL(:,idx);
                end

                if isempty(ERRH)
                    fitInfo.optimFailed(1,2) = true;
                    fitInfo.peakErrorFromOptimization(1,2) = max(errAll);
                    fitInfo.xGivingErrMax = xHseed(:,1);
                else
                    fitInfo.optimFailed(1,2) = false;
                    [fitInfo.peakErrorFromOptimization(1,2) idx] = max(ERRH);
                    fitInfo.xGivingErrMax = XH(:,idx);
                end
                fitInfo.xGivingErrMin = coded2X(fitInfo.xGivingErrMin,domrng,varTransNum);
                fitInfo.xGivingErrMax = coded2X(fitInfo.xGivingErrMax,domrng,varTransNum);

            end

            model.fitInfo = fitInfo;
            %model.userData = [criticalRange; trainingRange];
        end
        
        function out = vertcat(varargin)
            %TODO: is this needed?
            out = builtin('vertcat',varargin{:});
        end
        
        function bool = isempty(obj)
            bool=false;
            if isempty(obj(1).model)
                bool = true;
            end
        end
        function [list,sz] = displayProps(obj)
            list = [];
            sz = sprintf('%d-by-%d',size(obj));
        end
        
    end %public methods
end %classdef

% Local functions

function xC = x2Coded(x,domrng,varTransNum)

xC = zeros(size(x));
for i1 = 1:size(domrng,1)
    if varTransNum(i1) == 2;
        xC(i1,:) = 2/log10(domrng(i1,2)/domrng(i1,1))*(log10(x(i1,:)) - log10(domrng(i1,1))) - 1;
    else
        xC(i1,:) = 2*(1/diff(domrng(i1,:)))*(x(i1,:) - domrng(i1,1)) - 1;
    end
end
end

function x = coded2X(xC,domrng,varTransNum)

x = zeros(size(xC));
for i1 = 1:size(domrng,1)
    if varTransNum(i1) == 2;
        x(i1,:) = 10.^( 0.5*log10(domrng(i1,2)/domrng(i1,1))*(xC(i1,:) + 1) + log10(domrng(i1,1)) );
    else
        x(i1,:) = 0.5*(diff(domrng(i1,:)))*(xC(i1,:) + 1) + domrng(i1,1);
    end
end
end

function [cval eval] = criticalRangeConstFun(x,fx,critRng)
y = fx(x);
cval = [y-critRng(2); critRng(1)-y];
eval = 0;
end

function [cval eval] = criticalRangeUConstFun(x,fx,critRngU)
y = fx(x);
cval = y-critRngU;
eval = 0;
end

function [cval eval] = criticalRangeLConstFun(x,fx,critRngL)
y = fx(x);
cval = critRngL-y;
eval = 0;
end