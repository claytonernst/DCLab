classdef DCSurface < DClab.DCObject
    %DCSURFACE Construct a DCSurface object
    %
    %   %%%% HELP IS OUTDATED %%%%
    %
    %   Only the first two inputs are required. You can supply [] for any of
    %   the optional inputs to indicate the constructor should use the default
    %   value of the corresponding input.
    %
    %   SURF = DCSURFACE creates an empty object
    %
    %   SURF = DCSURFACE(RESPMODEL,DOMAIN) creates a DCSurface object SURF from
    %   the ResponseModel object RESPMODEL over the domain defined by the nx1
    %   structure DOMAIN. DOMAIN should have fields .name and .range so that
    %   DOMAIN(i).name is the name of the ith parameter and DOMAIN(i).range is
    %   a 1x2 containing a lower and upper bound on the ith dimension of the
    %   rectangular domain the surface is to be fit over. The components of
    %   DOMAIN do not need to be (and generally aren't) in the same order as
    %   the parameters occur in RESPMODEL.
    %
    %   SURF = DCSURFACE(RESPMODEL,DOMAIN,RESPTRANS) allows you to specify a
    %   transformation that will be applied to the output of RESPMODEL before
    %   fitting the surface. Valid values are 'none' and 'log10'. The default
    %   is 'none'.
    %
    %   SURF = DCSURFACE(RESPMODEL,DOMAIN,RESPTRANS,VARTRANS) allows you to
    %   specify a transformation for each of the dataset parameters that will
    %   be applied before fitting the surface. The points the model is
    %   evaluated at will be distributed in the transformed coordinates.
    %   VARTRANS should be a nx1 cell array with each cell containing 'none' or
    %   'log10'. The transformations should be supplied in the order of the
    %   dataset parameters. The default is repmat({'none'},n,1).
    %
    %   SURF = DCSURFACE(RESPMODEL,DOMAIN,RESPTRANS,VARTRANS,MOPAIRIDX) allows
    %   you to indicate which element of the Datset's array of
    %   ModelAndObservationPair objects RESPMODEL is taken from. I.e.,
    %   DCDataset.ModelAndObservationPair(MOPAIRIDX).ResponseModel should equal
    %   RESPMODEL. Through this index in the only way a DCSurface knows from
    %   whenst it came.
    %
    %   SURF = DCSURFACE(RESPMODEL,DOMAIN,RESPTRANS,VARTRANS,MOPAIRIDX,OPTS)
    %   uses the fit settings from the DCOptions object OPTS. These are those
    %   given by the 2nd output of DCOptions/decompose.
    %
    %   SURF = DCSURFACE(INTEGER) initializes an INTEGER-by-1 empty object.
    %
    %   Note, the fitting error is lower and upper bounds on S(x) - M(x).
    %
    %   See also ResponseModel, DCSurface.buildSurface
    %
    %   Constructor reference page in Help browser
    %      <a href="matlab:web(fullfile(fileparts(which('DCsetup')), ...
    %      'docs', 'html', 'functions', 'DCSurface.html'), ...
    %      '-helpbrowser')">DCSurface constructor</a>
    %
    %   Object reference page in Help browser
    %      <a href="matlab:web(fullfile(fileparts(which('DCsetup')), ...
    %      'docs', 'html', 'obj_dset.html'), ...
    %      '-helpbrowser')">Broken link</a>
    
    %TODO: Help update
    
    properties
        type = '';
        responseTransformation = '';
        surrogateModel = [];
        surrogateFitInfo = [];
        trueModelRangeEstimate = [];
        MOPairIdx = [];
        isPredictor = [];
        status = -1;
    end %properties
    
    properties (Access=private)
        api;
    end
    
    properties (Dependent)
        activeParameterIndex
        activeParameterTransformation
        peakError
    end %dependent properties
    
    methods
        
        %constructor
        function Surf = DCSurface(varargin)
            
            ni = nargin;
            switch ni
                case 0
                    %defaults
                    return
                case 1
                    sz = varargin{1};
                    assert(isnumeric(sz) && isscalar(sz) && round(sz)==sz && sz>0,'Usage: When called with a single input, DCSurface must receive a positive integer');
                    [Surf(1:sz,1).type]=deal('');
                    return
                case 8
                    [domain,RespModel,critRng,trnRng,respTrans,varTrans,MOPairIdx,opt] = deal(varargin{:});
                    disableErrChk = false;
                case 9
                    [domain,RespModel,critRng,trnRng,respTrans,varTrans,MOPairIdx,opt,disableErrChk] = deal(varargin{:});
                otherwise
                    error('Usage: Incorrect number of inputs for DCSurface')
            end
            
            % Input error checking and initialization of empty inputs
            if ~disableErrChk
                [message respTrans varTrans opt] = ...
                    DClab.DCSurface.initialize(RespModel,domain,respTrans,varTrans,MOPairIdx,opt);
                if ~isempty(message)
                    error(message)
                end
            end
            
            %Create object here. Then use private methods to compute active parameters,
            %polynomials, fitting errors, and an estimate of the model range.
            Surf.type = opt.surfaceType;
            Surf.responseTransformation = respTrans;
            Surf.MOPairIdx = MOPairIdx;
            Surf.isPredictor = RespModel.isPredictor;
            
            [Surf.api idxStruct filesUsed] = DClab.DCSurface.findActive(Surf,RespModel,trnRng,domain,varTrans,opt);
            [Surf.surrogateModel info] = DClab.DCSurface.buildSurface(Surf,RespModel,critRng,trnRng,domain,varTrans,idxStruct,opt,filesUsed);
            
            % Estimate the range of the model function
            minResp = inf;
            maxResp = -inf;
            minResp = min([minResp; info.trueModelRangeEstimate(1)]);
            maxResp = max([maxResp; info.trueModelRangeEstimate(2)]);
            Surf.trueModelRangeEstimate = [minResp maxResp];
            
            info = rmfield(info,'trueModelRangeEstimate');
            Surf.surrogateFitInfo = info;
            Surf.status = 1;
            
        end %constructor
        
        %dependent properties
        %TODO: object array used to return contents in cell, need to fix
        %elsewhere
        function out = get.activeParameterIndex(obj)
            if isempty(obj.api)
                out = [];
            else
                out = obj.api(:,1);
            end
        end
        function out = get.activeParameterTransformation(obj)
            if isempty(obj.api)
                out = [];
            else
                out = obj.api(:,2);
            end
        end
        function out = get.peakError(obj)
            if isempty(obj.surrogateFitInfo)
                out = [];
            else
                if obj.isPredictor
                    SF = 0.5;
                else
                    SF = 0.1;
                end
                if isempty(obj.surrogateFitInfo.peakErrorOnCriticalSamplePoints4Fit)
                    peak1 = [inf -inf];
                else
                    peak1 = obj.surrogateFitInfo.peakErrorOnCriticalSamplePoints4Fit;
                end
                if isempty(obj.surrogateFitInfo.peakErrorOnSamplePoints4Validation)
                    peak2 = [inf -inf];
                else
                    peak2 = obj.surrogateFitInfo.peakErrorOnSamplePoints4Validation;
                end
                if isempty(obj.surrogateFitInfo.peakErrorFromOptimization)
                    peak3 = [inf -inf];
                else
                    peak3 = obj.surrogateFitInfo.peakErrorFromOptimization;
                end
                err(1) = min([peak1(1) peak2(1) peak3(1)]);
                err(2) = max([peak1(2) peak2(2) peak3(2)]);
                
                if isinf(err(1)) || isinf(err(2))
                    disp('Inf error found')
                    keyboard
                end
                out = (1+SF)*err;
            end
        end
        
        function obj = set.activeParameterIndex(obj,value)
            obj.api = value;
        end
        
        function out = nSurfaces(obj)
            out = sum([obj.status]==1);
        end
        
        %other methods
        function bool = isempty(obj)
            emptyIdx = vertcat(obj.status) == -1;
            if all(emptyIdx)
                bool = true;
            else
                bool = false;
            end
        end
        function [list,sz] = displayProps(obj)
            list = [];
            sz = sprintf('%d-surface',nSurfaces(obj));
        end
        
    end %methods
    
    methods (Static)
        function Q2 = composeQuadWithAffine(Q1,A,b)
            % Q2 = composeQuadWithAffine(Q1,A,b)
            %
            % Let Q(y) = [1 y']*Q1*[1; y] and let y = Aff(x) = b + Ax;
            %
            % This function computes a new quadratic form Q2 such that
            %
            % [1 x']*Q2*[1; x] = (Q \circ Aff)(x)
            
            % Description.
            %
            % Note [1; y] = [1 0 ... 0]*[1;x]
            %               [b    A   ]
            %
            % Call T = [1 0 ... 0]
            %          [b    A   ]
            %
            % Then Q2 = T'*Q1*T.
            
            %TODO add input checking to this function and give it decent help
            
            % Make sure we're dealing with a symmetric quadratic form.
            if any(any(Q1 - 0.5*(Q1+Q1') > 10*eps))
                Q1 = 0.5*(Q1+Q1');
            end
            T = [1 zeros(1,size(A,2)); b A];
            Q2 = T'*Q1*T;
            
            % Make sure the final output is symmetric.
            Q2 = 0.5*(Q2+Q2');
        end
    end
    
    methods (Static,Access='private')
        
        function [message respTransOUT varTransOUT optOUT] = ...
                initialize(RespModel,domain,respTrans,varTrans,MOPairIdx,opt)
            %TODO: doesn't look at all DCSurface inputs, needs updating
            
            % Input error checking and initialization
            
            message = '';
            respTransOUT = '';
            varTransOUT = {};
            optOUT = [];
            
            if ~isa(RespModel,'DClab.ResponseModel')
                message = 'Inputs: RESPMODEL must be a ResponseModel object';
                return
            end
            
            % the DCSurface.findActive method will verify that domain is a subset of RespModel's
            % domain.
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
            
            if isempty(respTrans)
                respTransOUT = 'none';
            elseif ~ischar(respTrans) || ~ismember(respTrans,{'none','log10'})
                message = 'Inputs: RESPTRANS must be ''none'' or ''log10''';
                return
            else
                respTransOUT = respTrans;
            end
            
            if isempty(varTrans);
                varTransOUT = repmat({'none'},size(domain,1),1);
            elseif ~iscellstr(varTrans) || ~isequal(size(varTrans),size(domain))
                message = 'Inputs: varTrans must be a column cell array of chars that has the same dimensions as DOMAIN.';
                return
            elseif length([strmatch('none',char(varTrans)); strmatch('log10',char(varTrans))]) ~= length(domain)
                message = 'Inputs: variable transformations must be ''none'' or ''log10''.';
                return
            else
                varTransOUT = varTrans;
            end
            
            if ~isnumeric(MOPairIdx)
                message = 'Inputs: MOPAIRIDX must be a positive scalar';
                return
            elseif isempty(MOPairIdx)
                % Do nothing
            elseif ~isnumeric(MOPairIdx) || MOPairIdx < 1 || round(MOPairIdx) ~= MOPairIdx
                message = 'Inputs: MOPAIRIDX must be a positive scalar';
                return
            else
                % Do nothing
            end
            
            if isempty(opt)
                optOUT = DClab.DCOptions;
            elseif ~isa(opt,'DClab.DCOptions')
                message = 'Inputs: OPT must be a DCOptions object';
            else
                optOUT = opt;
            end
            
        end
        
        function [poly, info] = buildSurface(Surf,RespModel,critRng,trnRng,domain,varTrans,idxStruct,opt,filesToSkip)
            %BUILDSURFACE: private DCSurface function that creates a surrogate model.
            %
            %   [POLY FITINFO] = BUILDSURFACE(SURF,RESPMODEL,DOMAIN,VARTRANS,IDXSTRUCT,OPT)
            %   creates the polynomial surrogate model POLY based on RESPMODEL. FITINFO
            %   is a structure that contains error information and other stats. These
            %   are placed in the fields .surrogateModel and .surrogateFitInfo of the
            %   final object. The variable transformations used in the fit are passed
            %   in as VARTRANS. The output or response transformation is present as a
            %   field of SURF.
            %
            %   Inputs:
            %      SURF is a partially complete DCSurface object, with
            %         activeParameterIndex specified.
            %      RESPMODEL is the ResponseModel the surface is for.
            %      DOMAIN: an nx1 structure with fields .name and .range.
            %      VARTRANS: an nx1 cell array containing 'none' or 'log10'.
            %      IDXSTRUCT is a structure whose fields contain useful indices. See
            %         DCSurface.getIdxVectors
            %      OPT is a DCOptions object or structure with the required fields
            %
            %   Outputs:
            %      POLY
            %      INFO
            %
            %   See also DCSurface, DCSurface.findActive
            
            %old call with metamodels
            % function [poly, info] =
            % buildSurface(Surf,RespModel,Metamodel,critRng,trnRng,domain,varTrans,idxS
            % truct,opt,filesToSkip)
            
            
            %TODO we need to take advantage of the fileCreated stuff so we don't repeat
            %these in validataion.
            
            % To improve accuracy, increase the number of fitting iterations if we're
            % performing a prediction on this response
            if RespModel.isPredictor
                opt.minFitIter = opt.minFitIter+1;
            end
            
            % Define some useful indexing vectors
            dsetNames2ModelNames = idxStruct.dsetNames2ModelNames;
            dsetNames2SharedNames = idxStruct.dsetNames2SharedNames;
            modelNames2SharedNames = idxStruct.modelNames2SharedNames;
            
            guihand = opt.guiHandle; %tell us to display to screen or gui.
            
            %Don't do any work if quadratic or linear.
            isLin = strcmp(RespModel.type,'linearModel');
            isQuad = strcmp(RespModel.type,'quadraticModel');
            linOk = ismember(opt.surfaceType,{'LinSurf';'QuadSurf';'QuadOverQuad'});
            quadOk = ismember(opt.surfaceType,{'QuadSurf';'QuadOverQuad'});
            
            %TODO, in principle, we only need the variable transformations to
            %correspond for the active parameters.
            varTransTheSame = isequal(RespModel.variableTransformations,varTrans(dsetNames2ModelNames));
            respTransTheSame = strcmp(RespModel.responseTransformation,Surf.responseTransformation);
            
            if ( (isLin && linOk) || (isQuad && quadOk) ) && varTransTheSame && respTransTheSame
                algform = RespModel.model;
                if strcmp(RespModel.type,'quadraticModel')
                    poly = algform([1; modelNames2SharedNames+1],[1; modelNames2SharedNames+1]);
                else
                    %TODO we represent linear models in quadratic form. Someday we should
                    %fix this.
                    algform = algform([1; modelNames2SharedNames+1]);
                    poly = zeros(length(algform));
                    poly(1:end,1) = 0.5*algform;
                    poly(1,1:end) = poly(1,1:end) + 0.5*algform';
                end
                
                % Eliminate any nonactive parameters from poly. We need to "sum" the
                % effect of these terms into the constant term and blow them away.
                sharedRanges = vertcat(domain(dsetNames2SharedNames).range);
                sh2nonact = find(diff(sharedRanges,[],2) == 0);
                sh2actLOGICAL = diff(sharedRanges,[],2) ~= 0;
                
                % Check to ensure nothing funky is going on
                if ~isequal(dsetNames2SharedNames(sh2actLOGICAL),Surf.activeParameterIndex)
                    error('An inconsistency was discovered in the active parameters')
                end
                
                % If any of the parameters are constant we need to plug in the constant
                % value of the parameters, add the result to the constant term, then
                % squeeze the corresponding dimensions out of the matrix.
                
                transformedSharedRanges = zeros(size(sharedRanges));
                for i1 = 1:size(sharedRanges,1)
                    if strcmp(varTrans{dsetNames2SharedNames(i1)},'log10')
                        transformedSharedRanges(i1,:) = log10(sharedRanges(i1,:));
                    else
                        transformedSharedRanges(i1,:) = sharedRanges(i1,:);
                    end
                end
                
                n = size(poly,1);
                constantx = transformedSharedRanges(sh2nonact,1);
                poly(:,1) = poly(:,1) + sum(repmat(constantx',n,1).*poly(:,sh2nonact+1),2);
                poly(:,sh2nonact+1) = [];
                n = size(poly,2);
                poly(1,:) = poly(1,:) + sum(repmat(constantx,1,n).*poly(sh2nonact+1,:),1);
                poly(sh2nonact+1,:) = [];
                
                % Here we scale the surface to be defined on +/- one.
                % y = [1 rho']*poly*[1;rho];
                % rho = 0.5*A*(x+1) + b
                
                A = diag(diff(transformedSharedRanges(sh2actLOGICAL,:),[],2));
                b = transformedSharedRanges(sh2actLOGICAL,1);
                
                AA = 0.5*A;
                bb = b+0.5*A*ones(size(A,1),1);
                poly = DClab.DCSurface.composeQuadWithAffine(poly,AA,bb);
                
                info.maxIterReached = false;
                info.maxPntsReached = false;
                info.nSamplePoints4Fit = 0;
                info.nSamplePoints4Validation = 0;
                info.errorType = 'absolute';
                info.averageErrorOnSamplePoints4Fit = 0;
                info.averageErrorOnCriticalSamplePoints4Fit = 0;
                info.averageErrorOnSamplePoints4Validation = 0;
                info.peakErrorOnSamplePoints4Fit = [0 0];
                info.peakErrorOnCriticalSamplePoints4Fit = [0 0];
                info.peakErrorOnSamplePoints4Validation = [0 0];
                info.peakErrorFromOptimization = [0 0];
                info.xInModelOrderGivingErrMin = mean(vertcat(domain(dsetNames2ModelNames).range),2);
                info.xInModelOrderGivingErrMax = mean(vertcat(domain(dsetNames2ModelNames).range),2);
                %just look at the linear part to estimate the range;
                %TODO: need to fix this as it gives us neither an inner estimate (which
                %is consistent with what the iterative fitting/validate produces) or an
                %outer estimate.
                if strcmp(Surf.responseTransformation,'log10')
                    if length(poly) < 2
                        info.trueModelRangeEstimate = 10.^[poly(1) poly(1)];
                    else
                        info.trueModelRangeEstimate = 10.^(poly(1)+[-1 1]*sum(abs(poly(2:end,1))));
                    end
                else
                    if length(poly) < 2
                        info.trueModelRangeEstimate = [poly(1) poly(1)];
                    else
                        info.trueModelRangeEstimate = (poly(1)+[-1 1]*sum(abs(poly(2:end,1))));
                    end
                end
                
            else
                %An algebraic surrogate is not already available. Save several variables
                %in a structure and call the surface fitting algorithm.
                t.opt = opt;
                t.guihand = guihand;
                t.RespModel = RespModel;
                t.Metamodel = DClab.DCMetamodel;
                t.critRng = critRng;
                t.trnRng = trnRng;
                t.domain = domain;
                t.varTrans = varTrans;
                t.dsetNames2ModelNames = dsetNames2ModelNames;
                t.dsetNames2SharedNames = dsetNames2SharedNames;
                t.modelNames2SharedNames = modelNames2SharedNames;
                t.Surf = Surf;
                t.filesToSkip = filesToSkip;
                
                % This is a little weird.
                % Fit error comes out in the transformed y-coordinated, but response is
                % in the original coordinates.
                if strcmp(opt.surfaceFittingMode,'one-shot')
                    [poly fitErr paramValues response N filesToSkip yMin yMax subspaceInfo] = DClab.DCSurface.oneShotFit(t);
                    maxIterReached = false;
                    maxPntsReached = false;
                else
                    [poly fitErr maxIterReached maxPntsReached paramValues response N filesToSkip] = DClab.DCSurface.iterativeFittingAlgorithm(t);
                    yMin = min(response);
                    yMax = max(response);
                    subspaceInfo = [];
                end
                t.filesToSkip = filesToSkip;
                
                critIdx = find(~isnan(response) & response >= critRng(1) & response <= critRng(2));
                critFitErr = fitErr(critIdx);
                
                % Begin to fill out surrogate fit info structure
                info.maxIterReached = maxIterReached;
                info.maxPntsReached = maxPntsReached;
                info.nSamplePoints4Fit = N;
                info.nSamplePoints4Validation = []; %filled in later by validate
                info.errorType = 'absolute';
                info.averageErrorOnSamplePoints4Fit = mean(abs(fitErr));
                if isempty(critFitErr)
                    info.averageErrorOnCriticalSamplePoints4Fit =NaN;
                else
                    info.averageErrorOnCriticalSamplePoints4Fit = mean(abs(critFitErr));
                end
                info.averageErrorOnSamplePoints4Validation = []; %filled in later
                
                info.peakErrorOnSamplePoints4Fit = [min(fitErr) max(fitErr)];
                info.subspaceInfo = subspaceInfo;
                
                % Location and max error
                [maxE idX] = max(critFitErr);
                [minE idx] = min(critFitErr);
                info.peakErrorOnCriticalSamplePoints4Fit = [minE maxE];
                
                info.peakErrorOnSamplePoints4Validation = []; %filled in later
                info.peakErrorFromOptimization = []; %filled in later
                
                % These are in the model parameter order. validate will update these
                % values.
                info.xInModelOrderGivingErrMin = paramValues(critIdx(idx),:)';
                info.xInModelOrderGivingErrMax = paramValues(critIdx(idX),:)';
                info.trueModelRangeEstimate = [yMin yMax];
                
                % Add some more info to the structure and validate the fit. We fill in
                % .surrogateModel here for convenience. The DCSurface constructure makes
                % this addition in the final object.
                t.Surf.surrogateModel = poly;
                t.info = info;
                info = DClab.DCSurface.validate(t);
                
            end
        end
        
        function [poly,fitError] = constructLinFit(x,xDomainBnds,response,opt)
            %function [poly,peakError,avgError,allError] = ...
            %    constructFit(x,activeSurfBnds,response,fitType);
            %
            %  function that recieves a matrix of paramValues, the domain
            %  overwhich these live (for normalization purposes), a vector of
            %  responses, and computes a quadratic fit of type specifed by
            %  fitType in the normalized x's to the responses.
            %
            %  Inputs:
            %    x: a Nsamples x Ndims matrix
            %    xDomainBnds: a Ndims x 2 matrix defining the domain of the fit.
            %    response: a column of models outputs computed at each row of x
            %    fitType: can be 2 or inf
            
            ni = nargin;
            if ni==3
                opt = DClab.DCOptions;
            end
            fitType = opt.fitNorm;
            
            n = size(xDomainBnds,1);
            N = length(response);
            
            % determine and delete any degenerate dimensions for which the
            % interval of x's is just a point.
            idx = find(diff(xDomainBnds,1,2) == 0);
            goodidx = setdiff([1:n]',idx);
            x = x(:,goodidx);
            xDomainBnds = xDomainBnds(goodidx,:);
            normx = 2*(x - repmat(xDomainBnds(:,1)',N,1))*diag(1./diff(xDomainBnds,1,2)) - 1;
            %X = DClab.DCSurface.buildLinX(normx);
            
            %now do fit
            if isinf(fitType)
                [coeff, fitError, EXITFLAG] = DClab.DCSurface.infNormRegression(normx,response);
            elseif fitType == 2
                [coeff, fitError, EXITFLAG] = DClab.twoNormRegression(normx,response);
                %elseif strcmp(fitType,'infb')
                %  [coeff, peakError, avgError, allError, EXITFLAG] = ...
                %      infNormRegression2(normx,response);
                %elseif strcmp(fitType,'infc')
                %  [coeff, peakError, avgError, allError, EXITFLAG] = ...
                %      infNormRegression3(fitX,fitY);
                
            else
                error('Usage: Improper fit type specified')
            end
            
            % fill any degenerate dimension back into the poly.
            poly = zeros(n+1);
            poly(1,1) = coeff(1);
            
            poly(1,goodidx+1) = 0.5*coeff(2:end);
            poly(goodidx+1,1) = 0.5*coeff(2:end);
        end
        
        function [poly,fitError] = constructQuadFit(x,xDomainBnds,response,opt)
            %function [poly,peakError,avgError,allError] = ...
            %    constructQuadFit(x,activeSurfBnds,response,fitType);
            %
            %  function that recieves a matrix of paramValues, the domain
            %  overwhich these live (for normalization purposes), a vector of
            %  responses, and computes a quadratic fit of type specifed by
            %  fitType in the normalized x's to the responses.
            %
            %  Inputs:
            %    x: a Nsamples x Ndims matrix
            %    xDomainBnds: a Ndims x 2 matrix defining the domain of the fit.
            %    response: a column of models outputs computed at each row of x
            %    fitType: can be 2 or inf
            
            ni = nargin;
            if ni==3
                opt = DClab.DCOptions;
            end
            fitType = opt.fitNorm;
            
            n = size(xDomainBnds,1);
            N = length(response);
            
            % determine and delete any degenerate dimensions for which the
            % interval of x's is just a point.
            idx = find(diff(xDomainBnds,1,2) == 0);
            goodidx = setdiff([1:n]',idx);
            x = x(:,goodidx);
            xDomainBnds = xDomainBnds(goodidx,:);
            
            if isempty(xDomainBnds)
                normx = zeros(N,0);
            else
                normx = 2*(x - repmat(xDomainBnds(:,1)',N,1))*diag(1./diff(xDomainBnds,1,2)) - 1;
            end
            X = DClab.buildQuadX(normx);
            
            %now do fit
            if isinf(fitType)
                [coeff, fitError, EXITFLAG] = ...
                    DClab.DCSurface.infNormRegression(X,response);
            elseif fitType == 2
                [coeff, fitError, EXITFLAG] = ...
                    DClab.twoNormRegression(X,response);
                
                % elseif strcmp(fitType,'infb')
                %     [coeff, peakError, avgError, allError, EXITFLAG] = ...
                %         infNormRegression2(X,response);
                % elseif strcmp(fitType,'infc')
                %     [coeff, peakError, avgError, allError, EXITFLAG] = ...
                %         infNormRegression3(X,response);
                
            else
                error('Usage: Improper fit type specified')
            end
            
            temppoly = DClab.coeff2quadform(coeff,size(xDomainBnds,1));
            
            % fill any degenerate dimension back into the poly.
            poly = zeros(n+1);
            poly([1; goodidx+1],[1; goodidx+1]) = temppoly;
            
            
            
        end
        
        function [activeIdx idxStruct allFiles] = findActive(Surf,RespModel,trnRng,domain,varTrans,opt)
            %FINDACTIVE is a private method to determine active parameters of a DCSurface.
            %
            %   ACTIVEIDX = FINDACTIVE(SURF,RESPMODEL,METAMODEL,DOMAIN,VARTRANS,OPT)
            %   returns a nactivex2 array. The first column is a vector of indices into
            %   Surf.domain indicating the active parameters. The second column will be
            %   contain 1's and 2's. 1 indicates that no transformation was applied
            %   to the corresponding active parameter during the computation. 2
            %   indicates a log10 transformation was used. The input varTrans specifies
            %   the transformations.
            %
            %   [ACTIVEIDX IDXSTRUCT] = FINDACTIVE(...) additionally returns a
            %   structure containing a bunch of indexing vectors. These can be passed
            %   to buildSurface and validateSurface so we don't need to recalculate
            %   them in these function.
            %
            %   Inputs:
            %      SURF: a partially complete DCSurface object
            %      RESPMODEL: the ResponseModel the surface will be build for
            %      DOMAIN: an nx1 structure with fields .name and .range.
            %      VARTRANS: an nx1 cell array containing 'none' or 'log10'.
            %      OPT: a DCOptions object or structure having the required fields.
            %
            %   See also DCSurface, DCSurface/buildSurface
            
            %old meta model call
            %function [activeIdx idxStruct allFiles] =
            %findActive(Surf,RespModel,Metamodel,trnRng,domain,varTrans,opt)
            
            % Define some useful indexing vectors. See comments in DCSurface.getIdxVectors for
            % details. In short
            %  dsetParamNames(dsetNames2ModelNames) == modelParamNames)
            %  sharedParamNames == dsetParamNames(dsetNames2SharedNames)
            %  sharedParamNames == modelParamNames(modelNames2SharedNames)
            
            allFiles = [];
            
            dsetParamNames = {domain.name}';
            modelParamNames = RespModel.parameterList;
            [dsetNames2ModelNames dsetNames2SharedNames modelNames2SharedNames] = DClab.DCSurface.getIdxVectors(dsetParamNames,modelParamNames);
            
            % Cache these as the second output
            idxStruct.dsetNames2ModelNames = dsetNames2ModelNames;
            idxStruct.dsetNames2SharedNames = dsetNames2SharedNames;
            idxStruct.modelNames2SharedNames = modelNames2SharedNames;
            
            % Make sure we're working on a subset of the domain of RespModel
            RMdomain = RespModel.domain;
            RMdomrng = vertcat(RMdomain.range);
            RMdomrng(isinf(RMdomrng(:,1)),1) = -1e100;
            RMdomrng(isinf(RMdomrng(:,2)),2) = 1e100;
            
            if ~DClab.issubset(vertcat(domain(dsetNames2ModelNames).range),RMdomrng,1e-10,1e-10)
                error('DOMAIN is not a subset of the domain of RESPMODEL')
            end
            
            % The variable sharedRanges contains the ranges of the shared parameters.
            %
            % It is possible that some of these shared variables will have a singleton
            % range. candidateActiveIdx is an index into dsetParamNames that returns the
            % names of the parameters that are considered as we search for active
            % parameters. These will be a subset of the shared parameters.
            
            isLin = strcmp(RespModel.type,'linearModel');
            isQuad = strcmp(RespModel.type,'quadraticModel');
            
            sharedRanges = vertcat(domain(dsetNames2SharedNames).range);
            doesVaryIdx = diff(sharedRanges,[],2) > 0;
            
            % If ResponseModel is of type dcModel and we're doing metamodel-based
            % analysis, eliminate any parameters that are not active in the
            % metamodel from the candidateActiveIdx.
            % if ~isLin && ~isQuad && ~strcmp(opt.analysisMode,'original')
            %   modelActive = Metamodel.variableTransformations~=0;
            %   sharedActive = modelActive(modelNames2SharedNames);
            %   candidateActiveIdx = dsetNames2SharedNames(doesVaryIdx & sharedActive);
            %   modelNames2candidateActiveNames = modelNames2SharedNames(doesVaryIdx & sharedActive);
            % else
            candidateActiveIdx = dsetNames2SharedNames(doesVaryIdx);
            modelNames2candidateActiveNames = modelNames2SharedNames(doesVaryIdx);
            % end
            
            guihand = opt.guiHandle; %tell us to display to screen or gui.
            
            linOk = ismember(opt.surfaceType,{'LinSurf';'QuadSurf';'QuadOverQuad'});
            quadOk = ismember(opt.surfaceType,{'QuadSurf';'QuadOverQuad'});
            varTransTheSame = isequal(RespModel.variableTransformations,varTrans(dsetNames2ModelNames));
            respTransTheSame = strcmp(RespModel.responseTransformation,Surf.responseTransformation);
            if ( (isLin && linOk) || (isQuad && quadOk) ) && varTransTheSame && respTransTheSame
                % Don't do any work if quadratic or linear
                if ismember(opt.display,{'all';'ALL'})
                    str = [blanks(2) 'Algebraic Model: all parameters of ' RespModel.name ' considered active'];
                    DClab.dcdispstr(str,guihand,false)
                end
                activeIdx = [candidateActiveIdx ones(size(candidateActiveIdx))];
                log10Trans = strmatch('log10',char(varTrans(candidateActiveIdx)),'exact');
                activeIdx(log10Trans,2) = 2; %Second column of activeIdx indicates the transformation
                
            elseif opt.activeParamSelCutOff == 0
                % Don't do any work if all will be considered active
                if ismember(opt.display,{'all';'ALL'})
                    str = [blanks(2) 'opt.activeParamSelCutOff == 0: all parameters of ' RespModel.name ' considered active'];
                    DClab.dcdispstr(str,guihand,false)
                end
                
                activeIdx = [candidateActiveIdx ones(size(candidateActiveIdx))];
                log10Trans = strmatch('log10',char(varTrans(candidateActiveIdx)),'exact');
                activeIdx(log10Trans,2) = 2;
                
            else
                % Compute the active parameters
                if ismember(opt.display,{'all';'ALL'})
                    str = [blanks(2) 'Determining active parameters for ' RespModel.name];
                    DClab.dcdispstr(str,guihand,false)
                end
                
                cutOff = opt.activeParamSelCutOff;
                modelOrderRanges = vertcat(domain(dsetNames2ModelNames).range);
                
                candidateRanges = vertcat(domain(candidateActiveIdx).range);
                n = size(candidateRanges,1);
                Nsamples = opt.nPntsPerParam4ActiveParamSel*n;
                
                % If we're not doing metamodel based analysis or the ResponseModels are
                % algebraic, obtain evaluations from the true models.
                %   if isLin || isQuad
                
                % Load any existing model evaluations, regardless of whether the points
                % were distributed in lin or log space
                [paramValues response filesQueried] = loadSavedEvaluations(RespModel,modelOrderRanges,[],Nsamples);
                paramValues = paramValues';
                response = response';
                
                % Perform further model evaluations if necessary
                if length(response) < Nsamples
                    Ns = Nsamples - length(response);
                    try
                        [newParamValues newResponse fileCreated] = randomEvalWithSave(RespModel,modelOrderRanges,varTrans,Ns,opt.nComputer);
                        newParamValues = newParamValues';
                        newResponse = newResponse';
                    catch
                        str = 'randomEvalWithSave call failed in DCSurface/private/findActive';
                        DClab.dcdispstr(str,guihand,true)
                    end
                    
                    paramValues = [paramValues; newParamValues];
                    response = [response; newResponse];
                else
                    fileCreated = [];
                end
                allFiles = [filesQueried; fileCreated];
                %   else
                %     [paramValues, response] = randomEval(Metamodel,modelOrderRanges,varTrans,Nsamples);
                %     paramValues = paramValues';
                %     response = response';
                %   end
                
                %   % Rejection step:
                %   if any(isnan(response)) && any(isinf(trnRng))
                %     error('Response Model produced NaN values, but training range was infinite')
                %   end
                %
                %   junkIdx = innan(response) | response < trnRng(1) | response > trnRng(2);
                %   response = response(~junkIdx);
                %   paramValues = paramValues(~junkIdx,:);
                %
                %   if length(response) <= 0.1*opt.nPntsPerCoeff4ActiveParamSel*n
                %     disp('In DCSurface/findActive: Fewer than 10%% of the points lie in the training range.')
                %   end
                
                % Shuffle the parameter values and eliminate unneeded ones so that they
                % are in the order of the candidate active parameters.
                paramValues = paramValues(:,modelNames2candidateActiveNames);
                N = length(response);
                
                if isempty(candidateRanges)
                    normx = [];
                else
                    normx = zeros(size(paramValues));
                    log10Trans = strmatch('log10',char(varTrans(candidateActiveIdx)),'exact');
                    noTrans = strmatch('none',char(varTrans(candidateActiveIdx)),'exact');
                    scale = diag(1./diff(candidateRanges(noTrans,:),1,2));
                    normx(:,noTrans) = 2*(paramValues(:,noTrans) - repmat(candidateRanges(noTrans,1)',N,1))*scale - 1;
                    scale = diag(1./log10(candidateRanges(log10Trans,2)./candidateRanges(log10Trans,1)));
                    normx(:,log10Trans) = 2*log10(paramValues(:,log10Trans)./repmat(candidateRanges(log10Trans,1)',N,1))*scale - 1;
                end
                XX = [ones(Nsamples,1) normx];
                
                %TODO we need to describe the math behind this somewhere...
                
                %now do qr on XX and fit coeffs to the orthogonal part, since we're just
                %ranking, no need to actually do the inv(R)*coeff bit, but let multiply by
                %diag(1./diag(R)) to get scaling somewhat right.
                
                [Q,R] = qr(XX,0); % 2nd input specified "economy size" R
                
                if strcmp(Surf.responseTransformation,'none')
                    coeff = Q\response;
                elseif strcmp(Surf.responseTransformation,'log10')
                    coeff = Q\log10(response);
                else
                    error('Encountered invalid response transformation')
                end
                
                coeff = diag(1./diag(R))*coeff; %we really shouldn't do the drop
                %parameters without refitting thing unless we verify that R is pretty
                %diagonal. otherwise we can make mistakes... if we had orthogonal design
                %points we'd be ok...
                
                %note, the constant term is delibrately not included here.
                den = sum( (1/3)*(2^n)*coeff(2:end).^2 ) + eps;
                
                currentError = 0;
                candidatesToKeep = [1:n]'; %#ok
                while currentError <= cutOff
                    tempCurrentError = repmat(inf,n,1);
                    oldCandidatesToKeep = candidatesToKeep;
                    for i2 = candidatesToKeep'
                        tempIdx = setdiff(candidatesToKeep,i2);
                        %    tempNormx = normx(:,tempIdx);
                        %    tempXX = [ones(length(normx),1) tempNormx];
                        %    tempCoeff = coeff([1;1+tempIdx]);
                        %    trytempCoeff = Q(:,[1; tempIdx+1])\response;
                        
                        %take the 2 norms of the difference between the fit in all and the fit
                        %in fewer. we can do this by just zeroing out components
                        c1 = coeff;
                        c1([1; tempIdx+1]) = 0;
                        tempnum = sum( (1/3)*(2^n)*c1.^2 );
                        tempCurrentError(i2,1) = tempnum/den;
                    end
                    
                    [currentError, drop] = min(tempCurrentError);
                    candidatesToKeep = setdiff(candidatesToKeep,drop);
                end
                
                activeIdx = candidateActiveIdx(oldCandidatesToKeep);
                log10Trans = strmatch('log10',char(varTrans(activeIdx)),'exact');
                activeIdx = [activeIdx ones(size(activeIdx))];
                activeIdx(log10Trans,2) = 2; %Second column of activeIdx indicates the transformation
                
                if ismember(opt.display,{'all';'ALL'})
                    str = [blanks(4) num2str(size(activeIdx,1)) ' of ' num2str(length(modelParamNames)) ' are active'];
                    DClab.dcdispstr(str,guihand,false)
                end
                
                if strcmp(opt.display,'ALL')
                    str = [blanks(4) 'The active parameters are:'];
                    DClab.dcdispstr(str,opt.guiHandle,false)
                    str = {domain(activeIdx(:,1)).name};
                    DClab.dcdispstr(str,guihand,false)
                end
                
            end
            
            
            
        end
        
        function [ratio, a] = isConverged(poly1,poly2,subCubeBnds)
            
            n = size(subCubeBnds,1);
            
            %===code that uses the symbolic toolbox to integrate===
            % x = 1;
            % for i1 = 1:n
            %   syms(['x' num2str(i1)],'real');
            %   eval(['x = [x; x' num2str(i1) ' ];']);
            % end
            %
            % integrand1 = (x'*poly1*x - x'*poly2*x)^2;
            %
            % syms a real
            %
            % integrand2 = (x'*poly2*x - a)^2;
            %
            % try
            %
            %   for i1 = 1:n
            %     intOut = ['x' num2str(i1)];
            %
            %     integrand1 = int(integrand1,intOut,-1,1);
            %     integrand2 = int(integrand2,intOut,-1,1);
            %   end
            % catch
            %   disp('error in integration steps of QuadSurf')
            %   keyboard
            % end
            %
            % num = double(integrand1);
            % a = solve(diff(integrand2,'a'),'a');
            % den = double(eval(integrand2));
            %
            % ratio = num/den;
            % a = double(a);
            
            %===code that does the integration by hand===
            if isnumeric(poly1)
                
                
                S1 = poly1-poly2;
                S2 = poly2; %really want to integrate poly2 - a where a is an
                %unknown scalar
                
                quadprod1 = 0; %initialize
                quadprod2 = 0;
                if n > 0
                    
                    affine1 = S1(1,1);
                    linear1 = 2*S1(2:end,1);
                    quad1 = diag(S1);
                    quad1 = quad1(2:end);
                    interact1 = S1 - diag(diag(S1));%these are still divided by two!
                    interact1 = 2*triu(interact1);
                    interact1 = interact1([2:end],[2:end]);
                    
                    affine2 = S2(1,1);
                    linear2 = 2*S2(2:end,1);
                    quad2 = diag(S2);
                    quad2 = quad2(2:end);
                    interact2 = S2 - diag(diag(S2));%these are still divided by two!
                    interact2 = 2*triu(interact2);
                    interact2 = interact2([2:end],[2:end]);
                    
                    ii = 1;
                    for i2 = 1:length(quad1);
                        for i3 = i2+1:length(quad1)
                            quadprod1(ii,1) = quad1(i2)*quad1(i3);
                            quadprod2(ii,1) = quad2(i2)*quad2(i3);
                            ii = ii+1;
                        end
                    end
                    
                else
                    affine1 = 0;
                    quad1 = 0;
                    linear1 = 0;
                    interact1 = 0;
                    quadprod1 = 0;
                    
                    affine2 = 0;
                    quad2 = 0;
                    linear2 = 0;
                    interact2 = 0;
                    quadprod2 = 0;
                    
                end
                
                integrand1 = 2^n*(affine1^2 + (2/3)*affine1*sum(quad1) + inv(3)* ...
                    sum(linear1.^2) + inv(9)*sum(sum(interact1.^2)) + (2/9)*sum(quadprod1) + inv(5)* ...
                    sum(quad1.^2));
                
                %we don't need to really integrate S2, since many terms drop out
                %when we differentiate. we need to differentiate (affine2-a)^2 +
                %(2/3)*(affine2-a)*sum(quad2) + CONST and solve for a
                A1 = 2;
                A0 = -2*affine2 - (2/3)*sum(quad2);
                
                a = -A0/A1;
                
                affine2 = affine2-a;
                integrand2 = 2^n*(affine2^2 + (2/3)*affine2*sum(quad2) + inv(3)* ...
                    sum(linear2.^2) + inv(9)*sum(sum(interact2.^2)) + (2/9)*sum(quadprod2) + inv(5)* ...
                    sum(quad2.^2));
                
                ratio = integrand1/(integrand2+eps);
                
            else
                x = 2*rand(10000,n)-1;
                xx = [ones(10000,1) x];
                num1 = poly1.num;
                den1 = poly1.den;
                num2 = poly2.num;
                den2 = poly2.den;
                
                int1 = mean( (sum(xx*num1.*xx,2)./sum(xx*den1.*xx,2) - sum(xx*num2.*xx,2)./sum(xx*den2.*xx,2)).^2 );
                a = mean( sum(xx*num2.*xx,2)./sum(xx*den2.*xx,2) );
                int2 = mean( (sum(xx*num2.*xx,2)./sum(xx*den2.*xx,2)).^2 ) - a^2;
                
                ratio = int1/int2;
            end
            
            %Ngrid = 25;
            %[X1 X2] = meshgrid(linspace(-1,1,Ngrid));
            %x1 = X1(:);
            %x2 = X2(:);
            %xx = [ones(size(x1)) x1 x2];
            %
            %if isstruct(poly1);
            %  yold = sum(xx*num1.*xx,2)./sum(xx*den1.*xx,2);
            %  ynew = sum(xx*num2.*xx,2)./sum(xx*den2.*xx,2);
            %else
            %  yold = sum(xx*poly1.*xx,2);
            %  ynew = sum(xx*poly2.*xx,2);
            %end
            %
            % h1 = figure;
            % surf(X1,X2,reshape(yold,size(X1)))
            % title('prev surface')
            % a1 = axis;
            % h2 = figure;
            % surf(X1,X2,reshape(ynew,size(X1)))
            % title('current surface')
            % a2 = axis;
            % a3 = a2;
            % a3(5) = min(a1(5),a2(5));
            % a3(6) = max(a1(6),a2(6));
            % axis(a3);
            % figure(h1)
            % axis(a3);
            % h3 = figure;
            % surf(X1,X2,reshape(ynew-yold,size(X1)))
            % title('deviation surface')
            %
            % keyboard
            %
            % close([h1 h2 h3])
        end
        
        function [poly fitErr maxIterReached maxPntsReached paramValues response N filesToSkip yMin yMax] = iterativeFittingAlgorithm(t)
            
            opt = t.opt;
            guihand = t.guihand;
            RespModel = t.RespModel;
            % Metamodel = t.Metamodel;
            critRng = t.critRng;
            domain = t.domain;
            varTrans = t.varTrans;
            dsetNames2ModelNames = t.dsetNames2ModelNames;
            dsetNames2SharedNames = t.dsetNames2SharedNames;
            modelNames2SharedNames = t.modelNames2SharedNames;
            Surf = t.Surf;
            filesToSkip = t.filesToSkip;
            
            if strcmp(opt.analysisMode,'original') || ~strcmp(RespModel.type,'dcModel')
                metaBased = false;
            else
                metaBased = true;
            end
            
            if ismember(opt.display,{'all';'ALL'})
                str = [blanks(2) 'Generating surrogate for ' RespModel.name];
                DClab.dcdispstr(str,guihand,false)
            end
            
            %decide which fit constructor to use:
            switch opt.surfaceType
                case 'LinSurf'
                    fitterHandle = @DClab.DCSurface.constructLinFit;
                case 'QuadSurf'
                    fitterHandle = @DClab.DCSurface.constructQuadFit;
                    %case 'QuadOverQuad'
                    % fitterHan = @qratfit;
                    %fitterHan = @constructQOQFit;
                otherwise
                    error('DCSurface does not support that surfaceType')
            end
            
            % Load any existing model evaluations
            modelOrderRanges = vertcat(domain(dsetNames2ModelNames).range);
            modelOrderTrans = varTrans(dsetNames2ModelNames);
            
            %keeps track of filesQueried and filesUsed so we don't duplicate
            %the same old samples when we do validation.
            %setminus(filesQueried,filesUsed) will just be files that didn't
            %contain any useful points.
            
            %This block seems to be redundant... all these variables are just
            %overwritten.
            % if metaBased
            %   % No saved points.
            %   paramValues = [];
            %   response = [];
            %   Nsamples = 0;
            % else
            %
            %   [paramValues, response] = ...
            %       loadSavedEvaluations(RespModel,modelOrderRanges,modelOrderTrans,filesToSkip);
            %
            %   paramValues = paramValues';
            %   response = response';
            %   Nsamples = length(response);
            % end
            
            
            
            % Determine number of coefficients in fitting polynomial
            n = size(Surf.activeParameterIndex,1);
            switch opt.surfaceType
                case 'QuadSurf'
                    NCoeff = (n+1)*(n+2)/2; %returns how many coeffs required to fit a
                    %quad in n variables.
                case 'LinSurf'
                    NCoeff = n+1;
                case 'QuadOverQuad'
                    NCoeff = (n+1)*(n+2);
                otherwise
                    error('Surface Type not supported by DCSurface class, check your opt.surfaceType property')
            end
            
            % Number of points to use for the initial fit
            N = 2*NCoeff;
            
            if strcmp(opt.display,'ALL')
                str = [blanks(4) 'Gathering model evalutations'];
                DClab.dcdispstr(str,guihand,false);
            end
            
            if metaBased
                [paramValues, response] = randomEval(Metamodel,modelOrderRanges,modelOrderTrans,N);
                paramValues = paramValues';
                response = response';
                Nsamples = length(response);
            else
                
                [paramValues, response, newFiles] = ...
                    loadSavedEvaluations(RespModel,modelOrderRanges,modelOrderTrans,N,filesToSkip);
                
                
                filesToSkip = [filesToSkip; newFiles];
                
                paramValues = paramValues';
                response = response';
                Nsamples = length(response);
                
                if Nsamples < N
                    Ns = N - Nsamples;
                    
                    try
                        [newParamValues newResponse newFile] = ...
                            randomEvalWithSave(RespModel,modelOrderRanges,modelOrderTrans,Ns,opt.nComputer);
                    catch
                        str = 'randomEvalWithSave call failed in DCSurface.buildSurface';
                        DClab.dcdispstr(str,guihand,true)
                    end
                    filesToSkip = [filesToSkip; newFile];
                    
                    if any(isnan(newResponse))
                        error('Encounted NaN response values.  Check model')
                    end
                    newParamValues = newParamValues';
                    newResponse = newResponse';
                    paramValues = [paramValues; newParamValues];
                    response = [response; newResponse];
                    clear newParamValues;
                    clear newResponse; %save memory
                end
                Nsamples = length(response);
            end
            
            %TODO for speedup, remove some of these next 10 lines when tested
            
            %Both dsetNames2SharedNames and Surf.activeParameterIdx are already sorted.
            [trash sharedNames2ActiveNames] = intersect(dsetNames2SharedNames,Surf.activeParameterIndex);
            if ~isequal(trash,Surf.activeParameterIndex)
                error('Inconsisteny: the shared parameters should be a superset of active')
            end
            modelNames2ActiveNames = modelNames2SharedNames(sharedNames2ActiveNames);
            
            modelNames = RespModel.parameterList;
            if ~isequal({domain(Surf.activeParameterIndex).name}',modelNames(modelNames2ActiveNames))
                error('Problem with indexing')
            end
            
            % Convert newParamValues to the order of the active parameters and get rid
            % of nonactive ones.
            transActParamValues = paramValues(:,modelNames2ActiveNames);
            transActRanges = vertcat(domain(Surf.activeParameterIndex).range);
            %Transform these
            activeIsLog = Surf.activeParameterTransformation == 2;
            transActParamValues(:,activeIsLog) = log10(transActParamValues(:,activeIsLog));
            transActRanges(activeIsLog,:) = log10(transActRanges(activeIsLog,:));
            
            if strcmp(Surf.responseTransformation,'log10')
                transResponse = log10(response);
            else
                transResponse = response;
            end
            
            if opt.maxFitIter == 1
                
                if opt.useAllPnts4Fit == 1
                    % Do one fit with all available points and get out of here
                    [poly, fitErr] = feval(fitterHandle,transActParamValues,transActRanges,transResponse,opt);
                    critIdx = ~isnan(response) & response >= critRng(1) & response <= critRng(2);
                else
                    randIdx = randPerm(length(response));
                    [poly, fitErr] = feval(fitterHandle,transActParamValues(randIdx(1:N),:),transActRanges,transResponse(randIdx(1:N)),opt);
                    critIdx = ~isnan(response(randIdx(1:N))) & response(randIdx(1:N)) >= critRng(1) & response(randIdx(1:N)) <= critRng(2);
                end
                
                if ~any(critIdx)
                    disp('no critcal points found')
                    keyboard
                end
                
                critFitErr = fitErr(critIdx);
                
                if ismember(opt.display,{'all';'ALL'})
                    str = [blanks(4) 'Peak approximation error on critical samples: [' num2str(min(critFitErr)) ', ' num2str(max(critFitErr)) ']'];
                    DClab.dcdispstr(str,guihand,false)
                    str = [blanks(4) 'Maximum number of approximation tries reached'];
                    DClab.dcdispstr(str,guihand,false)
                end
                maxIterReached = true;
                maxPntsReached = false;
            else
                % Perform the iterative fitting algorithm. We are guarenteed that
                % maxFitIter is > 1
                
                iterIdx = 1;
                notDone = true;
                
                % Construct initial fit
                [polyOld fitErrOld] = feval(fitterHandle,transActParamValues(1:N,:),transActRanges,transResponse(1:N),opt);
                
                if strcmp(opt.display,'ALL')
                    str = [blanks(4) 'Peak error with approx0: [' num2str(min(fitErrOld)) ', ' num2str(max(fitErrOld)) ']'];
                    DClab.dcdispstr(str,guihand,false)
                end
                
                while notDone
                    
                    N = N*2; %TODO we're only doubling the number of points each iteration???
                    
                    if strcmp(opt.display,'ALL')
                        str = [blanks(4) 'Gathering model evalutations'];
                        DClab.dcdispstr(str,guihand,false)
                    end
                    
                    [newParamValues, newResponse, newFiles] = ...
                        loadSavedEvaluations(RespModel,modelOrderRanges,modelOrderTrans,N-Nsamples,filesToSkip);
                    
                    
                    filesToSkip = [filesToSkip; newFiles];
                    
                    newParamValues = newParamValues';
                    newResponse = newResponse';
                    
                    Nsamples = Nsamples+length(newResponse); %update number of available samples
                    
                    if Nsamples < N
                        Ns = N - Nsamples;
                        
                        if metaBased
                            [newParamValues2, newResponse2] = randomEval(Metamodel,modelOrderRanges,modelOrderTrans,Ns);
                        else
                            try
                                [newParamValues2 newResponse2 newFile] = ...
                                    randomEvalWithSave(RespModel,modelOrderRanges,modelOrderTrans,Ns,opt.nComputer);
                                
                            catch
                                str = 'distribute_evaluations call failed in DCSurface.buildSurface';
                                DClab.dcdispstr(str,guihand,true)
                            end
                            filesToSkip = [filesToSkip; newFile];
                        end
                        newParamValues2 = newParamValues2';
                        newResponse2 = newResponse2';
                    else
                        newParamValues2=[];
                        newResponse2 = [];
                    end
                    
                    newParamValues = [newParamValues; newParamValues2];
                    newResponse = [newResponse; newResponse2];
                    
                    paramValues = [paramValues; newParamValues;]; %save these for later
                    newTransActParamValues = newParamValues(:,modelNames2ActiveNames);
                    %Transform these
                    newTransActParamValues(:,activeIsLog) = log10(newTransActParamValues(:,activeIsLog));
                    
                    transActParamValues = [transActParamValues; newTransActParamValues];
                    response = [response; newResponse;];
                    if strcmp(Surf.responseTransformation,'log10')
                        transResponse = [transResponse; log10(newResponse)];
                    else
                        transResponse = [transResponse; newResponse];
                    end
                    Nsamples = length(transResponse); %update number of available samples
                    
                    clear newParamValues newResponse newParamValues2 newResponse2 newTransActParamValues
                    
                    % If we are displaying everything, compute and display the error of the
                    % old fit on the old and new points. Otherwise we don't care about it.
                    if strcmp(opt.display,'ALL')
                        xx = [ones(N,1) 2*(transActParamValues(1:N,:)-repmat(transActRanges(:,1)',N,1))*inv(diag(diff(transActRanges,[],2)))-1];
                        fitErrOld = sum(xx*polyOld.*xx,2)-transResponse(1:N);
                        str = [blanks(4) 'Peak error with approx' num2str(iterIdx-1) ': [' num2str(min(fitErrOld)) ', ' num2str(max(fitErrOld)) '] (on new points)'];
                        DClab.dcdispstr(str,guihand,false)
                    end
                    
                    % Compute the next fit
                    [polyNew fitErrNew] = feval(fitterHandle,transActParamValues(1:N,:),...
                        transActRanges,transResponse(1:N),opt);
                    
                    % Compute the surface deviation
                    [dev(iterIdx,1), a(iterIdx,1)] = DClab.DCSurface.isConverged(polyOld,polyNew,transActRanges); %#ok
                    
                    if strcmp(opt.display,'ALL')
                        str = [blanks(4) 'Peak error with approx' num2str(iterIdx) ': [' num2str(min(fitErrNew)) ', ' num2str(max(fitErrNew)) ']'];
                        DClab.dcdispstr(str,guihand,false)
                    end
                    if any(strmatch(opt.display,{'all';'ALL'}))
                        str = [blanks(4) 'Normalized 2-norm between approx' num2str(iterIdx-1) ' and approx' num2str(iterIdx) ': ' num2str(dev(iterIdx))];
                        DClab.dcdispstr(str,guihand,false)
                    end
                    
                    if strcmp(opt.plotFitProgress,'on')
                        figure(1);
                        title(['DCSurface fitting.  Transformation: ' trans])
                        subplot(3,1,1);
                        plot(dev);
                        ylabel('Error')
                        xlabel('iteration')
                        subplot(3,1,2);
                        plot(a)
                        % TODO: These plots need some labels/titles...
                        %ylabel()
                        xlabel('iteration')
                        subplot(3,1,3);
                        hist(transResponse(1:N),floor(N/2));
                        ylabel('# of samples')
                        xlabel('surface response')
                        pause(0.01) % Here we could slow it down to the user can absorb the data
                        % TODO should the pause time be an option?
                        drawnow;
                    end
                    
                    % Determine if we have done enough to consider exiting.
                    readyToCheckTol = iterIdx >= max(opt.minFitIter,opt.nSuccessfulFitIter);
                    
                    if  readyToCheckTol && all(dev(iterIdx-opt.nSuccessfulFitIter+1:end,1) <= opt.fitConvergenceTol)
                        
                        critIdx = ~isnan(response(1:N)) & response(1:N) >= critRng(1) & response(1:N) <= critRng(2);
                        if ~any(critIdx)
                            disp('No critcal points found in training data for quadratic fit')
                        end
                        notDone = false;
                        poly = polyNew;
                        fitErr = fitErrNew;
                        critFitErr = fitErr(critIdx);
                        maxIterReached = false;
                        maxPntsReached = false;
                        
                        if any(strmatch(opt.display,{'all';'ALL'}))
                            str1 = [blanks(4) 'The approximations have converged'];
                            str2 = [blanks(4) 'Peak approximation error on critical samples: [' num2str(min(critFitErr)) ', ' num2str(max(critFitErr)) ']'];
                            DClab.dcdispstr(str1,guihand,false)
                            DClab.dcdispstr(str2,guihand,false)
                        end
                        
                    elseif iterIdx == opt.maxFitIter
                        % Maximum iterations reached
                        
                        critIdx = ~isnan(response(1:N)) & response(1:N) >= critRng(1) & response(1:N) <= critRng(2);
                        if ~any(critIdx)
                            disp('No critcal points found in training data for quadratic fit')
                        end
                        notDone=false;
                        poly=polyNew;
                        fitErr = fitErrNew;
                        critFitErr = fitErr(critIdx);
                        maxIterReached = true;
                        maxPntsReached = false;
                        
                        if any(strmatch(opt.display,{'all';'ALL'}))
                            str1 = [blanks(4) 'Maximum number of approximation tries reached'];
                            str2 = [blanks(4) 'Peak approximation error on critical samples: [' num2str(min(critFitErr)) ', ' num2str(max(critFitErr)) ']'];
                            DClab.dcdispstr(str1,guihand,false)
                            DClab.dcdispstr(str2,guihand,false)
                        end
                        
                    elseif N > opt.maxPnts4Fit
                        % Maximum number of points reached
                        
                        critIdx = ~isnan(response(1:N)) & response(1:N) >= critRng(1) & response(1:N) <= critRng(2);
                        if ~any(critIdx)
                            disp('No critcal points found in training data for quadratic fit')
                        end
                        maxIterReached = false;
                        maxPntsReached = true;
                        notDone = false;
                        poly = polyNew;
                        fitErr = fitErrNew;
                        critFitErr = fitErr(critIdx);
                        
                        
                        if ismember(opt.display,{'all';'ALL'})
                            str1 = [blanks(4) 'Maximum number of model evaluations reached'];
                            str2 = [blanks(4) 'Peak approximation error on critical samples: [' num2str(min(critFitErr)) ', ' num2str(max(critFitErr)) ']'];
                            DClab.dcdispstr(str1,guihand,false)
                            DClab.dcdispstr(str2,guihand,false)
                        end
                        
                    else
                        polyOld = polyNew;
                        iterIdx = iterIdx+1;
                    end
                    
                end % while notDone
                clear polyOld polyNew a dev err;
                
                %If we're allowed to and if more samples are available, do a super
                %good fit.
                if opt.useAllPnts4Fit == 1 && Nsamples > N
                    [poly fitErr] = feval(fitterHandle,transActParamValues,transActRanges,transResponse,opt);
                    N = Nsamples;
                    critIdx = ~isnan(response) & response >= critRng(1) & response <= critRng(2);
                    if ~any(critIdx)
                        disp('No critcal points found in training data for quadratic fit')
                    end
                    critFitErr = fitErr(critIdx);
                    
                    if ismember(opt.display,{'all';'ALL'})
                        str1 = [blanks(4) 'Generating final approximation with all available samples'];
                        str2 = [blanks(4) 'Peak approximation error on critical samples: [' num2str(min(critFitErr)) ', ' num2str(max(critFitErr)) ']'];
                        DClab.dcdispstr(str1)
                        DClab.dcdispstr(str2)
                    end
                end
                
            end %if opt.maxFitIter == 1
            
            if strcmp(Surf.responseTransformation,'log10')
                response = 10.^(transResponse);
            else
                response = transResponse;
            end
        end
        
        function [poly fitErr paramValues response N filesToSkip yMin yMax subspaceInfo] = oneShotFit(t)
            
            opt = t.opt;
            guihand = t.guihand;
            RespModel = t.RespModel;
            Metamodel = t.Metamodel;
            trnRng = t.trnRng;
            critRng = t.critRng;
            domain = t.domain;
            varTrans = t.varTrans;
            dsetNames2ModelNames = t.dsetNames2ModelNames;
            dsetNames2SharedNames = t.dsetNames2SharedNames;
            modelNames2SharedNames = t.modelNames2SharedNames;
            Surf = t.Surf;
            filesToSkip = t.filesToSkip;
            
            if strcmp(opt.analysisMode,'original') || ~strcmp(RespModel.type,'dcModel')
                metaBased = false;
            else
                metaBased = true;
                if strcmp(Surf.responseTransformation,'log10')
                    trnRng = 10.^(log10(trnRng)+Metamodel.fitInfo.peakErrorFromOptimization);
                else
                    trnRng = trnRng+Metamodel.fitInfo.peakErrorFromOptimization;
                end
            end
            
            if ismember(opt.display,{'all';'ALL'})
                str = [blanks(2) 'Generating surrogate for ' RespModel.name];
                DClab.dcdispstr(str,guihand,false)
            end
            
            %decide which fit constructor to use:
            switch opt.surfaceType
                case 'LinSurf'
                    fitterHandle = @DClab.DCSurface.constructLinFit;
                case 'QuadSurf'
                    fitterHandle = @DClab.DCSurface.constructQuadFit;
                    %case 'QuadOverQuad'
                    % fitterHan = @qratfit;
                    %fitterHan = @constructQOQFit;
                otherwise
                    error('DCSurface does not support that surfaceType')
            end
            
            
            
            if strcmpi(opt.subspaceDiscovery,'off')
                %no subspace discovery
                
                modelOrderRanges = vertcat(domain(dsetNames2ModelNames).range);
                modelOrderTrans = varTrans(dsetNames2ModelNames);
                
                % Determine number of coefficients in fitting polynomial
                n = size(Surf.activeParameterIndex,1);
                switch opt.surfaceType
                    case 'QuadSurf'
                        NCoeff = (n+1)*(n+2)/2; %returns how many coeffs required to fit a
                        %quad in n variables.
                    case 'LinSurf'
                        NCoeff = n+1;
                    case 'QuadOverQuad'
                        NCoeff = (n+1)*(n+2);
                    otherwise
                        error('Surface Type not supported by DCSurface class, check your opt.surfaceType property')
                end
                
                % Number of points to use for the fit
                N = opt.nPntsPerCoeff4OneShot*NCoeff;
                
                if strcmp(opt.display,'ALL')
                    str = [blanks(4) 'Gathering model evalutations'];
                    DClab.dcdispstr(str,guihand,false);
                end
                
                if metaBased
                    [paramValues, response] = randomEval(Metamodel,modelOrderRanges,modelOrderTrans,N);
                    paramValues = paramValues';
                    response = response';
                else
                    % Load any existing model evaluations
                    [paramValues, response, newFiles] = ...
                        loadSavedEvaluations(RespModel,modelOrderRanges,modelOrderTrans,N,filesToSkip);
                    
                    filesToSkip = [filesToSkip; newFiles];
                    
                    paramValues = paramValues';
                    response = response';
                    Nsamples = length(response);
                    
                    if Nsamples < N
                        Ns = N - Nsamples;
                        
                        try
                            [newParamValues newResponse newFile] = randomEvalWithSave(RespModel,modelOrderRanges,modelOrderTrans,Ns,opt.nComputer);
                        catch
                            str = 'randomEvalWithSave call failed in DCSurface.buildSurface';
                            DClab.dcdispstr(str,guihand,true)
                        end
                        filesToSkip = [filesToSkip; newFile];
                        
                        if any(isnan(newResponse))
                            error('Encounted NaN response values. Use Metamodel-based analysis')
                        end
                        newParamValues = newParamValues';
                        newResponse = newResponse';
                        paramValues = [paramValues; newParamValues];
                        response = [response; newResponse];
                        clear newParamValues;
                        clear newResponse; %save memory
                    end
                end
                
                %TODO in the iteractive fitting algorithm, we do not reject training data
                
                subspaceInfo = [];
                
                
                % toss out nontraining range points
                junk = isnan(response) | response <= trnRng(1) | response >= trnRng(2);
                yT = response(~junk);
                xT = paramValues(~junk,:);
                if length(yT) < 0.05*N
                    disp('rejected more than 95% of the training data, retrying with enlarged ranges')
                    trnRng = mean(trnRng)+1.3*(mean(trnRng)-trnRng(1))*[-1 1];
                    junk = isnan(response) | response <= trnRng(1) | response >= trnRng(2);
                    yT = response(~junk);
                    xT = paramValues(~junk,:);
                    if length(yT) < 0.05*N
                        disp('again rejected more than 95% of the training data, I give up')
                    end
                    
                end
                
                %Note, we should only have N points to fit with. We want to save any others
                %for validation.
                
                %TODO for speedup, remove some of these next 10 lines when tested
                
                %Both dsetNames2SharedNames and Surf.activeParameterIdx are already sorted.
                [trash sharedNames2ActiveNames] = intersect(dsetNames2SharedNames,Surf.activeParameterIndex);
                if ~isequal(trash,Surf.activeParameterIndex)
                    error('Inconsisteny: the shared parameters should be a superset of active')
                end
                modelNames2ActiveNames = modelNames2SharedNames(sharedNames2ActiveNames);
                
                modelNames = RespModel.parameterList;
                if ~isequal({domain(Surf.activeParameterIndex).name}',modelNames(modelNames2ActiveNames))
                    error('Problem with indexing')
                end
                
                % Convert newParamValues to the order of the active parameters and get rid
                % of nonactive ones.
                transActParamValues = xT(:,modelNames2ActiveNames);
                transActRanges = vertcat(domain(Surf.activeParameterIndex).range);
                %Transform these
                activeIsLog = Surf.activeParameterTransformation == 2;
                transActParamValues(:,activeIsLog) = log10(transActParamValues(:,activeIsLog));
                transActRanges(activeIsLog,:) = log10(transActRanges(activeIsLog,:));
                
            else
                %subspace discovery methods:
                
                assert(~metaBased,'Metabased methods do not work with Subspace methods')
                
                %find subspace then sample it
                [S,sigs] = DClab.DCSurface.discoverSubspace(t);
                [paramValues,activeParamValues,response,filesToSkip] = DClab.DCSurface.sampleSubspace(t,S);
                
                %store info
                subspaceInfo.S = S;
                subspaceInfo.singularValues = sigs;
                %size(S,2)
                
                %setup variables for fitting
                N = size(paramValues,1);
                transActParamValues = activeParamValues*S; %project down
                transActRanges = [min(transActParamValues)' max(transActParamValues)'];
                yT = response;
                
                
                transActRangesFull = vertcat(domain(Surf.activeParameterIndex).range);
                %Transform these
                activeIsLog = Surf.activeParameterTransformation == 2;
                transActRangesFull(activeIsLog,:) = log10(transActRangesFull(activeIsLog,:));
            end
            
            
            
            
            
            % Do one fit
            if strcmp(Surf.responseTransformation,'log10')
                [poly, fitErr] = feval(fitterHandle,transActParamValues,transActRanges,log10(yT),opt);
            else
                [poly, fitErr] = feval(fitterHandle,transActParamValues,transActRanges,yT,opt);
            end
            
            %transform poly to full dimension if fit was on subspace
            if strcmpi(opt.subspaceDiscovery,'on')
                %adjust the model for the subspace (including several transforms)
                az = 2*diag(1./diff(transActRanges,[],2));
                bz = -az*transActRanges(:,1)-1;
                axi = diag(diff(transActRangesFull,[],2))/2;
                bx = -2*diag(1./diff(transActRangesFull,[],2))*transActRangesFull(:,1)-1;
                r = size(transActRanges,1);
                n = size(transActRangesFull,1);
                transform = [1 zeros(1,r); bz az]*blkdiag(1,S')*[1 zeros(1,n); -axi*bx axi];
                %poly = blkdiag(1,S)*poly*blkdiag(1,S');
                poly = transform'*poly*transform;
                poly = 0.5*(poly+poly'); %just to be sure
            end
            
            
            critIdx = ~isnan(yT) & yT >= critRng(1) & yT <= critRng(2);
            
            if ~any(critIdx)
                disp('No critical points found in training data for quadratic fit')
            end
            
            critFitErr = fitErr(critIdx);
            
            if ismember(opt.display,{'all';'ALL'})
                str = [blanks(4) 'Peak approximation error on critical samples: [' num2str(min(critFitErr)) ', ' num2str(max(critFitErr)) ']'];
                DClab.dcdispstr(str,guihand,false)
                str = [blanks(4) 'In one-shot surface fitting mode'];
                DClab.dcdispstr(str,guihand,false)
            end
            
            yMin = min(response);
            yMax = max(response);
            response = yT;
        end
        
        function info = validate(t)
            %VALIDATE: private DCSurface function that creates estimates the error of a surrogate model.
            %
            %   [INFO] = VALIDATE(T)
            %   estimates the error of the surrogate model SURF.surrogateModel for the
            %   true model RESPMODEL.model.
            %
            %   Inputs:
            %      T.SURF is a partially complete DCSurface object
            %      T.RESPMODEL is the ResponseModel the surface is for.
            %      T.DOMAIN: an nx1 structure with fields .name and .range.
            %      T.VARTRANS: an nx1 cell array containing 'none' or 'log10'.
            %      T.OPT is a DCOptions object or structure with the required fields
            %
            %   Outputs:
            %      INFO is a structure array that is placed in the field
            %      .surrogateFitInfo of the final DCSurface object.
            %
            %   See also DCSurface, DCSurface.buildSurface
            opt = t.opt;
            guihand = t.guihand;
            RespModel = t.RespModel;
            % Metamodel = t.Metamodel;
            critRng = t.critRng;
            domain = t.domain;
            varTrans = t.varTrans;
            dsetNames2ModelNames = t.dsetNames2ModelNames;
            dsetNames2SharedNames = t.dsetNames2SharedNames;
            modelNames2SharedNames = t.modelNames2SharedNames;
            Surf = t.Surf;
            info = t.info;
            filesToSkip = t.filesToSkip;
            
            [trash sharedNames2ActiveNames] = intersect(dsetNames2SharedNames,Surf.activeParameterIndex);
            modelNames2ActiveNames = modelNames2SharedNames(sharedNames2ActiveNames);
            
            % if strcmp(opt.analysisMode,'metamodelBasedB') && strcmp(RespModel.type,'dcModel')
            %   metaBased = true;
            % else
            %   metaBased = false;
            % end
            
            if ismember(opt.display,{'all';'ALL'})
                str = [blanks(2) 'Validating approximation for ' RespModel.name];
                DClab.dcdispstr(str,guihand,false)
            end
            
            %determine number of coefficients in fitting polynomial
            n = size(Surf.activeParameterIndex,1);
            switch opt.surfaceType
                case 'QuadSurf'
                    NCoeff = (n+1)*(n+2)/2; %returns how many coeffs required to fit a
                    %quad in n variables.
                case 'LinSurf'
                    NCoeff = n+1;
                case 'QuadOverQuad'
                    NCoeff = (n+1)*(n+2);
                otherwise
                    error('Surface Type not supported by DCSurface class')
            end
            N = NCoeff*opt.nPntsPerCoeff4Validation;
            
            modelOrderRanges = vertcat(domain(dsetNames2ModelNames).range);
            modelOrderTrans = varTrans(dsetNames2ModelNames);
            
            if N > 0
                
                % To improve accuracy, increase the number of samples points to check if we're
                % performing a prediction on this response
                if RespModel.isPredictor
                    N = N*3;
                end
                
                % Evaluate the model
                %   if metaBased
                %    [paramValues, response] = randomEval(Metamodel,modelOrderRanges,modelOrderTrans,N);
                %   else
                % Load all points that are left
                [paramValues, response] = ...
                    loadSavedEvaluations(RespModel,modelOrderRanges,modelOrderTrans,inf,filesToSkip);
                Nsamples = length(response);
                if Nsamples < N
                    Ns = N - Nsamples;
                    
                    try
                        [newParamValues newResponse] = ...
                            randomEvalWithSave(RespModel,modelOrderRanges,modelOrderTrans,Ns,opt.nComputer);
                    catch
                        str = 'randomEvalWithSave call failed in DCSurface.validate';
                        DClab.dcdispstr(str,guihand,true)
                    end
                    %      if any(isnan(newResponse))
                    %        error('Encounted NaN response values. Use Metamodel-based analysis')
                    %      end
                    paramValues = [paramValues newParamValues];
                    response = [response newResponse];
                    clear newParamValues;
                    clear newResponse; %save memory
                end
                %   end
                
                %TODO in the first pass where we load saved evaluations, do we want to
                %eliminate NaNs so that we create more new points, or should we keep them
                %in to maintain the distribution of nan/not nans.
                
                % Get rid of NaN;
                junk = isnan(response);
                response = response(~junk);
                paramValues = paramValues(:,~junk);
                N = length(response);
                
                if isempty(vertcat(domain(Surf.activeParameterIndex).range))
                    normx = zeros(0,N);
                else
                    % Convert paramValues to the order of the active parameters and get rid
                    % of nonactive ones.
                    transActParamValues = paramValues(modelNames2ActiveNames,:);
                    transActRanges = vertcat(domain(Surf.activeParameterIndex).range);
                    %Transform these
                    for i1 = 1:size(Surf.activeParameterIndex,1)
                        if Surf.activeParameterTransformation(i1) == 2
                            transActParamValues(i1,:) = log10(transActParamValues(i1,:));
                            transActRanges(i1,:) = log10(transActRanges(i1,:));
                        end
                    end
                    if strcmp(Surf.responseTransformation,'log10')
                        transResponse = log10(response);
                    else
                        transResponse = response;
                    end
                    normx = 2*diag(1./diff(transActRanges,[],2))*(transActParamValues - repmat(transActRanges(:,1),1,N)) - 1;
                end
                xx = [ones(1,N); normx];
                
                if isstruct(Surf.surrogateModel)
                    yS = sum(xx.*(Surf.surrogateModel.num*xx),1)./sum(xx.*(Surf.surrogateModel.den*xx),1);
                else
                    yS = sum(xx.*(Surf.surrogateModel*xx),1);
                end
                
                err = yS - transResponse;
                
                % Reject noncritical points.
                if ~all(isinf(critRng))
                    
                    critIdx = ~isnan(response(1:N)) & response(1:N) >= critRng(1) & response(1:N) <= critRng(2);
                    if ~any(critIdx)
                        disp('No critcal points found. Temporarily enlarging critical range to find seed pts.')
                        %      keyboard
                    end
                    errCrit = err(critIdx);
                    paramValuesCrit = paramValues(:,critIdx);
                    
                    % Make sure we have enough points to initialize the local search.
                    tmp = critRng;
                    while sum(critIdx) <= 10
                        tmp = mean(tmp) + 1.25*(tmp-mean(tmp));
                        critIdx = ~isnan(response(1:N)) & response(1:N) >= tmp(1) & response(1:N) <= tmp(2);
                    end
                    
                    errSeed = err(critIdx);
                    paramValuesSeed = paramValues(:,critIdx);
                else
                    paramValuesSeed = paramValues;
                    paramValuesCrit = paramValues;
                    errSeed = err;
                    errCrit = err;
                end
                clear response
                
                % Find the vectors x that gave the top four errors
                xsamp_min = DClab.findExtremeX(paramValuesSeed,errSeed,4,'min');
                xsamp_max = DClab.findExtremeX(paramValuesSeed,errSeed,4,'max');
                
                [sampErrMax idxH] = max(errCrit);
                [sampErrMin idxL] = min(errCrit);
                
                % Make sure the first seed point is the best among the truely critical
                % points at the particular extrema.
                
                if ~isempty(errCrit) && ~isempty(info.xInModelOrderGivingErrMax)
                    if sampErrMax < info.peakErrorOnCriticalSamplePoints4Fit(2)
                        xsamp_max(:,1) = info.xInModelOrderGivingErrMax;
                    else
                        xsamp_max(:,1) = paramValuesCrit(:,idxH);
                    end
                elseif ~isempty(errCrit)
                    xsamp_max(:,1) = paramValuesCrit(:,idxH);
                elseif ~isempty(info.xInModelOrderGivingErrMax)
                    xsamp_max(:,1) = info.xInModelOrderGivingErrMax;
                else
                    disp('No feasible seed point')
                end
                
                if ~isempty(errCrit) && ~isempty(info.xInModelOrderGivingErrMin)
                    if sampErrMin > info.peakErrorOnCriticalSamplePoints4Fit(1)
                        xsamp_min(:,1) = info.xInModelOrderGivingErrMin;
                    else
                        xsamp_min(:,1) = paramValuesCrit(:,idxL);
                    end
                elseif ~isempty(errCrit)
                    xsamp_min(:,1) = paramValuesCrit(:,idxL);
                elseif ~isempty(info.xInModelOrderGivingErrMin)
                    xsamp_min(:,1) = info.xInModelOrderGivingErrMin;
                else
                    % do nothing. we've already seen the display
                end
                
                if ismember(opt.display,{'all';'ALL'})
                    str = [blanks(4) 'Peak approximation error on validation samples: [' num2str(sampErrMin) ', ' num2str(sampErrMax) ']'];
                    DClab.dcdispstr(str,guihand,false)
                end
                
                % Do twice as many restarts if if we're performing a prediction on this
                % ResponseModel
                if RespModel.isPredictor
                    NRE = 2*opt.nLocalValidationSearches;
                else
                    NRE = opt.nLocalValidationSearches;
                end
                
                optimizedErrMin = zeros(NRE,1);
                optimizedErrMax = zeros(NRE,1);
                
                %TODO: make these options configurable or at least decide on good defaults.
                options = optimset('fmincon');
                
                %   if metaBased
                %     options = optimset(options,'GradConstr','off','GradObj','off',...
                %           'LargeScale','off','MaxIter',15,'DiffMinChange',1e-8,'MaxFunEval',1000, ...
                %           'TolFun',1e-5,'MaxSQPIter',10000,'Display','off','TolCon',1e-6,'algorithm','active-set');
                %   else
                % Need to be a bit more careful with operating directly on an ODE model
                % We implement DiffMinChange and RelLineSrchBnd
                options = optimset(options,'GradConstr','off','GradObj','off',...
                    'display','off','RelLineSrchBnd',0.05,...
                    'RelLineSrchBndDuration',20,'DiffMinChange',1e-4,...
                    'MaxFunEval',1000,'largescale','off',...
                    'TolFun',1e-5,'MaxSQPIter',10000,'TolCon',1e-6,'algorithm','active-set');
                %   end
                
                xopt_min = zeros(size(xsamp_min,1),NRE);
                xopt_max = zeros(size(xsamp_min,1),NRE);
                for i1 = 1:NRE
                    %Use random seed for some cases, sample points giving big errors for
                    %others.
                    
                    try
                        if i1 <= min(ceil(NRE/2),4)
                            [optimizedErrMin(i1) xopt_min(:,i1) optimizedErrMax(i1) xopt_max(:,i1)] = optimError(Surf,RespModel,critRng,modelOrderRanges,modelNames2ActiveNames,xsamp_min(:,i1),xsamp_max(:,i1),options);
                        else
                            x = modelOrderRanges(:,1) + rand(size(modelOrderRanges,1),1).*diff(modelOrderRanges,[],2);
                            [optimizedErrMin(i1) xopt_min(:,i1) optimizedErrMax(i1) xopt_max(:,i1)] = optimError(Surf,RespModel,critRng,modelOrderRanges,modelNames2ActiveNames,x,x,options);
                        end
                    catch
                        disp('bombed')
                        keyboard
                    end
                    
                    if strcmp(opt.display,'ALL')
                        str = [blanks(4) 'Peak approximation error from local search' num2str(i1) ': [' num2str(optimizedErrMin(i1)) ', ' num2str(optimizedErrMax(i1)) ']'];
                        DClab.dcdispstr(str,guihand,false)
                    end
                end
                
                % Take the best
                [errLow idx] = min(optimizedErrMin);
                xopt_min = xopt_min(:,idx);
                [errHigh idx] = max(optimizedErrMax);
                xopt_max = xopt_max(:,idx);
                
                % Try one more time (allowing more iterations) with the best candidate x
                % value.
                options = optimset(options,'MaxIter',20,'MaxFunEval',2000,'TolFun',1e-6);
                if strmatch(opt.display,'ALL')
                    options = optimset(options,'Display','iter');
                else
                    options = optimset(options,'Display','off');
                end
                
                [optimizedErrMin xopt_min optimizedErrMax xopt_max] = optimError(Surf,RespModel,critRng,modelOrderRanges,modelNames2ActiveNames,xopt_min,xopt_max,options);
                
                if isinf(optimizedErrMin) || isinf(optimizedErrMax)
                    disp('no point found in critical set')
                    keyboard
                end
                
                if ismember(opt.display,{'all','ALL'})
                    str = [blanks(4) 'Peak approximation error from tight tol local search: [' num2str(optimizedErrMin) ', ' num2str(optimizedErrMax) ']'];
                    DClab.dcdispstr(str,guihand,false)
                end
                
                %   if metaBased
                %     Merr = Metamodel.fitInfo.peakErrorFromOptimization;
                %     if (Merr(1) < optimizedErrMin) || (Merr(2) < optimizedErrMin)
                %       disp('Metamodel refine indicated')
                %   %    keyboard
                %     end
                %
                %     optimizedErrMin = optimizedErrMin+Merr(1);
                %     optimizedErrMax = optimizedErrMax+Merr(2);
                %     if ismember(opt.display,{'all','ALL'})
                %       str = [blanks(4) 'Final approximation error accounting for Metamodel error: [' num2str(optimizedErrMin) ', ' num2str(optimizedErrMax) ']'];
                %       DClab.dcdispstr(str,guihand,false)
                %     end
                %   end
            elseif opt.nLocalValidationSearches > 0
                % No validation points. If we are still supposed to perform a local
                % search, start it from all random points.
                
                % dummy values
                errCrit = NaN;
                sampErrMin = NaN;
                sampErrMax = NaN;
                transResponse = NaN;
                
                if ismember(opt.display,{'all';'ALL'})
                    str = [blanks(4) 'No validation points requested'];
                    DClab.dcdispstr(str,guihand,false)
                end
                
                % Do twice as many restarts if if we're performing a prediction on this
                % ResponseModel
                if RespModel.isPredictor
                    NRE = 2*opt.nLocalValidationSearches;
                else
                    NRE = opt.nLocalValidationSearches;
                end
                
                optimizedErrMin = zeros(NRE,1);
                optimizedErrMax = zeros(NRE,1);
                
                %TODO: make these options configurable or at least decide on good defaults.
                options = optimset('fmincon');
                
                %   if metaBased
                %     options = optimset(options,'GradConstr','off','GradObj','off',...
                %           'LargeScale','off','MaxIter',15,'DiffMinChange',1e-8,'MaxFunEval',1000, ...
                %           'TolFun',1e-5,'MaxSQPIter',10000,'Display','off','TolCon',1e-6);
                %   else
                % Need to be a bit more careful with operating directly on an ODE model
                % We implement DiffMinChange and RelLineSrchBnd
                options = optimset(options,'GradConstr','off','GradObj','off',...
                    'display','off','RelLineSrchBnd',0.05,...
                    'RelLineSrchBndDuration',20,'DiffMinChange',1e-4,...
                    'MaxFunEval',1000,'largescale','off','MaxIter',20,...
                    'TolFun',1e-5,'MaxSQPIter',10000,'TolCon',1e-6);
                %   end
                
                ntotal = length(modelOrderRanges(:,1));
                
                xopt_min = zeros(ntotal,NRE);
                xopt_max = zeros(ntotal,NRE);
                
                xsamp_min = info.xInModelOrderGivingErrMin;
                xsamp_max = info.xInModelOrderGivingErrMax;
                
                %Use random seed for all except first search. Then use the points among the
                %training data with the worst errors.
                for i1 = 1:NRE
                    if i1 == 1
                        if isempty(xsamp_min)
                            xL = modelOrderRanges(:,1) + rand(size(modelOrderRanges,1),1).*diff(modelOrderRanges,[],2);
                        else
                            xL = xsamp_min;
                        end
                        if isempty(xsamp_max)
                            xH = modelOrderRanges(:,1) + rand(size(modelOrderRanges,1),1).*diff(modelOrderRanges,[],2);
                        else
                            xH = xsamp_max;
                        end
                        [optimizedErrMin(i1) xopt_min(:,i1) optimizedErrMax(i1) xopt_max(:,i1)] = optimError(Surf,RespModel,critRng,modelOrderRanges,modelNames2ActiveNames,xL,xH,options);
                    else
                        x = modelOrderRanges(:,1) + rand(size(modelOrderRanges,1),1).*diff(modelOrderRanges,[],2);
                        [optimizedErrMin(i1) xopt_min(:,i1) optimizedErrMax(i1) xopt_max(:,i1)] = optimError(Surf,RespModel,critRng,modelOrderRanges,modelNames2ActiveNames,x,x,options);
                    end
                    
                    if strcmp(opt.display,'ALL')
                        str = [blanks(4) 'Peak approximation error from local search' num2str(i1) ': [' num2str(optimizedErrMin(i1)) ', ' num2str(optimizedErrMax(i1)) ']'];
                        DClab.dcdispstr(str,guihand,false)
                    end
                end
                
                % Take the best
                [errLow idx] = min(optimizedErrMin);
                xopt_min = xopt_min(:,idx);
                [errHigh idx] = max(optimizedErrMax);
                xopt_max = xopt_max(:,idx);
                
                % Try one more time (allowing more iterations) with the best candidate x
                % value.
                options = optimset(options,'MaxIter',25,'MaxFunEval',2000,'TolFun',1e-6);
                if strmatch(opt.display,'ALL')
                    options = optimset(options,'Display','iter');
                else
                    options = optimset(options,'Display','off');
                end
                
                [optimizedErrMin xopt_min optimizedErrMax xopt_max] = optimError(Surf,RespModel,critRng,modelOrderRanges,modelNames2ActiveNames,xopt_min,xopt_max,options);
                
                if isinf(optimizedErrMin) || isinf(optimizedErrMax)
                    disp('no point found in critical set')
                    keyboard
                end
                
                if ismember(opt.display,{'all','ALL'})
                    str = [blanks(4) 'Peak approximation error from tight tol local search: [' num2str(optimizedErrMin) ', ' num2str(optimizedErrMax) ']'];
                    DClab.dcdispstr(str,guihand,false)
                end
                
                %   if metaBased
                %     Merr = Metamodel.fitInfo.peakErrorFromOptimization;
                %     if (Merr(1) < optimizedErrMin) || (Merr(2) < optimizedErrMin)
                %       disp('Metamodel refine indicated')
                %   %    keyboard
                %     end
                %
                %     optimizedErrMin = optimizedErrMin+Merr(1);
                %     optimizedErrMax = optimizedErrMax+Merr(2);
                %     if ismember(opt.display,{'all','ALL'})
                %       str = [blanks(4) 'Final approximation error accounting for Metamodel error: [' num2str(optimizedErrMin) ', ' num2str(optimizedErrMax) ']'];
                %       DClab.dcdispstr(str,guihand,false)
                %     end
                %   end
                
            else
                % dummy values
                errCrit = NaN;
                sampErrMin = NaN;
                sampErrMax = NaN;
                transResponse = NaN;
                optimizedErrMin = NaN;
                optimizedErrMax = NaN;
                xopt_min = [];
                xopt_max = [];
                %xsamp_min = info.xInModelOrderGivingErrMin;
                %xsamp_max = info.xInModelOrderGivingErrMax;
            end
            
            % Finish filing out info structure
            info.nSamplePoints4Validation = N;
            if isempty(errCrit)
                info.averageErrorOnSamplePoints4Validation = NaN;
            else
                info.averageErrorOnSamplePoints4Validation = mean(abs(errCrit));
            end
            info.peakErrorOnSamplePoints4Validation = [sampErrMin sampErrMax];
            info.peakErrorFromOptimization = [optimizedErrMin optimizedErrMax];
            %if optimizedErrMin < min(sampErrMin,info.peakErrorOnCriticalSamplePoints4Fit(1))
            info.xInModelOrderGivingErrMin = xopt_min;
            %else
            %  info.xInModelOrderGivingErrMin = xsamp_min(:,1);
            %end
            %if optimizedErrMax > max(sampErrMax,info.peakErrorOnCriticalSamplePoints4Fit(2))
            info.xInModelOrderGivingErrMax = xopt_max;
            %else
            %  info.xInModelOrderGivingErrMax = xsamp_max(:,1);
            %end
            
            if strcmp(Surf.responseTransformation,'log10')
                info.trueModelRangeEstimate(1) = min([10.^transResponse info.trueModelRangeEstimate(1)]);
                info.trueModelRangeEstimate(2) = max([10.^transResponse info.trueModelRangeEstimate(2)]);
            else
                info.trueModelRangeEstimate(1) = min([transResponse info.trueModelRangeEstimate(1)]);
                info.trueModelRangeEstimate(2) = max([transResponse info.trueModelRangeEstimate(2)]);
            end
            
            function [minErr xopt_min maxErr xopt_max] = optimError(Surf,RespModel,critRng,modelOrderRanges,modelNames2ActiveNames,xseed_min,xseed_max,options)
                % Function to minimize and maximize S(x)-M(x), where S(x) is the
                % surrogateModel and M(x) is ResponseModel.model. (M(x) will possibly
                % include an output transformation as well).
                %
                % Inputs:
                %   Surf: partially complete DCSurface object
                %   RespModel: Model we're building the surrogate for
                %   modelOrderRanges: nx2 vector containing the ranges of the parameters,
                %     with the parameter order dictated by RespModel.model
                %   modelNames2ActiveNames: a vector of indices such that if modelParams =
                %     RespModel.paramList, modelParams(modelNames2ActiveNames) will be the
                %     active parameters in the order they appear in the dataset.
                %   xseed_min: nx1 vector of parameters (in the order of M(x)) used to seed
                %     the minimization.
                %   xseed_max: nx1 vector of parameters (in the order of M(x)) used to seed
                %     the maximization.
                %   options: an options structure to be used by fmincon
                %
                % Outputs:
                %   minErr: local minimum of S(x)-M(x)
                %   xopt_min: nx1 vector of parameters (in the order of M(x)) that produce minErr
                %   maxErr: local maximum of S(x)-M(x)
                %   xopt_max: nx1 vector of parameters (in the order of M(x)) that produce maxErr
                
                % Accommodate situations where some of the model parameters are frozen. We
                % don't want to optimize these variables.
                
                % We optimize over all model parameters that have a nonsingleton range.
                % doesVaryBool is true if the corresponding element of the model parameter
                % list (in its model order) has a nonsingleton range.
                % logVaryBool XOR linVaryBool should equal doesVaryBool
                
                doesVaryBool = diff(modelOrderRanges,[],2) ~= 0;
                
                % if isempty(Metamodel)
                linVaryBool = true(size(doesVaryBool)); %default all true. any nonactive is assumed to be lin transformed
                logVaryBool = false(size(doesVaryBool));
                
                linVaryBool(modelNames2ActiveNames,1) = Surf.activeParameterTransformation==1;
                linVaryBool = linVaryBool & doesVaryBool;
                logVaryBool(modelNames2ActiveNames,1) = Surf.activeParameterTransformation==2;
                logVaryBool = logVaryBool & doesVaryBool;
                % else
                %   linVaryBool = (Metamodel.variableTransformations==1) & doesVaryBool;
                %   logVaryBool = (Metamodel.variableTransformations==2) & doesVaryBool;
                % end
                
                
                
                var2opt = find(doesVaryBool);
                lin2opt = find(linVaryBool);
                log2opt = find(logVaryBool);
                
                % var2opt, lin2opt, log2opt will always be sorted.
                [trash normxlin] = intersect(var2opt,lin2opt);
                [trash normxlog] = intersect(var2opt,log2opt);
                
                % Linearly transform coordinates of seed point that vary to be between +/-1
                normxmin = zeros(size(var2opt));
                normxmax = zeros(size(var2opt));
                
                normxmin(normxlin) = 2*(xseed_min(linVaryBool) - modelOrderRanges(linVaryBool,1))./diff(modelOrderRanges(linVaryBool,:),[],2)-1;
                normxmin(normxlog) = 2*(log10(xseed_min(logVaryBool)./modelOrderRanges(logVaryBool,1)))./log10(modelOrderRanges(logVaryBool,2)./modelOrderRanges(logVaryBool,1))-1;
                
                normxmax(normxlin) = 2*(xseed_max(linVaryBool) - modelOrderRanges(linVaryBool,1))./diff(modelOrderRanges(linVaryBool,:),[],2)-1;
                normxmax(normxlog) = 2*(log10(xseed_max(logVaryBool)./modelOrderRanges(logVaryBool,1)))./log10(modelOrderRanges(logVaryBool,2)./modelOrderRanges(logVaryBool,1))-1;
                opt2modx = @(x) optx2modelx(x,linVaryBool,logVaryBool,normxlin,normxlog,modelOrderRanges);
                
                % % To avoid a conditional statement in the local function that evaluates the
                % % error at a given x, determine the output transformation here.
                % switch Surf.responseTransformation
                %  case 'none'
                %    if metaBased
                %      hand1 = @minErr_meta_noneY;
                %    else
                %      hand1 = @minErr_noneY;
                %    end
                %  case 'log10'
                %    if metaBased
                %      hand1 = @minErr_meta_logY;
                %    else
                %      hand1 = @minErr_logY;
                %    end
                %  otherwise
                %   error('Internal inconsistency: condition should never occur')
                % end
                
                % TODO: when we have refined Metamodels, we need to do more work.
                
                % Create constraint function for critical range.
                % if metaBased
                %
                %   % Err objective functions
                %   switch Surf.responseTransformation
                %     case 'none'
                %       objfunH = @(x) -minErr_meta_noneY(x,opt2modx,Metamodel,modelOrderRanges,Surf,modelNames2ActiveNames);
                %       objfunL = @(x) minErr_meta_noneY(x,opt2modx,Metamodel,modelOrderRanges,Surf,modelNames2ActiveNames);
                %     case 'log10'
                %       objfunH = @(x) -minErr_meta_logY(x,opt2modx,Metamodel,modelOrderRanges,Surf,modelNames2ActiveNames);
                %       objfunL = @(x) minErr_meta_logY(x,opt2modx,Metamodel,modelOrderRanges,Surf,modelNames2ActiveNames);
                %     otherwise
                %       error('Internal inconsistency: condition should never occur')
                %   end
                %
                %   % Constraint function
                %   Merr = Metamodel.fitInfo.peakErrorFromOptimization;
                %   if all(isinf(critRng))
                %     constHand = [];
                %   elseif isinf(critRng(1))
                %     constHand = @(x) mmConstU(x,opt2modx,Metamodel,critRng(2)+Merr(2));
                %   elseif isinf(critRng(2))
                %     constHand = @(x) mmConstL(x,opt2modx,Metamodel,critRng(1)+Merr(1));
                %   else
                %     constHand = @(x) mmConst(x,opt2modx,Metamodel,critRng+Merr);
                %   end
                %
                % else
                % Err objective functions
                switch Surf.responseTransformation
                    case 'none'
                        objfunH = @(x) -minErr_noneY(x,opt2modx,RespModel,modelOrderRanges,Surf,modelNames2ActiveNames);
                        objfunL = @(x) minErr_noneY(x,opt2modx,RespModel,modelOrderRanges,Surf,modelNames2ActiveNames);
                    case 'log10'
                        objfunH = @(x) -minErr_logY(x,opt2modx,RespModel,modelOrderRanges,Surf,modelNames2ActiveNames);
                        objfunL = @(x) minErr_logY(x,opt2modx,RespModel,modelOrderRanges,Surf,modelNames2ActiveNames);
                    otherwise
                        error('Internal inconsistency: condition should never occur')
                end
                % Constraint function
                if all(isinf(critRng))
                    constHand = [];
                elseif isinf(critRng(1))
                    constHand = @(x) rmConstU(x,opt2modx,RespModel,critRng(2));
                elseif isinf(critRng(2))
                    constHand = @(x) rmConstL(x,opt2modx,RespModel,critRng(1));
                else
                    constHand = @(x) rmConst(x,opt2modx,RespModel,critRng);
                end
                
                % end
                
                nVary = sum(doesVaryBool);
                
                %Minimize the error S(x)-M(x)
                [x,fval] = fmincon(objfunL,normxmin,[],[],[],[],-1*ones(nVary,1),ones(nVary,1),constHand,options);
                
                if ~isempty(constHand) && max(constHand(x)) > 5e-3
                    % try again
                    options = optimset(options,'MaxIter',2*options.MaxIter);
                    [x,fval] = fmincon(objfunL,normxmin,[],[],[],[],-1*ones(nVary,1),ones(nVary,1),constHand,options);
                    
                    if max(constHand(x)) > 5e-3
                        % give up
                        fval = inf;
                    end
                end
                
                minErr = fval;
                
                %Convert xopt back to the model coordinates
                xopt_min = opt2modx(x);
                %xopt_min(doesVaryBool) = 0.5*(x+1).*diff(modelOrderRanges(doesVaryBool,:),[],2) + modelOrderRanges(doesVaryBool,1);
                
                %Maximize the error S(x)-M(x)
                [x,fval] = fmincon(objfunH,normxmax,[],[],[],[],-1*ones(nVary,1),ones(nVary,1),constHand,options);
                
                if ~isempty(constHand) && max(constHand(x)) > 5e-3
                    %try again
                    options = optimset(options,'MaxIter',2*options.MaxIter,'RelLineSrchBnd',options.RelLineSrchBnd/5);
                    [x,fval] = fmincon(objfunH,normxmax,[],[],[],[],-1*ones(nVary,1),ones(nVary,1),constHand,options);
                    
                    if max(constHand(x)) > 5e-3
                        % give up
                        fval = inf;
                    end
                end
                maxErr = -fval;
                
                %Convert xopt back to the model coordinates
                %xopt_max = mean(modelOrderRanges,2);
                %xopt_max(doesVaryBool) = 0.5*(x+1).*diff(modelOrderRanges(doesVaryBool,:),[],2) + modelOrderRanges(doesVaryBool,1);
                xopt_max = opt2modx(x);
            end
            
            function err = minErr_noneY(x,optx2modx,RespModel,modelOrderRanges,Surf,modelNames2ActiveNames)
                % Inputs:
                %   x: linearly transformed to be between +/- 1. These are in model order,
                %      with parameters having a singleton range eliminated.
                %   doesVary: nx1 logical array. Component will be true if the
                %      corresponding model parameter has a nonsingleton range
                %   RespModel:
                %   modelOrderRanges: nx2 matrix containing the ranges of all model
                %      parameters, in the model order
                %   Surf:
                %   modelNames2ActiveNames: a vector of indices such that if modelParams =
                %     RespModel.paramList, modelParams(modelNames2ActiveNames) will be the
                %     active parameters in the order they appear in the dataset.
                
                
                % Transform x to be in the model coordinates, fill in singletons where
                % needed.
                
                %optim2realx
                %optim2activex
                realx = optx2modx(x);
                %realx(doesVary) = 0.5*(x+1).*diff(modelOrderRanges(doesVary,:),[],2) +
                %modelOrderRanges(doesVary,1);
                
                % Transform realx to be in the active parameter coordinates.
                transActParamValues = realx(modelNames2ActiveNames);
                transActRanges = modelOrderRanges(modelNames2ActiveNames,:);
                %Transform these
                isLog = Surf.activeParameterTransformation == 2;
                transActParamValues(isLog) = log10(transActParamValues(isLog));
                transActRanges(isLog,:) = log10(transActRanges(isLog,:));
                
                normx = 2*diag(1./diff(transActRanges,[],2))*(transActParamValues-transActRanges(:,1))-1;
                err = [1 normx']*Surf.surrogateModel*[1;normx] - rapideval(RespModel,realx);
            end
            
            function err = minErr_logY(x,optx2modx,RespModel,modelOrderRanges,Surf,modelNames2ActiveNames)
                % Inputs: see minErr_noneY
                
                % Transform x to be in the model coordinates, fill in singletons where
                % needed.
                %realx = mean(modelOrderRanges,2);
                %realx(doesVary) = 0.5*(x+1).*diff(modelOrderRanges(doesVary,:),[],2) + modelOrderRanges(doesVary,1);
                
                realx = optx2modx(x);
                
                % Transform realx to be in the active parameter coordinates.
                transActParamValues = realx(modelNames2ActiveNames);
                transActRanges = modelOrderRanges(modelNames2ActiveNames,:);
                %Transform these
                isLog = Surf.activeParameterTransformation == 2;
                transActParamValues(isLog) = log10(transActParamValues(isLog));
                transActRanges(isLog,:) = log10(transActRanges(isLog,:));
                
                normx = 2*diag(1./diff(transActRanges,[],2))*(transActParamValues-transActRanges(:,1)) - 1;
                
                err = [1 normx']*Surf.surrogateModel*[1;normx] - log10(rapideval(RespModel,realx));
                
                
            end
            
            % function err = minErr_meta_noneY(x,optx2modx,Metamodel,modelOrderRanges,Surf,modelNames2ActiveNames)
            % % Inputs:
            % %   x: linearly transformed to be between +/- 1. These are in model order,
            % %      with parameters having a singleton range eliminated.
            % %   doesVary: nx1 logical array. Component will be true if the
            % %      corresponding model parameter has a nonsingleton range
            % %   RespModel:
            % %   modelOrderRanges: nx2 matrix containing the ranges of all model
            % %      parameters, in the model order
            % %   Surf:
            % %   modelNames2ActiveNames: a vector of indices such that if modelParams =
            % %     RespModel.paramList, modelParams(modelNames2ActiveNames) will be the
            % %     active parameters in the order they appear in the dataset.
            %
            %
            % % Transform x to be in the model coordinates, fill in singletons where
            % % needed.
            % %realx = mean(modelOrderRanges,2);
            % %realx(doesVary) = 0.5*(x+1).*diff(modelOrderRanges(doesVary,:),[],2) + modelOrderRanges(doesVary,1);
            %
            % realx = optx2modx(x);
            %
            % % Transform realx to be in the active parameter coordinates.
            % transActParamValues = realx(modelNames2ActiveNames);
            % transActRanges = modelOrderRanges(modelNames2ActiveNames,:);
            % %Transform these
            % isLog = Surf.activeParameterTransformation == 2;
            % transActParamValues(isLog) = log10(transActParamValues(isLog));
            % transActRanges(isLog,:) = log10(transActRanges(isLog,:));
            %
            % normx = 2*diag(1./diff(transActRanges,[],2))*(transActParamValues-transActRanges(:,1))-1;
            % err = [1 normx']*Surf.surrogateModel*[1;normx] - eval(Metamodel,realx);
            % end
            
            % function err = minErr_meta_logY(x,optx2modx,Metamodel,modelOrderRanges,Surf,modelNames2ActiveNames)
            % % Inputs: see minErr_noneY
            %
            % % Transform x to be in the model coordinates, fill in singletons where
            % % needed.
            % %realx = mean(modelOrderRanges,2);
            % %realx(doesVary) = 0.5*(x+1).*diff(modelOrderRanges(doesVary,:),[],2) + modelOrderRanges(doesVary,1);
            %
            % realx = optx2modx(x);
            %
            % % Transform realx to be in the active parameter coordinates.
            % transActParamValues = realx(modelNames2ActiveNames);
            % transActRanges = modelOrderRanges(modelNames2ActiveNames,:);
            % %Transform these
            % isLog = Surf.activeParameterTransformation == 2;
            % transActParamValues(isLog) = log10(transActParamValues(isLog));
            % transActRanges(isLog,:) = log10(transActRanges(isLog,:));
            %
            % normx = 2*diag(1./diff(transActRanges,[],2))*(transActParamValues-transActRanges(:,1)) - 1;
            %
            % err = [1 normx']*Surf.surrogateModel*[1;normx] - log10(eval(Metamodel,realx));
            
            % if isstruct(poly)
            %   err = ( [1 normx']*poly.num*[1;normx] )/ ...
            %     ( [1 normx']*poly.den*[1;normx] ) ...
            %     - log10(eval(Surf.ModelAssn,xx));
            % else
            %   err = [1 normx']*poly*[1;normx] - log10(eval(Surf.ModelAssn,xx));
            % end
            % end
            
            % function [cVal eVal] = mmConst(x,optx2modx,Metamodel,critRng)
            % % Inputs: see minErr_noneY
            %
            % % Transform x to be in the model coordinates, fill in singletons where
            % % needed.
            % %realx = mean(modelOrderRanges,2);
            % %realx(doesVary) = 0.5*(x+1).*diff(modelOrderRanges(doesVary,:),[],2) + modelOrderRanges(doesVary,1);
            %
            % y = eval(Metamodel,optx2modx(x));
            % cVal = [y-critRng(2); critRng(1)-y];
            % eVal = 0;
            % end
            
            % function [cVal eVal] = mmConstL(x,optx2modx,Metamodel,critRngL)
            % % Inputs: see minErr_noneY
            %
            % % Transform x to be in the model coordinates, fill in singletons where
            % % needed.
            %
            % y = eval(Metamodel,optx2modx(x));
            % cVal = critRngL-y;
            % eVal = 0;
            % end
            
            % function [cVal eVal] = mmConstU(x,optx2modx,Metamodel,critRngU)
            % % Inputs: see minErr_noneY
            %
            % % Transform x to be in the model coordinates, fill in singletons where
            % % needed.
            % %realx = mean(modelOrderRanges,2);
            % %realx(doesVary) = 0.5*(x+1).*diff(modelOrderRanges(doesVary,:),[],2) + modelOrderRanges(doesVary,1);
            %
            % y = eval(Metamodel,optx2modx(x));
            % cVal = y-critRngU;
            % eVal = 0;
            % end
            
            function [cval eval] = rmConst(x,optx2modx,RespModel,critRng)
                % Inputs: see minErr_noneY
                
                y = rapideval(RespModel,optx2modx(x));
                cval = [y-critRng(2); critRng(1)-y];
                eval = 0;
            end
            function [cval eval] = rmConstL(x,optx2modx,RespModel,critRngL)
                % Inputs: see minErr_noneY
                
                y = rapideval(RespModel,optx2modx(x));
                cval = critRngL-y;
                eval = 0;
            end
            function [cval eval] = rmConstU(x,optx2modx,RespModel,critRngU)
                % Inputs: see minErr_noneY
                
                y = rapideval(RespModel,optx2modx(x));
                cval = y-critRngU;
                eval = 0;
            end
            function modelx = optx2modelx(x,linVary,logVary,normxlin,normxlog,modelOrderRanges)
                
                modelx = mean(modelOrderRanges,2);
                modelx(linVary) = 0.5*(x(normxlin)+1).*diff(modelOrderRanges(linVary,:),[],2) + modelOrderRanges(linVary,1);
                modelx(logVary) = modelOrderRanges(logVary,1).*10.^(0.5*(x(normxlog)+1).*log10(modelOrderRanges(logVary,2)./modelOrderRanges(logVary,1)));
            end
        end
        
        function [dsetNames2ModelNames,dsetNames2SharedNames,modelNames2SharedNames] =  getIdxVectors(dsetParamNames,modelParamNames)
            %GETIDXVECTORS returns several useful indexing vectors
            %
            %   IDX1 = GETIDXVECTORS(DSETPARAMNAMES,MODELPARAMNAMES) returns IDX1 such
            %   that DSETPARAMNAMES(IDX1) == MODELPARAMNAMES. It is assumed that
            %   MODELPARAMNAMES is a subset of DSETPARAMNAMES. This condition is not
            %   verified in this function. Both inputs should be column cell arrays.
            %
            %   [IDX1 IDX2] = GETIDXVECTORS(DSETPARAMNAMES,MODELPARAMNAMES)
            %   additionally returns IDX2 such that DSETPARAMNAMES(IDX2) is the shared
            %   parameters that are common to both inputs in the dataset parameter
            %   order.
            %
            %   [IDX1 IDX2 IDX3] = GETIDXVECTORS(DSETPARAMNAMES,MODELPARAMNAMES)
            %   additionally returns IDX3 such that MODELPARAMNAMES(IDX3) =
            %   sharedParamNames, where sharedParamNames are those common to both
            %   inputs in the dataset parameter order.
            %
            %   IDX1 = dsetNames2ModelNames
            %   IDX2 = dsetNames2SharedNames
            %   IDX3 = modelNames2SharedNames
            %
            
            [trash dsetNames2sortedSharedNames modelNames2sortedSharedNames] = intersect(dsetParamNames,modelParamNames);
            dsetNames2ModelNames(modelNames2sortedSharedNames,1) = dsetNames2sortedSharedNames;
            
            if nargout > 1
                [dsetNames2SharedNames sortedSharedNames2SharedNames] = sort(dsetNames2sortedSharedNames);
                modelNames2SharedNames = modelNames2sortedSharedNames(sortedSharedNames2SharedNames);
            end
            
            %These index vectors are incredibly cryptic. Suppose the dataset parameters
            %are
            %b a g h c k d (=dsetParamNames)
            %and the model parameters are
            %d g k a c (=modelParamNames)
            %
            %then the sorted shared parameters are a c d g k so
            %dsetNames2sortedSharedNames is [2 5 7 3 6] and
            %modelNames2sortedSharedNames is [4 5 1 2 3]
            %
            %Thus dsetNames2ModelNames is [3 6 2 5 7] and
            %modelParamNames = dsetParamNames(dsetNames2ModelNames)
            %
            %The shared dataset parameters (in the dataset order) are a g c k d. These
            %are given by dsetParamNames([2 3 5 6 7]) and thus
            %dsetNames2SharedNames = [2 3 5 6 7]
            %
            %Additionally sortedSharedNames2SharedNames = [1 4 2 5 3]
            %
            %To get the modelNames in the order of the sharedNames we need
            %modelNames([4 2 5 3 1]) thus modelNames2SharedNames = [4 2 5 3 1]
            
            
        end
        
        function [coeff,fitError,EXITFLAG] = infNormRegression(x,pred)
            % function to solve the following linear program:
            %
            %  min_{theta,gamma} gamma
            %  s.t. |A*theta - b| <= gamma
            %
            %  the inputs are pred = b, A = x
            %
            %  we use this to fit a response surface to input output data.
            
            
            [n, m] = size(x);
            [o, p] = size(pred);
            
            if n ~=o || p ~= 1
                error('you have dimension issues')
            end
            
            f = zeros(m+1,1);
            f(1) = 1;
            
            %x*coeff-pred<= gamma
            A1 = [-ones(n,1) x];
            B1 = pred;
            
            %-gamma <= x*coeff-pred
            A2 = [-ones(n,1) -x];
            B2 = -pred;
            
            A = [A1; A2];
            B = [B1; B2];
            
            %if ispc
            options = optimset('linprog');
            options.MaxIter = 5000;
            options.Diagnostics = 'on';
            options.LargeScale = 'off';
            options.Display = 'off';
            options.Diagnostics = 'off';
            try
            [xopt, fval, EXITFLAG] = linprog(f,A,B,[],[],[],[],[],options);
            catch
                fprintf('Memory error?\nSize of A: %0.2gx%0.2g\n\n',size(A))
            end
            %else
            %  [xopt, fval, EXITFLAG] = mylinprog(f,A,B);
            %end
            
            coeff = xopt(2:end);
            
            fitError = x*coeff - pred;
        end
        
        function [S,sigs,nPnts] = discoverSubspace(t)
            
            %inputs
            opt = t.opt;
            guihand = t.guihand;
            RespModel = t.RespModel;
            domain = t.domain;
            varTrans = t.varTrans;
            dsetNames2ModelNames = t.dsetNames2ModelNames;
            dsetNames2SharedNames = t.dsetNames2SharedNames;
            modelNames2SharedNames = t.modelNames2SharedNames;
            Surf = t.Surf;
            
            nNewPntsEvald = 0;
            
            %we'll want to do the gradient calculations in a normalized space
            
            %reorder for the model
            modelOrderRanges = vertcat(domain(dsetNames2ModelNames).range);
            modelOrderTrans = varTrans(dsetNames2ModelNames);
            
            % Determine number of coefficients in fitting polynomial
            n = size(Surf.activeParameterIndex,1);
            
            %active variable indexing:
            %Both dsetNames2SharedNames and Surf.activeParameterIdx are already sorted.
            [trash sharedNames2ActiveNames] = intersect(dsetNames2SharedNames,Surf.activeParameterIndex);
            if ~isequal(trash,Surf.activeParameterIndex)
                error('Inconsisteny: the shared parameters should be a superset of active')
            end
            modelNames2ActiveNames = modelNames2SharedNames(sharedNames2ActiveNames);
            
            modelNames = RespModel.parameterList;
            if ~isequal({domain(Surf.activeParameterIndex).name}',modelNames(modelNames2ActiveNames))
                error('Problem with indexing')
            end
            activeIsLog = Surf.activeParameterTransformation == 2;
            
            activeOrderRanges = modelOrderRanges(modelNames2ActiveNames,:);
            activeOrderRanges(activeIsLog,:) = log10(activeOrderRanges(activeIsLog,:));
            
            % Number of points to use for the fit
            N = (n+1)*n;
            
            if strcmp(opt.display,'ALL')
                str = [blanks(4) 'Gathering model evalutations'];
                DClab.dcdispstr(str,guihand,false);
            end
            
            %scaling transformations:
            A = diag(2./diff(activeOrderRanges,[],2));
            bmat = repmat(sum(activeOrderRanges,2)'./diff(activeOrderRanges,[],2)',n+1,1);
            
            % Load any existing model evaluations
            [paramValues, response, newFiles4Grad] = ...
                loadSavedEvaluations4Grad(RespModel,modelOrderRanges,modelOrderTrans,N,[]);
            
            filesToSkip4Grad = newFiles4Grad;
            
            paramValues = paramValues';
            response = response';
            Nsamples = length(response);
            
            if Nsamples < 2*(n+1)
                Ns = (2*(n+1) - Nsamples)/(n+1);
                
                try
                    [newParamValues newResponse newFile4Grad] = randomEvalWithSave4Grad(RespModel,modelOrderRanges,modelOrderTrans,Ns,opt.derivRange,opt.nComputer);
                catch
                    keyboard
                    str = 'randomEvalWithSave4Grad call failed in DCSurface.buildSurface';
                    DClab.dcdispstr(str,guihand,true)
                end
                filesToSkip4Grad = [filesToSkip4Grad; newFile4Grad];
                
                if any(isnan(newResponse))
                    error('Encounted NaN response values. Use Metamodel-based analysis')
                end
                newParamValues = newParamValues';
                newResponse = newResponse';
                paramValues = [paramValues; newParamValues];
                response = [response; newResponse];
                clear newParamValues;
                clear newResponse; %save memory
            end
            
            % Convert newParamValues to the order of the active parameters and get rid
            % of nonactive ones.
            paramValues = paramValues(:,modelNames2ActiveNames);
            paramValues(:,activeIsLog) = log10(paramValues(:,activeIsLog));
            
            %compute gradients
            gradient = zeros(n,n);
            N = floor(size(response,1)/(n+1));
            for i1=1:N
                idx = (i1-1)*(n+1)+1;
                vec = [ones(n+1,1) -bmat + paramValues(idx:(idx+n),:)*A]\response(idx:(idx+n));
                gradient(i1,:) = vec(2:end)';
            end
            
            %when do we stop looking for a subspace?
            searchLimit = sqrt(2*n^2+4*n+9/4) -n -3/2;
            
            %begin iteration of derivatives
            while LOCALgetRank(gradient(1:N,:),opt)==N && N<searchLimit
                
                N = N+1;
                
                %add another gradient
                try
                    [newParamValues newResponse newFile4Grad] = randomEvalWithSave4Grad(RespModel,modelOrderRanges,modelOrderTrans,1,opt.derivRange,opt.nComputer);
                catch
                    str = 'randomEvalWithSave4Grad call failed in DCSurface.buildSurface';
                    DClab.dcdispstr(str,guihand,true)
                end
                filesToSkip4Grad = [filesToSkip4Grad; newFile4Grad];
                
                if any(isnan(newResponse))
                    error('Encounted NaN response values. Use Metamodel-based analysis')
                end
                
                newParamValues = newParamValues(modelNames2ActiveNames,:);
                newParamValues(activeIsLog,:) = log10(newParamValues(activeIsLog,:));
                newSz = size(newParamValues,2);
                nNewPntsEvald = nNewPntsEvald + newSz;
                bmat = repmat(sum(activeOrderRanges,2)'./diff(activeOrderRanges,[],2)',newSz,1);
                vec = [ones(newSz,1) -bmat + newParamValues'*A]\newResponse';
                gradient(N,:) = vec(2:end)';
                
            end
            
            [r,S,sigs] = LOCALgetRank(gradient(1:N,:),opt);
            S = A*S; %scale back to original units
            nPnts = nNewPntsEvald;
            
            
            
            function [r,S,sigs] = LOCALgetRank(gradient,opt)
                
                [u,s,v] = svd(gradient);
                sigs = diag(s);
                
                for i2=2:length(sigs)
                    if opt.subspaceThreshold*sigs(1)>sigs(i2) %opt.subspaceThreshold*sigs(i1-1)>sigs(i1)
                        r = i1-1;
                        break;
                    else
                        r = i1;
                    end
                end
                
                S = v(:,1:r);
                
            end
        end
        
        function [design,activeDesign,response,filesToSkip] = sampleSubspace(t,S)
            
            
            opt = t.opt;
            RespModel = t.RespModel;
            domrng = vertcat(t.domain.range);
            varTrans = t.varTrans;
            dsetNames2ModelNames = t.dsetNames2ModelNames;
            dsetNames2SharedNames = t.dsetNames2SharedNames;
            modelNames2SharedNames = t.modelNames2SharedNames;
            filesToSkip = t.filesToSkip;
            domain = t.domain;
            Surf = t.Surf;
            
            
            domrng = vertcat(domrng(dsetNames2ModelNames,:));
            modelOrderTrans = varTrans(dsetNames2ModelNames);
            
            %Both dsetNames2SharedNames and Surf.activeParameterIdx are already sorted.
            [trash sharedNames2ActiveNames] = intersect(dsetNames2SharedNames,Surf.activeParameterIndex);
            if ~isequal(trash,Surf.activeParameterIndex)
                error('Inconsisteny: the shared parameters should be a superset of active')
            end
            modelNames2ActiveNames = modelNames2SharedNames(sharedNames2ActiveNames);
            
            modelNames = RespModel.parameterList;
            if ~isequal({domain(Surf.activeParameterIndex).name}',modelNames(modelNames2ActiveNames))
                error('Problem with indexing')
            end
            
            activeIsLog = Surf.activeParameterTransformation == 2;
            
            
            [n,r] = size(S);
            
            %Now we need to sample the polytope.
            switch opt.surfaceType
                case 'QuadSurf'
                    NCoeff = (r+1)*(r+2)/2; %returns how many coeffs required to fit a
                    %quad in n variables.
                case 'LinSurf'
                    NCoeff = r+1;
                case 'QuadOverQuad'
                    NCoeff = (r+1)*(r+2);
                otherwise
                    error('Surface Type not supported by DCSurface class, check your opt.surfaceType property')
            end
            % Number of points to use for the fit
            N = opt.nPntsPerCoeff4OneShot*NCoeff;
            
            
            %First we load saved evals, and see how many we can use:
            [design, response, newFiles] = ...
                loadSavedEvaluations(RespModel,domrng,modelOrderTrans,N,filesToSkip);
            filesToSkip = [filesToSkip; newFiles];
            activeDesign = design(modelNames2ActiveNames,:);
            activeDesign(activeIsLog,:) = log10(activeDesign(activeIsLog,:));
            %and pull design from there
            if ~isempty(design) && strcmp(opt.surfaceFittingMode,'one-shot') && size(design,2)>=NCoeff
                designIdx = LOCALgetDesign(S'*activeDesign,opt);
            else
                designIdx = [];
            end
            
            designNew = [];
            
            
            
            %heavily sample H
            doesVaryIdx = find(diff(domrng,[],2) > 0);
            nvary = length(doesVaryIdx);
            tmp = DClab.lhsdesign(N*n,nvary,'smooth','off','criterion','none');
            %log10Trans = strmatch('log10',char(varTrans(doesVaryIdx)),'exact');
            %noTrans = strmatch('none',char(varTrans(doesVaryIdx)),'exact');
            
            scale = diag(diff(domrng(doesVaryIdx,:),1,2));
            x(doesVaryIdx,:) = repmat(domrng(doesVaryIdx,1),1,N*n) + scale*tmp';
            xInActiveOrder = x;
            xInActiveOrder = xInActiveOrder(modelNames2ActiveNames,:);
            xInActiveOrder(activeIsLog,:) = log10(xInActiveOrder(activeIsLog,:));
            
            %scale = diag(diff(domrng(doesVaryIdx(noTrans),:),1,2));
            %x(doesVaryIdx(noTrans),:) = repmat(domrng(doesVaryIdx(noTrans),1),1,N*n) + scale*tmp(:,noTrans)';
            
            %scale = diag(log10(domrng(doesVaryIdx(log10Trans),2)./domrng(doesVaryIdx(log10Trans),1)));
            %x(doesVaryIdx(log10Trans),:) = 10.^(repmat(log10(domrng(doesVaryIdx(log10Trans),1)),1,N*n) + scale*tmp(:,log10Trans)');
            
            
            
            if length(designIdx)<NCoeff && strcmp(opt.surfaceFittingMode,'one-shot') %design didn't give enough... (maybe loaded point weren't any good)
                
                %ok now have "dense" sample x (really? how dense?)
                %solve SDP to get design
                designIdxNew = LOCALgetDesign(S'*xInActiveOrder,opt);
                activeDesignNew = xInActiveOrder(:,designIdxNew);
                designNew = x(:,designIdxNew);
            else
                designIdxNew = [];
            end
            
            
            if size(design,2) + length(designIdxNew)<N
                more = N-length(designIdx);
                idxLeft = setdiff(1:N,designIdx);
                idxOfLeft = ceil(length(idxLeft)*rand(1,more));
                moreIdx = idxLeft(idxOfLeft);
                activeDesignNew = [activeDesignNew xInActiveOrder(:,moreIdx)];
                designNew = [designNew x(:,moreIdx)];
                design = [design designNew];
                activeDesign = [activeDesign activeDesignNew];
            end
            
            %ok, now actually you need the ys
            if ~isempty(designNew)
                [responseNew newFile] = evalWithSave(RespModel,designNew,domrng,modelOrderTrans,opt.nComputer);
                filesToSkip = [filesToSkip; newFile];
                response = [response responseNew];
            end
            
            response = response';
            design = design';
            activeDesign = activeDesign';
            
            
            
            
            
            
            % sample H
            % project down to subspace
            % do variance SDP problem...
            % then pick randomly from remaining points to get desired number
            
            % evaluate points, and send back...
            
            
            
            %maybe mode two that only does random? (for iterative version)
            
            function idx = LOCALgetDesign(x,opt)
                
                [n,N] = size(x);
                
                idx = ceil(N*rand((n+1)*(n+2)/2,1));
                return
                
                
                
                
                
                %%%%%%%%%%%%%%%%%%%%%%
                %THIS WHOLE PROBLEM SEEMS TO RUN OUT OF MEMORY.... SO WON'T DO IT FOR NOW
                %%%%%%%%%%%%%%%%%%%%%%
                
                %SHOULD ALSO NORMALIZE VARIABLES
                
                %x is n-by-nevals
                [n,N] = size(x);
                
                %subtract off mean to center points
                mx = mean(x,2);
                x0 = x-repmat(mx,1,N);
                M = max(abs(x0),[],2);
                x0 = x0./repmat(M,1,N);
                
                %convert to regression basis
                switch opt.surfaceType
                    case 'LinSurf'
                        v = [ones(1,N);x0];
                    case 'QuadSurf'
                        v = [ones(1,N);x0;zeros(n*(n+1)/2,N)];
                        offset = 0;
                        for i1=1:n
                            for i2=1:i1
                                v(n+1+offset+i2,:) = x0(i1,:).*x0(i2,:);
                            end
                            offset = offset + i1;
                        end
                    case 'QuadOverQuad'
                        % What the heck should these be?  Same as Quad?
                        v = [ones(1,N);x0;zeros(n*(n+1)/2,N)];
                        offset = 0;
                        for i1=1:n
                            for i2=1:i1
                                v(n+1+offset+i2,:) = x0(i1,:).*x0(i2,:);
                            end
                            offset = offset + i1;
                        end
                    otherwise
                        error('Invalid fit type')
                end
                
                [n,N] = size(v);  %redefine n
                
                
                %Setup SDP
                % y = [lambda1, lambda2, ..., lambdaN, t];
                % Will solve prob:
                %  max b'*y
                %  s.t.  c-At*y >=_K 0
                
                %objective
                b = spalloc(N+1,1,1);
                b(end) = 1;
                
                %linear constraints
                K.l = N+1;
                clin = spalloc(N+1,1,1);
                clin(end) = 1;
                Alin = spalloc(N+1,N+1,2*N+1);
                Alin(1:N,1:N) = -speye(N);
                Alin(N+1,1:N) = ones(1,N);
                
                %lmi constraints
                K.s = n;
                clmi = spalloc(n^2,1,0);
                Almi = zeros(n^2,N+1);
                for i1=1:N
                    Almi(:,i1) = -reshape(v(:,i1)*v(:,i1)',n^2,1);
                end
                Almi(:,end) = reshape(eye(n),n^2,1);
                Almi = sparse(Almi);
                
                c = [clin;clmi];
                At = [Alin;Almi];
                
                switch opt.display
                    case {'all','ALL'}
                        pars.fid = 1;
                    otherwise
                        pars.fid = 0;
                end
                
                
                [x,y,info] = sedumi(At',b,c,K,pars);
                
                if info.pinf==1 || info.dinf==1
                    error('The design SDP was infeasible')
                end
                
                [sy,si] = sort(y(1:N));
                idx = si(cumsum(sy)<(0.999));
                % design = x(:,idx);
            end
        end
    end %private,static methods
    
    
end %classdef