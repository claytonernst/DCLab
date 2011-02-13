classdef DCOptions < DClab.DCObject
    %OPTIONS
    %
    %   object = DCOptions creates a default object
    %
    %   object = DCOptions(OPTION,VALUE,OPTION,VALUE,...)
    %
    %   object = DCOptions('gui') launches a GUI used to interactively specify
    %   the options.
    %
    %   object = DCOptions(structin) creates an object with the property names
    %   present in the structure structin set to the corresponding values. The
    %   remaining properties are set to the default values.
    %
    %   The constructor for an Options object to be used with many of
    %   the functions in the data collaboration toolbox.  The options
    %   can be set by object.OPTION=VALUE or
    %   set(object,'OPTION',VALUE), and similiarly options values can
    %   be extracted from the object.  The options generally below to four
    %   categories: display, fitting, optimization, and diagnostics. See
    %   DCOptions/decompose to see how each option is classified.
    %
    %   The following is a list of the options, their default values,
    %   and a description.
    %
    %     display: ['iter']  Level of display. Valid values are
    %        'off', no output will be displayed
    %        'iter', progess of branch&bound iterations will be displayed
    %        'final', the final results of a calculation
    %        'notify', display will only occur if weird stuff happens
    %        'all', all progress is displayed
    %        'ALL', additionally display results of canned optimization calls:
    %               sedumi, fmincon iterations, etc.
    %     guiHandle: [] Handle to listbox in a gui to optionally send messages
    %        to a gui instead of the command window
    %     plotFitProgress: ['off']  Turns on plots generated during the
    %        construction of DCSurface objects. These display the convergence
    %        of the fits.  (other option is 'on')
    %
    %     surfaceType: ['QuadSurf']  Type of surface to use to fit the
    %        models with in the Cube object.  Default is 'QuadSurf', a
    %        quadratic surface object.  (other options are 'LinSurf' and
    %        'QuadOverQuad')
    %     subspaceDiscovery: ['off']  When set to 'on', fitting algorithm
    %        involves discovering a subspace on which the models vary.
    %        Uses subspaceThreshold to decide where to make a cutoff of
    %        principle subspace directions
    %     subspaceThreshold: [0.1]  The subspace discovery algorithm stops
    %        when eigenvalue (from the subspace discovery algorithm) is
    %        suspaceThreshold of the previous value or less.
    %     derivRange: [0.01] When evaluating points for derivative
    %        calculations, the spread will be derivRange of the parameter
    %        ranges around the point of interested
    %     surfaceTransformation:  [{'linXlinY'}]  This is a cell array of
    %        strings listing  the different transformations to be used when
    %        making a fit.  This must be a 'superset' of constraints.  Valid
    %        values are any cell array containing one or more of the
    %        following: 'linXlinY', 'logXlinY', 'linXlogY', and 'logXlogY'.
    %     fitNorm: [2]  This is the type of norm to use when creating
    %        the model fits (the surrogate models).  The fit tries to
    %        reduce the error of the fitted points based on this norm.
    %        The default is inf, for the infinity of sup norm. The other
    %        options is 2, for the 2-norm.
    %     minFitIter: [3]  Similar to maxFitIter this is the minimum
    %        number of times the algorithm will try to improve the fit.
    %     maxFitIter: [7]  When making a surrogate fit to the true
    %        models the surface object makes many fits until the
    %        fitting error converges.  This is a maximum on the number
    %        of time it will try to improve the fit.
    %     nSuccessfulFitIter: [2]  The number of times successive iterations
    %        a fit must meet the tolerance on the error convergence
    %        described below, before it is considered a good fit.
    %     maxPnts4Fit: [inf]  This option is turned on and off internally
    %        and is not meant to be used by the user.  Please leave it
    %        with the default value of inf.
    %     useAllPnts4Fit: [1]  This flag is used when making fits.  If
    %        useAllPnts is 1 then once the fit has converged, if more
    %        evaluations of the real models have been computed than
    %        were necessary for the fit, a 'better' fit is made using
    %        all available points.  If many more points are stored in
    %        the evaluations directory than are needed (possibly from
    %        previous runs) then if all the points are used, fitting may
    %        take a much longer time. (other option is 0)
    %     fitConvergenceTol: [0.05]  A tolerance used for assessing convergence of a
    %        model approximation to the actual function. A normalized
    %        version of either the two or infinity norm of successive
    %        approximations must be below tol
    %     nPntsPerParam4ActiveParamSel: [25]  Determination of active parameters
    %        will consider the (number of parameters) x (nPntsPerCoeff4ActiveParamSel)
    %        model evaluations.
    %     activeParamSelCutOff: [0.05]  A tolerance used for truncating
    %        the active parameter list. Setting this value to zero will
    %        insure all model parameters are used in the corresponding
    %        surrogate model. 100*activeParamSelCutOff roughly corresponds to the
    %        percent error introduced by eliminating nonactive parameters. We
    %        suggest never exceeding 0.2.
    %     nPntsPerCoeff4Validation: [10]  Determination of validation error
    %        will consider the (number of parameters) x (nPntsPerCoeff4Validation)
    %        model evaluations.
    %     nLocalValidationSearches: [3]  Number of times to estimate the
    %        maximum fitting error of a quadratic surface using a local search.
    %        Ceil(nLocalValidationSearchs/2) attempts are made using an ititial
    %        guess from the top sample points, while the remaining attempts are
    %        initiallized from random points.
    %     nLocalValidationSearches4Metamodel: [3]  Number of times to estimate
    %        the maximum fitting error of a metamodel using a local search.
    %        Ceil(nLocalValidationSearchs/2) attempts are made using an ititial
    %        guess from the top sample points, while the remaining attempts are
    %        initiallized from random points.
    %     maxDenSwing: [5] When fitting with ratios of quadratics, this
    %        limits the denominator in some way. How?
    %     nComputer: [0]  The number of computers to be used in the
    %        parallel distribution scheme.
    %
    %     omitInnerBound: [0]  For functions that create bounds on an answer
    %        (e.g. prediction, consistency) setting this to 1 (true) will cause
    %        the function to not caculate the inner bounds.
    %     omitOuterBound: [0]  Same as above but for outer bounds.
    %     maxBranchBoundIter: [1] Maximum iterations the branch and bound
    %        algorithm will execute.  This needs to be an integer between 1 and
    %        inf.
    %     branchBoundTermTol: [0.02] Termination tolerance for a branch and bound
    %        routine. The code will halt when the upper and lower bounds on the
    %        optimum value differ by less than this amount.
    %     branchBoundQGapTol: [0.01] Tolerance at which the algorithm switches
    %        from trying to reduce the fitting error to trying to reduce the
    %        duality gap. If the estimated q-gap (from local search) is less
    %        then branchBoundQGapTol, the algorithm subdivides so as to reduce
    %        the duality gap.
    %     nRestart: [2]  Number of times to restart any local search routines
    %        performed by fmincon
    %     tolFun: [1e-5] Tolerance passed to fmincon as TolFun using optimset.
    %     tolCon: [1e-5] Tolerance passed to fmincon as TolCon and TolX (most of the time...we're not consistent).
    %     sedumiParEps: [1e-9] Tolerance passed to sedumi.
    %     constraints: [1 0 0 0]  defines which transformations to
    %        perform the optimization on.  The first column is linXlinY, the
    %        second is logXlinY, the third is linXlogY, and the fourth column
    %        is logXlogY, where X is the parameters and Y is the model outputs.
    %        1 means include the indicated transformation in the optimization,
    %        0 means exclude them.
    %
    %     fileName2Save: ['']  If nonempty, the main algorithms will
    %        periodically save data about the progress of the algorithm to a
    %        mat-file bearing this name. This way the programmers may look at
    %        the progress without stopping a process.  This is mainly designed
    %        for the programmers as the information that is saved may be
    %        cryptic.
    %     stopGui: ['off']  Turns on a small GUI that allows the user to
    %        stop algorithm without loosing the information.  See
    %        specific help files to see if this is implemented in a
    %        particular function. (other option is 'on'). Currently
    %        implementation is spotty.
    %
    %
    % See also PolyDataset, ResponsePrediction, ConsistencyTest, ParameterOptimization
    %
    %   Constructor reference page in Help browser
    %      <a href="matlab:web(fullfile(fileparts(which('DCsetup')), ...
    %      'docs', 'html', 'functions', 'DatasetUnit.html'), ...
    %      '-helpbrowser')">Broken link</a>
    %
    %   Object reference page in Help browser
    %      <a href="matlab:web(fullfile(fileparts(which('DCsetup')), ...
    %      'docs', 'html', 'obj_dcopt.html'), ...
    %      '-helpbrowser')">DCOptions object</a>

    %TODO the structure and help of this function needs work.
    
    
    properties
        
        %options for display
        display = 'iter';
        guiHandle = [];
        plotFitProgress = 'off';

        %options for fitting
        analysisMode = 'original'; %deprecate this?  It's metamodel stuff,
        %right?
        surrogateType = 'hybridSVM';
        subspaceDiscovery = 'off';
        subspaceThreshold = 0.1;
        derivRange = 0.01;
        surfaceType = 'QuadSurf';
        surfaceFittingMode = 'iterative';
        surfaceTransformation = 'linXlinY';
        fitNorm = 2;
        minFitIter = 3;
        maxFitIter = 7;
        nSuccessfulFitIter = 2;
        maxPnts4Fit = inf;
        useAllPnts4Fit = true;
        fitConvergenceTol = 0.05;
        nPntsPerCoeff4OneShot = 20;
        nPntsPerParam4ActiveParamSel = 25;
        activeParamSelCutOff = 0.05;
        nPntsPerCoeff4Validation = 10;
        nPntsPerCoeff4Metamodel = 250;
        nPntsPerCoeff4MetamodelValidation = 250;
        nLocalValidationSearches = 3;
        nLocalValidationSearches4Metamodel = 3;
        maxDenSwing = 5;
        nComputer = 0;
        findIOTransforms = true;

        %options for optimiztion
        omitInnerBound = false;
        omitOuterBound = false;
        maxBranchBoundIter = 1;
        branchBoundTermTol = 0.02;
        branchBoundQGapTol = 0.01;
        nRestart = 2;
        tolFun = 1e-5;
        tolCon = 1e-5;
        sedumiParEps = 1e-9;
        fullBlockHRelaxation = false;
        constraints = [1 0 0 0];
        paramOptimObjectiveFctnNorm = 'two';

        %miscellaneous options
        fileName2Save = '';
        
    end %properties
    
    
    methods 
        
        %constructor
        function obj = DCOptions(varargin)
            ni = nargin;
            assert(~mod(ni,2),'Incorrect number of inputs');
            if ni>0
                for i1=1:ni/2
                    obj.(varargin{2*i1-1}) = varargin{2*i1};
                end
            end
        end
        
        %special sets
        function obj = set.display(obj,value)
            if ~(strcmp(value,'iter') ||...
                    strcmp(value,'off') ||...
                    strcmp(value,'final') ||...
                    strcmp(value,'notify') ||...
                    strcmp(value,'all') ||...
                    strcmp(value,'ALL'))
                error('display option must be ''iter'', ''off'', ''final'', ''notify'', ''all'', or ''ALL''')
            end
            obj.display = value;
        end
        function obj = set.plotFitProgress(obj,value)
            if ~(strcmp(value,'off') ||...
                    strcmp(value,'on'))
                error('plotFitProgress option must be ''off'' or ''on''')
            end
            obj.plotFitProgress=value;
        end
        function obj = set.subspaceDiscovery(obj,value)
            assert((strcmp(value,'off') || strcmp(value,'on')),'subspaceDisovery option must be ''off'' or ''on''');
            obj.subspaceDiscovery = value;
        end
        function obj = set.subspaceThreshold(obj,value)
            assert((value>0 && value<=1),'subspaceThreshold option takes values in the interval (0,1]')
            obj.subspaceThreshold = value;
        end
        function obj = set.omitInnerBound(obj,value)
            if ~(all(value==0) || all(value==1))
                error('omitInnerBound option must be true or false')
            end
            obj.omitInnerBound = logical(value);
        end
        function obj = set.omitOuterBound(obj,value)
            if ~(all(value==0) || all(value==1))
                error('omitOuterBound option must be true or false')
            end
            obj.omitOuterBound = logical(value);
        end
        function obj = set.fullBlockHRelaxation(obj,value)
            if ~(all(value==0) || all(value==1))
                error('fullBlockHRelaxation option must be true or false')
            end
            obj.fullBlockHRelaxation = logical(value);
        end
        function obj = set.paramOptimObjectiveFctnNorm(obj,value)
            if ~(strcmp(value,'two') ||...
                    strcmp(value,'one') ||...
                    strcmp(value,'inf'))
                error('paramOptimObjectiveFctnNorm option must be ''two'', ''one'', or ''inf''')
            end
            obj.paramOptimObjectiveFctnNorm=value;
        end
        function obj = set.analysisMode(obj,value)
            if ~(strcmp(value,'original') ||...
                    strcmp(value,'metamodelBasedA') ||...
                    strcmp(value,'metamodelBasedB'))
                error('analysisMode option must be ''original'', ''metamodelBasedA'', or ''metamodelBasedB''')
            end
            obj.analysisMode=value;
        end
        function obj = set.surrogateType(obj,value)
            if ~(strcmp(value,'hybridSVM') ||...
                    strcmp(value,'2NormQuad') ||...
                    strcmp(value,'infNormQuad') ||...
                    strcmp(value,'infNormRatQuad') ||...
                    strcmp(value,'kriging') ||...
                    strcmp(value,'eSVM'))
                error('surrogateType option must be ''hybridSVM'', ''2NormQuad'', ''infNormQuad'', ''infNormRatQuad'', ''kriging'', or ''eSVM''')
            end
            obj.surrogateType=value;
        end
        function obj = set.surfaceType(obj,value)
            if ~(strcmp(value,'QuadSurf') ||...
                    strcmp(value,'LinSurf') ||...
                    strcmp(value,'QuadOverQuad'))
                error('surfaceType option must be ''QUadSurf'', ''LinSurf'', or ''QuadOverQuad''')
            end
            obj.surfaceType=value;
        end
        function obj = set.surfaceFittingMode(obj,value)
            if ~(strcmp(value,'iterative') ||...
                    strcmp(value,'one-shot'))
                error('surfaceFittingMode option must be ''iterative'', or ''one-shot''')
            end
            obj.surfaceFittingMode=value;
        end
        function obj = set.fitNorm(obj,value)
            if ~(isequal(value,inf)  || isequal(value,2))
                error('fiNorm option must be inf or 2')
            end
            obj.fitNorm = value;
        end
        function obj = set.useAllPnts4Fit(obj,value)
            if ~(all(value==0) || all(value==1))
                error('useAllPnts4Fit option must be true or false')
            end
            obj.useAllPnts4Fit = logical(value);
        end
        function obj = set.surfaceTransformation(obj,value)
            obj.constraints = zeros(1,4);
            if strmatch('linXlinY',value);
                obj.constraints(1) = 1;
            end
            if strmatch('logXlinY',value);
                obj.constraints(2) = 1;
            end
            if strmatch('linXlogY',value);
                obj.constraints(3) = 1;
            end
            if strmatch('logXlogY',value);
                obj.constraints(4) = 1;
            end
            obj.surfaceTransformation = value;
        end
        function obj = set.findIOTransforms(obj,value)
            if (~all(value==1) && ~all(value==0)) || ~isscalar(value)
                error('Must be scalar logical, 1, or 0')
            end
            obj.findIOTransforms = value;
        end
        
        function [list,sz] = displayProps(obj)
            props=properties(obj);
            width = max(cellfun(@length,props));
            nprops = length(props);
            list = cell(nprops,1);
            for i1=1:nprops
                list{i1} = sprintf(['%' width 's: %s'],props{i1},obj.(props{i1}));
            end
            list = [];
            sz = '';
        end

    end %methods
    


end %classdef