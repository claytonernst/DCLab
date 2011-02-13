classdef ResponseModel < DClab.DCObject
    %RESPONSEMODEL: constructor for object that describes a model for the value of a physical quantity
    %
    %   RM = ResponseModel(@DCMODEL)
    %   RM = ResponseModel(@DCMODEL,ADDLINPUT1,...,ADDLINPUTn)
    %   RM = ResponseModel(COEFFMATRIX,DOMAIN)
    %   RM = ResponseModel(MODELSTRUCT,DOMAIN)
    %   RM = ResponseModel(COEFFMATRIX,DOMAIN,OUPUTUNC)
    %   RM = ResponseModel(COEFFMATRIX,DOMAIN,OUPUTUNC,GUIDISPCALLBACK)
    %
    %   Creation when using a user supplied dcModel function mfile. The
    %   calling syntax this mfile must support is described in the html
    %   file linked below. See /DClabV1p1/fileTemplates/dcModel.m for an
    %   example.:
    %
    %      RM = ResponseModel(@DCMODEL) creates an object with all properties
    %      taken from the mfile DCMODEL.m upon construction. During object
    %      construction DCMODEL.m is set to read-only and afterwards it is only
    %      used to evaluate the model.
    %
    %      RM = ResponseModel(@DCMODEL,ADDLINPUT1,...,ADDLINPUTn) creates an
    %      object where the syntax to evaluate the model is at x is
    %      y = DCMODEL('simulate',x,ADDLINPUT1,..,ADDLINPUTn). These additional
    %      inputs are automatically passed to DCMODEL when it is invoked by
    %      ResponseModel/eval.
    %
    %   Creation when using a linear or quadratic model
    %
    %      M = ResponseModel(COEFFMATRIX,DOMAIN) creates an object whose
    %      algebraic model equation is drawn from COEFFMATRIX. The structure
    %      array DOMAIN specifies the domain of the model.
    %
    %      M = ResponseModel(MODELSTRUCT,DOMAIN) is similar to the previous
    %      syntax, except here the algebraic model equation expresses a
    %      relationship between transformed versions of y and x. See the input
    %      description regarding MODELSTRUCT below.
    %
    %      M = ResponseModel(COEFFMATRIX,DOMAIN,OUTPUTUNC)
    %      M = ResponseModel(MODELSTRUCT,DOMAIN,OUTPUTUNC) accepts an
    %      additional input that specifies any uncertainty in the output of the
    %      algebraic model. OUTPUTUNC should be provided if the model output is
    %      potentially biased. See the input description regarding OUTPUTUNC
    %      below.
    %
    %      M = ResponseModel(COEFFMATRIX,DOMAIN,OUTPUTUNC,@GUIDISPCALLBACK)
    %      creates an object possessing a callback function that is invoked by
    %      GUIs when the user requested additional information about the object.
    %
    %   ===Inputs when using a dcModel function===:
    %
    %   @DCMODEL is a handle to a function mfile on the path that must conform
    %      to specific syntax. More information is provided in the html
    %      documentation, or by examining the template file
    %      /DClabV1p1/fileTemplates/dcModel.m
    %   ADDLINPUTi: these are passed as additional inputs to the dcModel
    %      function. They may be used, for example, specify the input
    %      magnitude, initial conditions, or to differentiate between multiple
    %      responses that are derived from the same simulation. If model
    %      evaluations are to be stored (see the saveEnabled property below)
    %      these inputs should be (or be convertible to) single line strings.
    %
    %   ===Inputs when using a linear or quadratic model===
    %
    %   COEFFMATRIX: When specifying a linear model, this should be a
    %      (n+1)-by-1 vector, where n is the number of model parameters, such
    %      that the model equation is y = COEFFMATRIX'*[1; x1; ... xn]. When
    %      specifying a quadratic model, this should be a (n+1)-by-(n+1)
    %      symmetric matrix such that the model equation is
    %      y = [1 x1 ... xn]*COEFFMATRIX*[1; x1; ... xn].
    %   MODELSTRUCT: a structure with fields .value, .responseTransformation,
    %      and .variableTransformations. .value should contain a (n+1)-by-(n+1)
    %      symmetric matrix or (n+1)-by-1 vector). .responseTransformation must
    %      be 'none' or 'log10'. .variableTransformations must an nx1 cell
    %      array where each entry is 'none' or 'log10'. For the linear model
    %      case, the model equation is
    %      respTrans(y) = value'*[1; varTrans1(x1); ... varTransn(xn)]. The
    %      equation for a quadratic model is analogous.
    %   DOMAIN: an nx1 struct array with fields .name and .range specifying the
    %      domain of the model. The first component of the structure should
    %      correspond to x1 in the above expressions, etc. The range is
    %      typically a bounded interval (e.g., [-1, 1]) but can be infinite if
    %      the algebraic model expression is thought to be valid for all values
    %      of a particular variable.
    %   OUTPUTUNC[optional, default=0]: This input is used to supply bounds on
    %      the uncertainty (bias) in the model output. For example, this would
    %      be the peak fitting error when the algebraic model is a low-order
    %      approximation to a more detailed model. outputUncertainty should be
    %      the scalar if the uncertainty in the model output is symmetric. If
    %      it is asymmetric this should be a 1x2 vector, with
    %      OUTPUTUNC(1) <= eta(x) - M(x) <= OUTPUTUNC(2), where M(x) is the
    %      algebraic model and eta(x) is the true model.
    %   OUTPUTUNC[non absolute interval-type uncertainty] If the output
    %      uncertainty is expressed as a percent uncertainty or in transformed
    %      coordinates, OUTPUTUNC should be a structure with fields .value,
    %      .type, and .transformation. .value should be a scalar or a 1x2,
    %      .type can be 'absolute' or 'relative', while transformation can be
    %      'none' or 'log10'. If either .type or .transformation are omitted,
    %      their respective defaults of 'absolute' and 'none' will be assumed.
    %      For the relative case, value should be a relative fraction, i.e. 0.1
    %      instead of 10%. Let
    %         case 1 denote .type = 'absolute' and .transformation = 'none'.
    %         case 2 denote .type = 'relative' and .transformation = 'none'.
    %         case 3 denote .type = 'absolute' and .transformation = 'log10'.
    %         case 4 denote .type = 'relative' and .transformation = 'log10'
    %      The uncertainty expressions for the four cases are
    %         symmetric uncertainty (.value = a scalar)
    %            case 1: |eta(x) - M(x)| <= .value
    %            case 2: |eta(x) - M(x)|/M(x) <= .value
    %            case 3:  |log10(eta(x)/M(x))| <= .value
    %            case 4: |log10(eta(x)/M(x))|/log10(M(x)) <= .value
    %         asymmetric uncertainty (.value = 1x2, 1st component negative)
    %            case 1: .value(1) <= eta(x) - M(x) <= .value(2)
    %            case 2: .value(1) <= (eta(x)-M(x))/M(x) <= .value(2)
    %            case 3: .value(1) <= log10(eta(x)/M(x)) <= .value(2)
    %            case 4: .value(1) <= log10(eta(x)/M(x))/log10(M(x)) <= .value(2)
    %      Limitations:
    %         case 1: none.
    %         case 2: M(x) must be postive for all x in the domain.
    %         case 3: both eta(x) and M(x) must be positive for all x in the domain.
    %         case 4: eta(x) must be positive and M(x) must be greater than 1
    %            for all x in the domain.
    %   @GUIDISPCALLBACK[optional]: a function handle to an mfile that will be
    %      called with the syntax GUIDISPCALLBACK(RM) that displays further
    %      info about the object RM. Typically such an mfile would make use of
    %      the .name or .userdata properties of the object in order to identify
    %      it and launch an appropriate graphical display (e.g., webpage).
    %
    %   The created object has the following properties (fields):
    %
    %      type: either 'dcModel', 'quadraticModel', or 'linearModel'
    %      name: a single line string providing the name of the response. The
    %         constructor generates this automatically, but it may be changed
    %         using SET.
    %      model: either the handle to a dcModel mfile or a matrix
    %         representation of an algebraic model
    %      domain: a structure array with fields .name and .range
    %      outputUncertainty: a scalar or 1x2 indicating bounds on the uncertainty
    %         in the model output
    %      outputUncertaintyType: either 'absolute' or 'relative'
    %      outputUncertaintyTransformation: either 'none' or 'log10'
    %      uncertaintyCase: either 1,2,3, or 4, dictated by the cases described
    %         above.
    %      guiDisplayCallback: a function_handle to a display function
    %      responseTransformation: either 'linear' or 'log10'
    %      variableTransformations: nx1 cell array with each cell either
    %         'linear', or 'log10'
    %      additionalInputs: additional inputs that are passed to a dcModel file
    %         when the 'simulate' flag is invoked.
    %      saveEnabled: a boolean indicating whether evaluations are to be
    %         stored. This functionality is only supported when the
    %         ResponseModel object is constructed from a DCModel mfile and
    %         any additional inputs are convertible to single line strings.
    %      multipleResponsesEnabled: a boolean indicating whether multiple
    %         responses (usually features of a simulated time-trace) will be
    %         computed from a single dcModel evaluation.
    %     isPredictor: used internally by the software to determine if a
    %         ResponsePrediction is being performed on this ResponseModel
    %
    %   See also ResponseObservation, ModelAndObservationPair
    %
    %   Constructor reference page in Help browser
    %      <a href="matlab:web(fullfile(fileparts(which('DCsetup')), ...
    %      'docs', 'html', 'functions', 'ResponseModel.html'), ...
    %      '-helpbrowser')">ResponseModel constructor</a>
    %
    %   Object reference page in Help browser
    %      <a href="matlab:web(fullfile(fileparts(which('DCsetup')), ...
    %      'docs', 'html', 'obj_modresp.html'), ...
    %      '-helpbrowser')">ResponseModel object</a>
    %
    %   dcModel mfile description in Help browser
    %      <a href="matlab:web(fullfile(fileparts(which('DCsetup')), ...
    %      'docs', 'html', 'misc', 'dcModel.html'), ...
    %      '-helpbrowser')">ResponseModel object</a>

    %TODO outputUncertainty in mfile should be able to be a vector or matrix when
    %responseList is present

    %TODO should we drop the 'output' that precedes uncertainty? That would
    %make this consistent with the ResponseObservation and FreeParameter
    %objects.

    properties
        type;
        name = '';
        model;
        domain = struct([]);
        guiDisplayCallback;
        outputUncertainty;
        outputUncertaintyType = 'absolute';
        outputUncertaintyTransformation = 'none';
        responseTransformation = '';
        variableTransformations = {};
        additionalInputs = {};
        saveEnabled;
        multipleResponsesEnabled;
        responseSimulationSegment;
        isPredictor;
    end
    
    properties (Dependent)
        parameterList;
        outputUncertaintyPlusMinus;
        responseList;
        getResponseSimulationSegment;
    end

    properties (Dependent, SetAccess=private)
        uncertaintyCase;
    end

    methods
        function Obj = ResponseModel(varargin)


            ni = nargin;
            no = nargout;
            error( nargoutchk(0,1,no) );

            if ni > 0
                switch(class(varargin{1}))
                    case 'function_handle'
                        Obj = fromDCModel(Obj,varargin{:});
                    case {'double','struct'}
                        assert(ni >= 2,'Usage: M = RESPONSEMODEL(ALGFORM,MODELDOMAIN,...)');
                        Obj = fromAlg(Obj,varargin{:});
                    otherwise
                        error('Usage: First input to RESPONSEMODEL should be a MODEL_HANDLE, DOUBLE, or STRUCT');
                end



                Obj.isPredictor = false;

            end
        end
        
        function bool=isempty(obj)
            bool=true;
            if ~isempty(obj.type)
                bool=false;
            end
        end

        function out=get.parameterList(obj)
            if ~isempty(obj.domain)
                out = {obj.domain.name}';
            else
                out = {};
            end
        end
        function out=get.outputUncertaintyPlusMinus(obj)
            if isscalar(obj.outputUncertainty)
                out = [-obj.outputUncertainty obj.outputUncertainty];
            else
                out = obj.outputUncertainty;
            end
        end
        function out=get.responseList(obj)
            try
                out = obj.model('getResponseList');
            catch ME %#ok
                out = {};
            end
        end
        function out=get.getResponseSimulationSegment(obj)
            try
                out = obj.model('getResponseSimulationSegment');
            catch ME %#ok
                out = [];
            end
        end
        function value=get.uncertaintyCase(obj)
            if strcmp(obj.outputUncertaintyType,'absolute') && strcmp(obj.outputUncertaintyTransformation,'none')
                value = 1;
            elseif strcmp(obj.outputUncertaintyType,'relative') && strcmp(obj.outputUncertaintyTransformation,'none')
                value = 2;
            elseif strcmp(obj.outputUncertaintyType,'absolute') && strcmp(obj.outputUncertaintyTransformation,'log10')
                value = 3;
            elseif strcmp(obj.outputUncertaintyType,'relative') && strcmp(obj.outputUncertaintyTransformation,'log10')
                value = 4;
            else
                value = [];
            end
        end
        function [list,sz] = displayProps(obj)
            if ~isempty(obj.name)
                list{1} = sprintf('Model: %s',obj.name);
            else
                list = [];
            end
            sz = '';
        end
        
        function obj=set.outputUncertainty(obj,value)
            obj.outputUncertainty = value;
        end
        function obj=set.outputUncertaintyType(obj,value)
            obj.outputUncertaintyType = value;
        end
        function obj=set.outputUncertaintyTransformation(obj,value)
            obj.outputUncertaintyTransformation = value;
        end
    end %public methods

    methods (Access=private)
        % --------------------- PRIVATE METHODS ------------------------- %
        function Obj = fromDCModel(Obj,dcModel,varargin)

            % Make dcModel.m read-only
            try
                fileattrib([func2str(dcModel) '.m'],'-w');
            catch ME %#ok
                %that's ok, move on
            end

            % Create the object.
            Obj.type = 'dcModel';
            Obj.name = '';
            Obj.model = dcModel;

            % model domain is the only this this file is required to provide
            try
                Obj.domain = feval(dcModel,'getModelDomain',[],varargin{:});
            catch ME
                handleInfo = functions(dcModel); %TODO the 'functions' function is subject to change in later matlab releases
                if isempty(handleInfo.file)
                    error('Inputs: DCMODEL is not a handle to an m-file on the MATLAB path')
                else
                    error('Invalid DCMODEL m-file: feval(DCMODEL,''getModelDomain'') failed.')
                end
            end

            try
                Obj.outputUncertainty = feval(dcModel,'getOutputUncertainty',[],varargin{:});
            catch %#ok
                Obj.outputUncertainty = 0;
            end

            try
                Obj.outputUncertaintyType = feval(dcModel,'getOutputUncertaintyType',[],varargin{:});
            catch %#ok
                Obj.outputUncertaintyType = 'absolute';
            end

            try
                Obj.outputUncertaintyTransformation = feval(dcModel,'getOutputUncertaintyTransformation',[],varargin{:});
            catch %#ok
                Obj.outputUncertaintyTransformation = 'none';
            end


            try
                Obj.guiDisplayCallback = feval(dcModel,'getGuiDisplayCallback',[],varargin{:});
            catch %#ok
                Obj.guiDisplayCallback = [];
            end

            %If you define the model in an m-file, there is no reason for you to need
            %transformations
            Obj.responseTransformation = 'none';
            Obj.variableTransformations = repmat({'none'},size(Obj.domain));

            Obj.additionalInputs = varargin;

            try
                Obj.saveEnabled = feval(dcModel,'isSaveEnabled',[],varargin{:});
            catch %#ok
                Obj.saveEnabled = false;
            end

            try
                Obj.multipleResponsesEnabled = feval(dcModel,'isMultipleResponsesEnabled',[],varargin{:});
            catch %#ok
                Obj.multipleResponsesEnabled = false;
            end
            
            try
                Obj.responseSimulationSegment = feval(dcModel,'getResponseSimulationSegment',[],varargin{:});
            catch %#ok
                Obj.responseSimulationSegment = [];
            end
            
            try
                Obj.name = feval(dcModel,'getName',[],varargin{:});
            catch %#ok
                Obj.name = func2str(dcModel);
            end

        end
        function Obj = fromAlg(Obj,coeffMatrix,domain,varargin)
            %
            % CoeffMatrix will be a structure or a double.

            if isstruct(coeffMatrix)
                if isfield(coeffMatrix,'value') && isfield(coeffMatrix,'responseTransformation') && isfield(coeffMatrix,'variableTransformations')
                    respTransformation = coeffMatrix.responseTransformation;
                    varTransformations = coeffMatrix.variableTransformations;
                    coeffMatrix = coeffMatrix.value;
                else
                    error('Inputs: structure COEFFMATRIX contains insufficient fields')
                end
            else
                respTransformation = 'none'; %default is none
                varTransformations = repmat({'none'},length(coeffMatrix)-1,1); %default is all none
            end

            if size(coeffMatrix,1) == size(coeffMatrix,2) && numel(coeffMatrix) > 1
                Obj.type = 'quadraticModel';
                Obj.name = 'unnamedQuad';
                Obj.model = 0.5*(coeffMatrix+coeffMatrix');
            else
                Obj.type = 'linearModel';
                Obj.name = 'unnamedLin';
                Obj.model = coeffMatrix;
            end
            Obj.domain = domain;

            n = length(varargin);

            if n == 0
                Obj.outputUncertainty = 0;
                Obj.outputUncertaintyType = 'absolute';
                Obj.outputUncertaintyTransformation = 'none';
                Obj.guiDisplayCallback = [];
                Obj.responseTransformation = respTransformation;
                Obj.variableTransformations = varTransformations;
            else
                if isstruct(varargin{1})
                    if ~isempty(setdiff(fieldnames(varargin{1}),{'value','type','transformation'}))
                        warning('DCLAB:unusedInputs','structure OUTPUTUNC contains unused fields')
                    end
                    if isfield(varargin{1},'value')
                        Obj.outputUncertainty = varargin{1}.value;
                    else
                        error('Inputs: structure OUTPUTUNC contains insufficient fields')
                    end
                    if isfield(varargin{1},'type')
                        Obj.outputUncertaintyType = varargin{1}.type;
                    else
                        Obj.outputUncertaintyType = 'absolute';
                    end
                    if isfield(varargin{1},'transformation')
                        Obj.outputUncertaintyTransformation = varargin{1}.transformation;
                    else
                        Obj.outputUncertaintyTransformation = 'none';
                    end


                elseif isnumeric(varargin{1})
                    Obj.outputUncertainty = varargin{1};
                    Obj.outputUncertaintyType = 'absolute';
                    Obj.outputUncertaintyTransformation = 'none';
                else
                    error('Inputs: OUTPUTUNC of improper class')
                end

                Obj.guiDisplayCallback = [];
                Obj.responseTransformation = respTransformation;
                Obj.variableTransformations = varTransformations;

                if n == 2
                    Obj.guiDisplayCallback = varargin{2};
                end
                if n > 2
                    Obj.guiDisplayCallback = varargin{2};
                    warning('DCLAB:unusedInputs','RESPONSEMODEL received unused inputs')
                end
            end

            Obj.additionalInputs = {};
            Obj.saveEnabled = 0;
            Obj.multipleResponsesEnabled = 0;


        end

    end %private methods

end %classdef