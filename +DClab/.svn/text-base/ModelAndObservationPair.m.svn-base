classdef ModelAndObservationPair < DClab.DCObject
    % MOPair = ModelAndObservationPair(RespObs,RespModel,name,guiDispCallback)
    %
    % MODELANDOBSERVATIONPAIR: constructor for a ModelAndObservationPair
    % object. This object groups a ResponseObservation and a ResponseModel
    % (both describing the same response or physical quantity) to form a constraint
    % on the model's parameter vector.
    %
    % Syntax:
    %   MOPair = ModelAndObservationPair(RespObs,RespModel,name)
    %   MOPair = ModelAndObservationPair(RespObs,RespModel,name,critRange)
    %   MOPair = ModelAndObservationPair(RespObs,RespModel,name,critRng,trnRng)
    %   MOPair = ModelAndObservationPair(RespObs,RespModel,name,critRng,trnRng,guiDispCallback)
    %   MOPair = ModelAndObservationPair(integer) initializes a integerx1 empty object.
    %
    % Inputs:
    %   RespObs: a ResponseObservation object
    %   RespModel: a ResponseModel object
    %   name: a single row character array intended to be used as a
    %     unique identifier. All components of a multidimensional
    %     ModelAndObservationPair object (formed by vertcat) must have unique names.
    %   critRng[opt]: critical range for surrogate error assessment; a 1-by-2
    %     vector. The fitting error of a metamodel or quadratic surface will be
    %     assessed over the preimage under RespModel or this set.
    %   trnRng[opt]: training range for surrogate creation; a 1-by-2 vector.
    %     The metamodel and quadratic surfaces will be constructed from
    %     training data whose y-components line in [trnRng(1) trnRng(2)];
    %   guiDispCallback[opt]: function handle to an mfile that will be
    %     called with the syntax guiDispCallback(MOPair) that displays further
    %     info about the object MOPair. Typically such an mfile would make use
    %     of the .userdata property of the object in order to identify it
    %     and launch an appropriate graphical display (e.g., webpage).
    %
    % See also ResponseModel, ResponseObservation
    %
    %   Constructor reference page in Help browser
    %      <a href="matlab:web(fullfile(fileparts(which('DCsetup')), ...
    %      'docs', 'html', 'functions', 'ModelAndObservationPair.html'), ...
    %      '-helpbrowser')">ModelAndObservationPair constructor</a>
    %
    %   Object reference page in Help browser
    %      <a href="matlab:web(fullfile(fileparts(which('DCsetup')), ...
    %      'docs', 'html', 'obj_mopair.html'), ...
    %      '-helpbrowser')">ModelAndObservationPair object</a>

    properties
        ResponseObservation = DClab.ResponseObservation;
        ResponseModel = DClab.ResponseModel;
        name = '';
        criticalRange;
        trainingRange;
        guiDisplayCallback;
        status = -1;
    end

    properties (Dependent)
        observedValue;
        observationUncertaintyPlusMinus;
        observationUncertainty;
        uncertaintyCase
        constraintLowerLimit
        constraintUpperLimit
    end

    methods


        function Obj = ModelAndObservationPair(RespObs,RespModel,name,critRng,trnRng,guiDispCallback)

            %TODO, based on the output uncertainty transformations and sign of L,
            %we will disallow certain combinations.

            ni = nargin;
            no = nargout;
            error( nargoutchk(0,1,no) );

            switch ni
                case 0
                    return
                case 1   %initiallize an mx1 empty object

                    % Special case 2:
                    input1 = RespObs;
                    assert(isnumeric(input1) && isscalar(input1) && ceil(input1) == input1 && input1 > 0,'Usage: When called with a single input, ModelAndObservationPair must recieve a positive integer');
                    [Obj(1:input1,1).ResponseObservation] = deal(DClab.ResponseObservation);
                    [Obj.ResponseModel] = deal(DClab.ResponseModel);
                    [Obj.name] = deal('');
                    [Obj.criticalRange] = deal([]);
                    [Obj.trainingRange] = deal([]);
                    [Obj.guiDisplayCallback] = deal([]);
                    [Obj.status] = deal(-1);
                    return

                case 2
                    error('Usage: Incorrect number of input arguments')

                case 3
                    critRng = [];
                    trnRng = [];
                    guiDispCallback = [];
                case 4
                    trnRng = [];
                    guiDispCallback = [];
                case 5
                    guiDispCallback = [];
                otherwise
                    error( nargchk(0,6,ni) );
            end

            % Initialize empty inputs:
            if isempty(critRng)
                critRng = [-inf inf];
            end
            if isempty(trnRng)
                trnRng = [-inf inf];
            end

            % Create object
            Obj.ResponseObservation = RespObs;
            Obj.ResponseModel = RespModel;
            Obj.name = name;
            Obj.criticalRange = critRng;
            Obj.trainingRange = trnRng;
            Obj.guiDisplayCallback = guiDispCallback;
            if ~isempty(RespObs) && ~isempty(RespModel) && ~isempty(name)
                Obj.status = 1;
            elseif isempty(RespObs) && isempty(RespModel) && isempty(name)
                Obj.status = -1;
            else
                Obj.status = 0;
            end

        end

        function plist = neededParameters(obj)
            % function plist = neededParameters(obj)
            %
            % Method to obtain a sorted list of the names of the parameters required
            % by the collection of scalar pairs in a ModelAndObservationPair object.
            % The returned list plist  is sorted both out of convenience, and because
            % nothing else makes sense.
            %
            % See also ModelAndObservationPair, ModelAndObservationPair/getPrivate

            m = length(obj);
            plist = cell(m,1);
            for i1 = 1:m
                plist{i1} = obj(i1).ResponseModel.parameterList;
            end
            plist = unique(vertcat(plist{:}));
        end

        function [bool emptyIdx] = isempty(obj)
            %ISEMPTY Isempty for FreeParameter class
            %
            %   BOOL = ISEMPTY(OBJ) returns true if all componenents of a multidim OBJ
            %   are empty. A component is considered empty if it is missing a name,
            %   ResponseObservation, and ResponseModel.
            %
            %   [BOOL EMPTYIDX] = ISEMPTY(OBJ) additionally returns a vector of
            %   indicies of the empty components.
            %
            %
            %   See also ModelAndObservationPair, ModelAndObservationPair/status

            error(nargoutchk(0,2,nargout))
            error(nargchk(1,1,nargin))

            emptyIdx = vertcat(obj.status) == -1;
            if all(emptyIdx)
                bool = true;
            else
                bool = false;
            end

            if nargout == 2
                emptyIdx = find(vertcat(obj.complete) == -1);
            end
        end
        function [PairArray erased] = vertcat(varargin)
            %VERTCAT Vertical concatenation of ModelAndObservationPair objects.
            %
            %   PAIRARRAY = VERTCAT(OBJ1,OBJ2,...) or PAIRARRAY = [OBJ1; OBJ2; ...]
            %   contatenates the input objects to form the multidimension output
            %   PAIRARRAY. The MOPairs in PAIRARRAY will never share a common name. See
            %   "Behavior" below.
            %
            %   [PAIRARRAY ERASED] = VERTCAT(OBJ1,...,OBJn) additionally returns a vector
            %   of indices indicating which components of
            %   tmp = builtin('vertcat',varargin{:}), which will not strip out
            %   MOPairs with the same name,  were included in PAIRARRAY. This output is
            %   used to properly eliminate the corresponding DCSurfaces when vertcat is
            %   invoked by PolyDataset/merge.
            %
            %   Behavior:
            %
            %   If an "Incomplete" MOPair is encountered, an error with be thrown, as
            %   this condition can create difficulties/ambiguities with the
            %   requirement that the output of this method be an object with no shared
            %   names. "Empty" ModelAndObservationPair are not included in the final
            %   output.
            %
            %   Each input is potentially a multidimensional vector object
            %   consisting of several MOPairs. Consider all individual MOPairs present
            %   in the input list. If multiple MOPairs share an identical name,
            %   the code will check if they are identical in all respects. If so, only
            %   one copy will be present in the final output. If nonidentical MOPairs
            %   with a common name are encountered, a warning will be issued. Again
            %   only MOPair with the given name will be present in the final output,
            %   and it will correspond to the earliest MOPair appearing in the input
            %   list with that same name.
            %
            %   See also ModelAndObservationPair

            no = nargout;
            error(nargoutchk(0,2,no));

            %Remove any native empties from the varargin list
            for i1 = length(varargin):-1:1
                if isa(varargin{i1},'DClab.ModelAndObservationPair')
                    %do nothing, this is the typical case
                elseif builtin('isempty',varargin{i1})
                    varargin(i1) = [];
                else
                    error('Usage: concatentation of nonidentical object types')
                end
            end

            %At this point varargin is a cell array with each cell containing a
            %ModelAndObservationPair object. Combine these with MATLAB's struct vertcat
            %method, and then check for duplicate names.

            PairArray = builtin('vertcat',varargin{:});

            %if there are any incompletes, just bomb. We need to remove empties in
            %all cases. Then we need to check for duplicates. If they occur we
            %have some work to do...

            if sum(vertcat(PairArray.status)==0) > 0
                error('Inputs: an incomplete object was supplied to MODELANDOBSERVATIONPAIR/VERTCAT')
            end

            %Remove any empty MOPairs from the final output
            PairArray(vertcat(PairArray.status)==-1) = [];

            if isempty(PairArray)
                PairArray = DClab.ModelAndObservationPair;
            end

            %Check for duplicates (as indicated by shared names)
            ulist = {PairArray.name}';
            uulist = unique(ulist);
            erased = [];

            if length(uulist)~= length(ulist)
                %We have duplicates and we need to eliminate them.

                shouldWarn = false;
                for i1 = 1:length(uulist)
                    idx = strmatch(uulist(i1),ulist,'exact');
                    if length(idx) > 1
                        for i2 = 2:length(idx)
                            if ~isequal(PairArray(idx(1)),PairArray(idx(i2)))
                                shouldWarn = true;
                            end
                        end
                        erased = [erased; idx(2:end)];
                    end
                end
                PairArray(erased) = []; %eliminate the duplicates

                if shouldWarn
                    warning off backtrace

                    %Call warning with an empty as the 2nd input to force it to display a
                    %formatted string.
                    warning(['Concatenating dataset units with identical names \n '...
                        'but other differences, earliest input object will be used'],[]);

                    warning on backtrace
                end
            end


        end
        function [list,sz] = displayProps(obj)
            list = [];
            sz = sprintf('%d-by-%d',size(obj));
        end

        function out=nPairs(obj)
            out = sum([obj.status] == 1);
        end
        function out=get.observedValue(obj)
            out = obj.ResponseObservation.observedValue;
        end
        function out=get.observationUncertaintyPlusMinus(obj)
            out = obj.ResponseObservation.uncertaintyPlusMinus;
        end
        function out=get.observationUncertainty(obj)
            out = obj.ResponseObservation.uncertainty;
        end
        function out=get.uncertaintyCase(obj)
            tmp = reshape(1:16,4,4);

            %TODO if RMcase or ROcase is empty this will give 1x0 empties. 0x0
            %empty would be nicer.
            RMcase = obj.ResponseModel.uncertaintyCase;
            ROcase = obj.ResponseObservation.uncertaintyCase;
            out = tmp(ROcase,RMcase);
        end
        function out=get.constraintLowerLimit(obj)
            if obj.status ~= 1
                out = [];
                return
            end
            uncCase = obj.uncertaintyCase;
            ou = obj.ResponseModel.outputUncertaintyPlusMinus;
            u = obj.ResponseObservation.uncertaintyPlusMinus;
            d = obj.ResponseObservation.observedValue;
            switch uncCase
                case 1
                    out = d + u(1) - ou(2);
                case 2
                    out = d*(1 + u(1)) - ou(2);
                case 3
                    out = d*10^u(1) - ou(2);
                case 4
                    out = d^(1+u(1)) - ou(2);
                case 5
                    if ou(1) >= -1
                        out = (d+u(1))/(1+ou(2));
                    else
                        out = min((d+u(1))/(1+ou(2)),(d+u(2))/(1+ou(1)));
                    end
                case 6
                    if ou(1) >= -1
                        out = d*(1+u(1))/(1+ou(2));
                    else
                        out = min(d*(1+u(1))/(1+ou(2)),d*(1+u(2))/(1+ou(1)));
                    end
                case 7
                    if ou(1) >= -1
                        out = (d*10^u(1))/(1+ou(2));
                    else
                        out = min((d*10^u(1))/(1+ou(2)),(d*10^u(2))/(1+ou(1)));
                    end
                case 8
                    if ou(1) >= -1
                        out = (d*10^u(1))/(1+ou(2));
                    else
                        out = min((d*10^u(1))/(1+ou(2)),(d*10^u(2))/(1+ou(1)));
                    end
                case 9
                    out = (d+u(1))*10^-ou(2);
                case 10
                    out = d*(1+u(1))*10^-ou(2);
                case 11
                    out = d*10^(u(1)-ou(2));
                case 12
                    out = d^(1+u(1))*10^-ou(2);
                case 13
                    if ou(1) >= -1
                        out = (d+u(1))^(1/(1+ou(2)));
                    else
                        out = min((d+u(1))^(1/(1+ou(2))),(d+u(2))^(1/(1+ou(1))));
                    end
                case 14
                    if ou(1) >= -1
                        out = (d*(1+u(1)))^(1/(1+ou(2)));
                    else
                        out = min((d(1+u(1)))^(1/(1+ou(2))),(d*(1+u(2)))^(1/(1+ou(1))));
                    end
                case 15
                    if ou(1) >= -1
                        out = 10^((log10(d)+u(1))/(1+ou(2)));
                    else
                        out = min(10^((log10(d)+u(1))/(1+ou(2))),10^((log10(d)+u(2))/(1+ou(1))));
                    end
                case 16
                    if ou(1) >= -1
                        out = d^((1+u(1))/(1+ou(2)));
                    else
                        out = min(d^((1+u(1))/(1+ou(2))),d^((1+u(2))/(1+ou(1))));
                    end
                otherwise
                    error('Unexpected uncertaintyCase, condition should never occur')
            end
        end
        function out=get.constraintUpperLimit(obj)
            if obj.status ~= 1
                out = [];
                return
            end
            uncCase = obj.uncertaintyCase;
            ou = obj.ResponseModel.outputUncertaintyPlusMinus;
            u = obj.ResponseObservation.uncertaintyPlusMinus;
            d = obj.ResponseObservation.observedValue;
            switch uncCase
                case 1
                    out = d + u(2) - ou(1);
                case 2
                    out = d*(1 + u(2)) - ou(1);
                case 3
                    out = d*10^u(2) - ou(1);
                case 4
                    out = d^(1+u(2)) - ou(1);
                case 5
                    if ou(1) > -1
                        out = (d+u(2))/(1+ou(1));
                    else
                        out = Inf;
                    end
                case 6
                    if ou(1) > -1
                        out = d*(1+u(2))/(1+ou(1));
                    else
                        out = Inf;
                    end
                case 7
                    if ou(1) > -1
                        out = (d*10^u(2))/(1+ou(1));
                    else
                        out = Inf;
                    end
                case 8
                    if ou(1) > -1
                        out = (d*10^u(2))/(1+ou(1));
                    else
                        out = Inf;
                    end
                case 9
                    out = (d+u(2))*10^-ou(1);
                case 10
                    out = d*(1+u(2))*10^-ou(1);
                case 11
                    out = d*10^(u(2)-ou(1));
                case 12
                    out = d^(1+u(2))*10^-ou(1);
                case 13
                    if ou(1) > -1
                        out = (d+u(2))^(1/(1+ou(1)));
                    else
                        out = Inf;
                    end
                case 14
                    if ou(1) > -1
                        out = (d*(1+u(2)))^(1/(1+ou(1)));
                    else
                        out = Inf;
                    end
                case 15
                    if ou(1) > -1
                        out = 10^((log10(d)+u(2))/(1+ou(1)));
                    else
                        out = Inf;
                    end
                case 16
                    if ou(1) > -1
                        out = d^((1+u(2))/(1+ou(1)));
                    else
                        out = Inf;
                    end
                otherwise
                    error('Unexpected uncertaintyCase, condition should never occur')
            end
        end

    end %public methods

end % classdef
