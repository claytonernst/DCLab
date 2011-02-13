classdef FreeParameter < DClab.DCObject
    %FREEPARAMETER Construct FreeParameter object
    %
    %   FP = FreeParameter creates an empty object
    %
    %   FP = FreeParameter(NAME,RANGE)
    %
    %   FP = FreeParameter(NAME,NOMINAL,UNCERTAINTY) creates an object FP
    %   stating that the parameter NAME has a nominal or best estimate value
    %   of NOMINAL but can potentially takes any value in the interval
    %   [NOMINAL-UNCERTAINTY, NOMINAL+UNCERTAINTY]. Syntax for more expressive
    %   forms of parametric uncertainty is shown below.
    %
    %   FP = FreeParameter(NAME,NOMINAL,UNCERTAINTY,GUIDISPLAYCALLBACK)
    %   additionally adds the function_handle GUIDISPLAYCALLBACK to the
    %   object. This function will be invoked as GUIDISPLAYCALLBACK(FP) when
    %   a gui user requests more information about the object.
    %
    %   FP = FreeParameter(INTEGER) initializes an INTEGERx1 empty object.
    %
    %   Examples:
    %   >> FP = FreeParameter('p1',[8 12])
    %      creates an object stating the value of parameter 'p1' lies
    %      between 8 and 12.
    %   >> FP = FreeParameter('p1',10,2)
    %      creates an object stating the value of parameter 'p1' lies
    %      between 10 +/- 2.  This is the same as the previous example.
    %   >> unc.value = 0.05, unc.type = 'relative', unc.transformation = 'none';
    %   >> FP = FreeParameter('p1',10,unc)
    %      creates an object stating the value of parameter 'p1' is within 5% of 10
    %
    %   Supported forms of uncertainty:
    %
    %   The uncertainty can be symmetric or asymmetric about the nominal value,
    %   It can be relative or absolute, and it can be expressed in transformed
    %   coordinates (currently only log10 is supported). The most simple case
    %   it absolute uncertainty with no transformation. Here the input
    %   UNCERTAINTY can be a scalar if the uncertainty is symmetric.
    %   Alternatively it should be a 1x2 double array containing bounds on the
    %   uncertainty s.t.
    %            nom + unc(1) <= true value <= nom + unc(2).
    %   Note unc(1) should be negative.
    %
    %   To use other forms, the input UNCERTAINTY must be a structure with
    %   fields 'value', 'type', and 'transformation'. 'value' should be a
    %   scalar or 1x2 vector as above. 'type' must be 'absolute', or
    %   'relative'. If relative, value should be a relative fraction, i.e, 0.1
    %   rather than 10(%). 'transformation' must be 'none' or 'log10'. Let p*
    %   denote the true value of the parameter under consideration. Note a
    %   log10 transformation only makes sense if nom and p* have the same sign.
    %   These different types allow the uncertainty to be expressed as
    %     symmetric case (value = a scalar)
    %       absolute, no trans: |p* - nom| <= value
    %       absolute, log10: |log10(p*/nom)| <= value
    %       relative, no trans: |p* - nom|/|nom| <= value
    %       relative, log10: |log10(p*/nom)|/|log10(nom)| <= value
    %     asymmetric case (value = 1x2, 1st component negative)
    %       absolute, no trans: value(1) <= p* - nom <= value(2)
    %       absolute, log10: value(1) <= log10(p*/nom) <= value(2)
    %       relative, no trans: value(1) <= (p*-nom)/|nom| <= value(2)
    %       relative, log10: value(1) <= log10(p*/nom)/|log10(nom)| <= value(2)
    %   Although both transformation cases serve to bound the possible
    %   values of p*, because the coordinates in which this is expressed differ,
    %   sensitivity calculations with respect to the uncertainty in the
    %   object will differ.
    %
    %
    %   Object Properties:
    %       name: a unique name for the parameter. This should correspond to
    %           the names used by the ResponseModels.
    %       nominal: the nominal value for the parameter
    %       uncertainty: a scalar or 1x2 vector indicating the uncertainty
    %       uncertaintyType: either 'absolute' or 'relative'
    %       uncertaintyTransformation: either 'none' or 'log10'
    %       guiDisplayCallback: a function_handle to a display function
    %       description: a description of the parameter.
    %       complete: used internally to determine if sufficient properties are
    %           set for the object to be functional
    %
    %   Scalar FreeParameter objects can be vertically concatenated
    %   to form a larger dimensional FreeParameter object.
    %
    %   See also DCDataset, ResponseModel
    %
    %   Constructor reference page in Help browser
    %      <a href="matlab:web(fullfile(fileparts(which('DCsetup')), ...
    %      'docs', 'html', 'functions', 'FreeParameter.html'), ...
    %      '-helpbrowser')">FreeParameter constructor</a>
    %
    %   Object reference page in Help browser
    %      <a href="matlab:web(fullfile(fileparts(which('DCsetup')), ...
    %      'docs', 'html', 'obj_freeparam.html'), ...
    %      '-helpbrowser')">FreeParameter object</a>

    properties
        name = '';
        nominal
        uncertainty
        uncertaintyType
        uncertaintyTransformation
    end

    properties (Dependent=true, SetAccess=private)
        uncertaintyCase
        uncertaintyVector
        range
        status
    end



    methods

        %constructor
        function ParamObj = FreeParameter(name,nominal,uncertainty)

            ni = nargin;
            if ni>0
                switch ni
                    case 1 %initialize an empty name x 1 object
                        int = name;
                        assert(isnumeric(int) && isscalar(int) && ceil(int) == int && int > 0,'Usage: When called with a single input, FreeParameter must receive a positive integer');
                        ParamObj(int,1).name = '';
                    case 2
                        ParamObj.name = name;
                        ParamObj.nominal = mean(nominal);
                        ParamObj.uncertainty = diff(nominal)/2;
                        ParamObj.uncertaintyType = 'absolute';
                        ParamObj.uncertaintyTransformation = 'none';
                    case 3
                        ParamObj.name = name;
                        ParamObj.nominal = nominal;
                        if isa(uncertainty,'struct')
                            if isfield(uncertainty,'value') && isfield(uncertainty,'type') && isfield(uncertainty,'transformation')
                                ParamObj.uncertainty = uncertainty.value;
                                ParamObj.uncertaintyType = uncertainty.type;
                                ParamObj.uncertaintyTransformation = uncertainty.transformation;
                            else
                                error('Inputs: structure UNCERTAINTY contains insufficient fields')
                            end
                        else
                            ParamObj.uncertainty = uncertainty;
                            ParamObj.uncertaintyType = 'absolute';
                            ParamObj.uncertaintyTransformation = 'none';
                        end
                        %ParamObj.status = 1; %this is dependent now...so
                        %we don't set it?

                    otherwise
                        error('Usage: Incorrect number of inputs');
                end
            end

        end % end constructor

        function value=get.uncertaintyCase(ParamObj)
            if strcmp(ParamObj.uncertaintyType,'absolute') && strcmp(ParamObj.uncertaintyTransformation,'none')
                value = 1;
            elseif strcmp(ParamObj.uncertaintyType,'relative') && strcmp(ParamObj.uncertaintyTransformation,'none')
                value = 2;
            elseif strcmp(ParamObj.uncertaintyType,'absolute') && strcmp(ParamObj.uncertaintyTransformation,'log10')
                value = 3;
            elseif strcmp(ParamObj.uncertaintyType,'relative') && strcmp(ParamObj.uncertaintyTransformation,'log10')
                value = 4;
            else
                value = [];
            end
        end
        function value=get.uncertaintyVector(ParamObj)
            
            if ParamObj.status ~= 1
                error('Simultaneously getting the uncertaintyVector for all elements of multidim FreeParam requires each to be ''complete''')
            else
                if isscalar(ParamObj.uncertainty)
                    value = [-ParamObj.uncertainty ParamObj.uncertainty];
                else
                    value = ParamObj.uncertainty;
                end
            end
        end
        function range=get.range(ParamObj)
            if ParamObj.status ~= 1
                range = zeros(0,2);
                return
            else
                uncVect = ParamObj.uncertaintyVector;
                nom = ParamObj.nominal;
                
                % If relative, multiply by abs(nom) to make absolute
                if strcmp(ParamObj.uncertaintyType,'relative')
                    if strcmp(ParamObj.uncertaintyTransformation,'log10')
                        uncVect = uncVect*abs(log10(nom));
                    else
                        uncVect = uncVect*abs(nom);
                    end
                end

                if strcmp(ParamObj.uncertaintyTransformation,'log10')
                    range = nom*10.^uncVect;
                else
                    range = nom+uncVect;
                end
            end
        end
        function value=nParameters(ParamObj)
            value = sum([ParamObj.status]==1);
        end
        function out = get.status(ParamObj)
            if ~isempty(ParamObj.name) && ~isempty(ParamObj.nominal) && ~isempty(ParamObj.uncertainty)
                out = 1;
            elseif isempty(ParamObj.name) && isempty(ParamObj.nominal) && isempty(ParamObj.uncertainty)
                out = -1;
            else
                out = 0;
            end
        end
        
        function [list,sz] = displayProps(obj)
            list = [];
            sz = sprintf('%d-by-%d (%d complete parameters)',size(obj),nParameters(obj));
        end

        function idx = findIdx(FreeParam,names)
            % function idx = findIdx(FreeParam,names)
            %
            % This method locates the parameters indicated by the cell array names in
            % the multidimensional FreeParameter object FreeParam. idx(1)
            % corresponds to the location of names{1}, etc. This method only works on
            % "complete" objects, i.e., those for which all(status(obj)) = 1.
            %
            % Inputs
            %   FreeParam: A complete parameter assertion object
            %   names: a column or row cell array of unique strings
            % Outputs
            %   idx: A vector of the same dimension as names. If
            %     names{i} is not in ParamAssn, idx(i) will be NaN.
            %
            % See also FreeParameter, FreeParameter/status

            if ~all(vertcat(FreeParam.status)==1)
                error('Inputs: This method only works on "complete" objects, i.e., those for which all(status(obj)) = 1.')
            end
            if ~iscellstr(names) || length(names)~=length(unique(names)) || numel(names)~=length(names)
                error('Inputs: 2nd argument incorrect type, dimensions, or not unique')
            end

            idx = NaN(size(names));
            pnames = {FreeParam.name};
            [trash pnames2sortFoundNames inames2sortFoundNames] = intersect(pnames,names); 
            idx(inames2sortFoundNames) = pnames2sortFoundNames;


        end

        function [bool emptyIdx] = isempty(obj)
            %ISEMPTY Isempty for FreeParameter class
            %
            %   BOOL = ISEMPTY(OBJ) returns true if all componenents of a multidim OBJ
            %   are empty. A component is considered empty if it is missing a name,
            %   nominal, and uncertainty.
            %
            %   [BOOL EMPTYIDX] = ISEMPTY(OBJ) additionally returns a vector of
            %   indicies of the empty components.
            %
            %   See also FreeParameter
            
            error(nargoutchk(0,2,nargout))
            error(nargchk(1,1,nargin))

            emptyIdx = vertcat(obj.status) == -1;
            if all(emptyIdx)
                bool = true;
            else
                bool = false;
            end

            if nargout == 2
                emptyIdx = find(vertcat(obj.status) == -1);
            end

        end

        function FP = scaleUnc(FP,scaleFactor)
            %SCALEUNC multiplicatively alters the uncertainty in a FreeParameter object
            %
            %   NEWFP = SCALEUNC(FP,SCALEFACTOR) creates a new (possibly
            %   multidimensional) FreeParameter object NEWFP by multiplying the
            %   uncertainties in FP by the positive SCALEFACTOR.
            %
            %   This method exists primarily to be used by trustregion algorithms.
            %
            %   See also FreeParameter, FreeParameter/shiftNominal, CMeasTrust

            error(nargchk(2,2,nargin))
            error(nargoutchk(0,1,nargout))

            n = length(FP);

            if ~isnumeric(scaleFactor) || ~isequal(size(scaleFactor),[1 1]) || scaleFactor <= 0
                error('The scaleFactor supplied to scaleUnc must be a positive scalar')
            end

            for i1 = 1:n
                FP(i1).uncertainty = FP(i1).uncertainty*scaleFactor;
            end

        end

        function P = shiftNominal(P,newNominal)
            %SHIFTNOMINAL method of FreeParameter translates the rectangular domain
            %
            %   P = SHIFTNOMINAL(P,NEWNOMINAL) sets the nominal value of each component
            %   of the possibly multidimensional FreeParameter object P to the
            %   corresponding component of NEWNOMINAL. size(P) must equal
            %   size(NEWNOMINAL).
            %
            %   This method exists primarily to be used by trustregion algorithms.
            %
            %   See also FreeParameter, FreeParameter/scaleUnc, CMeasTrust

            error(nargchk(2,2,nargin))
            error(nargoutchk(0,1,nargout))

            n = length(P);
            if ~isnumeric(newNominal) || ~isequal(size(newNominal),[n 1])
                error('newNominal must be an nx1 double array.')
            end

            for i1 = 1:n
                if strcmp(P(i1).uncertaintyTransformation,'log10') && newNominal(i1) <= 0
                    error('For components with ''log10'' uncertainty type, nominal must be positive')
                end
                P(i1).nominal = newNominal(i1);
            end

        end


    end %end public methods

end %end classdef


