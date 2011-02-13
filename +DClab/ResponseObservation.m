classdef ResponseObservation < DClab.DCObject
    %RESPONSEOBSERVATION: constructor for object that describes the observed value of a physical quantity
    %
    %   ROBS = ResponseObservation(OBSVALUE,UNC)
    %   ROBS = ResponseObservation(OBSVALUE,UNC,@GUIDISPCALLBACK)
    %
    %   Examples:
    %   >> ROBS = ResponseObservation(10,2)
    %     creates an object where the measured value of y is 10 with an
    %     absolute uncertainty of 2, i.e. y = 10 +/- 2 .
    %   >> unc.value = 0.1, unc.type = 'relative', unc.transformation = 'none';
    %   >> ROBS = ResponseObservation(10,unc)
    %     creates an object where the measured value of y is 10 with a relative
    %     uncertainty of 0.1 (10%).
    %
    %   Input datatypes
    %
    %   OBSVALUE: a scalar containing the measured value.
    %   UNC: If the uncertainty in the observed value of y is symmetric, this
    %      should be a scalar. Alternatively it should be a 1x2 double array
    %      containing bounds on the measurement uncertainty s.t.
    %      OBSVALUE + UNC(1) <= y <= OBSVALUE + UNC(2).
    %   UNC[nonabsolute or transformation case]: If the uncertainty is
    %      expressed as a percent uncertainty or in transformed coordinates,
    %      OUTPUTUNC should be a structure with fields .value, .type, and
    %      .transformation. .value should be a scalar or a 1x2, .type can be
    %      'absolute' or 'relative', while transformation can be 'none' or
    %      'log10'. If either .type or .transformation are omitted, their
    %      respective defaults of 'absolute' and 'none' will be assumed. For
    %      the relative case, value should be a relative fraction, i.e. 0.1
    %      instead of 10%. Let
    %         case 1 denote .type = 'absolute' and .transformation = 'none'.
    %         case 2 denote .type = 'relative' and .transformation = 'none'.
    %         case 3 denote .type = 'absolute' and .transformation = 'log10'.
    %         case 4 denote .type = 'relative' and .transformation = 'log10'
    %      The uncertainty expressions for the four cases are
    %         symmetric uncertainty (.value = a scalar)
    %            case 1: |y - d| <= .value
    %            case 2: |y - d|/d <= .value
    %            case 3:  |log10(y/d)| <= .value
    %            case 4: |log10(y/d)|/log10(d) <= .value
    %         asymmetric uncertainty (.value = 1x2, 1st component negative)
    %            case 1: .value(1) <= y - d <= .value(2)
    %            case 2: .value(1) <= (y-d)/d <= .value(2)
    %            case 3: .value(1) <= log10(y/d) <= .value(2)
    %            case 4: .value(1) <= log10(y/d)/log10(d) <= .value(2)
    %      Limitations:
    %         case 1: none.
    %         case 2: d must be postive.
    %         case 3: both y and d must be positive.
    %         case 4: y must be positive and d must be greater than 1.
    %      Although each case serves to bound the possible values of y, because
    %      the coordinates in which this is expressed differ, sensitivity
    %      calculations with respect to the uncertainty in the object will
    %      differ. Additionally, the expressions for the consistency measure
    %      differ somewhat depending on how the uncertainty is expressed.
    %   @GUIDISPCALLBACK[optional]: function handle to an mfile that will be
    %      called with the syntax GUIDISPCALLBACK(ROBJ) that displays further
    %      info about the object ROBS. Typically such an mfile would make use
    %      of the .userdata property of the object in order to identify it and
    %      launch an appropriate graphical display (e.g., webpage).
    %
    %   See also ResponseModel, ModelAndObservationPair
    %
    %   Constructor reference page in Help browser
    %      <a href="matlab:web(fullfile(fileparts(which('DCsetup')), ...
    %      'docs', 'html', 'functions', 'ResponseObservation.html'), ...
    %      '-helpbrowser')">ResponseObservation constructor</a>
    %
    %   Object reference page in Help browser
    %      <a href="matlab:web(fullfile(fileparts(which('DCsetup')), ...
    %      'docs', 'html', 'obj_expass.html'), ...
    %      '-helpbrowser')">ResponseObservation object</a>
    
    properties
        observedValue;
        uncertainty;
        uncertaintyType = 'absolute';
        uncertaintyTransformation = 'none';
        guiDisplayCallback;
    end
    
    properties (Dependent)
        uncertaintyPlusMinus;
        range;
    end
    
    properties (Dependent, SetAccess=private)
        uncertaintyCase;
    end
    
    methods
        function Obj = ResponseObservation(obsValue,unc,guiDispCallback)
            
            
            error( nargoutchk(0,1,nargout) );
            ni = nargin;
            error( nargchk(0,3,ni) );
            
            if ni==1
                error('Usage: Incorrect number of input arguments')
            elseif ni>1
                Obj.observedValue = obsValue;
                if isa(unc,'struct')
                    if ~isempty(setdiff(fieldnames(unc),{'value','type','transformation'}))
                        warning('DCLAB:unusedInputs','structure UNC contains unused fields')
                    end
                    assert(isfield(unc,'value'),'Inputs: structure UNC contains insufficient fields')
                    Obj.uncertainty = unc.value;
                    if isfield(unc,'type')
                        Obj.uncertaintyType = unc.type;
                    else
                        Obj.uncertaintyType = 'absolute';
                    end
                    if isfield(unc,'transformation')
                        Obj.uncertaintyTransformation = unc.transformation;
                    else
                        Obj.uncertaintyTransformation = 'none';
                    end
                    
                    
                else
                    Obj.uncertainty = unc;
                    Obj.uncertaintyType = 'absolute';
                    Obj.uncertaintyTransformation = 'none';
                end
                Obj.guiDisplayCallback = [];
            end
            
            if ni == 3
                Obj.guiDisplayCallback = guiDispCallback;
            end
            
        end
        
        function bool = isempty(obj)
            
            bool=false;
            if isempty(obj.observedValue) || isempty(obj.uncertainty)
                bool = true;
            end
        end
        function [list,sz] = displayProps(obj)
            unc = obj.uncertainty;
            val = obj.observedValue;
            type = obj.uncertaintyType;
            trans = obj.uncertaintyTransformation;
            if length(unc)<2
                unc = [-unc unc];
            end
            if isempty(unc)
                list = [];
                sz = 'empty';
            else
                switch type
                    case 'absolute'
                        switch trans
                            case 'none'
                                list{1} = sprintf('%0.2g <= y - %0.2g <= %0.2g',unc(1),val,unc(2));
                            case 'log10'
                                list{1} = sprintf('%0.2g <= log10(y/%0.2g) <= %0.2g',unc(1),val,unc(2));
                        end
                    case 'relative'
                        switch trans
                            case 'none'
                                list{1} = sprintf('%0.2g <= (y - %0.2g)/(%0.2g) <= %0.2g',unc(1),val,val,unc(2));
                            case 'log10'
                                list{1} = sprintf('%0.2g <= log10(y/%0.2g)/log10(%0.2g) <= %0.2g',unc(1),val,val,unc(2));
                        end
                end
                sz = '';
            end
        end
        
        function out=get.uncertaintyPlusMinus(obj)
            if isscalar(obj.uncertainty)
                out = [-obj.uncertainty obj.uncertainty];
            else
                out = obj.uncertainty;
            end
        end
        function out=get.range(obj)
            if isscalar(obj.uncertainty)
                uncPM = [-obj.uncertainty obj.uncertainty];
            else
                uncPM = obj.uncertainty;
            end
            if isempty(obj.uncertainty) || isempty(obj.observedValue)
                out = [];
            elseif obj.uncertaintyCase == 1
                out = obj.observedValue+uncPM;
            elseif obj.uncertaintyCase == 2
                out = obj.observedValue*(1+uncPM);
            elseif obj.uncertaintyCase == 3
                out = obj.observedValue*10.^uncPM;
            elseif obj.uncertaintyCase == 4
                out = obj.observedValue.^(1+uncPM);
            else
                error('Internal inconsistency, condition should never occur')
            end
            
        end
        function out = get.uncertaintyCase(obj)
            if strcmp(obj.uncertaintyType,'absolute') && strcmp(obj.uncertaintyTransformation,'none')
                out = 1;
            elseif strcmp(obj.uncertaintyType,'relative') && strcmp(obj.uncertaintyTransformation,'none')
                out = 2;
            elseif strcmp(obj.uncertaintyType,'absolute') && strcmp(obj.uncertaintyTransformation,'log10')
                out = 3;
            else
                out = 4;
            end
        end
        
    end %methods
    
end %classdef
