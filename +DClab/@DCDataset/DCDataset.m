classdef DCDataset < DClab.DCObject
    %DCDATASET Construct Dataset object
    %
    %   D = DCDATASET(MOPAIRS,FREEPARAMS) creates a dataset from the array of
    %   ModelAndObservationPair objects MOPAIRS and the array of FreeParameter
    %   objects FREEPARAMS.
    %
    %   Limitations: To ensure everything has a unique name, MOPAIRS must be a
    %   scalar empty or a multidimensional object having only "complete"
    %   (status==1) elements. Similar restrictions apply to FREEPARAMS.
    %
    %   See also ModelAndObservationPair, FreeParameter
    %
    %   Constructor reference page in Help browser
    %      <a href="matlab:web(fullfile(fileparts(which('DCsetup')), ...
    %      'docs', 'html', 'functions', 'DCDataset.html'), ...
    %      '-helpbrowser')">DCDataset constructor</a>
    %
    %   Object reference page in Help browser
    %      <a href="matlab:web(fullfile(fileparts(which('DCsetup')), ...
    %      'docs', 'html', 'obj_dset.html'), ...
    %      '-helpbrowser')">DCDataset object</a>

    % Last modified 3/26/09 tmr

    properties
        ModelAndObservationPair = DClab.ModelAndObservationPair;
        FreeParameter = DClab.FreeParameter;
        name = '';
    end
    properties (Dependent=true,SetAccess=private)
        parameterList
        pairList
        nPairs
        nParameters
    end %dependent properties


    methods

        %dependent gets
        function pL = get.parameterList(obj)
            pL = {obj.FreeParameter.name}';
        end
        function pL = get.pairList(obj)
            pL = {obj.ModelAndObservationPair.name}';
        end
        function nP = get.nPairs(obj)
            nP = nPairs(obj.ModelAndObservationPair);
        end
        function nP = get.nParameters(obj)
            nP = nParameters(obj.FreeParameter);
        end

        %constructor
        function D = DCDataset(MOPairs,FreeParams)


            %TODO: I thought both the parameter assertions and dataset units had to be
            %"complete" in order to form a dataset from them -- rpf 6/27/07

            ni = nargin;
            no = nargout;
            error( nargoutchk(0,1,no) );

            switch ni
                case 0
                    %defaults
                case 1
                    assert(isa(MOPairs,'DClab.DCDataset'),'Usage: incorrect number of input arguments');
                    D.ModelAndObservationPair = MOPairs.ModelAndObservationPair;
                    D.FreeParameter = MOPairs.FreeParameter;
                    D.name = MOPairs.name;
                    D.userData = MOPairs.userData;
                case 2
                    if isempty(MOPairs)
                        D.ModelAndObservationPair = DClab.ModelAndObservationPair;
                    else
                        D.ModelAndObservationPair = MOPairs;
                    end
                    if isempty(FreeParams)
                        D.FreeParameter = DClab.FreeParameter;
                    else
                        D.FreeParameter = FreeParams;
                    end
                    D.name = '';
                otherwise
                    error( nargchk(0,2,ni,'struct') );
            end

        end %constructor
        
        function [list,sz] = displayProps(obj)
            list = {sprintf('%d parameters',obj.nParameters);...
                sprintf('%d model-and-observation pairs',obj.nPairs)};
            sz = '';
        end
        
        function deleteSavedEvaluations(obj)
           %Deletes all saved evaluations associated with the DCDataset,
           %i.e. all points associated with each of its ResponseModels
           m=obj.nPairs ;
           for i1=1:m
               deleteSavedEvaluations(obj.ModelAndObservationPair(i1).ResponseModel);
           end
        end

        
    end %methods


end %classdef