classdef PolyDataset < DClab.DCDataset
    %POLYDATASET Construct a PolyDataset object
    %
    %   PDSET = POLYDATASET creates an empty object
    %
    %   PDSET = POLYDATASET(DSET) creates a PolyDataset object from the DCDataset
    %   object DSET.
    %
    %   PDSET = POLYDATASET(DSET,OPTS) uses the display and fitting options
    %   from the DCOptions object OPTS.
    %
    %   Description:
    %   Users may create PolyDataset objects, but they are primarily created
    %   internally by analysis functions such as ConsistencyTest and
    %   ResponsePrediction.
    %
    %   The PolyDataset class inherits from DCDataset. Additional fields are
    %   .displaySettings: contains the display settings from a DCOptions object.
    %   .PiecewiseSurrogateModelTree: implements a binary tree that contains
    %      response surfaces developed for each ResponseModel in the
    %      ModelAndObservationPairs., together with the surrogate fitting
    %      settings from a DCOptions object. Even if the ensemble of
    %      ResponseModels require all n parameters, the domain of the surrogate
    %      models will be n-dimensional.
    %
    %   A PolyDataset is "partitionable". I.e., the response surfaces can be
    %   defined piecewise over rectangular subdomains of the rectangular domain
    %   specified by the FreeParameters. This behavior is performed
    %   automatically by the PARTITIONDOMAIN method of the PolyDataset class,
    %   and works by calling methods of PiecewiseSurrogateModelTree.
    %
    % See also DCDataset, PiecewiseSurrogateModelTree
    %
    %   Constructor reference page in Help browser
    %      <a href="matlab:web(fullfile(fileparts(which('DCsetup')), ...
    %      'docs', 'html', 'functions', 'PolyDataset.html'), ...
    %      '-helpbrowser')">PolyDataset constructor</a>
    %
    %   Object reference page in Help browser
    %      <a href="matlab:web(fullfile(fileparts(which('DCsetup')), ...
    %      'docs', 'html', 'obj_polydset.html'), ...
    %      '-helpbrowser')">PolyDataset object</a>

    properties
        displaySettings;
        PiecewiseSurrogateModelTree = DClab.PiecewiseSurrogateModelTree;
        %userData; %already defined in DCDataset
    end
    
    properties (Dependent)
        fittingSettings;
        nSubdivisions;
        leafNodes;
    end
    
    methods
        function PD = PolyDataset(Dset,opt)

            ni = nargin;
            error(nargoutchk(0,1,nargout))
            
            if ni==0
                Dset=DClab.DCDataset;
            end
            PD = PD@DClab.DCDataset(Dset);

            switch ni
                case 0
                    return
                case 1
                    opt = DClab.DCOptions;
                otherwise
                    error(nargchk(1,2,ni));
            end

            assert(isa(Dset,'DClab.DCDataset') && ~isempty(Dset),'The first input to POLYDATASET must be nonempty and of type DCDataset.')
            assert(isa(opt,'DClab.DCOptions'),'The second input to POLYDATASET must be of type DCOptions.')

            %construct the surrogateTree
            FP = Dset.FreeParameter;
            paramList = {FP(:).name}';
            rng= vertcat(FP(:).range);
            domain1 = DClab.createDomainStructure(paramList,rng); %actual domain

            for i1 = 1:nParameters(FP)
                FP(i1).uncertainty=1.05*FP(i1).uncertainty;
            end
            rng= vertcat(FP(:).range);
            domain2 = DClab.createDomainStructure(paramList,rng); %for training

            PD.displaySettings = decompose(opt);
            if Dset.nPairs == 0
                PD.PiecewiseSurrogateModelTree = DClab.PiecewiseSurrogateModelTree({},domain1,domain2,[],[],opt,true);
            else
                m = Dset.nPairs;
                RMCell = cell(m,1);
                critRng = zeros(m,2);
                trnRng = zeros(m,2);
                for i1 = 1:m
                    RMCell{i1} = Dset.ModelAndObservationPair(i1).ResponseModel;
                    critRng(i1,:) = Dset.ModelAndObservationPair(i1).criticalRange;
                    trnRng(i1,:) = Dset.ModelAndObservationPair(i1).trainingRange;
                end
                PD.PiecewiseSurrogateModelTree = DClab.PiecewiseSurrogateModelTree(RMCell,domain1,domain2,critRng,trnRng,opt,true);
            end
        end
        
        function out=get.fittingSettings(PD)
            out = {PD.PiecewiseSurrogateModelTree.fittingSettings}';
        end
        function out=get.nSubdivisions(PD)
            out = nSubdivisions(PD.PiecewiseSurrogateModelTree);
        end
        function out=get.leafNodes(PD)
            out = leafNodes(PD.PiecewiseSurrogateModelTree);
        end
        function [list,sz] = displayProps(obj)
            list = [];
            sz = '';
        end
        
        function bool = isPerfectFit(PD)
            LN = PD.leafNodes;
            bool = true;
            for i1=1:length(LN)
                if PD.PiecewiseSurrogateModelTree.DCSurface(LN(i1)).peakError>0
                    bool=false;
                    break;
                end
            end
        end
    end %public methods
    
end %classdef