classdef XSect < DClab.DCObject
    % xsectObj = XSect(Dset,param2use,unit2use,Npts,innerOrTrue,x)
    % xsectObj = XSect(PDset,param2use,unit2use,Npts,innerOrTrue,x)
    % 3DXSectObj = XSect(2DXSectObj,newParam)
    %
    % XSECT is the constructor for 2- or 3-d feasible set X-section class
    %  Inputs:
    %    Dset/PDset: a DCDataset or PolyDataset object
    %    param2use: a vector containing the parameter numbers or a cell
    %      array containing the parameter names to use for plotting. Must
    %      be of length 2 or 3 (the length determines the dimension).
    %    unit2use[opt]: a vector containing the dataset unit numbers or a
    %      cell array containing the Unit names to use for plotting.
    %    Npts[opt]: controls how finely the unit cube is gridded.
    %      Numbers between 25 and 500 are typical, default is 50.
    %    innerOrTrue[opt]: 0 means use inner approx, 1 means use true
    %      models. default=0.
    %    x[opt]: the parameter vector at which to freeze all params not in
    %      param2use. Must satisfy the FreeParameter bounds of PDset.
    %      Default is the FreeParameter nominal values.
    %
    % Can also be called with an already created two-dimensional object,
    % an a new parameter.  This will elevate the object to a 3D object
    % using the previous points and options.
    %
    % See also XSect/plot
    %

    properties
        PDset;
        units = '';
        params = '';
        paramBnds;
        x;
        Npts;
        innerOrTrue = '';
        combFeasSet;
        combOFeasSet;
        indvFeasSet = {};
        indvOFeasSet = {};
        ndims;
    end

    methods
        function xsect = XSect(varargin)

            error( nargoutchk(0,1,nargout) );
            ni = nargin;
            if ni>0

                % Input error checking and define any defaults
                [PDset param2use unit2use Npts innerOrTrue x] = getInputs(varargin{:});

                xsect.ndims = length(param2use);

                unc = vertcat(PDset.FreeParameter(:).range);
                if any( unc(:,1) > x | unc(:,2) < x)
                    error('Inputs: x does not satisfy the parameter assertions of the PolyDataset')
                end

                lb = unc(param2use,1);
                ub = unc(param2use,2);

                % Generate Npts points at the center of the boxes of an Npts+1 grid
                x1 = linspace(lb(1),ub(1),Npts+1);
                x1 = x1(2:end)-0.5*diff(x1(1:2));
                x2 = linspace(lb(2),ub(2),Npts+1);
                x2 = x2(2:end)-0.5*diff(x2(1:2));
                if xsect.ndims==3
                    x3 = linspace(lb(3),ub(3),Npts+1);
                    x3 = x3(2:end)-0.5*diff(x3(1:2));
                end

                if xsect.ndims==2
                    [X Y] = meshgrid(x1,x2);
                elseif xsect.ndims==3
                    [X Y Z] = meshgrid(x1,x2,x3);
                end

                %=========Generate feasible set data======================

                % Reshape parameter values to vectorize code
                N = numel(X);
                XX = reshape(X,1,N);
                YY = reshape(Y,1,N);
                if xsect.ndims==3
                    ZZ = reshape(Z,1,N);
                end
                clear X Y Z

                P = repmat(x,1,N);
                if xsect.ndims==2
                    P(param2use,:) = [XX;YY];
                else
                    P(param2use,:) = [XX;YY;ZZ];
                end

                % Initialize combined feasible set as all true
                combFeasSetVect = true(1,N);
                combOFeasSetVect = true(1,N);

                m = length(unit2use);
                feasSetVect = cell(m,1);
                ofeasSetVect = cell(m,1);
                %t = PDset.options.trans;
                for i1 = 1:m
                    unit = unit2use(i1);
                    range = PDset.ModelAndObservationPair(unit).ResponseObservation.range;

                    % Compute outer feasible set
                    if i1==1
                        [yApprox yBnds yInts node] = evalSurrogateModels(PDset,P,unit); %#ok
                    else
                        [yApprox yBnds yInts node] = evalSurrogateModels(PDset,P,unit,PDset.parameterList,node); %#ok
                        %[y, yInt1, yInt2, node] = evalSurrogateModels(PD,x,pairs,paramList,node)
                    end

                    ylb = yBnds{1}(1,:);
                    yub = yBnds{1}(2,:);
                    ofeasSetVect{i1} = (range(1)<=min(yub,[],1)) & (max(ylb,[],1) <= range(2));

                    % Compute either the true feasible set or the inner approximation
                    if strcmp(innerOrTrue,'true')
                        y = eval(PD.ModelAndObservationPair(unit).ResponseModel,P,PDset.parameterList);
                        feasSetVect{i1} = (range(1)<=y(unit).actual) & (y(unit).actual <= range(2));
                    else
                        feasSetVect{i1} = (range(1)<=max(ylb,[],1)) & (min(yub,[],1) <= range(2));
                    end
                    combFeasSetVect = combFeasSetVect & feasSetVect{i1};
                    combOFeasSetVect = combOFeasSetVect & ofeasSetVect{i1};
                end

                unitNames = {PDset.ModelAndObservationPair(:).name}';
                paramNames = {PDset.FreeParameter(:).name}';

                xsect.PDset = PDset;
                xsect.units = unitNames(unit2use);
                xsect.params = paramNames(param2use);
                xsect.paramBnds = [lb ub];
                xsect.x = x;
                xsect.Npts = Npts;
                xsect.innerOrTrue = innerOrTrue;
                xsect.combFeasSet = combFeasSetVect;
                xsect.combOFeasSet = combOFeasSetVect;
                xsect.indvFeasSet = feasSetVect;
                xsect.indvOFeasSet = ofeasSetVect;
            end
        end

        function [list,sz] = displayProps(obj)
            list = [];
            sz = '';
        end

    end %public methods

end %classdef

%LOCAL Function
function [PDset param2use unit2use Npts innerOrTrue x] = getInputs(PDset,param2use,unit2use,Npts,innerOrTrue,x)

ni = nargin;

%parse input list and insert defaults
switch ni
    case 1
        error('Incorrect number of input arguments');
    case 2
        unit2use = [];
        Npts = [];
        innerOrTrue = 'inner';
        x = [];
    case 3
        Npts = [];
        innerOrTrue = 'inner';
        x = [];
    case 4
        innerOrTrue = 'inner';
        x = [];
    case 5
        x = [];
    case 6
        %  do nothing
    otherwise
        error( nargchk(0,5,ni) );
end

if ~isa(PDset,'DClab.PolyDataset') && isa(PDset,'DClab.DCDataset')
    PDset = DClab.PolyDataset(PDset);
elseif ~isa(PDset,'DClab.DCDataset')
    error('First input must be either a DCDataset or PolyDataset object')
end

% Massage unit2use to get a column vector
if isempty(unit2use)
    unit2use = (1:PDset.nPairs)';
elseif isnumeric(unit2use)
    if length(unit2use)~=numel(unit2use) || max(unit2use) > PDset.nPairs
        error('Inputs: invalid vector unit2use')
    end
elseif iscell(unit2use)
    uList = {PDset.ModelAndObservationPair.name}'; %asked for unitList?
    if ~iscellstr('ischar',unit2use) || length(unit2use)~=numel(unit2use) || ...
            ~isempty(setdiff(unit2use,uList))
        error('Inputs: invalid cell array unit2use')
    end
    [trash unit2sort] = sort(unit2use); %#ok
    [trash sort2unit] = sort(unit2sort); %#ok
    [trash PDUnit2sort] = intersect(uList,sunit2use); %#ok
    unit2use = PDUnit2sort(sort2unit);
else
    error('Inputs: invalid unit2use')
end
unit2use = reshape(unit2use,length(unit2use),1);

% Make sure param2use is a valid vector of parameter numbers or cell
% array of parameter names
if isnumeric(param2use)
    if (numel(param2use)~=2 && numel(param2use)~=3) || max(param2use) > PDset.nParameters
        error('Inputs: invalid vector param2use')
    end
elseif iscell(param2use)
    pList = PDset.parameterList;
    if ~iscellstr(param2use) || (numel(param2use)~=2 && numel(param2use)) || ...
            ~isempty(setdiff(param2use,pList))
        error('Inputs: invalid cell array param2use')
    end
    [trash param2sort] = sort(param2use); %#ok
    [trash sort2param] = sort(param2sort); %#ok
    [trash PDParam2sort] = intersect(pList,param2use); %#ok
    param2use = PDParam2sort(sort2param);
else
    error('Inputs: invalid param2use')
end
param2use = reshape(param2use,length(param2use),1);

if isempty(Npts)
    Npts = 50;
end

if isempty(innerOrTrue)
    innerOrTrue = 'inner';
elseif ~any(strmatch(innerOrTrue,{'inner','true'}))
    error('Inputs: innerOr true must be ''inner'' or ''true''')
else
    %do nothing
end

if isempty(x)
    x = vertcat(PDset.FreeParameter(:).nominal);
elseif ~isnumeric(x) || length(x)~=numel(x) || length(x)~=PDset.nParameters
    error('Inputs: invalid dimensions for x')
else
    x = reshape(x,length(x),1);
end
end

