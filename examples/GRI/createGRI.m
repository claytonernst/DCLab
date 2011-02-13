function D = createGRI(flag1,flag2,flag3)
%CREATEGRI assembles the GRI-Mech Dataset.
%
%   D = createGRI returns the GRI dataset of 77 dataset units over a 102
%   dimensional parameter domain.
%
%   D = createGRI(modifiedFlag,realvarsFlag,versionFlag) creates the dataset
%   with a specified set of options.  Any number of these inputs (0 to 3)
%   can be included or left empty to use the defaults.
%       modifiedFlag:  can be 'newunc' (default), 'modified', or
%           'original'.  The modified data is using the corrected Stanford
%           data.  The original data is pre-correction.  The new
%           uncertainty data is the corrected Stanford data, and the new
%           uncertainties provided by Xiaoqing You.
%       newuncFlag:  can be 'transvars' (default) or 'realvars'.
%           Transformed variables have the variable bounds of [-1, 1], and
%           all sensitivities are in these units.  The real variables use
%           the original units.
%       versionFlag:  can be '3.0' (default) or '2.11'.  Changes the
%           nominal parameter vector to a specific GRI-Mech version.  Only
%           used if modifiedFlag is set to 'newunc'.

%RecursionLimit = get(0,'RecursionLimit');
%set(0,'RecursionLimit',500);

import DClab.*

ni = nargin;
switch ni
    case 0
        flag1 = 'newunc';
        flag2 = 'transVars';
        flag3 = '3.0';
    case 1
        flag2 = 'transVars';
        flag3 = '3.0';
        if isempty(flag1)
            flag1='newunc';
        end
    case 2
        flag3 = '3.0';
        if isempty(flag1)
            flag1 = 'newunc';
        end
        if isempty(flag2)
            flag2 = 'transVars';
        end
    case 3
        if isempty(flag1)
            flag1 = 'newunc';
        end
        if isempty(flag2)
            flag2 = 'transVars';
        end
        if isempty(flag3)
            flag3 = '3.0';
        end
    otherwise
        error(nargchk(ni,0,2,'struct'))
end

switch lower(flag1)
    case {'original','modified','newunc'}
        %cool
    otherwise
        error('Incorrect flag.  First flag must be ''original'', ''modified'', ''newunc'', or empty');
end
switch lower(flag2)
    case {'transvars','realvars'}
        %cool
    otherwise
        error('Incorrect flag.  Second flag must be ''transvars'', ''realVar'', or empty');
end
switch lower(flag3)
    case {'3.0','2.11'}
        %cool
    otherwise
        error('Incorrect flag.  Third flag must be ''3.0'', ''2.11'', or empty');
end

load GRIDataStruct
trg = GRIDataStruct; %this is essentially Michael's GRI target object

t = load('GRIspans.mat');

%load nominal values, names, and Xiaoqing's uncertainties
switch flag3
    case '3.0'
        load GRIuncVec3
    case '2.11'
        load GRIuncVec211
end


%Create cell arrays to contain the 77 model assertions, experiment
%assertions, and dataset units.
n = length(trg);
RO = cell(n,1);
RM = cell(n,1);
Pair = cell(n,1);

if strcmpi(flag1,'modified') || strcmpi(flag1,'newunc')
    % Do nothing, the GRIDataStruct has the modified values in it.
elseif strcmpi(flag1,'original')
    trg(57).value = 970;
    trg(58).value = 218;
end

allnums = unique([trg.var]);

for i1 = 1:n
    %measured data
    d = trg(i1).value;

    %only keeping the data to four decimal places may give some of our guis
    %an easier time
    if strcmpi(flag2,'transvars')
        d = round(10000*log10(d))/10000;
    else
        d = round(10000*d)/10000;
    end

    if strcmpi(flag1,'newunc')
        u.value = uncVec(i1,:);
    else
        u.value = 0.1;
    end
    u.type='absolute';
    if strcmpi(flag2,'transvars')
        u.transformation='none';
    else
        u.transformation='log10';
    end
    %Uncertainty model is absolute uncertainty on a log transformation:
    % log10(d) +/- value

    % Create the i1^th response observation. The constructor syntax is
    %   Obj = ResponseObservation(data,unc,dispCallback)
    RO{i1} = ResponseObservation(d,u);

    % Create the i1^th response model. The constructor syntax is
    %   Obj = ResponseModel(modelStruct,domain,outputUnc,dispCallback)
    % domain is an nx1 struct array with fields .name and .range
    
    [trash,ia,ib] = intersect(allnums,trg(i1).var);
    pNames = cell(length(trg(i1).var),1);
    pNames(ib) = allnames(ia); %strtrim(cellstr(num2str(trg(i1).var')));
    pRanges = {};
    [trash,idx1,idx2] = intersect(t.spans(:,1),trg(i1).var');
    spantmp = zeros(length(idx2),1);
    spantmp(idx2) = t.spans(idx1,2);

    if strcmpi(flag2,'realvars')
        for i2=1:length(pNames)
            pRanges(i2,1) = {[1/spantmp(i2) spantmp(i2)]};
            pRanges(i2,1) = {[-inf inf]};
        end
    else
        for i2=1:length(pNames)
            pRanges(i2,1) = {[-inf inf]};
        end
    end
        
    paramDomain = struct('name',pNames,'range',pRanges);

    if strcmpi(flag2,'realvars')
        %Adjust quadratic. trg(i1).quadform is a polynomial function of
        %normalized variables x\in [-1,1]. We'd like to transform these to
        %polynomials in log10(\rho) \in [-log10(span) log10(span)].
        %
        % (x+1)/2 = (log10(\rho)+log10(span))/(2*log10(span))

        % We thus want to compose the quadratics with the affine transformation
        %
        % x = (log10(\rho) + log10(span))/log10(span) - 1 =
        % log10(\rho)/log10(span)
        %
        B = zeros(length(idx2),1);
        A = diag(1./log10(spantmp));

        quadformtmp = DCSurface.composeQuadWithAffine(0.5*trg(i1).quadform,A,B);

        %  spanadjust = repmat([1; 1./log10(spantmp)],1,length(idx2)+1);
        %  quadformtmp = trg(i1).quadform.*spanadjust.*spanadjust'/2;
        quadform.value = quadformtmp;
        quadform.responseTransformation = 'log10';
        quadform.variableTransformations = repmat({'log10'},length(idx2),1);
        % Call the constructor
        RM{i1} = ResponseModel(quadform,paramDomain);
    else
        RM{i1} = ResponseModel(0.5*trg(i1).quadform,paramDomain);
    end
    %set the name property
    RM{i1} = set(RM{i1},'name',[trg(i1).type ' (' trg(i1).key '): ' trg(i1).comment]);

    % Create a dataset unit from the experiment and model
    % assertions. The constructor syntax is
    %   Obj = DatasetUnit(EA,MA,name,dispCallback)
    % the name value will be used in gui displays and the callback function
    % is called by certain gui actions (double clicks)
    %
    % let's make a anonymous function callbacks that will open the
    % GRI-MECH 3.0 webpage for each dataset unit
    key = trg(i1).key;

    name = [trg(i1).type ' (' trg(i1).key '): ' trg(i1).comment];
    Pair{i1,1} = ModelAndObservationPair(RO{i1},RM{i1},name);
end

% Create 102 parameter assertions. This object can be multidimensional
% and has constructor syntax
%   Obj = ParameterAssertion(name,range,nominal)
params = unique([trg.var]);
n = length(params);
P = cell(n,1);
if strcmpi(flag2,'realvars')
    for i1 = 1:n
        name = allnames{i1};
        %there isn't a webpage for parameter 999
        u.type = 'absolute'; u.transformation = 'log10';
        if strcmpi(flag1,'newunc')
            u.value = [-log10(t.spans(i1,2)) log10(t.spans(i1,2))]-p0(i1)*log10(t.spans(i1,2));
            P{i1} = FreeParameter(name,t.spans(i1,2)^p0(i1),u);
        else
            u.value = [-log10(t.spans(i1,2)) log10(t.spans(i1,2))];
            P{i1} = FreeParameter(name,1,u);
        end

    end
else
    for i1 = 1:n
        name = allnames{i1}; %num2str(params(i1));
        u.type = 'absolute'; u.transformation = 'none';
        if strcmpi(flag1,'newunc')
            u.value = [-1 1]-p0(i1);
            P{i1} = FreeParameter(name,p0(i1),u);
        else
            u.value = [-1 1]; 
            P{i1} = FreeParameter(name,0,u);
        end

    end
end
FreeParam = vertcat(P{:});

% Finally create the dataset. This object includes all created
% dataset units and the multidimensional parameter assertion
% object. It has constructor syntax
%  D = Dataset(DatasetUnitCellArray,ParamAss)
D = DCDataset(vertcat(Pair{:}),FreeParam);
D.name = 'GRI-Mech 3.0 Dataset';



