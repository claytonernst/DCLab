function gui(obj)
%how to call gui with object
%builds DC objects from DCContainer
%calls other gui6


%build FreeParameter object
fpLen = length(obj.parameterKey);
fp = cell(fpLen,1);
for i=1:fpLen
    d = mean(obj.parameterRange{i});
    unc.value = obj.parameterRange{i}-d;
    unc.type = 'absolute';
    unc.transformation = 'none';
    %fp{i} = FreeParameter(obj.parameterKey{i},d,unc);
    fp{i} = DClab.FreeParameter(obj.parameterID{i},d,unc);
    fp{i} = set(fp{i},'userData',obj.parameterKey{i});
end
fp = vertcat(fp{:});

%build MOP object
mopLen = length(obj.modelKey);
mop = cell(mopLen,1);
for i=1:mopLen
    d = mean(obj.uncRange{i});
    unc.value = obj.uncRange{i}-d;
    unc.type = 'absolute';
    unc.transformation = 'none';
    respObs = DClab.ResponseObservation(d,unc);
    modStruct.value = obj.modelCoeffs{i};
    modStruct.responseTransformation = 'none';
    m = length(obj.modelParamIDs{i});
    modStruct.variableTransformations = repmat({'none'},m,1);
    modelDomain = struct('name',cell(m,1),'range',cell(m,1));
    for j=1:m
        idx = find(strcmp(obj.parameterID,obj.modelParamIDs{i}{j}));
        %modelDomain(j).name = obj.parameterKey{idx};
        modelDomain(j).name = obj.parameterID{idx};
        modelDomain(j).range = obj.parameterRange{idx};
    end
    respMod = DClab.ResponseModel(modStruct,modelDomain);
    set(respMod,'name',obj.modelKey{i});
    mop{i} = DClab.ModelAndObservationPair(respObs,respMod,obj.modelKey{i});
end
mop = vertcat(mop{:});

%targetmodels
%targets = vertcat(get(mop,'ResponseModel'));
targets = {mop.ResponseModel}';

%open gui:
gui6(mop,fp,targets)

%ds = DCDataset(mop,fp);