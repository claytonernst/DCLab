function gui(obj)
%how to call gui with object
%builds DC objects from DCContainer
%calls other gui6


%build FreeParameter object
fpLen = length(obj.parameterKey);
fp = cell(fpLen,1);
for i=1:fpLen
    d = obj.paramInitialValue{i};
    unc.value = obj.parameterRange{i}-d;
    unc.type = 'absolute';
    unc.transformation = 'none';
%    cb = @(fpObj) obj.callback2show(obj.parameterID{i});
%    fp{i} = DClab.FreeParameter(obj.parameterID{i},d,unc,cb);
    fp{i} = DClab.FreeParameter(obj.parameterID{i},d,unc);    
    info.parameterKey=obj.parameterKey{i};
    info.parameterID=obj.parameterID{i};
    info.parameterRange=obj.parameterRange{i};
    info.paramInitialValue=obj.paramInitialValue{i};
    info.paramTrans=obj.paramTrans{i};
    info.centerSpan=obj.centerSpan{i};
    info.optVarUnits=obj.optVarUnits{i};
    info.parameterLinks=obj.parameterLinks{i};
    info.trialModelID=obj.trialModelID;
    info.callback2show=obj.callback2show;
    info.callback2show2=obj.callback2show2;
    fp{i}.userData = info;
end
fp = vertcat(fp{:});

%build MOP object
mopLen = length(obj.targetID); %length(obj.modelKey);
mop = cell(mopLen,1);
for i=1:mopLen

    %RO
    d = obj.targetValue{i};
    unc.value = obj.uncRange{i}-d;
    unc.type = 'absolute';
    unc.transformation = 'none';
 %   cb = @(roObj) obj.callback2show(obj.targetID{i});
 %   respObs = DClab.ResponseObservation(d,unc,cb);
    respObs = DClab.ResponseObservation(d,unc);
    info.targetKey=obj.targetKey{i};
    info.targetID=obj.targetID{i};
    info.targetLabel=obj.targetLabel{i};
    info.uncRange=obj.uncRange{i};
    info.targetValue=obj.targetValue{i};
    info.targetTrans=obj.targetTrans{i};
    info.targetUnits=obj.targetUnits{i};
    info.targetLinks=obj.targetLinks{i};
    info.trialModelID=obj.trialModelID;
    info.callback2show=obj.callback2show;
    info.callback2show2=obj.callback2show2;
    respObs.userData = info;


    %RM
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
 %   cb = @(rmObj) obj.callback2show(obj.modelID{i});
    respMod = DClab.ResponseModel(modStruct,modelDomain);
    info.modelKey=obj.modelKey{i};
    info.modelID=obj.modelID{i};
    info.modelCoeffs=obj.modelCoeffs{i};
    info.modelParamIDs=obj.modelParamIDs{i};
    info.trialModelID=obj.trialModelID;
    info.callback2show=obj.callback2show;
    info.callback2show2=obj.callback2show2;
    respMod.name = obj.modelKey{i};
    respMod.userData = info;

    %MOP
    mop{i} = DClab.ModelAndObservationPair(respObs,respMod,obj.targetKey{i});
end
mop = vertcat(mop{:});

if isempty(obj.predictionModelCoeffs)
    %targetmodels
    targets = {mop.ResponseModel}';

    %open gui:
    gui6(mop,fp,targets)
    
else
    nPredMods = length(obj.predictionModelCoeffs);
    targets = cell(nPredMods,1);
    for i=1:nPredMods
        info = [];
        info.targetUnits = 'µs';
        info.targetLabel = 't_ign';
        %RM
        modStruct.value = obj.predictionModelCoeffs{i}; %obj.modelCoeffs{i};
        modStruct.responseTransformation = 'none';
        m = length(obj.predictionModelParamIDs{i});  %length(obj.modelParamIDs{i});
        modStruct.variableTransformations = repmat({'none'},m,1);
        modelDomain = struct('name',cell(m,1),'range',cell(m,1));
        
        for j=1:m
            idx = find(strcmp(obj.parameterID,obj.predictionModelParamIDs{i}{j})); %obj.modelParamIDs{i}{j}
            modelDomain(j).name = obj.parameterID{idx};
            modelDomain(j).range = obj.parameterRange{idx};
        end
        respMod = DClab.ResponseModel(modStruct,modelDomain);
        info.modelKey=obj.predictionModelKey{i};
        info.modelID=obj.predictionModelID{i};
        info.modelCoeffs=obj.predictionModelCoeffs{i};
        info.modelParamIDs=obj.predictionModelParamIDs{i};
        %info.trialModelID=obj.trialModelID; %??
        info.callback2show=obj.callback2show;
        info.callback2show2=obj.callback2show2;
        respMod.name = obj.predictionModelKey{i};
        respMod.userData = info;
        targets{i} = respMod;
    end
    gui6(mop,fp,targets,true);
    
end

