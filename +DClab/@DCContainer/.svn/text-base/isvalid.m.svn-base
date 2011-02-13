function bool = isvalid(obj)

%right now I'm only really checking lengths

bool = true;

% obj.parameterNames = {};
% obj.parameterNominal = {};
% obj.parameterUnc = {};
% obj.parameterType = {};
% obj.parameterTrans = {};
% obj.modelNames = {};
% obj.models = {};
% obj.modelDomain = {};
% obj.modelOutputTrans = {};
% obj.modelInputTrans = {};
% obj.data = {};
% obj.unc = {};
% obj.uncType = {};
% obj.uncTrans = {};
% obj.targetFlag = [];

%all fields should be vector cell arrays
for i=fieldnames(obj)';
    if ~iscell(obj.(i{1})) || ~isvector(obj.(i{1}))
        if ~strcmp(i,'targetFlag')
            bool = false;
            warning('DC:invalid','all fields (except targetFlag) of DCContainer of vector cell arrays')
        end
    end
end

%parameter fields
pLen = length(obj.parameterNames);
if pLen~=length(obj.parameterNominal)
    bool = false;
    warning('DC:invalid','parameterNames and parameterNominal are different lengths');
end
if pLen~=length(obj.parameterUnc)
    bool = false;
    warning('DC:invalid','parameterName and parameterUnc are different lengths');
end
if pLen~=length(obj.parameterType)
    bool = false;
    warning('DC:invalid','parameterName and parameterType are different lengths');
end
if pLen~=length(obj.parameterTrans)
    bool = false;
    warning('DC:invalid','parameterName and parameterTrans are different lengths');
end

%models and experiments and targetFlag
mLen = length(obj.models);
if mLen~=length(obj.modelNames)
    bool = false;
    warning('DC:invalid','models and modelNames are different lengths');
end
if mLen~=length(obj.modelDomain)
    bool = false;
    warning('DC:invalid','models and modelDomain are different lengths');
end
if mLen~=length(obj.modelOutputTrans)
    bool = false;
    warning('DC:invalid','models and modelOutputTrans are different lengths');
end
if mLen~=length(obj.modelInputTrans)
    bool = false;
    warning('DC:invalid','models and modelInputTrans are different lengths');
end
if mLen~=length(obj.data)
    bool = false;
    warning('DC:invalid','models and data are different lengths');
end
if mLen~=length(obj.unc)
    bool = false;
    warning('DC:invalid','models and unc are different lengths');
end
if mLen~=length(obj.uncType)
    bool = false;
    warning('DC:invalid','models and uncType are different lengths');
end
if mLen~=length(obj.uncTrans)
    bool = false;
    warning('DC:invalid','models and uncTrans are different lengths');
end
if mLen~=length(obj.targetFlag)
    bool = false;
    warning('DC:invalid','models and targetsFlag are different lengths');
end