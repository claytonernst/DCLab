function obj=adjustDomains(obj)


for i=1:length(obj.modelDomain)
    
    for j=1:length(obj.modelDomain{i})
        
        idx = find(strcmp(obj.modelDomain{i}(j).name,obj.parameterNames));
        
        type = obj.parameterType{idx};
        nom = obj.parameterNominal{idx};
        unc = obj.parameterUnc{idx};
        trans = obj.parameterTrans{idx};
        
        if strcmp(type,'absolute') && strcmp(trans,'none')
            obj.modelDomain{i}(j).range = [unc(1)+nom, unc(2)+nom];
        elseif strcmp(type,'absolute') && strcmp(trans,'log10')
            obj.modelDomain{i}(j).range = [nom*exp(unc(1)), nom*exp(unc(2))];
        elseif strcmp(type,'relative') && strcmp(trans,'none')
            obj.modelDomain{i}(j).range = [abs(nom)*unc(1)+nom, abs(nom)*unc(2)+nom];
        elseif strcmp(type,'relative') && strcmp(trans,'log10')
            obj.modelDomain{i}(j).range = [nom*exp(abs(log10(nom))*unc(1)), nom*exp(abs(log10(nom)*unc(2)))];
        else
            error('Invalid parameter type/transformation pair')
        end
    end
end