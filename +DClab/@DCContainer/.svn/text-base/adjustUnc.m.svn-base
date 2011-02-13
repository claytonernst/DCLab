function obj = adjustUnc(obj)

lenM = length(obj.models);

%TODO: what if stuff is uncType is 'relative'?

for i=1:lenM
    
    switch obj.modelOutputTrans{i}
        
        case 'none'
            
            if strcmp(obj.uncTrans{i},'log10')
                d = obj.data{i};
                unc = obj.unc{i};
                obj.uncTrans{i}='none';
                dtmp = mean([d*10^unc(1),d*10^unc(2)]);
                obj.data{i} = dtmp;
                obj.unc{i} = [d*10^unc(1)-dtmp d*10^unc(2)-dtmp];
            end
            
        case 'log10'
            
            if strcmp(obj.uncTrans{i},'none')
                d = obj.data{i};
                unc = obj.unc{i};
                obj.uncTrans{i}='log10';
                dtmp = mean([log10(unc(1)+d),log10(unc(2)+d)]);
                obj.data{i} = 10^dtmp;
                obj.unc{i} = [log10(unc(1)+d)-dtmp,log10(unc(2)+d)-dtmp];
            end
            
    end
    
end