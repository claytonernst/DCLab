function fitTable(fitObj,type)

% Created 2/25/09 Xiaoqing You, University of California Berkeley
% Last modified 3/1/09 Xiaoqing You, University of California Berkeley
% Copyright 2009 Xiaoqing You. All rights reserved.

ss = get(0,'screensize');
posx = max([0,ss(3)-800])/2;
posy = max([0,ss(4)-600])/2;

%figure
figH = figure('Menubar', 'none', ...
              'NumberTitle', 'off', ...
              'Name', 'Optimization Results',...
              'ToolBar','none',...
              'resize','off',...
              'dockcontrols','off',...
              'position',[posx posy 800 600]); 

PDset = fitObj.PDset;

switch type
  case 'data'
    set(figH,'Name','Optimal model predictions and experimental observations');
    MOP = PDset.ModelAndObservationPair;
    RO = {MOP.ResponseObservation}';
    for i=1:length(RO)
        uD = RO{i}.userData;
        label{i,1} = uD.targetLabel;
        trans{i,1} = uD.targetTrans;
        units{i,1} = uD.targetUnits;
    end    
    eList = {MOP.name}';
    weights = fitObj.weights;
    D = cell2mat({MOP.observedValue}');
    U = cell2mat({MOP.observationUncertaintyPlusMinus}');
    wresid = fitObj.wresid;
    Mx = D + wresid;
    for i=1:length(D)
      T(i,1) = 10.^D(i,1);
      M(i,1) = 10.^Mx(i,1);
      rError(i,1) = 100.*(M(i,1)-T(i,1))./T(i,1);
      dat(i,:) = {char(eList(i)),char(label{i,1}),weights(i),D(i,1),U(i,1)+D(i,1),U(i,2)+D(i,1),Mx(i,1), ...
          num2str(T(i,1),'%0.4g'),num2str(M(i,1),'%0.3g'),num2str(rError(i,1),'%0.2g'),char(units{i,1})};
    end
    cnames = {'Name','Label','Weight','d','Lower bound', 'Upper bound','M(x)','target value','Optimal model prediction','Relative error(%)','units'};
    tb = uitable('Data',dat,'ColumnName',cnames,... 
            'Parent',figH,'Position',[0 0 800 600]);

  case 'param'
    set(figH,'Name','Optimal parameters');
    FP = PDset.FreeParameter;
    for i = 1:length(FP)
      uDstruct = FP(i).userData;
      pList{i,1} = uDstruct.parameterKey;
      pTrans{i,1} = uDstruct.paramTrans;
      centerSpan{i,1} = uDstruct.centerSpan;
      pUnits{i,1} = uDstruct.optVarUnits;
    end
    x0 = cell2mat({FP.nominal}');
    %range = cell2mat({PDset.FreeParameter.range}');
    bestx = fitObj.bestx;
    dx = bestx - x0;
    for i=1:length(x0)
      switch pTrans{i,1}
          case 'center*span^x'
              bestA(i,1) = centerSpan{i,1}(1).*centerSpan{i,1}(2).^bestx(i);
              %Arange1(i,1) = centerSpan{i,1}(1)./centerSpan{i,1}(2);
              %Arange2(i,1) = centerSpan{i,1}(1).*centerSpan{i,1}(2);
              %boundsKind{i,1} = 'relative';
          case 'center+span*x'
              bestA(i,1) = centerSpan{i,1}(1)+centerSpan{i,1}(2).*bestx(i);
              %Arange1(i,1) = centerSpan{i,1}(1)-centerSpan{i,1}(2);
              %Arange2(i,1) = centerSpan{i,1}(1)+centerSpan{i,1}(2);  
              %boundsKind{i,1} = 'absolute';
      end
      dat(i,:) = {char(pList(i,1)),x0(i,1),bestx(i,1),bestA(i,1)};
    end
    cnames = {'Name','x0','x(optimal)','A(optimal)/A0 or Hf298(optimal)-Hf298(0) '};
    tb = uitable('Data',dat,'ColumnName',cnames,... 
            'Parent',figH,'Position',[0 0 800 600]);
end




