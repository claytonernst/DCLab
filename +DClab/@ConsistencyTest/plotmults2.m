function handles = plotmults2(Obj)

%for exporting plots, not ideal on screen display, maximize, then export as
%emf

Dset = Obj.Dataset;
n = Dset.nParameters;
m = Dset.nPairs;

sens = Obj.upperBndSens;

sensalpha = sens.paraml;
sensbeta = sens.paramu;
sensl = sens.expl;
sensu = sens.expu;

axisFontSize = 12;
labelFontSize = 14;

figHand = figure;

plotLambAlphaHand = subplot(2,2,1);
lambalphaHand = bar(sensalpha);

plotLambBetaHand = subplot(2,2,2);
lambbetaHand = bar(sensbeta);

plotLambLHand = subplot(2,2,3);
lamblHand = bar(sensl);

plotLambUHand = subplot(2,2,4);
lambuHand = bar(sensu);

%get and set axes dimensions
xalpha = get(plotLambAlphaHand,'Xlim');
yalpha = get(plotLambAlphaHand,'Ylim');
xbeta = get(plotLambBetaHand,'Xlim');
ybeta = get(plotLambBetaHand,'Ylim');
y = max(abs(yalpha(1)),ybeta(2));
set(plotLambAlphaHand,'Xlim',[xalpha(1)-1 n+1],'Ylim',[-y 0]);
set(plotLambBetaHand,'Xlim',[xbeta(1)-1 n+1],'Ylim',[0 y]);

xl = get(plotLambLHand,'Xlim');
yl = get(plotLambLHand,'Ylim');
xu = get(plotLambUHand,'Xlim');
yu = get(plotLambUHand,'Ylim');
y = max(abs(yl(1)),yu(2));
set(plotLambLHand,'Xlim',[xl(1)-1 m+1],'Ylim',[-y 0]);
set(plotLambUHand,'Xlim',[xu(1)-1 m+1],'Ylim',[0 y]);

%tweak plots
tweakplot('lambalpha',plotLambAlphaHand,axisFontSize,labelFontSize);
tweakplot('lambbeta',plotLambBetaHand,axisFontSize,labelFontSize);
tweakplot('lambl',plotLambLHand,axisFontSize,labelFontSize);
tweakplot('lambu',plotLambUHand,axisFontSize,labelFontSize);

set(plotLambAlphaHand,'Position',[0.178 0.59 0.303 0.3439]);
set(plotLambBetaHand,'Position',[0.6424  0.59 0.3155 0.3439]);
set(plotLambLHand,'Position',[0.178 0.1  0.303 0.3439]);
set(plotLambUHand,'Position',[0.6424 0.1  0.3155 0.3439]);

set(figHand,'position',[1 1 620 370])

% store output data
handles.main = figHand;
handles.alpha.axis = plotLambAlphaHand;
handles.alpha.bars = lambalphaHand;
handles.beta.axis = plotLambBetaHand;
handles.beta.bars = lambbetaHand;
handles.l.axis = plotLambLHand;
handles.l.bars = lamblHand;
handles.u.axis = plotLambUHand;
handles.u.bars = lambuHand;

%note: set XTicks before changing XTickLabel



%========================== local functions ==========================

% Local function to tweak plot appearance
function [] = tweakplot(axisName,axisHandle,axisFontSize,labelFontSize)

switch axisName
  case 'lambalpha'
    xString = 'Parameter number';
    yString = '$-\frac{\partial}{\partial \alpha_i}\overline{\mathrm{C}}_\mathrm{D}$';
  case 'lambbeta'
    xString = 'Parameter number';
%    yString = '\lambda^{(\beta)}';
    yString = '$\frac{\partial}{\partial \beta_i}\overline{\mathrm{C}}_\mathrm{D}$';
 case 'lambl'
    xString = 'Attribute number';
%    yString = '-\lambda^{({\itl})}';
    yString = '$-\frac{\partial}{\partial l_a}\overline{\mathrm{C}}_\mathrm{D}$';
 case 'lambu'
    xString = 'Attribute number';
%    yString = '\lambda^{({\itu})}';
    yString = '$\frac{\partial}{\partial u_a}\overline{\mathrm{C}}_\mathrm{D}$';
 otherwise
    %do nothing or warn user of broken code
end

set(axisHandle,'FontName','times','FontSize',axisFontSize, ...
  'TickLength',[0.03 0.025]);

if length(get( get(axisHandle,'Children'),'Xdata')) > 20
  set(axisHandle,'XMinorTick','on')
end

xlabelHand = get(axisHandle,'xlabel');
ylabelHand = get(axisHandle,'ylabel');
set(xlabelHand,'interpreter','latex','String',xString,'Units','normalized');
set(ylabelHand,'interpreter','latex','String',yString,'Units','normalized');

ppx = [0.5, -0.14, 0];
if strcmp(axisName,'lambalpha')
  ppy = [-0.33,0.38,0];
elseif  strcmp(axisName,'lambbeta')
  ppy = [-0.29,0.38,0];
elseif strcmp(axisName','lambl')
  ppy = [-0.27,0.40,0];
else
  ppy = [-0.27,0.40,0];
end
  
set(xlabelHand,'FontSize',labelFontSize,'FontWeight','normal','FontName','times',...
    'Position',ppx);
set(ylabelHand,'FontSize',labelFontSize,'Rotation',0.001,'FontName','times', ...
               'Position',ppy,'FontWeight','normal');
% The ylabel rotation is set near zero and not zero because of a 
% matlab bug. If you save the figure and reopen it, the 
% rotation is reset to default(90) if it was set to zero.

%to do: check that twenty is about right for when to turn on
%         XMinorTicks
%       get this to work if parameters were not separated
