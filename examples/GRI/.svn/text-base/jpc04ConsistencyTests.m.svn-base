%set selector to 0 to generate and plot data, set it to 1 to load
%data from mat files and plot data

selector = 0;
switch selector
 case 0
  %load original dataset
  D = createGRI('original');
  opt = DCOptions('nRestart',0);
  
  Pairs = D.ModelAndObservationPair;
  %change all uncertainties to 0.08 from 0.1
  for i1 = 1:Pairs.nPairs
    Pairs(i1).ResponseObservation.uncertainty = 0.08;
  end
  D.ModelAndObservationPair = Pairs;

  nominal_test = ConsistencyTest(D,opt);

  %change data of experiment 57 to reflect the updated value
  Dupdate1 = D;
  Dupdate1.ModelAndObservationPair(57).ResponseObservation.observedValue = log10(700);

  oh1a_test = ConsistencyTest(Dupdate1,opt);

  %change data of experiment 58 to reflect the updated value
  Dupdate2 = D;
  Dupdate2.ModelAndObservationPair(58).ResponseObservation.observedValue = log10(255);

  oh1b_test = ConsistencyTest(Dupdate2,opt);

  %change data of both experiment 57 and 58 to reflect update.
  Dboth = Dupdate1;
  Dboth.ModelAndObservationPair(58).ResponseObservation.observedValue = log10(255);

  both_test = ConsistencyTest(Dboth,opt);

  %opt = DCOptions('nRestart',29,'display','all');
  %both_supertest = ConsistencyTest(Dboth,opt);

  save nominal_test nominal_test
  save oh1a_test oh1a_test
  save oh1b_test oh1b_test
  save both_test both_test
  %save both_supertest both_supertest
 case 1
  load nominal_test
  load oh1a_test
  load oh1b_test 
  load both_test
 otherwise
  error('  Selector must be 0 or 1')
end

%generate plots
h1 = plotmults2(nominal_test);
h2 = plotmults2(oh1a_test);
h3 = plotmults2(oh1b_test);
h4 = plotmults2(both_test);

%tweak plots layouts and axis labels

topheight = [-0.01 0.18]; %yaxis limits for the plots
bottomheight = [-1 8];
%barcolor = [0.5 0.5 0.5]; %color of the bars 
barcolor = 'flat';

topticks = [0 20 40 60 80 100]; %tick labels to top two plots

%nominal plot layout
set(h1.alpha.axis,'Ylim',topheight,'Xtick',topticks);
set(h1.alpha.bars,'facecolor',barcolor);
set(h1.beta.axis,'Ylim',topheight,'Xtick',topticks);
set(h1.beta.bars,'facecolor',barcolor);
set(h1.l.axis,'Ylim',bottomheight);
set(h1.l.bars,'facecolor',barcolor);
set(h1.u.axis,'Ylim',bottomheight);
set(h1.u.bars,'facecolor',barcolor);

%update to 57 layout
set(h2.alpha.axis,'Ylim',topheight,'Xtick',topticks);
set(h2.alpha.bars,'facecolor',barcolor);
set(h2.beta.axis,'Ylim',topheight,'Xtick',topticks);
set(h2.beta.bars,'facecolor',barcolor);
set(h2.l.axis,'Ylim',bottomheight);
set(h2.l.bars,'facecolor',barcolor);
set(h2.u.axis,'Ylim',bottomheight);
set(h2.u.bars,'facecolor',barcolor);

%update to 58 layout
set(h3.alpha.axis,'Ylim',topheight,'Xtick',topticks);
set(h3.alpha.bars,'facecolor',barcolor);
set(h3.beta.axis,'Ylim',topheight,'Xtick',topticks);
set(h3.beta.bars,'facecolor',barcolor);
set(h3.l.axis,'Ylim',bottomheight);
set(h3.l.bars,'facecolor',barcolor);
set(h3.u.axis,'Ylim',bottomheight);
set(h3.u.bars,'facecolor',barcolor);

%update to both 57 and 58 layout
set(h4.alpha.axis,'Ylim',topheight,'Xtick',topticks);
set(h4.alpha.bars,'facecolor',barcolor);
set(h4.beta.axis,'Ylim',topheight,'Xtick',topticks);
set(h4.beta.bars,'facecolor',barcolor);
set(h4.l.axis,'Ylim',bottomheight);
set(h4.l.bars,'facecolor',barcolor);
set(h4.u.axis,'Ylim',bottomheight);
set(h4.u.bars,'facecolor',barcolor);

topoffset = 0.015; %distance about top of bars to add labels on
                    %the top two plots

%distance about top of bars to add labels on the bottom two plots
bottomoffset = topoffset*(bottomheight(2) ...
                          - bottomheight(1))/(topheight(2) - topheight(1));

%horizontal offset to place labels above bars on the top two plots
%so the numbers are centered.
topslide = 2.5;

%horizontal offset to place labels above bars on the bottom two plots
%so the numbers are centered.
bottomslide = 2;

%cutoffs above which we add a label for the nominal plot
alpha1cut = 0.024;
beta1cut = 0.024;
l1cut = 0.2;
u1cut = 0.2;

%cutoffs above which we add a label for the updated 57 plot
alpha2cut = 0.032;
beta2cut = 0.032;
l2cut = 0.3;
u2cut = 0.3;

%cutoffs above which we add a label for the updated 58 plot
alpha3cut = 0.0198;
beta3cut = 0.0198;
l3cut = 0.215;
u3cut = 0.215;

%cutoffs above which we add a label for the both updated plot
alpha4cut = 0.0206;
beta4cut = 0.0206;
l4cut = 0.15;
u4cut = 0.15;

labsz = 10;

% we still need to tweak the labels some for a nice display. Remove
% labels that are too close together and change the horizontal
% offset of labels that contain only one digit. We do this and
% create the labels here

% for the nominal plot
sens = get(h1.alpha.bars,'YData');
idx = find(sens > alpha1cut);
axes(h1.alpha.axis);
for i1 = idx
  th = text(i1-topslide,sens(i1)+topoffset,num2str(i1));
  set(th,'fontweight','bold','fontname','times','fontsize',labsz);
end

sens = get(h1.beta.bars,'YData');
idx = find(sens > beta1cut);
axes(h1.beta.axis);
for i1 = idx
  if i1==18
    th = text(i1-topslide+2.5,sens(i1)+topoffset,num2str(i1));
  else
    th = text(i1-topslide,sens(i1)+topoffset,num2str(i1));
  end
  set(th,'fontweight','bold','fontname','times','fontsize',labsz);
end

sens = get(h1.l.bars,'YData');
idx = find(sens > l1cut);
axes(h1.l.axis);
for i1 = idx
  th = text(i1-bottomslide,sens(i1)+bottomoffset,num2str(i1));
  set(th,'fontweight','bold','fontname','times','fontsize',labsz);
end

sens = get(h1.u.bars,'YData');
idx = find(sens > u1cut);
axes(h1.u.axis);
for i1 = idx
  th = text(i1-bottomslide,sens(i1)+bottomoffset,num2str(i1));
  set(th,'fontweight','bold','fontname','times','fontsize',labsz);
end

% for the updated oh1a plot
sens = get(h2.alpha.bars,'YData');
idx = find(sens > alpha2cut);
idx = setdiff(idx,59); %dump, crowding
axes(h2.alpha.axis);
for i1 = idx
  if i1 == 57
    th = text(i1-topslide-3,sens(i1)+topoffset,num2str(i1));
  elseif i1 == 3
    th = text(i1-topslide+1.5,sens(i1)+topoffset,num2str(i1));
  else
    th = text(i1-topslide,sens(i1)+topoffset,num2str(i1));
  end
  set(th,'fontweight','bold','fontname','times','fontsize',labsz);
end

sens = get(h2.beta.bars,'YData');
idx = find(sens > beta2cut);
idx = setdiff(idx,[2 41]);
axes(h2.beta.axis);
for i1 = idx
  if i1 == 1
    th = text(i1-topslide+2,sens(i1)+topoffset,num2str(i1));
  elseif i1 == 49
    th = text(i1-topslide-2,sens(i1)+topoffset,num2str(i1));
  else
    th = text(i1-topslide,sens(i1)+topoffset,num2str(i1));
  end
  set(th,'fontweight','bold','fontname','times','fontsize',labsz);
end

sens = get(h2.l.bars,'YData');
idx = find(sens > l2cut);
axes(h2.l.axis);
for i1 = idx
  th = text(i1-bottomslide,sens(i1)+bottomoffset,num2str(i1));
  set(th,'fontweight','bold','fontname','times','fontsize',labsz);
end

sens = get(h2.u.bars,'YData');
idx = find(sens > u2cut);
axes(h2.u.axis);
for i1 = idx
  th = text(i1-bottomslide,sens(i1)+bottomoffset,num2str(i1));
  set(th,'fontweight','bold','fontname','times','fontsize',labsz);
end

% for the updated oh1b plot
sens = get(h3.alpha.bars,'YData');
idx = find(sens > alpha3cut);
axes(h3.alpha.axis);
for i1 = idx
  th = text(i1-topslide,sens(i1)+topoffset,num2str(i1));
  set(th,'fontweight','bold','fontname','times','fontsize',labsz);
end

sens = get(h3.beta.bars,'YData');
idx = find(sens > beta3cut);
idx = setdiff(idx,50);
axes(h3.beta.axis);
for i1 = idx
  th = text(i1-topslide,sens(i1)+topoffset,num2str(i1));
  set(th,'fontweight','bold','fontname','times','fontsize',labsz);
end

sens = get(h3.l.bars,'YData');
idx = find(sens > l3cut);
axes(h3.l.axis);
for i1 = idx
  th = text(i1-bottomslide,sens(i1)+bottomoffset,num2str(i1));
  set(th,'fontweight','bold','fontname','times','fontsize',labsz);
end

sens = get(h3.u.bars,'YData');
idx = find(sens > u3cut);
axes(h3.u.axis);
for i1 = idx
  if i1==2
    th = text(i1-bottomslide+1,sens(i1)+bottomoffset,num2str(i1));
  else
    th = text(i1-bottomslide,sens(i1)+bottomoffset,num2str(i1));
  end
  set(th,'fontweight','bold','fontname','times','fontsize',labsz);
end

% for the both updated plot
sens = get(h4.alpha.bars,'YData');
idx = find(sens > alpha4cut);
axes(h4.alpha.axis);
for i1 = idx
  if i1==57 
    th = text(i1-topslide-3.5,sens(i1)+topoffset,num2str(i1));
  else
    th = text(i1-topslide,sens(i1)+topoffset,num2str(i1));
  end
  set(th,'fontweight','bold','fontname','times','fontsize',labsz);
end

sens = get(h4.beta.bars,'YData');
idx = find(sens > beta4cut);
idx = setdiff(idx,50);
axes(h4.beta.axis);
for i1 = idx
  if i1==34
    th = text(i1-topslide-1.5,sens(i1)+topoffset,num2str(i1));
  elseif i1 == 1
    th = text(i1-topslide+1.5,sens(i1)+topoffset,num2str(i1));
  else
    th = text(i1-topslide,sens(i1)+topoffset,num2str(i1));
  end
  set(th,'fontweight','bold','fontname','times','fontsize',labsz);
end

sens = get(h4.l.bars,'YData');
idx = find(sens > l4cut);
axes(h4.l.axis);
for i1 = idx
  th = text(i1-bottomslide,sens(i1)+bottomoffset,num2str(i1));
  set(th,'fontweight','bold','fontname','times','fontsize',labsz);
end

sens = get(h4.u.bars,'YData');
idx = find(sens > u4cut);
axes(h4.u.axis);
for i1 = idx
  if i1 == 2
    th = text(i1-bottomslide+1,sens(i1)+bottomoffset,num2str(i1));
  else
    th = text(i1-bottomslide,sens(i1)+bottomoffset,num2str(i1));
  end
  set(th,'fontweight','bold','fontname','times','fontsize',labsz);
end
