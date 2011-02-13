function loc = report(ctest,format)
% ConsistencyTest report.
%
% loc = report(ctestObj) creates an HTML report on the results and
% computation of the ConsistencyTest object ctestObj.  This will create
% an m-file in the current directory as well as a subdirectory "html" where
% the html file will be placed. The filesystem location of the generated
% html file will be returned as a string in the output loc.
%
% loc = report(ctestObj,format) allows for various output formats.  format
% is one of the following strings
%    'html' (default): creates HTML file to viewed with a web browser
%    'doc':  creates a Microsoft Word document.  Word must be installed on
%            your machine for this option
%    'ppt':  creates a Microsoft Powerpoint document.  Powerpoint much be
%            installed on your machine for this option
%    'xml':  creates an XML file
%    'latex': creates a latex file
% With every option the subdirectory created is still called "html"

if nargin<2
    format = 'html';
end

%create file
c = clock;
d = date;
d(strfind(d,'-')) = [];
tag = sprintf('%s%d%d',d,c(4:5));
filename = ['ConsistencyTestReport' tag '.m'];
fid = fopen(filename,'w');
if ~exist('html','dir')
    mkdir html
end

%save which are old figs so we can delete new ones later
oldFigs = get(0,'child');

%Title
fprintf(fid,'%%%% ConsistencyTest Report\n');

%names
dname = ctest.DatasetName;
if isempty(dname)
    dname = 'Not specified';
end
fprintf(fid,'%% Using Dataset: %s\n\n',dname);

%consistency test range
LB = ctest.LB;
UB = ctest.UB;
if isnan(LB)
    LB = 'Not calculated';
else
    LB = sprintf('%0.4g',LB);
end
if isnan(UB)
    UBo = 'Not calculated';
else
    UB = sprintf('%0.4g',UB);
end
fprintf(fid,'%%%% Consistency Measure Bounds: [%s, %s]\n',LB,UB);
if ctest.LB>0
    fprintf(fid,'%% The dataset is consistent.\n');
elseif ctest.UB<0
    fprintf(fid,'%% The dataset is inconsistent.\n');
else
    fprintf(fid,'%% The consistency test was inconclusive.  Haven''t proven consistency one way or the other.  Try increasing the number of iterations.\n')
end

%sensitivities
sens = ctest.upperBndSens;
figure;
%upperboundsens
a1 = subplot(2,2,1);
bar(sens.paramu);
ylabel('Sens. to u.b.')
a2 = subplot(2,2,2);
bar(sens.expu);
a3 = subplot(2,2,3);
bar(sens.paraml);
ylabel('Sens. to l.b.')
xlabel('Parameter')
a4 = subplot(2,2,4);
bar(sens.expl);
xlabel('Experiment')
M = 1.1*max(abs([sens.paramu; sens.paraml; sens.expu; sens.expl]));
m = ctest.iter(end).PDset.nPairs;
n = ctest.iter(end).PDset.nParameters;
set(a1,'ylim',[0 M],'xlim',[0.5 n+0.5],'xticklabel',[])
set(a2,'ylim',[0 M],'xlim',[0.5 m+0.5],'xticklabel',[],'yticklabel',[])
set(a3,'ylim',[-M 0],'xlim',[0.5 n+0.5])
set(a4,'ylim',[-M 0],'xlim',[0.5 m+0.5],'yticklabel',[])
LOCALmassageFig(gcf,5,[],2,2,1);
pos1 = get(a1,'position');
pos2 = get(a2,'position');
pos3 = get(a3,'position');
pos4 = get(a4,'position');
a3w = n*(pos3(3)+pos4(3))/(n+m);
a4w = m*(pos3(3)+pos4(3))/(n+m);
pos4(1) = pos4(1)-(pos3(3)-a3w);
pos3(3) = a3w;
pos4(3) = a4w;
pos1(3) = pos3(3);
pos2(3) = pos4(3);
pos2(1) = pos4(1);
set(a1,'position',pos1)
set(a2,'position',pos2)
set(a3,'position',pos3)
set(a4,'position',pos4)
print(['html/upperbndSens-' tag],'-dpng','-r96')


fprintf(fid,'%%%% Sensitivities\n');
fprintf(fid,'%%\n%% Sensitivities to the Upper bound of the consistency test:\n');
fprintf(fid,'%%\n%% <<upperbndSens-%s.png>>\n%%\n',tag);




%b&b
nIter = length(ctest.iter);
opt = ctest.iter(end).optimOpts;
UB = zeros(nIter,1);
LB = zeros(nIter,1);
for i=1:nIter
    if opt.omitInnerBound
        LB(i) = min(ctest.iter(i).upper.bnd);
    else
        LB(i) = min(ctest.iter(i).lower.bnd);
    end
    if opt.omitOuterBound
        UB(i) = max(ctest.iter(i).lower.bnd);
    else
        UB(i) = max(ctest.iter(i).upper.bnd);
    end
end
figure
if ~opt.omitOuterBound
    plot(1:nIter,UB,'r--')
    hold on
end
if ~opt.omitInnerBound
    plot(1:nIter,LB,'g--')
end
LOCALfuzzyInterval([LB LB UB UB],[],[0.5 4]);
set(gca,'xtick',1:nIter)
ylabel('Consistency Measure')
xlabel('Iteration')
hold off
LOCALmassageFig(gcf,min([5 nIter]),min([5 nIter])*10/nIter,1,1,1)
legend('Upper bound','Lower bound')
print(['html/iteration-' tag],'-dpng','-r96')

if opt.omitInnerBound || opt.omitOuterBound
    stopCrit = 'Need both inner and outer bounds for iteration.  See OuterboundSubDivide object.';
elseif nIter<opt.maxBranchBoundIter
    stopCrit = 'Branch and bound tolerance met';
else
    stopCrit = 'Reached maximum branch and bound iterations';
end
totTime = num2str(sum([ctest.iter.runtime]));

fprintf(fid,'%%%% Branch and Bound Iterations\n');
fprintf(fid,'%% General information:\n');
fprintf(fid,'%% \n%% * Total iterations: %d\n%% * Stopping criterion: %s\n%% * Total computation time: %s\n%% \n',nIter,stopCrit,totTime);
if nIter>1
    fprintf(fid,'%% Consistency test bounds during iteration:\n');
    fprintf(fid,'%% \n%% <<iteration-%s.png>> \n%% \n',tag);
end



%Fits on final iteration
fprintf(fid,'%%%% Final Iteration Surrogate Fits\n');
PDset = ctest.iter(end).PDset;
leaves = PDset.leafNodes;
for i=1:m
    
    fprintf(fid,'%% \n%% *%s*\n',PDset.ModelAndObservationPair(i).name);
    if strcmp(PDset.ModelAndObservationPair(i).ResponseModel.type,'linearModel')
        fprintf(fid,'%% \n%% * Linear model: no fitting error\n%% \n');
    elseif strcmp(PDset.ModelAndObservationPair(i).ResponseModel.type,'quadraticModel')
        fprintf(fid,'%% \n%% * Quadratic model: no fitting error\n%% \n');
    else
        
        idx = cell(nIter,1);
        nPnts4Fits = 0;
        nPnts4Valid = 0;
        avgErr = zeros(nIter,1);
        peakVErr = zeros(nIter,1);
        peakOptErr = zeros(nIter,1);
        for j=1:nIter
            surf = PDset.PiecewiseSurrogateModelTree(leaves(j)).DCSurface(i);
            idx{i} = surf.activeParameterIndex;
            nPnts4Fits = nPnts4Fits + surf.surrogateFitInfo.nSamplePoints4Fit;
            nPnts4Valid = nPnts4Valid + surf.surrogateFitInfo.nSamplePoints4Validation;
            avgErr(i) = surf.surrogateFitInfo.averageErrorOnSamplePoints4Validation;
            peakVErr(i) = max(abs(surf.surrogateFitInfo.peakErrorOnSamplePoints4Validation));
            peakOptErr(i) = max(abs(surf.surrogateFitInfo.peakErrorFromOptimization));
        end
        fprintf(fid,'%% \n%% * %d of %d parameters are active\n',length(unique(vertcat(idx{:}))),length(PDset.ModelAndObservationPair(i).ResponseModel.parameterList));
        fprintf(fid,'%% * Num. evals for fit: %d\n',nPnts4Fits);
        fprintf(fid,'%% * Num. evals for validation: %d\n',nPnts4Valid);
        fprintf(fid,'%% * Average validation error: %0.2g\n',mean(avgErr));
        fprintf(fid,'%% * Peak validation error: %0.2g\n',max(peakVErr));
        fprintf(fid,'%% * Peak error from local search: %0.2g\n%% \n',max(peakOptErr));
        
        if nIter>1
            avgErrPlot = zeros(nIter,1);
            peakVErrPlot = zeros(nIter,1);
            peakOptErrPlot = zeros(nIter,1);
            for j=1:nIter
                leaves = ctest.iter(j).leafNodes;
                
                nSurf = length(leaves);
                avgErr = zeros(nSurf,1);
                peakVErr = zeros(nSurf,1);
                peakOptErr = zeros(nSurf,1);
                for k=1:length(leaves)
                    surf = PDset.PiecewiseSurrogateModelTree(leaves(k)).DCSurface(i);
                    avgErr(i) = surf.surrogateFitInfo.averageErrorOnSamplePoints4Validation;
                    peakVErr(i) = max(abs(surf.surrogateFitInfo.peakErrorOnSamplePoints4Validation));
                    peakOptErr(i) = max(abs(surf.surrogateFitInfo.peakErrorFromOptimization));
                end
                avgErrPlot(j) = mean(avgErr);
                peakVErrPlot(j) = max(peakVErr);
                peakOptErrPlot(j) = max(peakOptErr);
            end
            figH = figure;
            plot([avgErrPlot peakVErrPlot peakOptErrPlot])
            xlabel('Iteration')
            ylabel('Fit Error')
            LOCALmassageFig(gcf,5,[],1,1,1);
            legend('Average Valid. Error','Peak Valid. Error','Peak Optimized Error') 
            print(['html/Fiterror_' num2str(i) '-' tag],'-dpng','-r96');
            delete(figH);
            fprintf(fid,'%% Fit Error v. Iteration\n');
            fprintf(fid,'%% \n%% <<Fiterror_%d-%s.png>>\n%% \n',i,tag);
            
        end
                
        
    end
    
    
end
    



fclose(fid);

newFigs = setdiff(get(0,'child'),oldFigs);
delete(newFigs)

options.format = format;
loc = publish(filename,options);





function LOCALmassageFig(figH,w,h,nr,nc,tightFlag)
%massageFig(figH,w,h,nr,nc,tightFlag)

ni = nargin;
phi = (sqrt(5)-1)/2;
if ni==2
    h = phi*w;
elseif ni==3
    if isempty(h) && isempty(w)
        w = 5;
        h = 5*phi;
    elseif isempty(h)
        h = w*phi;
    elseif isempty(w)
        w = h/phi;
    end
end
if ni>=4
    if ni<5
        nc=[];
    end
    childs = get(figH,'child');
    axesHansTmp = childs(strcmp(get(get(figH,'child'),'type'),'axes'));
    nAxes = length(axesHansTmp);
    if isempty(nr) && isempty(nc)
        nr=sqrt(nAxes);
        nc=nr;
    elseif isempty(nr)
        nr = nAxes/nc;
    else
        nc = nAxes/nr;
    end
    assert(round(nr)==nr && round(nc)==nc,'number of axes doesn''t match number of rows and columns')
    if isempty(h) && isempty(w)
        w = 5;
        h = 5*nr*phi/nc;
    elseif isempty(h)
        h = w*nr*phi/nc;
    elseif isempty(w)
        w = h*nc/phi/nr;
    end
end
if ni<6
    tightFlag = false;
end


%set figure units to inches and change size
set(figH,'units','inches','position',[0.5 0.5 w h]);

if tightFlag
    
    %axesHansTmp = get(figH,'child');
    nAxes = length(axesHansTmp);
    axesHans = zeros(nAxes,1);
    for k=1:nAxes  %change the order?
        p = get(axesHansTmp(k),'position');
        j=ceil((nc+0.5)*p(1));
        i=nr-ceil((nr+0.5)*p(2))+1;
        axesHans((i-1)*nc+j) = axesHansTmp(k);
    end
    
    %let's get some measurements
    ow = zeros(nr,2);
    bh = zeros(2*nr,1);
    for i=1:nr
        ti = get(axesHans((i-1)*nc+1),'tightinset');
        ti2 = get(axesHans(i*nc),'tightinset');
        ow(i,1) = ti(1);
        ow(i,2) = ti2(3);
        bh(2*i-1) = ti(4);
        bh(2*i) = ti(2);
    end
    owL = max(ow(:,1));
    owR = max(ow(:,2));
    bh = 2*max(bh(1:end-1));
    oh = zeros(nc,2);
    bw = zeros(2*nc,1);
    for i=1:nc
        ti = get(axesHans((nr-1)*nc+i),'tightinset');
        ti2 = get(axesHans(i),'tightinset');
        oh(i,1) = ti(2);
        oh(i,2) = ti2(4);
        bw(2*i-1) = ti(1);
        bw(2*i) = ti(3);
    end
    ohB = max(oh(:,1));
    ohT = max(oh(:,2));
    bw = 2*max(bw(2:end));
    
    
    %axes width/height
    ah = (1-ohT-ohB-(nr-1)*bh)/nr;
    aw = (1-owL-owR-(nc-1)*bw)/nc;
    
    %position arrays
    p1 = zeros(nc,1);
    for i=1:nc
        p1(i) = owL + (i-1)*(aw+bw);
    end
    p2 = zeros(nr,1);
    for i=1:nr
        p2(i) = ohB + (nr-i)*(ah+bh);
    end
    
    %fix sizes
    for i=1:nc
        for j=1:nr
            set(axesHans((j-1)*nc+i),'position',[p1(i) p2(j) aw ah])
        end
    end
    
end
        
function handles = LOCALfuzzyInterval(varargin)
%
% handles = fuzzyInterval(dat,loc,w)
% handles = fuzzyInterval(ahan,dat,loc,w)
%
% ahan: axes handle
% dat: n-by-4 double, [LBo LBi UBi UBo]
% loc: n-by-1 double, x-axis center points of bars
% w: width of bars

ni = nargin;
error(nargchk(1,4,ni));


switch ni
    case 1
        if ishandle(varargin{1})
            error('Data must be provided')
        end
        ahan = gca;
        dat = varargin{1};
        loc = [];
        w = [];
    case 2
        if ishandle(varargin{1})
            ahan = varargin{1};
            dat = varargin{2};
            loc = [];
        else
            ahan = gca;
            dat = varargin{1};
            loc = varargin{2};
        end
        w = [];
    case 3
        if ishandle(varargin{1})
            ahan = varargin{1};
            dat = varargin{2};
            loc = varargin{3};
            w = [];
        else
            ahan = gca;
            dat = varargin{1};
            loc = varargin{2};
            w = varargin{3};
        end
    case 4
        if ishandle(varargin{1})
            ahan = varargin{1};
            dat = varargin{2};
            loc = varargin{3};
            w = varargin{4};
        else
            error('Improper inputs')
        end
end
            
      
sz = size(dat);
if sz(2)~=4
    error('Size of dat needs to be n-by-4')
end
if isempty(loc)
    loc=1:sz(1);
end
if length(loc)~= sz(1)
    error('Location vector needs to have the same number of points')
end
if ~any(size(loc)==1)
    error('Location vector needs to be a vector')
end
if isempty(w)
    w = [0.75 4];
end
if all(size(w)==1)
    w = [w 2];
elseif any(size(w)~=[1 2])
    error('Improper width.  Should be scalar or 1-by-2')
end
if w(1)>1 || w(1)<=0
    error('w should be be in the interval [0,1)')
end


handles = struct('patches',cell(sz(1),1),'line',cell(sz(1),1));
axes(ahan)
np = get(ahan,'nextplot');
for i1=1:sz(1)
    
    handles(i1).patches.low = patch(loc(i1)+w(1)*[-1 1 1 -1]/2,dat(i1,[1 1 2 2]),'k','parent',ahan);
    hold on
    handles(i1).patches.high = patch(loc(i1)+w(1)*[-1 1 1 -1]/2,dat(i1,[3 3 4 4]),'k','parent',ahan);
    handles(i1).line = plot(ahan,loc([i1 i1]),dat(i1,[2 3]),'k','linewidth',w(2));
    
end
set(ahan,'nextplot',np)