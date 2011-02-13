function objOut = guiLite(obj,modal)
%
%open DCOptions GUI.  Build/edit a portion of DCOptions
%
% One exception, the object passed out of the function has the same
% guiHandle property as the object passed into it, even if the user resorts
% to defaults. This facilitates interfacing with the main gui.

nout = nargout;

ss = get(0,'screensize');
posx = max([0,ss(3)-670])/2;
posy = max([0,ss(4)-450])/2;

%figure
figH = figure('MenuBar','none',...
              'NumberTitle','off',...
              'name','DCOptions',...
              'ToolBar','none',...
              'resize','off',...
              'dockcontrols','off',...
              'closerequestfcn','uiresume(gcf)',...
              'position',[posx posy 600 260]);
if nargin==2 && strcmp(modal,'modal')
    set(figH,'WindowStyle','modal')
end

%Buttons
uicontrol(figH,'style','pushbutton',...
    'string','Resort to Defaults',...
    'callback',@local_defaults,...
    'position',[15 10 125 20],...
    'selected','off');
uicontrol(figH,'style','pushbutton',...
    'string','Help',...
    'callback','doc(''dcopt'')',...
    'position',[155 10 50 20],...
    'selected','off');
uicontrol(figH,'style','pushbutton',...
    'string','OK',...
    'callback','uiresume(gcf)',...
    'position',[470 10 40 20 ],...
    'selected','off');
uicontrol(figH,'style','pushbutton',...
    'string','Cancel',...
    'callback','setappdata(gcf,''opt'',getappdata(gcf,''original''));uiresume(gcf)',...
    'position',[525 10 60 20],...
    'selected','off');

%title
uicontrol(figH,'style','text',...
    'string','Options',...
    'fontsize',16,...
    'horizontalalignment','left',...
    'backgroundcolor',[0.8 0.8 0.8],...
    'position',[20 225 80 25]);

%panel
mainpan = uipanel(figH,'title','',...
    'units','pixels',...
    'position',[10 40 580 175]);



%left column
omitInnerBound = uicontrol(mainpan,'style','checkbox',...
    'string',' Calculate Inner Bound(s)',...
    'horizontalalignment','left',...
    'position',[20 130 200 20],...
    'callback',{@local_checkbox_inverse,'omitInnerBound'});
omitOuterBound = uicontrol(mainpan,'style','checkbox',...
    'horizontalalignment','left',...
    'string',' Calculate Outer Bound(s)',...
    'position',[20 92 200 20],...
    'callback',{@local_checkbox_inverse,'omitOuterBound'});
uicontrol(mainpan,'style','text',...
    'horizontalalignment','left',...
    'position',[20,50,250,20],...
    'string','Maximum Branch & Bound iterations:');
maxBranchBoundIter = uicontrol(mainpan,'style','edit',...
    'string',' Perform a Branch and Bound Algorithm',...
    'position',[260 55 40 20],...
    'callback',{@local_generic_set,'maxBranchBoundIter'});
uicontrol(mainpan,'style','text',...
    'string','Branch and Bound Convergence Tolerence:',...
    'horizontalalignment','left',...
    'position',[20 7 180 30]);
branchBoundTermTol = uicontrol(mainpan,'style','edit',...
    'position',[180 10 80 20],...
    'callback',{@local_generic_set,'branchBoundTermTol'});

%right column
uicontrol(mainpan,'style','text',...
    'string','Command Line Display:',...
    'Horizontalalignment','left',...
    'position',[320 130 220 20]);
display = uicontrol(mainpan,'style','popupmenu',...
    'position',[480 135 70 20],...
    'string',{' off',' final',' notify',' iter',' all',' ALL'},...
    'callback',{@local_pulldowns,'display'});
uicontrol(mainpan,'style','text',...
    'string','Objective Function Convergence Tolerence:',...
    'horizontalalignment','left',...
    'position',[320 90 180 30]);
tolFun = uicontrol(mainpan,'style','edit',...
    'position',[480 88 80 20],...
    'callback',{@local_generic_set,'tolFun'});
uicontrol(mainpan,'style','text',...
    'string','Number of restarts:',...
    'horizontalalignment','left',...
    'position',[320 50 150 20]);
nRestart = uicontrol(mainpan,'style','edit',...
    'position',[480 53 40 20],...
    'callback',{@local_generic_set,'nRestart'});
uicontrol(mainpan,'style','text',...
    'string','Parameter Optimization Objective Norm:',...
    'horizontalalignment','left',...
    'position',[320 7 150 30]);
paramOptimObjectiveFctnNorm = uicontrol(mainpan,'style','popupmenu',...
    'position',[480 10 50 20],...
    'string',{'two','one','inf'},...
    'callback',{@local_pulldowns,'paramOptimObjectiveFctnNorm'});





%--------save handles-------------
handles.fig = figH;
handles.display = display;
handles.nRestart = nRestart;
handles.omitInnerBound = omitInnerBound;
handles.omitOuterBound = omitOuterBound;
handles.maxBranchBoundIter = maxBranchBoundIter;
handles.tolFun = tolFun;
handles.branchBoundTermTol = branchBoundTermTol;
handles.paramOptimObjectiveFctnNorm = paramOptimObjectiveFctnNorm;

setappdata(figH,'handles',handles);
setappdata(figH,'nout',nout);
setappdata(figH,'isSaved',0);
setappdata(figH,'original',obj);

local_setForm([],[],obj);


uiwait(figH);
objOut = getappdata(figH,'opt');
objOut.guiHandle = obj.guiHandle; %DCOptions/guiLite should not modify this property
delete(figH);


%-------------------------------------------
%-----------import data---------------------
function local_setForm(obj,eventD,opt)

handles = getappdata(gcf,'handles');

%display
switch opt.display
    case 'off'
        set(handles.display,'value',1);
    case 'final'
        set(handles.display,'value',2);
    case 'notify'
        set(handles.display,'value',3);
    case 'iter'
        set(handles.display,'value',4);
    case 'all'
        set(handles.display,'value',5);
    case 'ALL'
        set(handles.display,'value',6);
end

%omitInnerBound
set(handles.omitInnerBound,'value',~opt.omitInnerBound);
%omitOuterBound
set(handles.omitOuterBound,'value',~opt.omitOuterBound);
%maxBranchBoundIter
set(handles.maxBranchBoundIter,'string',num2str(opt.maxBranchBoundIter));

%nRestart
set(handles.nRestart,'string',num2str(opt.nRestart));

%magObjectFunTol
set(handles.tolFun,'string',num2str(opt.tolFun));

%branchBoundTermTol
set(handles.branchBoundTermTol,'string',num2str(opt.branchBoundTermTol));

%paramOptimObjectiveFctnNorm
switch opt.paramOptimObjectiveFctnNorm
    case 'two'
        set(handles.paramOptimObjectiveFctnNorm,'value',1);
    case 'one'
        set(handles.paramOptimObjectiveFctnNorm,'value',2);
    case 'inf'
        set(handles.paramOptimObjectiveFctnNorm,'value',3);
end

setappdata(gcf,'opt',opt);

%-----------------------------------------------
%----------BEGIN SET CALLBACKS------------------

%----generic----
function local_generic_set(obj,eventD,prop)

opt = getappdata(gcf,'opt');

try
    try
        a = eval(get(obj,'string'));
    catch
        a = get(obj,'string');
    end
    opt.(prop)=a;
    setappdata(gcf,'isSaved',1);
catch
    msgbox('Invalid Entry: Option was not modified','Error','Error','modal');
    set(obj,'ForegroundColor',[1 0 0]);
    return
end

%change color to black:
set(obj,'ForegroundColor',[0 0 0]);
setappdata(gcf,'opt',opt);



%----checkbox inverse----
function local_checkbox_inverse(obj,eventD,prop)

opt = getappdata(gcf,'opt');
handles = getappdata(gcf,'handles');
switch prop
    case 'omitInnerBound'
        if ~opt.omitInnerBound
            opt.omitOuterBound=0;
            set(handles.omitOuterBound,'value',1);
        end
    case 'omitOuterBound'
        if ~opt.omitOuterBound
            opt.omitInnerBound=0;
            set(handles.omitInnerBound,'value',1);
        end
end
    
opt.(prop)=~get(obj,'value');

setappdata(gcf,'opt',opt);
setappdata(gcf,'isSaved',1);


%----pull-downs----
function local_pulldowns(obj,eventD,prop)

opt = getappdata(gcf,'opt');

switch prop
    case 'display'
        types = {'off','final','notify','iter','all','ALL'};
    case 'paramOptimObjectiveFctnNorm'
        types = {'two','one','inf'};
end

opt.(prop)=types{get(obj,'value')};

setappdata(gcf,'opt',opt);
setappdata(gcf,'isSaved',1);

%----reset all options?----
function local_defaults(obj,eventD)

opt = getappdata(gcf,'opt');
optDef = DClab.DCOptions;

fn = fieldnames(opt);

for i1 = fn'
    if strcmp(i1{1},'guiHandle')
        advOptMod = false; %we don't want to replace the guiHandle property with the default
    elseif strcmp(i1{1},'userData')
        advOptMod = ~isempty(opt.userData);
    else
        advOptMod = ~isequal(opt.(i1{1}),optDef.(i1{1}));
    end
  
    if advOptMod 
        
        if ~(strcmp(i1{1},{'maxBranchBoundIter','branchBoundTermTol','omitInnerBound','omitOuterBound','display','tolFun','nRestart','paramOptimObjectiveFctnNorm'}))
            answer = questdlg('There are advanced options that are not set to the defaults.  Would you like to revert these options to default as well?','Revert to all Defaults','No');
            switch answer
                case 'Cancel'
                    return
                case 'Yes'
                    break
                case 'No'
                    tmp = DClab.DCOptions;
                    optDef = opt;
                    optDef.maxBranchBoundIter = tmp.maxBranchBoundIter;
                    optDef.branchBoundTermTol = tmp.branchBoundTermTol;
                    optDef.omitInnerBound = tmp.omitInnerBound;
                    optDef.omitOuterBound = tmp.omitOuterBound;
                    optDef.display = tmp.display;
                    optDef.tolFun = tmp.tolFun;
                    break
            end
        end
    end
    
end

local_setForm([],[],optDef);
                    




