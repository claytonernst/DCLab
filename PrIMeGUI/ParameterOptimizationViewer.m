function varargout = ParameterOptimizationViewer(varargin)
% PARAMETEROPTIMIZATIONVIEWER M-file for ParameterOptimizationViewer.fig
%      PARAMETEROPTIMIZATIONVIEWER, by itself, creates a new PARAMETEROPTIMIZATIONVIEWER or raises the existing
%      singleton*.
%
%      H = PARAMETEROPTIMIZATIONVIEWER returns the handle to a new PARAMETEROPTIMIZATIONVIEWER or the handle to
%      the existing singleton*.
%
%      PARAMETEROPTIMIZATIONVIEWER(cellarray,elem) receives a cell array of
%      BestFitX objects, and opens the gui with the info from the
%      element elem.
%
%      PARAMETEROPTIMIZATIONVIEWER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PARAMETEROPTIMIZATIONVIEWER.M with the given input arguments.
%
%      PARAMETEROPTIMIZATIONVIEWER('Property','Value',...) creates a new PARAMETEROPTIMIZATIONVIEWER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ParameterOptimizationViewer_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ParameterOptimizationViewer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ParameterOptimizationViewer

% Last Modified by GUIDE v2.5 26-Feb-2009 15:36:00

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ParameterOptimizationViewer_OpeningFcn, ...
                   'gui_OutputFcn',  @ParameterOptimizationViewer_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else 
  gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before ParameterOptimizationViewer is made visible.
function ParameterOptimizationViewer_OpeningFcn(hObject, eventdata, handles, varargin) %#ok
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ParameterOptimizationViewer (see VARARGIN)

%blow away any existing timers saved in the handles structure
if isfield(handles,'timer')
  try
    stop(handles.timer);
    %delete(handles.timer);
    handles = rmfield(handles,'timer');
  catch
    %do nothing
  end
end 

% Choose default command line output for ParameterOptimizationViewer
handles.output = hObject;

%Add needed fields to the handles structure
if length(varargin) == 2
  handles.fitObjCell = varargin{1}; %cell array of ParameterOptimization objects
  handles.elem = varargin{2}; %element of cell we will display
  fitObj = handles.fitObjCell{handles.elem};
  if isempty(fitObj)
    handles.unitArray = DClab.ModelAndObservationPair;
    handles.paramAss = DClab.FreeParameter;
    handles.wresid = 0;
    handles.bestx = 0;
    handles.normBestX = 0;
    handles.weights = 1;
    handles.param.highlighted = 1;
    handles.param.order = 1;
    handles.exp.highlighted = 1;
    handles.exp.order = 1;
    handles.ExpList = {' '};
    handles.ParamList = {' '};
    tstr = 'Empty ParameterOptimization Object';
  else
    PDset = fitObj.PDset;
    handles.unitArray = PDset.ModelAndObservationPair;
    handles.paramAss = PDset.FreeParameter;
  
    %stuff for experiment display
    handles.wresid = fitObj.wresid;
    handles.weights = fitObj.weights;
    handles.exp.highlighted = 1;
    handles.exp.order = (1:length(handles.wresid))';
    handles.ExpList = {handles.unitArray.name}';
    
    %stuff for parameter display
    handles.bestx = fitObj.bestx;
    nom = vertcat(handles.paramAss.nominal);
    %bnds = vertcat(handles.paramAss.range);
    %maxRad = max(bnds(:,2)-nom,nom-bnds(:,1)); %interval radius relative to nominal 
    handles.normBestx = (handles.bestx - nom); %./maxRad; %scale to be between +/- 1
    %handles.ParamList = {handles.paramAss.name}';
    handles.ParamList = {};
     for i = 1:length(handles.paramAss)
           uDstruct = handles.paramAss(i).userData;
             pKey = uDstruct.parameterKey;
             handles.ParamList = [handles.ParamList,pKey];
     end
        handles.ParamList = handles.ParamList';  
 
    handles.param.highlighted = 1;
    handles.param.order = (1:length(handles.bestx))';
    
    DsetName = fitObj.DsetName;
    if isempty(DsetName)
      name = 'Unnamed Dataset';
    else
      name = ['Dataset ' DsetName];
    end
    OptimOptns = fitObj.OptimOptns;
    costUB = fitObj.costUB;
    switch OptimOptns.method
        case 1
            tstr = [name ': the value of LS-H objective function is ' num2str(costUB, '%0.3d')];
        case 2
            tstr = [name ': the value of LS-F objective function is ' num2str(costUB, '%0.3d')];
        case 3
            tstr = [name ': the value of 1N-F objective function is ' num2str(costUB, '%0.3d')];
        case 4
            tstr = [name ': the values of MO-F objective functions are ' num2str(costUB(1), '%0.3d') ', ' num2str(costUB(2), '%0.3d')];
        case 6
            tstr = [name ': the value of LD-H objective function is ' num2str(costUB, '%0.3d')];
    end
  end
  
  % Keep track of the old list box top, the old order, and the old
  % highlighted rows. These will be compared to the current values by a
  % timer callback function, which scrolls and highlights the bar graphs
  handles.exp.prevLBTop = 1;
  handles.exp.prevOrder = (1:length(handles.wresid))';
  handles.exp.prevLBval = 1;
  
  handles.param.prevLBTop = 1;
  handles.param.prevOrder = (1:length(handles.bestx))';
  handles.param.prevLBval = 1;  
else
  disp('==weird inputs==')
  %keyboard
end

%=======Modify the UI appearance=======
%Enable states of navigation buttons
nobj = size(handles.fitObjCell,1);
if handles.elem == nobj
  set(handles.NextPB,'enable','off')
else
  set(handles.NextPB,'enable','on')
end
if handles.elem == 1
  set(handles.PrevPB,'enable','off')
else
  set(handles.PrevPB,'enable','on')
end

%Start with experiment and parameter lists in Dataset order
set(handles.SortByNum,'Value',1);
set(handles.SortAscending,'Value',0);
set(handles.SortDescending,'Value',0);

%Position and content of listboxes
set(handles.ExperimentList,'String',handles.ExpList,'Value',1);
% If the names are short, we need to shorten up the listbox by trimming
% some off the bottom
tmp = size(char(handles.ExpList),2);
pos = get(handles.ExperimentList,'Position');
pos(2) = 7.1; %set to default
pos(4) = 11.1;
if tmp < 45 %this needs to get tweaked
  pos(2) = pos(2)+1;
  pos(4) = pos(4)-1;
end
set(handles.ExperimentList,'Position',pos)

set(handles.ParameterList,'String',handles.ParamList,'Value',1);
% If the names are short, we need to shorten up the listbox by trimming
% some off the bottom
tmp = size(char(handles.ParamList),2);
pos = get(handles.ParameterList,'Position');
pos(2) = 5.3; %set to default
pos(4) = 11.1;
if tmp < 45
  pos(2) = pos(2)+1;
  pos(4) = pos(4)-1;
end
set(handles.ParameterList,'Position',pos)

%plot and color the bars
handles = plotEBars(handles,1);
handles = plotPBars(handles,1);
colorEBars(handles,1);
colorPBars(handles,1);

%The title string
set(handles.titleStr,'string',tstr)
 
%=================end modify UI appearence

% Update handles structure
guidata(hObject, handles);

% There is no callback for the scroll bar on the listbox, so to move 
% the bar graphs with the scroll bar, we create a timer that periodically 
% call the the timer callback. This callback determines the current 
% position of the scrollbar and updates the graphs accordingly. The period
% is 0.1 sec. This doesn't seem to hog the processor and provides a
% reasonable refresh rate. Timer.UserData is a graphics handle. With
% this available, the timer callback can access and edit the handles structure.
handles.timer = timer('executionmode','fixedrate',...
  'TimerFcn',@timer_Callback,'Period',0.1,'startdelay',0.25,'UserData',hObject);
start(handles.timer);

%Call the listbox callbacks to display the textual stuff
ExperimentList_Callback(handles.ExperimentList,eventdata,handles);
ParameterList_Callback(handles.ParameterList,eventdata,handles);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ParameterOptimizationViewer wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = ParameterOptimizationViewer_OutputFcn(hObject, eventdata, handles) %#ok
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

%===Create functions. These Execute during object creation, after setting all properties.===
%edit and listbox controls usually have a white background on Windows.
function ExperimentList_CreateFcn(hObject, eventdata, handles) %#ok
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ParameterList_CreateFcn(hObject, eventdata, handles) %#ok
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ExpDataEdit_CreateFcn(hObject, eventdata, handles) %#ok
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ExpUncEdit_CreateFcn(hObject, eventdata, handles) %#ok
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ExpWeightEdit_CreateFcn(hObject, eventdata, handles) %#ok
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ParamNomEdit_CreateFcn(hObject, eventdata, handles) %#ok
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ParamRangeEdit_CreateFcn(hObject, eventdata, handles) %#ok
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function OptimValEdit_CreateFcn(hObject, eventdata, handles) %#ok
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ExpResidEdit_CreateFcn(hObject, eventdata, handles) %#ok
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%===Callbacks===

% --- Executes on button press in SortByNum.
function SortByNum_Callback(hObject, eventdata, handles) %#ok
% This function returns the order of the experiment list and the multiplier
% bars to the initial 1:m order.

% Clear the other checkboxes
set(handles.SortAscending,'Value',0);
set(handles.SortDescending,'Value',0);

% We will make it impossible to unclick a checkbox by clicking on it. To
% unclick it, the user must click SortAscending or SortDescending
if get(hObject,'Value')
  handles.exp.order = (1:length(handles.wresid))';
  set(handles.ExperimentList,'String',handles.ExpList);

  handles.param.order = (1:length(handles.bestx))';
  set(handles.ParameterList,'String',handles.ParamList);
  
  % If an experiment is highlighted, we need to
  % a) change the value of the experiment listbox
  % b) change which bar is highlighted
  hExp = handles.exp.highlighted;
  if ~isempty(hExp)
    LBval = find(handles.exp.order == hExp);
    set(handles.ExperimentList,'Value',LBval)
    
    % Refresh the display of the experimentList listbox so the selected 
    % experiment is in dark blue.
    uicontrol(handles.ExperimentList)
  end

  % If a parameter is highlighted, we need to
  % a) change the value of the parameter listbox
  % b) change which bar is highlighted
  hParam = handles.param.highlighted;
  if ~isempty(hParam)
    LBval = find(handles.param.order == hParam);
    set(handles.ParameterList,'Value',LBval)
    
    % Refresh the display of the ParameterList listbox so the selected 
    % parameter is in dark blue.
    uicontrol(handles.ParameterList)
  end

  guidata(hObject,handles);
else
  set(hObject,'Value',1);
end
  
% --- Executes on button press in SortAscending.
function SortAscending_Callback(hObject, eventdata, handles) %#ok
% This function orders the experiments by the ascending order of their 
% multipliers with the top experiment having the largest?? multiplier.

% Clear the other checkboxes
set(handles.SortByNum,'Value',0);
set(handles.SortDescending,'Value',0);

% We will make it impossible to unclick a checkbox by clicking on it. To
% unclick it, the user must click SortByNum or SortDescending
if get(hObject,'Value')
  % Set the order variable and display the names in this order. The timer 
  % callback will notice of the order was changed and update the plots. 
  [trash idx] = sort(abs(handles.wresid),1,'ascend'); %#ok
  handles.exp.order = idx;
  set(handles.ExperimentList,'String',handles.ExpList(idx,:));

  [trash idx] = sort(abs(handles.normBestx),1,'ascend'); %#ok
  handles.param.order = idx;
  set(handles.ParameterList,'String',handles.ParamList(idx,:));
  
  % If an experiment is highlighted, we need to
  % a) change the value of the experiment listbox
  % b) change which bar is highlighted
  hExp = handles.exp.highlighted;
  if ~isempty(hExp)
    LBval = find(handles.exp.order == hExp);
    set(handles.ExperimentList,'Value',LBval)

    % Refresh the display of the experimentList listbox so the selected 
    % experiment is in dark blue.
    uicontrol(handles.ExperimentList)
  end

  hParam = handles.param.highlighted;
  if ~isempty(hParam)
    LBval = find(handles.param.order == hParam);
    set(handles.ParameterList,'Value',LBval)
    
    % Refresh the display of the ParameterList listbox so the selected 
    % parameter is in dark blue.
    uicontrol(handles.ParameterList)
  end
  
  guidata(hObject,handles);
else
  set(hObject,'Value',1);
end

% --- Executes on button press in SortDescending.
function SortDescending_Callback(hObject, eventdata, handles) %#ok
% This function orders the experiments by the descending order of their 
% multipliers with the top experiment having the smallest?? multiplier.

% Clear the other checkboxes
set(handles.SortByNum,'Value',0);
set(handles.SortAscending,'Value',0);

% We will make it impossible to unclick a checkbox by clicking on it. To
% unclick it, the user must click SortAscending or SortByNum
if get(hObject,'Value')
  % Set the order variable and display the names in this order. The timer 
  % callback will notice of the order was changed and update the plots. 
  [temp idx] = sort(abs(handles.wresid),1,'descend'); %#ok
  handles.exp.order = idx;
  set(handles.ExperimentList,'String',handles.ExpList(idx,:));

  [trash idx] = sort(abs(handles.normBestx),1,'descend'); %#ok
  handles.param.order = idx;
  set(handles.ParameterList,'String',handles.ParamList(idx,:));
  
  % If an experiment is highlighted, we need to
  % a) change the value of the experiment listbox
  % b) change which bar is highlighted
  hExp = handles.exp.highlighted;
  if ~isempty(hExp)
    LBval = find(handles.exp.order == hExp);
    set(handles.ExperimentList,'Value',LBval)
    
    % Refresh the display of the experimentList listbox so the selected 
    % experiment is in dark blue.
    uicontrol(handles.ExperimentList)
  end
  
  hParam = handles.param.highlighted;
  if ~isempty(hParam)
    LBval = find(handles.param.order == hParam);
    set(handles.ParameterList,'Value',LBval)
    
    % Refresh the display of the ParameterList listbox so the selected 
    % parameter is in dark blue.
    uicontrol(handles.ParameterList)
  end
  
  guidata(hObject,handles);
else
  set(hObject,'Value',1);
end

% --- Executes on selection change in ExperimentList.
function ExperimentList_Callback(hObject, eventdata, handles) %#ok

% If the user selects one the the listbox elements
% a) determine which experiment is was
% b) update handles.highlightedExp
% c) refresh the textual display
% The timer callback will take care of everything else

LBval = get(hObject,'Value');
expnum = handles.exp.order(LBval);
handles.exp.highlighted = expnum;

% Update the textual display
d = vertcat(handles.unitArray(expnum).observedValue);
u = {handles.unitArray(expnum).observationUncertainty}';
w = handles.weights(expnum);
set(handles.ExpDataEdit,'String',num2str(d,'%0.4g'))
if numel(u{1}) == 1
  set(handles.ExpUncEdit,'String',num2str(u{1},'%0.4g'))
elseif numel(u{1}) == 2  
  set(handles.ExpUncEdit,'String',['[' num2str(u{1}(1),'%7.4f') ',' num2str(u{1}(2),'%7.4f') ']'])
else
  error('Internal inconsistency, condition should never occur')
end
set(handles.ExpWeightEdit,'String',num2str(w,'%0.4g'))
set(handles.ExpResidEdit,'String',num2str(handles.wresid(expnum),'%0.4g'))
%{
if strcmp(get(handles.figure1,'SelectionType'),'open')
  hand = get(handles.unitArray(expnum),'dispCallback');
  if ~isempty(hand)
    hand(handles.unitArray(expnum));
  end
end

guidata(hObject,handles)
%}

function eBar_downfcn(obj, eventdata) %#ok

handles = guidata(obj);
% Determine the axis containing the patch and grab all its children
ax = get(obj, 'parent');  % axes containing the patch
chi = get(obj,'children');

% the ith column of xdata is the x-range of the ith bar on the graph
xdata = get(chi, 'xdata');
xdata = xdata([1,2],:);
% the ith column of ydata is the y-range of the ith bar on the graph
ydata = get(chi, 'ydata');
ydata = ydata([1, 3], :);

cpoint = get(ax, 'currentpoint');
cpoint = cpoint([1, 3]);    % (x,y) of current point

% Determine which bar the current point is within. The product is negative
% id cpoint(1) is less then xdata(1,i) and greater and xdata(2,i)
barnumber = find(prod(cpoint(1) - xdata) <= 0 & prod(cpoint(2) - ydata) <= 0);

% Occationally the buttom down is called but the coordinates are not
% contained within any bar. In such case, just exit
if isempty(barnumber)
    return
end

maxShow = 10; %this is the maximum number of bars we can fit in the display
nShow = length(handles.wresid); %this is the number of experiments
nDisp = min(maxShow,nShow); %this is the number of bars we will draw

% As defined, barnumber is the number counted from the bottom of the graph,
% we want the number counted from the top of the graph.
barnumber = find((nDisp:-1:1) == barnumber);

%the graph we clicked on is only as current as prev listbox top, so use
%prev listbox top
top = handles.exp.prevLBTop;
LBval = barnumber+top-1;

set(handles.ExperimentList,'Value',LBval)
uicontrol(handles.ExperimentList)
ExperimentList_Callback(handles.ExperimentList, eventdata, handles)

function pBar_downfcn(obj, eventdata) %#ok

if ~ishandle(obj)
  disp('no hand')
  return
end
  
handles = guidata(obj);
% Determine the axis containing the patch and grab all its children
ax = get(obj, 'parent');  % axes containing the patch
chi = get(obj,'children');

% the ith column of xdata is the x-range of the ith bar on the graph
xdata = get(chi, 'xdata');
xdata = xdata([1,2],:);
% the ith column of ydata is the y-range of the ith bar on the graph
ydata = get(chi, 'ydata');
ydata = ydata([1, 3], :);

cpoint = get(ax, 'currentpoint');
cpoint = cpoint([1, 3]);    % (x,y) of current point

% Determine which bar the current point is within. The product is negative
% id cpoint(1) is less then xdata(1,i) and greater and xdata(2,i)
barnumber = find(prod(cpoint(1) - xdata) <= 0 & prod(cpoint(2) - ydata) <= 0);

% Occationally the buttom down is called but the coordinates are not
% contained within any bar. In such case, just exit
if isempty(barnumber)
    return
end

maxShow = 10; %this is the maximum number of bars we can fit in the display
nShow = length(handles.bestx); %this is the number of experiments
nDisp = min(maxShow,nShow); %this is the number of bars we will draw

% As defined, barnumber is the number counted from the bottom of the graph,
% we want the number counted from the top of the graph.
barnumber = find((nDisp:-1:1) == barnumber);

%the graph we clicked on is only as current as prev listbox top, so use
%prev listbox top
top = handles.param.prevLBTop;
LBval = barnumber+top-1;

set(handles.ParameterList,'Value',LBval)
uicontrol(handles.ParameterList)
ParameterList_Callback(handles.ParameterList, eventdata, handles)

function handles = plotEBars(handles,eTop)
% This function draws the bar graphs. The second input is the number of the
% exp. listbox element that is at the top of the currently displayed elements.
% I.e., if the top element is the list box is visible, eTop = 1. The third
% input is the analogous number for the parameter listbox

maxShow = 10; %this is the maximum number of bars we can fit in the display
nShow = length(handles.wresid); %this is the number of experiments
nDisp = min(maxShow,nShow); %this is the number of bars we will draw

% Order the multipliers in the order they are to be displayed then just
% grab the ones corresponding to the visible listbox elements.
order = handles.exp.order;

lamb = handles.wresid(order);
lambToPlot = lamb(eTop:nDisp+eTop-1);

%axes(handles.ResidAxis)
%cla(handles.ResidAxis)
%text(0.3,1.1,0,'$w_e |M_e(x)-d_e|$','Units','Normalized','interpreter','latex','fontsize',12);
%set(handles.ResidAxis,'NextPlot','add')
h = barh(handles.ResidAxis,lambToPlot(end:-1:1));
TextInfo.Units              = 'pixels'   ;   
TextInfo.FontSize           = 12;
TextInfo.FontWeight         = get(handles.figure1,'DefaultTextFontWeight');
TextInfo.HorizontalAlignment= 'center'     ;
TextInfo.HandleVisibility   = 'callback' ;
TextInfo.VerticalAlignment  = 'bottom';
TextInfo.Color              = get(0,'FactoryUIControlForegroundColor');
text('Parent',handles.ResidAxis,TextInfo,'Position',[0.5 1.05],'String','$M(x)-d$','Units','Normalized','interpreter','latex','fontsize',12);

hh = get(h,'children');

set(h, 'buttondownfcn', {@eBar_downfcn})
set(hh,'cdatamapping', 'direct','facevertexcdata', ones(nDisp,1)) %all the same color
handles.eBar = h;
handles.ePatch = hh;     

set(handles.ResidAxis,'YTickLabel',{});

% Some funky code to determing the x-axis limits. First determine the 
% largest mult (among ALL so the scale is constant)
mx = max(abs(lamb)); %mx = max(lamb)
if mx==0
  mx = 1;
end
pwr = log10(mx); %detmine its order of mag
if pwr <= 1;
  pwr = floor(pwr);
else
  pwr = ceil(pwr);
end
coeff = mx/10^pwr; %now mx = coeff x 10^pwr and coeff \in [1, 10]

% Let's get the right limit within a 25% from mx. Since we know the range 
% of coeff, the rounding of ceil is of the order we want.
x2 = ceil(coeff*4)/4 * 10^pwr; %since 
x1 = -x2; %x1 = 0;
xlim = [x1 x2];

% The top is nDisp, the bottom is nDisp-maxShow (since the listbox 
% is tall enough for max show and the bars to line up with it).
ylim = [nDisp-maxShow+0.5 nDisp+0.5]; 
set(handles.ResidAxis,'XLim',xlim,'YLim',ylim);

function handles = plotPBars(handles,pTop)
% This function draws the bar graphs. The second input is the number of the
% exp. listbox element that is at the top of the currently displayed elements.
% I.e., if the top element is the list box is visible, eTop = 1. The third
% input is the analogous number for the parameter listbox

maxShow = 10; %this is the maximum number of bars we can fit in the display
nShow = length(handles.bestx); %this is the number of experiments
nDisp = min(maxShow,nShow); %this is the number of bars we will draw

% Order the multipliers in the order they are to be displayed then just
% grab the ones corresponding to the visible listbox elements.
order = handles.param.order;

lamb = handles.normBestx(order);
lambToPlot = lamb(pTop:nDisp+pTop-1);

h = barh(handles.OptParamAxis,lambToPlot(end:-1:1));
text('Parent',handles.OptParamAxis,'Position',[0.4 1.1],'String','$x-x_0$','Units','Normalized','interpreter','latex','fontsize',12);
hh = get(h,'children');
set(h, 'buttondownfcn', {@pBar_downfcn})
try
  set(hh,'cdatamapping', 'direct','facevertexcdata', ones(nDisp,1)) %all the same color
catch
  disp('asdf')
  %keyboard
end
handles.pBar = h;
handles.pPatch = hh;     

set(handles.OptParamAxis,'YTickLabel',{});

% Some funky code to determing the x-axis limits. First determine the 
% largest mult (among ALL so the scale is constant)
mx = max(abs(lamb)); %mx = max(lamb)
if mx==0
  mx = 1;
end
pwr = log10(mx); %detmine its order of mag
if pwr <= 1;
  pwr = floor(pwr);
else
  pwr = ceil(pwr);
end
coeff = mx/10^pwr; %now mx = coeff x 10^pwr and coeff \in [1, 10]

% Let's get the right limit within a 25% from mx. Since we know the range 
% of coeff, the rounding of ceil is of the order we want.
x2 = ceil(coeff*4)/4 * 10^pwr; %since 
x1 = -x2; %x1 = 0;
xlim = [x1 x2];
%xlim = [-1.1 1.1];

% The top is nDisp, the bottom is nDisp-maxShow (since the listbox 
% is tall enough for max show and the bars to line up with it).
ylim = [nDisp-maxShow+0.5 nDisp+0.5]; 
set(handles.OptParamAxis,'XLim',xlim,'YLim',ylim);

function colorEBars(handles,barnumber)
% barnumber is the number of bars from the top of the visible listbox 
% display, with barnumber = 1 meaning the first element is visible

maxShow = 10; %this is the maximum number of bars we can fit in the display
nShow = length(handles.wresid); %this is the number of experiments
nDisp = min(maxShow,nShow); %this is the number of bars we will draw

if barnumber < 1 || barnumber > maxShow
  return
end

% Set the faceCData, with the first element being the color of the top bar on
% the display
faceCData = ones(nDisp,1);
faceCData(barnumber) = 64;

% Invert the order of faceCData, since matlab treats the first element as
% the color of the bottom bar.
faceCData = faceCData(end:-1:1);

%Update colors of both bar plots
%Update colors of both bar plots
try
  set(handles.ePatch, 'facevertexcdata', faceCData);
catch
  disp('problem with colorbars')
  %keyboard
end

function colorPBars(handles,barnumber)
% barnumber is the number of bars from the top of the visible listbox 
% display, with barnumber = 1 meaning the first element is visible

maxShow = 10; %this is the maximum number of bars we can fit in the display
nShow = length(handles.bestx); %this is the number of experiments
nDisp = min(maxShow,nShow); %this is the number of bars we will draw

if barnumber < 1 || barnumber > maxShow
  return
end

% Set the faceCData, with the first element being the color of the top bar on
% the display
faceCData = ones(nDisp,1);
faceCData(barnumber) = 64;

% Invert the order of faceCData, since matlab treats the first element as
% the color of the bottom bar.
faceCData = faceCData(end:-1:1);

%Update colors of both bar plots
try
  set(handles.pPatch, 'facevertexcdata', faceCData);
catch
  disp('problem with colorbars')
  %keyboard
end


% --- Executes on selection change in ParameterList.
function ParameterList_Callback(hObject, eventdata, handles) %#ok

% If the user selects one the the listbox elements
% a) determine which parameter it was
% b) update handles.param.highlighted
% c) refresh the textual display
% The timer callback will take care of everything else

LBval = get(hObject,'Value');
paramnum = handles.param.order(LBval);
handles.param.highlighted = paramnum;

% Update the textual display
%lamb = handles.param.lamb(paramnum);
%set(handles.pLamb,'String',num2str(lamb,'%0.4g'))
% Update the textual display
nom = handles.paramAss(paramnum).nominal;
range = handles.paramAss(paramnum).range;
opt = handles.bestx(paramnum);
set(handles.ParamNomEdit,'String',num2str(nom,'%0.4g'))
set(handles.ParamRangeEdit,'String',['[' num2str(range(1),'%0.4g') ',' num2str(range(2),'%0.4g') ']'])
set(handles.OptimValEdit,'String',num2str(opt,'%0.4g'))
%{
if strcmp(get(handles.figure1,'SelectionType'),'open')
  hand = get(handles.paramAss(paramnum),'dispCallback');
  if ~isempty(hand)
    hand(handles.paramAss(paramnum));
  end
end

guidata(hObject,handles)
%}
function timer_Callback(varargin)
% This is the callback of the timer. Everytime it executes it checks if the
% listbox has been moved since the previous execution (by comparing the
% current listboxtop to the saved prev. value).

timer = varargin{1};
try
  handles = guidata(timer.UserData); %.UserData is a graphics handle
catch
  return
end
  
top = get(handles.ExperimentList,'ListBoxTop');
order = handles.exp.order;
LBval = get(handles.ExperimentList,'value');

if top ~= handles.exp.prevLBTop || ~isequal(order,handles.exp.prevOrder) || LBval ~= handles.exp.prevLBval  
  handles = plotEBars(handles,top);
  hExp = handles.exp.highlighted;
  if ~isempty(hExp) 
    % Determine the bar that should be highlighted and color it. In the
    % call to colorbars, we need to subtract off top-1, since the 1:top-1 are
    % not drawn.
    colorEBars(handles,LBval-(top-1));
  end
  handles.exp.prevLBTop = top; %update the PrevLBTop, etc
  handles.exp.prevOrder = order;
  handles.exp.prevLBval = LBval;
end

top = get(handles.ParameterList,'ListBoxTop');
order = handles.param.order;
LBval = get(handles.ParameterList,'value');

if top ~= handles.param.prevLBTop || ~isequal(order,handles.param.prevOrder) || LBval ~= handles.param.prevLBval  
  handles = plotPBars(handles,top);
  hParam = handles.param.highlighted;
  if ~isempty(hParam) 
    % Determine the bar that should be highlighted and color it. In the
    % call to colorbars, we need to subtract off top-1, since the 1:top-1 are
    % not drawn.
    colorPBars(handles,LBval-(top-1));
  end
  handles.param.prevLBTop = top; %update the PrevLBTop, etc
  handles.param.prevOrder = order;
  handles.param.prevLBval = LBval;
end

guidata(timer.UserData,handles);


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles) %#ok
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
if isfield(handles,'timer')
  stop(handles.timer);
  %delete(handles.timer);
end

delete(hObject);

function PrevPB_Callback(hObject, eventdata, handles) %#ok
%stop the timer
if isfield(handles,'timer')
  try
    stop(handles.timer);
    %delete(handles.timer);
  catch
    %do nothing
  end
end

%Call with the previous element
ParameterOptimizationViewer(handles.fitObjCell,handles.elem-1);

% --- Executes on button press in NextPB.
function NextPB_Callback(hObject, eventdata, handles) %#ok
%stop the timer
if isfield(handles,'timer')
  try
    stop(handles.timer);
    %delete(handles.timer);
  catch
    %do nothing
  end
end

%Call with the next element
ParameterOptimizationViewer(handles.fitObjCell,handles.elem+1);


% --------------------------------------------------------------------
function parameterID_menu_Callback(hObject, eventdata, handles)
% hObject    handle to ShowXML_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
LBval = get(handles.ParameterList,'value');

uDstruct = handles.paramAss(LBval).userData;
pid = uDstruct.parameterID;
if iscell(pid) && isscalar(pid)
    uDstruct.callback2show(pid{1});
elseif iscell(pid)
    uDstruct.callback2show(pid{1});
elseif ischar(pid)
    uDstruct.callback2show(pid);
end


% --------------------------------------------------------------------
function Parameter_Links_Menu_Callback(hObject, eventdata, handles)
% hObject    handle to Parameter_Links_Menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
LBval = get(handles.ParameterList,'value');

uDstruct = handles.paramAss(LBval).userData;
pid = uDstruct.parameterLinks;
if iscell(pid) && isscalar(pid)
    uDstruct.callback2show(pid{1});
elseif iscell(pid)
    uDstruct.callback2show2(pid{1},pid{2});
elseif ischar(pid)
    uDstruct.callback2show(pid);
end


% --------------------------------------------------------------------
function trial_model_menu_Callback(hObject, eventdata, handles)
% hObject    handle to trial_model_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
LBval = get(handles.ParameterList,'value');

uDstruct = handles.paramAss(LBval).userData;
pid = uDstruct.trialModelID;
if iscell(pid) && isscalar(pid)
    uDstruct.callback2show(pid{1});
elseif iscell(pid)
    uDstruct.callback2show(pid{1});
elseif ischar(pid)
    uDstruct.callback2show(pid);
end

% --------------------------------------------------------------------
function target_menu_Callback(hObject, eventdata, handles)
% hObject    handle to ShowXML_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
LBval = get(handles.ExperimentList,'value');

uDstruct = handles.unitArray(LBval).ResponseObservation.userData;
pid = uDstruct.targetID;
if iscell(pid) && isscalar(pid)
    uDstruct.callback2show(pid{1});
elseif iscell(pid)
    uDstruct.callback2show(pid{1});
elseif ischar(pid)
    uDstruct.callback2show(pid);
end


% --------------------------------------------------------------------
function Target_Links_menu_Callback(hObject, eventdata, handles)
% hObject    handle to Target_Links_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
LBval = get(handles.ExperimentList,'value');

uDstruct = handles.unitArray(LBval).ResponseObservation.userData;
pid = uDstruct.targetLinks;
if iscell(pid) && isscalar(pid)
    uDstruct.callback2show(pid{1});
elseif iscell(pid)
    uDstruct.callback2show(pid{1});
elseif ischar(pid)
    uDstruct.callback2show(pid);
end


% --------------------------------------------------------------------
function Model_menu_Callback(hObject, eventdata, handles)
% hObject    handle to Model_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
LBval = get(handles.ExperimentList,'value');

uDstruct = handles.unitArray(LBval).ResponseModel.userData;
pid = uDstruct.modelID;
if iscell(pid) && isscalar(pid)
    uDstruct.callback2show(pid{1});
elseif iscell(pid)
    uDstruct.callback2show2(pid{1},pid{2});
elseif ischar(pid)
    uDstruct.callback2show(pid);
end


% --------------------------------------------------------------------
function Exp_menu_Callback(hObject, eventdata, handles)
% hObject    handle to Exp_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Param_Menu_Callback(hObject, eventdata, handles)
% hObject    handle to Param_Menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function ShowXML_menu_Callback(hObject, eventdata, handles)
% hObject    handle to ShowXML_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
LBval = get(handles.ParameterList,'value');
uDstruct = handles.paramAss(LBval).userData;
mid = uDstruct.trialModelID;

%y = gate2primeData('get',{mid});
%mname = getPreferredKey(y);
mname = gate2primeData('getPreferredKey',{'primeID',mid});

set(handles.trial_model_menu, 'Label', char(mname));
pid = uDstruct.parameterLinks;
recMark = pid{2}(1:length(pid{2})-8);
      switch recMark
         case 'rk'
            y = 'Rate Coefficient';
         case 'thp'
            y = 'Thermodynamic Polynomial';
      end
set(handles. Parameter_Links_Menu, 'Label', y);


% --- Executes on button press in DataTablePB.
function DataTablePB_Callback(hObject, eventdata, handles)
% hObject    handle to DataTablePB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

fitObj = handles.fitObjCell{handles.elem};
fitTable(fitObj,'data');


% --- Executes on button press in ParamTablePB.
function ParamTablePB_Callback(hObject, eventdata, handles)
% hObject    handle to ParamTablePB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

fitObj = handles.fitObjCell{handles.elem};
fitTable(fitObj,'param');

