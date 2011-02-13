function varargout = ConsistTestViewer(varargin)
% CONSISTTESTVIEWER M-file for ConsistTestViewer.fig
%      CONSISTTESTVIEWER, by itself, creates a new CONSISTTESTVIEWER or raises the existing
%      singleton*.
%
%      H = CONSISTTESTVIEWER returns the handle to a new CONSISTTESTVIEWER or the handle to
%      the existing singleton*.
%
%      CONSISTTESTVIEWER(cellarray,elem) receives a cell array of
%      ConsistTest objects, and opens the gui with the info from the
%      element elem.
%
%      CONSISTTESTVIEWER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CONSISTTESTVIEWER.M with the given input arguments.
%
%      CONSISTTESTVIEWER('Property','Value',...) creates a new CONSISTTESTVIEWER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ConsistTestViewer_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ConsistTestViewer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ConsistTestViewer

% Last Modified by GUIDE v2.5 22-Feb-2010 11:19:54

%pragmas
%#function ConsistTestViewer colorEBars colorPBars
%#function ConsistTestViewer_OpeningFcn ConsistTestViewer_OutputFcn
%#function dsetorder_radio_callback eBar_downfcn eLamb_CreateFcn
%#function ExperimentList_Callback ExperimentList_CreateFcn export_Callback
%#function figure1_CloseRequestFcn File_Callback LBEdit_CreateFcn
%#function NextPB_Callback ParameterList_Callback ParameterList_CreateFcn
%#function pBar_downfcn pLamb_CreateFcn plotEBars plotPBars PrevPB_Callback
%#function sortAscd_radio_callback SortAscending_Callback
%#function SortByNum_Callback SortDescending_Callback
%#function sortDscd_radio_callback timer_Callback UBEdit_CreateFcn

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ConsistTestViewer_OpeningFcn, ...
                   'gui_OutputFcn',  @ConsistTestViewer_OutputFcn, ...
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

% --- Executes just before ConsistTestViewer is made visible.
function ConsistTestViewer_OpeningFcn(hObject, eventdata, handles, varargin) %#ok
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ConsistTestViewer (see VARARGIN)

%blow away any existing timers saved in the handles structure
if isfield(handles,'timer')
  try
    stop(handles.timer);
    %delete(handles.timer);
    handles = rmfield(handles,'timer');
  catch
    disp([lasterr ': Noncritical error deleting timer methods'])
    %do nothing
  end
end 

% Choose default command line output for ConsistTestViewer
handles.output = hObject;

%Add needed fields to the handles structure
if length(varargin) == 2
  handles.testObjCell = varargin{1};
  handles.elem = varargin{2};
  testObj = handles.testObjCell{handles.elem};
  if isempty(testObj)
    handles.LB = -1;
    handles.UB = 1;
    handles.unitArray = DClab.ModelAndObservationPair;
    handles.paramAss = DClab.FreeParameter;
    handles.exp.lamb = 0;
    handles.exp.lambR = 0;
    handles.param.lamb = 0;
    handles.param.lambR = 0;
    handles.exp.highlighted = 1;
    handles.exp.order = 1;
    hanldes.exp.scale = 'linear';
    handles.param.highlighted = 1;
    handles.param.order = 1;
    handles.ExpList = {' '};
    handles.ParamList = {' '};
    tstr = 'Empty ConsistTest Object';
  else
    
    handles.LB = testObj.LB;
    handles.UB = testObj.UB;
    PDset = testObj.PolyDataset;
    handles.unitArray = PDset.ModelAndObservationPair;
    handles.paramAss = PDset.FreeParameter;
    handles.exp.scale = 'linear';

    mults = testObj.upperBndMults;
    if isempty(mults)
      handles.exp.lamb = 0;
      handles.exp.lambR = 0;
      handles.ExpList = {'Multipliers unavailable: upper bound was not computed'};
      handles.param.lamb = 0;
      handles.param.lambR = 0;
      handles.ParamList = {'Multipliers unavailable: upper bound was not computed'};
    else
      %handles.exp.lamb = max([mults.expl mults.expu],[],2);
      handles.exp.lamb = max(mults.expl,0);
      handles.exp.lambR = max(mults.expu,0);
      handles.ExpList = {handles.unitArray.name}';
      %handles.param.lamb = max(mults.param,[],2);
      %handles.param.lamb = max([mults.paraml mults.paramu],[],2);
      handles.param.lamb = max(mults.paraml,0);
      handles.param.lambR = max(mults.paramu,0);
      handles.ParamList = {};
     for i = 1:length(handles.paramAss)
             uDstruct = handles.paramAss(i).userData;
             pKey = uDstruct.parameterKey;
             handles.ParamList = [handles.ParamList,pKey];
     end
        handles.ParamList = handles.ParamList';  
    end
    handles.exp.highlighted = 1;
    handles.exp.order = (1:length(handles.exp.lamb))';
    
    %for i1 = 1:length(handles.exp.lamb)
    %  u = handles.unitArray(i1).ResponseObservation.uncVect);
    %  u = u(2);
    %  handles.exp.lamb(i1) = handles.exp.lamb(i1)/u;
    %end
    handles.param.highlighted = 1;
    handles.param.order = (1:length(handles.param.lamb))';

    DsetName = testObj.DatasetName;
    if isempty(DsetName)
      name = 'The Dataset';
    else
      name = ['Dataset ' DsetName];
    end
    if handles.LB >= 0
      tstr = [name ' is Consistent'];
    elseif handles.UB < 0
      tstr = [name ' is Inconsistent'];
    elseif isnan(handles.LB(1))
      tstr = ['Consistency of ' name ' is Inconclusive: Lower Bound was not Computed'];
    elseif isnan(handles.UB(1))
      tstr = ['Consistency of ' name ' is Inconclusive: Upper Bound was not Computed'];
    else
      tstr = ['Consistency of ' name ' is Inconclusive: Try Increasing the number of Branch & Bound Iterations'];
    end
  end
  % Keep track of the old list box top, the old order, and the old
  % highlighted rows. These will be compared to the current values by a
  % timer callback function, which scrolls and highlights the bar graphs
  handles.exp.prevLBTop = 1;
  handles.exp.prevOrder = (1:length(handles.exp.lamb))';
  handles.exp.prevLBval = 1;

  handles.param.prevLBTop = 1;
  handles.param.prevOrder = (1:length(handles.param.lamb))';
  handles.param.prevLBval = 1;
  
else
  error('==weird inputs==')
end

%===Tweak the UI appearance===
%Navegation enable
nobj = size(handles.testObjCell,1);
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

%Start with experiment and parameter lists in DCDataset order
set(handles.dsetorder_radio,'Value',1);

%Position and content of listboxes
set(handles.ExperimentList,'String',handles.ExpList,'Value',1);
% If the names are short, we need to shorten up the listbox by trimming
% some off the bottom
tmp = size(char(handles.ExpList),2);
pos = get(handles.ExperimentList,'Position');
pos(2) = 26.35; % pos(2) = 27.154; %set to default
pos(4) = 11.1; % pos(4) = 11.3077;
if tmp < 45
  pos(2) = pos(2)+1;
  pos(4) = pos(4)-1;
end 
set(handles.ExperimentList,'Position',pos)

set(handles.ParameterList,'String',handles.ParamList,'Value',1);
% If the names are short, we need to shorten up the listbox by trimming
% some off the bottom
tmp = size(char(handles.ParamList),2);
pos = get(handles.ParameterList,'Position');
pos(2) = 11.38; %pos(2) = 11; %set to default
pos(4) = 11.1; % pos(4) = 11.3077;
if tmp < 45
  pos(2) = pos(2)+1;
  pos(4) = pos(4)-1;
end
set(handles.ParameterList,'Position',pos)

%fix weird text size
pos = get(handles.text24,'position');
pos(3) = 100;
set(handles.text24,'position',pos);

% Draw the Cmeas interval. 
axes(handles.Interval);
cla(handles.Interval)
if isnan(handles.LB(1)) && isnan(handles.UB(1))
  error('Both bounds NaN: Check your code')
elseif isnan(handles.LB(1))
  line([0.2 0.3],ones(1,2)*handles.UB(1),'LineWidth',3);
  line([0 1],[0 0],'Color',[0 0 0]);
  
  top = max(0,handles.UB(1));
  bottom = min(0,handles.UB(1));
  a(1) = 0; a(2) = 0.5;
  a(3) = bottom-0.05*(top-bottom);
  a(4) = top+0.05*(top-bottom);
  axis(a);
  set(handles.CmeasPanel,'Title',{'Consistency Measure','Upper Bound'})
 
  set([handles.UBEdit handles.Upperbnd_text],'Visible','on')
  set([handles.LBEdit handles.Lowerbnd_text],'Visible','off')
  set(handles.UBEdit,'String',num2str(handles.UB(1),'%0.5g'))
 
elseif isnan(handles.UB(1))
  line([0.2 0.3],ones(1,2)*handles.LB(1),'LineWidth',3);
  line([0 1],[0 0],'Color',[0 0 0]);
  
  top = max(0,handles.LB(1));
  bottom = min(0,handles.LB(1));
  a(1) = 0; a(2) = 0.5;
  a(3) = bottom-0.05*(top-bottom);
  a(4) = top+0.05*(top-bottom);
  axis(a);
  set(handles.CmeasPanel,'Title',{'Consistency Measure','Lower Bound'})
  
  set([handles.UBEdit handles.Upperbnd_text],'Visible','off')
  set([handles.LBEdit handles.Lowerbnd_text],'Visible','on')
  set(handles.LBEdit,'String',num2str(handles.LB(1),'%0.5g'))
else
  line([0.2 0.3],ones(1,2)*handles.LB(1),'LineWidth',2);
  line([0.2 0.3],ones(1,2)*handles.UB(1),'LineWidth',2);
  line([0.25 0.25],[handles.LB(1) handles.UB(1)],'LineWidth',3);
  line([0 1],[0 0],'Color',[0 0 0]);

  top = max(0,handles.UB(1));
  bottom = min(0,handles.LB(1));
  a(1) = 0; a(2) = 0.5;
  a(3) = bottom-0.05*(top-bottom);
  a(4) = top+0.05*(top-bottom);
  axis(a);
  set(handles.CmeasPanel,'Title',{'Consistency Measure'})
  
  set([handles.UBEdit handles.Upperbnd_text],'Visible','on')
  set([handles.LBEdit handles.Lowerbnd_text],'Visible','on')
  set(handles.LBEdit,'String',num2str(handles.LB(1),'%0.5g'))
  set(handles.UBEdit,'String',num2str(handles.UB(1),'%0.5g'))
end
set(handles.Interval,'XTickLabel',{});

%plot and color the bars
handles = plotEBars(handles,1);
handles = plotPBars(handles,1);
colorEBars(handles,1);
colorPBars(handles,1);

%title string
set(handles.titleStr,'string',tstr)

%parameter and target list titles
set(handles.abcTargets,'string',[DsetName ' targets']);
set(handles.abcParameters,'string',[DsetName ' parameters']);

%===end tweak UI===

% Update handles structure
guidata(hObject, handles);

% There is no callback for the scroll bar on the listbox, so to move 
% the bar graphs with the scroll bar, we create a timer that periodically 
% call the timer callback. This callback determines the current 
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

% UIWAIT makes ConsistTestViewer wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = ConsistTestViewer_OutputFcn(hObject, eventdata, handles) %#ok
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

%===Create functions===
% These executes during object creation, after setting all properties.
% Listbox and edit controls usually have a white background on Windows.
function ExperimentList_CreateFcn(hObject, eventdata, handles) %#ok
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ParameterList_CreateFcn(hObject, eventdata, handles) %#ok
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function eLamb_CreateFcn(hObject, eventdata, handles) %#ok
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function pLamb_CreateFcn(hObject, eventdata, handles) %#ok
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function UBEdit_CreateFcn(hObject, eventdata, handles) %#ok
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function LBEdit_CreateFcn(hObject, eventdata, handles) %#ok
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%===Callback functions===

% --- Executes on button press in SortByNum.
function SortByNum_Callback(hObject, eventdata, handles) %#ok
% This function returns the order of the experiment list and the multiplier
% bars to the initial 1:m order.

return

% Clear the other checkboxes
set(handles.SortAscending,'Value',0);
set(handles.SortDescending,'Value',0);

% We will make it impossible to unclick a checkbox by clicking on it. To
% unclick it, the user must click SortAscending or SortDescending
if get(hObject,'Value')
  handles.exp.order = (1:length(handles.exp.lamb))';
  set(handles.ExperimentList,'String',handles.ExpList);

  handles.param.order = (1:length(handles.param.lamb))';
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

  % If an experiment is highlighted, we need to
  % a) change the value of the experiment listbox
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

return

% Clear the other checkboxes
set(handles.SortByNum,'Value',0);
set(handles.SortDescending,'Value',0);

% We will make it impossible to unclick a checkbox by clicking on it. To
% unclick it, the user must click SortByNum or SortDescending
if get(hObject,'Value')
  % Set the order variable and display the names in this order. The timer 
  % callback will notice of the order was changed and update the plots. 
  [trash idx] = sort(handles.exp.lamb,1,'ascend'); %#ok
  handles.exp.order = idx;
  set(handles.ExperimentList,'String',handles.ExpList(idx,:));

  [trash idx] = sort(handles.param.lamb,1,'ascend'); %#ok
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

return

% Clear the other checkboxes
set(handles.SortByNum,'Value',0);
set(handles.SortAscending,'Value',0);

% We will make it impossible to unclick a checkbox by clicking on it. To
% unclick it, the user must click SortAscending or SortByNum
if get(hObject,'Value')
  % Set the order variable and display the names in this order. The timer 
  % callback will notice of the order was changed and update the plots. 
  [temp idx] = sort(handles.exp.lamb,1,'descend'); %#ok
  handles.exp.order = idx;
  set(handles.ExperimentList,'String',handles.ExpList(idx,:));

  [trash idx] = sort(handles.param.lamb,1,'descend'); %#ok
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
% a) determine which experiment it was
% b) update handles.highlightedExp
% c) refresh the textual display
% The timer callback will take care of everything else

LBval = get(hObject,'Value');
expnum = handles.exp.order(LBval);
handles.exp.highlighted = expnum;

% Update the textual display
lamb = handles.exp.lamb(expnum);
set(handles.eLamb,'String',num2str(lamb,'%0.4g'))
set(handles.eLambR,'String',num2str(handles.exp.lambR(expnum),'%0.4g'));

%udRM = get(get(handles.unitArray{expnum},'ResponseModel'),'userData');
%udRO = get(get(handles.unitArray{expnum},'ResponseObservation'),'userData');

%set(handles.ShowXML_menu,'Label',['Reaction: ' ud.parameterID]);

% if strcmp(get(handles.figure1,'SelectionType'),'alternate')
%   hand = get(handles.unitArray(expnum),'guiDisplayCallback');
%   if ~isempty(hand)
%     hand(handles.unitArray(expnum));
%   end
% end

guidata(hObject,handles)

function eBar_downfcn(obj, eventdata) %#ok

handles = guidata(gcbo);

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
nShow = length(handles.exp.lamb); %this is the number of experiments
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

handles = guidata(gcbo);

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
nShow = length(handles.param.lamb); %this is the number of experiments
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
nShow = length(handles.exp.lamb); %this is the number of experiments
nDisp = min(maxShow,nShow); %this is the number of bars we will draw

% Order the multipliers in the order they are to be displayed then just
% grab the ones corresponding to the visible listbox elements.
order = handles.exp.order;

if strcmp(handles.exp.scale,'linear')
    lamb = handles.exp.lamb(order);
    lambR = handles.exp.lambR(order);
else
    lamb = log10(handles.exp.lamb(order));
    lambR = log10(handles.exp.lambR(order));
end

% Some funky code to determing the x-axis limits. First determine the 
% largest mult (among ALL so the scale is constant)
mx = max([lamb; lambR]); 
mn = min([lamb; lambR]);
lamb1=lamb;
lamb2=lamb;
if ceil(mx)-floor(mn)>10
    lamb1((mx-lamb)>10) = mx-10;
    lamb2((mx-lamb)>10) = mx-10;
    mx = max([lamb1; lamb2]);
    mn = min([lamb1; lamb2]); 
end
lambToPlot = lamb(eTop:nDisp+eTop-1);
lambToPlotR = lambR(eTop:nDisp+eTop-1);

if isempty(mx)
    xlim = [0 1];
else
    % Let's get the right limit within a 25% from mx. Since we know the range
    % of coeff, the rounding of ceil is of the order we want.
    switch handles.exp.scale
        case 'linear'
            xlim = [0 1.25*mx];
        case 'log'
            x1 = floor(mn);
            x2 = ceil(mx);
            xlim = [0 x2-x1];
            lambToPlot = lambToPlot-x1;
            lambToPlotR = lambToPlotR-x1;
    end
end



h = barh(handles.TopPanel,lambToPlot(end:-1:1));
hR = barh(handles.TopPanelR,lambToPlotR(end:-1:1));
hh = get(h,'children');
hhR = get(hR,'children');
set([h hR], 'buttondownfcn', {@eBar_downfcn})
set([hh hhR],'cdatamapping', 'direct','facevertexcdata', ones(nDisp,1)) %all the same color
handles.eBar = h;
handles.eBarR = hR;
handles.ePatch = hh;
handles.ePatchR = hhR;

set([handles.TopPanel handles.TopPanelR],'YTickLabel',{});

switch handles.exp.scale
    case 'linear'
        xlabel(handles.TopPanel,'\lambda')
        xlabel(handles.TopPanelR,'\lambda')
    case 'log'
        xlabel(handles.TopPanel,'log10(\lambda)');
        xlabel(handles.TopPanelR,'log10(\lambda)');
end


% The top is nDisp, the bottom is nDisp-maxShow (since the listbox 
% is tall enough for max show and the bars to line up with it).
ylim = [nDisp-maxShow+0.5 nDisp+0.5]; 
set(handles.TopPanel,'XLim',xlim,'YLim',ylim);
set(handles.TopPanelR,'XLim',xlim,'YLim',ylim);
if strcmp(handles.exp.scale,'log')
    %now relabel ticks to correct numbers
    xticks = get(handles.TopPanel,'xtick');
    xticksR = get(handles.TopPanelR,'xtick');
    set(handles.TopPanel,'xticklabel',xticks+x1);
    set(handles.TopPanelR,'xticklabel',xticksR+x1);
end

function handles = plotPBars(handles,pTop)
% This function draws the bar graphs. The second input is the number of the
% exp. listbox element that is at the top of the currently displayed elements.
% I.e., if the top element is the list box is visible, eTop = 1. The third
% input is the analogous number for the parameter listbox

maxShow = 10; %this is the maximum number of bars we can fit in the display
nShow = length(handles.param.lamb); %this is the number of experiments
nDisp = min(maxShow,nShow); %this is the number of bars we will draw

% Order the multipliers in the order they are to be displayed then just
% grab the ones corresponding to the visible listbox elements.
order = handles.param.order;
switch handles.exp.scale
    case 'linear'
        lamb = handles.param.lamb(order);
        lambR = handles.param.lambR(order);
    case 'log'
        lamb = log10(abs(handles.param.lamb(order)));
        lambR = log10(abs(handles.param.lambR(order)));
end

% Some funky code to determing the x-axis limits. First determine the 
% largest mult (among ALL so the scale is constant)
mx = max([lamb; lambR]); 
mn = min([lamb; lambR]);
lamb1=lamb;
lamb2 = lambR;
if ceil(mx)-floor(mn)>10
    lamb1((mx-lamb)>10) = mx-10;
    lamb2((mx-lamb)>10) = mx-10;
    mx = max([lamb1; lamb2]);
    mn = min([lamb1; lamb2]); 
end
lambToPlot = lamb(pTop:nDisp+pTop-1);
lambToPlotR = lambR(pTop:nDisp+pTop-1);

if isempty(mx)
    xlim = [0 1];
else
    % Let's get the right limit within a 25% from mx. Since we know the range
    % of coeff, the rounding of ceil is of the order we want.
    switch handles.exp.scale
        case 'linear'
            xlim = [0 1.25*mx];
        case 'log'
            x1 = floor(mn);
            x2 = ceil(mx);
            xlim = [0 x2-x1];
  %          if x2==x1
  %              xlim = [0 1];
  %           end
            lambToPlot = lambToPlot-x1;
            lambToPlotR = lambToPlotR-x1;
    end
end

h = barh(handles.BottomPanel,lambToPlot(end:-1:1));
hR = barh(handles.BottomPanelR,lambToPlotR(end:-1:1));
hh = get(h,'children');
hhR = get(hR,'children');
set([h hR], 'buttondownfcn', {@pBar_downfcn})
set([hh hhR],'cdatamapping', 'direct','facevertexcdata', ones(nDisp,1)) %all the same color
handles.pBar = h;
handles.pBarR = hR;
handles.pPatch = hh;
handles.pPatchR = hhR;

set(handles.BottomPanel,'YTickLabel',{});
set(handles.BottomPanelR,'YTickLabel',{});

switch handles.exp.scale
    case 'linear'
        xlabel(handles.BottomPanel,'\lambda')
        xlabel(handles.BottomPanelR,'\lambda')
    case 'log'
        xlabel(handles.BottomPanel,'log10(\lambda)');
        xlabel(handles.BottomPanelR,'log10(\lambda)');
end

% The top is nDisp, the bottom is nDisp-maxShow (since the listbox 
% is tall enough for max show and the bars to line up with it).
ylim = [nDisp-maxShow+0.5 nDisp+0.5]; 
set(handles.BottomPanel,'XLim',xlim,'YLim',ylim);
set(handles.BottomPanelR,'XLim',xlim,'YLim',ylim);
if strcmp(handles.exp.scale,'log')
    %now relabel ticks to correct numbers
    xticks = get(handles.BottomPanel,'xtick');
    set(handles.BottomPanel,'xticklabel',xticks+x1);
    xticksR = get(handles.BottomPanelR,'xtick');
    set(handles.BottomPanelR,'xticklabel',xticksR+x1);
end


function colorEBars(handles,barnumber)
% barnumber is the number of bars from the top of the visible listbox 
% display, with barnumber = 1 meaning the first element is visible

maxShow = 10; %this is the maximum number of bars we can fit in the display
nShow = length(handles.exp.lamb); %this is the number of experiments
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
set(handles.ePatch, 'facevertexcdata', faceCData);
set(handles.ePatchR, 'facevertexcdata', faceCData);

function colorPBars(handles,barnumber)
% barnumber is the number of bars from the top of the visible listbox 
% display, with barnumber = 1 meaning the first element is visible

maxShow = 10; %this is the maximum number of bars we can fit in the display
nShow = length(handles.param.lamb); %this is the number of experiments
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
set(handles.pPatch, 'facevertexcdata', faceCData);
set(handles.pPatchR, 'facevertexcdata', faceCData);

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
lamb = handles.param.lamb(paramnum);
set(handles.pLamb,'String',num2str(lamb,'%0.4g'))
set(handles.pLambR,'String',num2str(handles.param.lambR(paramnum),'%0.4g'))

%ud = get(handles.paramAss(paramnum),'userData');

%set(handles.ShowXML_menu,'Label',['Reaction: ' ud.parameterID]);

% if strcmp(get(handles.figure1,'SelectionType'),'alternate')
%   hand = get(handles.paramAss(paramnum),'guiDisplayCallback');
%   if ~isempty(hand)
%     hand(handles.paramAss(paramnum));
%   end
% end

guidata(hObject,handles)


function timer_Callback(varargin)
% This is the callback of the timer. Everytime it executes it checks of the
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
  try
    stop(handles.timer);
    %delete(handles.timer);
  catch
    disp([lasterr ': noncritical problem with timer'])
  end
end

delete(hObject);

% --- Executes on button press in PrevPB.
function PrevPB_Callback(hObject, eventdata, handles) %#ok

%stop the timer
if isfield(handles,'timer')
  try
    stop(handles.timer);
    %delete(handles.timer);
    handles = rmfield(handles,'timer');
  catch
    disp([lasterr ': noncritical problem with timer'])
  end
end
guidata(hObject,handles)

%Call with the previous element
ConsistTestViewer(handles.testObjCell,handles.elem-1);

% --- Executes on button press in NextPB.
function NextPB_Callback(hObject, eventdata, handles) %#ok

%stop the timer
if isfield(handles,'timer')
  try
    stop(handles.timer);
    %delete(handles.timer);
    handles = rmfield(handles,'timer');
  catch
    disp([lasterr ': noncritical problem with timer'])
  end
end
guidata(hObject,handles)

%Call with the next element
ConsistTestViewer(handles.testObjCell,handles.elem+1);

% --------------------------------------------------------------------
function File_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function export_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to export (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
msgbox('Feature currently unavailable','We''re sorry','warn','modal')






% --- Executes on button press in dsetorder_radio.
function dsetorder_radio_Callback(hObject, eventdata, handles)
% hObject    handle to dsetorder_radio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of dsetorder_radio

% We will make it impossible to unclick a radio button by clicking on it. To
% unclick it, the user must click SortAscending or SortDescending
if get(hObject,'Value')
  handles.exp.order = (1:length(handles.exp.lamb))';
  set(handles.ExperimentList,'String',handles.ExpList);

  handles.param.order = (1:length(handles.param.lamb))';
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

  % If an experiment is highlighted, we need to
  % a) change the value of the experiment listbox
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

% --- Executes on button press in sortAscd_radio.
function sortAscd_radio_Callback(hObject, eventdata, handles)
% hObject    handle to sortAscd_radio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of sortAscd_radio

% We will make it impossible to unclick a radio by clicking on it. To
% unclick it, the user must click SortByNum or SortDescending
if get(hObject,'Value')
  % Set the order variable and display the names in this order. The timer 
  % callback will notice of the order was changed and update the plots. 
  [trash idx] = sort(handles.exp.lamb+handles.exp.lambR,1,'ascend'); %#ok
  handles.exp.order = idx;
  set(handles.ExperimentList,'String',handles.ExpList(idx,:));
  
  [trash idx] = sort(handles.param.lamb+handles.param.lambR,1,'ascend'); %#ok
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

% --- Executes on button press in sortDscd_radio.
function sortDscd_radio_Callback(hObject, eventdata, handles)
% hObject    handle to sortDscd_radio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of sortDscd_radio

% We will make it impossible to unclick a radio by clicking on it. To
% unclick it, the user must click SortAscending or SortByNum
if get(hObject,'Value')
  % Set the order variable and display the names in this order. The timer 
  % callback will notice of the order was changed and update the plots. 
  [temp idx] = sort(handles.exp.lamb+handles.exp.lambR,1,'descend'); %#ok
  handles.exp.order = idx;
  set(handles.ExperimentList,'String',handles.ExpList(idx,:));

  [trash idx] = sort(handles.param.lamb+handles.param.lambR,1,'descend'); %#ok
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


% --- Executes when selected object is changed in uipanel11.
function uipanel11_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel11 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

handles.exp.scale = lower(get(hObject,'string'));
topE = get(handles.ExperimentList,'ListBoxTop');
topP = get(handles.ParameterList,'ListBoxTop');

handles = plotEBars(handles,topE);
colorEBars(handles,handles.exp.highlighted-(topE-1));
handles = plotPBars(handles,topP);
colorPBars(handles,handles.param.highlighted-(topP-1));

guidata(hObject,handles);


% --------------------------------------------------------------------
function parameterID_menu_Callback(hObject, eventdata, handles)
% hObject    handle to ShowXML_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
paramnum = handles.param.highlighted;

uDstruct = handles.paramAss(paramnum).userData;
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

paramnum = handles.param.highlighted;

uDstruct = handles.paramAss(paramnum).userData;
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

paramnum = handles.param.highlighted;

uDstruct = handles.paramAss(paramnum).userData;
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
expnum = handles.exp.highlighted;

uDstruct = handles.unitArray(expnum).ResponseObservation.userData;
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
expnum = handles.exp.highlighted;
uDstruct = handles.unitArray(expnum).ResponseObservation.userData;
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
expnum = handles.exp.highlighted;

uDstruct = handles.unitArray(expnum).ResponseModel.userData;
pid = uDstruct.modelID;
if iscell(pid) && isscalar(pid)
    uDstruct.callback2show(pid{1});
elseif iscell(pid)
    uDstruct.callback2show2(pid{1},pid{2});
elseif ischar(pid)
    uDstruct.callback2show(pid);
end


% --------------------------------------------------------------------
function Exp_Menu_Callback(hObject, eventdata, handles)
% hObject    handle to Exp_Menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Param_Menu_Callback(hObject, eventdata, handles)
% hObject    handle to Param_Menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function ShowXML_menu_Callback(hObject, eventdata, handles)
% hObject    handle to ShowXML_menu_Callback_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
paramnum = handles.param.highlighted;
uDstruct = handles.paramAss(paramnum).userData;

mid = uDstruct.trialModelID;
% y = gate2primeData('get',{mid});
% mname = getPreferredKey(y);
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


% --- Executes on button press in help.
function help_Callback(hObject, eventdata, handles)
% hObject    handle to help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h = msgbox(['Example 1: If consistency measure is [0.3 0.4], this means the dataset is consistent. '...
    'All parameter and experiment uncertainties can be decreased by at least 30% but not more than 40% and remain consistent.' 10 10,...
    'Example 2: If consistency measure is [-0.2 -0.15], then the dataset is inconsistant. '...
    'It is sufficient to increase all parameter and experiment uncertainties by 20% (and perhaps as little as 15%) to achieve consistency.' 10 10,...
    'Example 3: If consistency measure is [-0.1 0.2], then the consitency test is inconclusive. '...
    'It is sufficient to increase all parameter and experiment uncertainties by 10% to achieve consistency.'],'Consistency Measure', 'modal');




function pLambR_Callback(hObject, eventdata, handles)
% hObject    handle to pLambR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pLambR as text
%        str2double(get(hObject,'String')) returns contents of pLambR as a double


% --- Executes during object creation, after setting all properties.
function pLambR_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pLambR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function eLambR_Callback(hObject, eventdata, handles)
% hObject    handle to eLambR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of eLambR as text
%        str2double(get(hObject,'String')) returns contents of eLambR as a double


% --- Executes during object creation, after setting all properties.
function eLambR_CreateFcn(hObject, eventdata, handles)
% hObject    handle to eLambR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
