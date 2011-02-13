function varargout = PredictionViewer(varargin)
% PREDICTIONVIEWER M-file for PredictionViewer.fig
%      PREDICTIONVIEWER, by itself, creates a new PREDICTIONVIEWER or raises the existing
%      singleton*.
%
%      H = PREDICTIONVIEWER returns the handle to a new PREDICTIONVIEWER or the handle to
%      the existing singleton*.
%
%      PREDICTIONVIEWER(cellarray,elem) receives a cell array of
%      Prediction objects, and opens the gui with the info from the
%      element elem.
%
%
%      PREDICTIONVIEWER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PREDICTIONVIEWER.M with the given input arguments.
%
%      PREDICTIONVIEWER('Property','Value',...) creates a new PREDICTIONVIEWER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before PredictionViewer_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to PredictionViewer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help PredictionViewer

% Last Modified by GUIDE v2.5 22-Feb-2010 12:31:28

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @PredictionViewer_OpeningFcn, ...
    'gui_OutputFcn',  @PredictionViewer_OutputFcn, ...
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

% --- Executes just before PredictionViewer is made visible.
function PredictionViewer_OpeningFcn(hObject, eventdata, handles, varargin) %#ok
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to PredictionViewer (see VARARGIN)

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

% Choose default command line output for PredictionViewer
handles.output = hObject;

%Add needed fields to the handles structure
if length(varargin) == 2
    handles.predObjCell = varargin{1};
    handles.elem = varargin{2};
    obj = handles.predObjCell{handles.elem};
    if isempty(obj)
        handles.LBo = -1;
        handles.LBi = -1;
        handles.UBi = 1;
        handles.UBo = 1;
        handles.unitArray = DClab.ModelAndObservationPair;
        handles.paramAss = DClab.FreeParameter;
        handles.exp.lambUB = 0;
        handles.exp.lambLB = 0;
        handles.exp.lambUBR = 0;
        handles.exp.lambLBR = 0;
        handles.ExpList = {' '};
        handles.param.lambUB = 0;
        handles.param.lambLB = 0;
        handles.param.lambUBR = 0;
        handles.param.lambLBR = 0;
        handles.ParamList = {' '};
        handles.exp.highlighted = 1;
        handles.exp.order = 1;
        handles.exp.scale = 'linear';
        handles.param.highlighted = 1;
        handles.param.order = 1;
        tstr = 'Empty Prediction Object';
    else
        handles.LBo = obj.LBo;
        handles.LBi = obj.LBi;
        handles.UBi = obj.UBi;
        handles.UBo = obj.UBo;
        handles.exp.scale = 'linear';

        PDset = obj.PolyDataset;
        handles.unitArray = PDset.ModelAndObservationPair;
        handles.paramAss = PDset.FreeParameter;

        cutoff = 0;%1e-8; %We will zero out multipliers < cutoff (why?)

        %Grab the multipliers corresponding to the upper and lower bounds on
        %the prediction interval. Insert dummy stuff if the outer bounds were
        %not computed.
        mults = obj.outerBndMults;
        if isempty(mults)
            handles.exp.lambUB = 0;
            handles.exp.lambLB = 0;
            handles.exp.lambUBR = 0;
            handles.exp.lambLBR = 0;
            handles.ExpList = {'Multipliers unavailable: no outer bounds were computed'};
            handles.param.lambUB = 0;
            handles.param.lambLB = 0;
            handles.param.lambUBR = 0;
            handles.param.lambLBR = 0;
            handles.ParamList = {'Multipliers unavailable: no outer bounds were computed'};
        else
            if PDset.nPairs == 0 %no experiment multipliers to display
                handles.exp.lambUB = 0;
                handles.exp.lambLB = 0;
                handles.exp.lambUBR = 0;
                handles.exp.lambLBR = 0;
                handles.ExpList = {' '};
            else
                handles.exp.lambLB = mults.lower.expl;
                handles.exp.lambLB(handles.exp.lambLB < cutoff) = 0;
                handles.exp.lambLBR = mults.lower.expu;
                handles.exp.lambLBR(handles.exp.lambLBR < cutoff) = 0;
                handles.exp.lambUB = mults.upper.expl;
                handles.exp.lambUB(handles.exp.lambUB < cutoff) = 0;
                handles.exp.lambUBR = mults.upper.expu;
                handles.exp.lambUBR(handles.exp.lambUBR < cutoff) = 0;
                handles.ExpList = {handles.unitArray.name}';
            end
            handles.param.lambUB = mults.upper.paraml;
            handles.param.lambUB(handles.param.lambUB < cutoff) = 0;
            handles.param.lambUBR = mults.upper.paramu;
            handles.param.lambUBR(handles.param.lambUBR < cutoff) = 0;
            handles.param.lambLB = mults.lower.paraml;
            handles.param.lambLB(handles.param.lambLB < cutoff) = 0;
            handles.param.lambLBR = mults.lower.paramu;
            handles.param.lambLBR(handles.param.lambLBR < cutoff) = 0;
            handles.ParamList = {};
            for i = 1:length(handles.paramAss)
            uDstruct = handles.paramAss(i).userData;
             pKey = uDstruct.parameterKey;
             handles.ParamList = [handles.ParamList,pKey];
            end
            handles.ParamList = handles.ParamList';  
        end


        % Suppose the experiments in the dataset are numbered 1 to m. Then
        % .highlightedExp contains the experiment number of the highlighted
        % experiment and .order contains the order they are currently displayed.
        % (with order(1) as the top experiment in the list)
        handles.exp.highlighted = 1; %initially highlight the first (top)
        handles.exp.order = (1:length(handles.exp.lambUB))';

        handles.param.highlighted = 1;
        handles.param.order = (1:length(handles.param.lambUB))';

        %Assemble the title string
        RM = obj.ResponseModel;
        Dname = obj.DatasetName;
        if isempty(Dname)
            Dname = 'Unnamed';
        end
        if isempty(RM.name)
            tstr = ['Interval Prediction for Unnamed Model from Dataset ' Dname];
        else
            tstr = ['Interval Prediction for Surrogate Model ' RM.name ' from Dataset ' Dname];
        end
    end

    %   for i1 = 1:length(handles.exp.lambUB)
    %     u = handles.unitArray(i1).ExpAss.uncVect;
    %     u = u(2);
    %     handles.exp.lambUB(i1) = handles.exp.lambUB(i1)/u;
    %     handles.exp.lambLB(i1) = handles.exp.lambLB(i1)/u;
    %   end

    % Keep track of the old list box top, the old order, and the old
    % highlighted rows. These will be compared to the current values by a
    % timer callback function, which scrolls and highlights the bar graphs
    handles.exp.prevLBTop = 1;
    handles.exp.prevOrder = (1:length(handles.exp.lambUB))';
    handles.exp.prevLBval = 1;

    handles.param.prevLBTop = 1;
    handles.param.prevOrder = (1:length(handles.param.lambUB))';
    handles.param.prevLBval = 1;

else
    disp('==weird inputs==')
    %keyboard
end

%===Tweak the UI appearance===
%Navegation enable
nobj = size(handles.predObjCell,1);
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
set(handles.dsetorder_radio,'Value',1);

%Position and content of listboxes
set(handles.ExperimentList,'String',handles.ExpList,'Value',1);
% If the names are short, we need to shorten up the listbox by trimming
% some off the bottom
tmp = size(char(handles.ExpList),2);
pos = get(handles.ExperimentList,'Position');
pos(2) = 27.5; % pos(2) = 27.5384; %set to default
pos(4) = 12.1; % pos(4) = 12.2308;
if tmp < 25
    pos(2) = pos(2)+1;
    pos(4) = pos(4)-1;
end
set(handles.ExperimentList,'Position',pos)

set(handles.ParameterList,'String',handles.ParamList,'Value',1);
% If the names are short, we need to shorten up the listbox by trimming
% some off the bottom
tmp = size(char(handles.ParamList),2);
pos = get(handles.ParameterList,'Position');
pos(2) = 11.2; %pos(2) = 11.1538; %set to default
pos(4) = 12.1; % pos(4) = 12.2308;
if tmp < 25
    pos(2) = pos(2)+1;
    pos(4) = pos(4)-1;
end
set(handles.ParameterList,'Position',pos)

% Draw the prediction interval.
axes(handles.Interval);
cla(handles.Interval)
dstr = {'Feasible parameter values are';...
    'found at which the surrogate model ';...
    'produces the endpoints of the inner ';...
    'interval.';'';...
    'Constrained by the dataset, the';...
    'surrogate model does not take values ';...
    'outside the outer interval.'};
uD = handles.predObjCell{handles.elem}.ResponseModel.userData;

if isnan(handles.LBo) && isnan(handles.LBi)
    %keyboard
    error('Both lower bounds NaN: Check your code')
elseif isnan(handles.LBo) %just inner bounds
    %inner bound interval
    line([0.45 0.55],ones(1,2)*handles.LBi,'LineWidth',2);
    line([0.45 0.55],ones(1,2)*handles.UBi,'LineWidth',2);
    line([0.5 0.5],[handles.LBi handles.UBi],'LineWidth',3);

    a = axis(handles.Interval);
    top = handles.UBi;
    bottom = handles.LBi;
    a(1) = 0; a(2) = 1;
    a(3) = bottom-0.05*(top-bottom);
    a(4) = top+0.1*(top-bottom);
    axis(a);

    text(0.4,top+0.05*(top-bottom),'inner')
    set(handles.intervalDesc,'string',dstr(1:3))

    %pos1 = get(handles.InnerIntEdit,'Position');
    %pos2 = get(handles.OuterIntEdit,'Position');
    %shift = 0.5*(pos1(1)-pos2(1));
    %pos1(1) = pos1(1)-shift;
    %set(handles.InnerIntEdit,'Position',pos1)
    set([handles.InnerIntEdit handles.InnerInt_text],'Visible','on')
    set([handles.OuterIntEdit handles.OuterInt_text],'Visible','off')
    set(handles.InnerIntEdit,'String',['[' num2str(10^handles.LBi,'%0.4g') ', ' num2str(10^handles.UBi,'%0.4g') '] ' uD.targetUnits])

elseif isnan(handles.LBi) %just outer bounds
    %outer bound interval
    line([0.45 0.55],ones(1,2)*handles.LBo,'LineWidth',2);
    line([0.45 0.55],ones(1,2)*handles.UBo,'LineWidth',2);
    line([0.5 0.5],[handles.LBo handles.UBo],'LineWidth',3);

    a = axis(gca);
    top = handles.UBo;
    bottom = handles.LBo;
    a(1) = 0; a(2) = 1;
    a(3) = bottom-0.05*(top-bottom);
    a(4) = top+0.1*(top-bottom);
    axis(a);

    text(0.4,top+0.05*(top-bottom),'outer')
    set(handles.intervalDesc,'string',dstr(5:end))

    set([handles.InnerIntEdit handles.InnerInt_text],'Visible','off')
    set([handles.OuterIntEdit handles.OuterInt_text],'Visible','on')
    set(handles.OuterIntEdit,'String',['[' num2str(10^handles.LBo,'%0.4g') ', ' num2str(10^handles.UBo,'%0.4g') '] ' uD.targetUnits])

else
    %outer bound interval
    line([0.2 0.3],ones(1,2)*handles.LBo,'LineWidth',2);
    line([0.2 0.3],ones(1,2)*handles.UBo,'LineWidth',2);
    line([0.25 0.25],[handles.LBo handles.UBo],'LineWidth',3);

    %inner bound interval
    line([0.7 0.8],ones(1,2)*handles.LBi,'LineWidth',2);
    line([0.7 0.8],ones(1,2)*handles.UBi,'LineWidth',2);
    line([0.75 0.75],[handles.LBi handles.UBi],'LineWidth',3);

    a = axis(gca);
    top = handles.UBo;
    bottom = handles.LBo;
    a(1) = 0; a(2) = 1;
    a(3) = bottom-0.05*(top-bottom);
    a(4) = top+0.1*(top-bottom);
    axis(a);

    text(0.13,top+0.05*(top-bottom),'outer')
    text(0.63,top+0.05*(top-bottom),'inner')
    set(handles.intervalDesc,'string',dstr)

    set([handles.InnerIntEdit handles.InnerInt_text],'Visible','on')
    set([handles.OuterIntEdit handles.OuterInt_text],'Visible','on')
    set(handles.InnerIntEdit,'String',['[' num2str(10^handles.LBi,'%0.3g') ', ' num2str(10^handles.UBi,'%0.3g') '] ' uD.targetUnits])
    set(handles.OuterIntEdit,'String',['[' num2str(10^handles.LBo,'%0.3g') ', ' num2str(10^handles.UBo,'%0.3g') '] ' uD.targetUnits])
end
set(handles.Interval,'XTickLabel',{}); %clear the XTickLabels

if isstruct(uD) && isfield(uD,'targetLabel')
    ylabel(['log10( ' texlabel(uD.targetLabel,'literal') ' )'])
end

%plot and color the bars (this action modifies the handles structure)
handles = plotEBars(handles,1);
handles = plotPBars(handles,1);
colorEBars(handles,1);
colorPBars(handles,1);

%title string
set(handles.titleStr,'string',tstr)
if length(tstr) > 100
    set(handles.titleStr,'tooltipstring',tstr)
else
    set(handles.titleStr,'tooltipstring','')
end

set(handles.abcTargets,'string',[Dname ' targets'])
set(handles.abcParameters,'string',[Dname ' parameters'])

%===end tweak UI===

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
    'TimerFcn',@timer_Callback,'Period',0.15,'startdelay',0.25,'UserData',hObject);
start(handles.timer);

%Call the listbox callbacks to display the textual stuff
ExperimentList_Callback(handles.ExperimentList,eventdata,handles);
ParameterList_Callback(handles.ParameterList,eventdata,handles);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes PredictionViewer wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = PredictionViewer_OutputFcn(hObject, eventdata, handles) %#ok
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

%===Create functions===
% These executes during object creation, after setting all properties.
% Listbox and edit controls usually have a white background on Windows.
% --- Executes during object creation, after setting all properties.
function edit_CreateFcn(hObject, eventdata, handles) %#ok
%create function used for the four edit boxes used for display
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ExperimentList_CreateFcn(hObject, eventdata, handles) %#ok
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ParameterList_CreateFcn(hObject, eventdata, handles) %#ok
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%===Callback functions===

% --- Executes on button press in SortByNum.
function SortByNum_Callback(hObject, eventdata, handles) %#ok
% This function returns the order of the experiment list and the multiplier
% bars to the initial 1:m order.

% Clear the other checkboxes
set(handles.SortRight,'Value',0);
set(handles.SortLeft,'Value',0);

% We will make it impossible to unclick a checkbox by clicking on it. To
% unclick it, the user must click SortLeft or SortRight
if get(hObject,'Value')
    % Reset the order variable to 1:m and display the names in this order.
    % The timer callback will notice of the order was changed and update the
    % plots.
    handles.exp.order = (1:length(handles.exp.lambUB))';
    set(handles.ExperimentList,'String',handles.ExpList);

    handles.param.order = (1:length(handles.param.lambUB))';
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

    hParam = handles.param.highlighted;
    if ~isempty(hParam)
        LBval = find(handles.param.order == hParam);
        set(handles.ParameterList,'Value',LBval)
        uicontrol(handles.ParameterList)
    end

    guidata(hObject,handles);
else
    set(hObject,'Value',1);
end

% --- Executes on button press in SortRight.
function SortRight_Callback(hObject, eventdata, handles) %#ok
% This function orders the experiments by their multipliers corresponding
% to the lower bound on the prediction interval, with the top experiment
% having the largest multiplier.

% Clear the other checkboxes
set(handles.SortByNum,'Value',0);
set(handles.SortLeft,'Value',0);

% We will make it impossible to unclick a checkbox by clicking on it. To
% unclick it, the user must click SortByNum or SortLeft
if get(hObject,'Value')
    % Set the order variable and display the names in this order. The timer
    % callback will notice of the order was changed and update the plots.
    [temp idx] = sort(handles.exp.lambUB,1,'descend'); %#ok
    handles.exp.order = idx;
    set(handles.ExperimentList,'String',handles.ExpList(idx,:));

    [temp idx] = sort(handles.param.lambUB,1,'descend'); %#ok
    handles.param.order = idx;
    set(handles.ParameterList,'String',handles.ParamList(idx,:));


    % The timer callback will notice of the order was changed and update the
    % plots.

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
        uicontrol(handles.ParameterList)
    end

    guidata(hObject,handles);
else
    set(hObject,'Value',1);
end

% --- Executes on button press in SortLeft.
function SortLeft_Callback(hObject, eventdata, handles) %#ok
% This function orders the experiments by their multipliers corresponding
% to the upper bound on the prediction interval, with the top experiment
% having the largest multiplier.

% Clear the other checkboxes
set(handles.SortByNum,'Value',0);
set(handles.SortRight,'Value',0);

% We will make it impossible to unclick a checkbox by clicking on it. To
% unclick it, the user must click SortByNum or SortRight
if get(hObject,'Value')
    % Set the order variable, display the names in this order, and
    % reorder the bar graphs. The timmer callback will notice of the order
    % was changed and update the plots.
    [temp idx] = sort(handles.exp.lambLB,1,'descend'); %#ok
    handles.exp.order = idx;
    set(handles.ExperimentList,'String',handles.ExpList(idx,:));

    [temp idx] = sort(handles.param.lambLB,1,'descend'); %#ok
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

switch get(handles.maxradio,'Value')
    case 0 %minimization
        % Update the textual display
        lambUB = handles.exp.lambLBR(expnum);
        set(handles.eUB,'String',num2str(lambUB,'%0.4g'))
        lambLB = handles.exp.lambLB(expnum);
        set(handles.eLB,'String',num2str(lambLB,'%0.4g'))
    case 1 %maximization
        % Update the textual display
        lambUB = handles.exp.lambUBR(expnum);
        set(handles.eUB,'String',num2str(lambUB,'%0.4g'))
        lambLB = handles.exp.lambUB(expnum);
        set(handles.eLB,'String',num2str(lambLB,'%0.4g'))
end

%Call the display Callback if available
% if strcmp(get(handles.figure1,'SelectionType'),'alternate')
%     hand = get(handles.unitArray(expnum),'guiDisplayCallback');
%     if ~isempty(hand)
%         hand(handles.unitArray(expnum));
%     end
% end
guidata(hObject,handles)

% --- Executes on selection change in ParameterList.
function ParameterList_Callback(hObject, eventdata, handles) %#ok

% If the user selects one the the listbox elements
% a) determine which experiment is was
% b) update handles.highlightedExp
% c) refresh the textual display
% The timer callback will take care of everything else

LBval = get(hObject,'Value');
paramnum = handles.param.order(LBval);
handles.param.highlighted = paramnum;

switch get(handles.maxradio,'Value')
    case 0 %minimization
        % Update the textual display
        lambUB = handles.param.lambLBR(paramnum);
        set(handles.pUB,'String',num2str(lambUB,'%0.4g'))
        lambLB = handles.param.lambLB(paramnum);
        set(handles.pLB,'String',num2str(lambLB,'%0.4g'))
    case 1 %maximization
        % Update the textual display
        lambUB = handles.param.lambUBR(paramnum);
        set(handles.pUB,'String',num2str(lambUB,'%0.4g'))
        lambLB = handles.param.lambUB(paramnum);
        set(handles.pLB,'String',num2str(lambLB,'%0.4g'))
end

%Call the display Callback if available
% if strcmp(get(handles.figure1,'SelectionType'),'alternate')
%     hand = get(handles.paramAss(paramnum),'dispCallback');
%     if ~isempty(hand)
%         hand(handles.paramAss(paramnum));
%     end
% end
guidata(hObject,handles)


function eBar_downfcn(obj, eventdata)

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

maxShow = 11; %this is the maximum number of bars we can fit in the display
nShow = length(handles.exp.lambUB); %this is the number of experiments
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

%guidata(gcbo,handles);

function pBar_downfcn(obj, eventdata)

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

% Occationally the button down is called but the coordinates are not
% contained within any bar. In such case, just exit
if isempty(barnumber)
    return
end

maxShow = 11; %this is the maximum number of bars we can fit in the display
nShow = length(handles.param.lambUB); %this is the number of experiments
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

function handles = plotEBars(handles,top)
% This function draws the bar graphs. The second input is the number of the
% listbox element that is at the top of the currently displayed elements.
% I.e., if the top element is the list box is visible, top = 1.

maxShow = 11; %this is the maximum number of bars we can fit in the display
nShow = length(handles.exp.lambUB); %this is the number of experiments
nDisp = min(maxShow,nShow); %this is the number of bars we will draw

% Order the multipliers in the order they are to be displayed then just
% grab the ones corresponding to the visible listbox elements.
order = handles.exp.order;
switch get(handles.maxradio,'Value')
    case 1
        if strcmp(handles.exp.scale,'log')
            lambUB = log10(handles.exp.lambUBR(order));
            lambLB = log10(handles.exp.lambUB(order));
        else
            lambUB = handles.exp.lambUBR(order);
            lambLB = handles.exp.lambUB(order);
        end
    case 0
        if strcmp(handles.exp.scale,'log')
            lambUB = log10(handles.exp.lambLBR(order));
            lambLB = log10(handles.exp.lambLB(order));
        else
            lambUB = handles.exp.lambLBR(order);
            lambLB = handles.exp.lambLB(order);
        end
end

% Some funky code to determing the x-axis limits. First determine the
% largest mult (among ALL so the scale is constant)
mx = max(max([lambUB lambLB]));
mn = min(min([lambUB lambLB]));
lambUB1=lambUB; lambLB1=lambLB;
if ceil(mx)-floor(mn)>10
    lambUB1((mx-lambUB)>10) = mx-10;
    lambLB1((mx-lambLB)>10) = mx-10;
    mx = max(max([lambUB1 lambLB1]));
    mn = min(min([lambUB1 lambLB1]));
end
lambUBToPlot = lambUB(top:nDisp+top-1);
lambLBToPlot = lambLB(top:nDisp+top-1);

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
            lambUBToPlot = lambUBToPlot-x1;
            lambLBToPlot = lambLBToPlot-x1;
    end
end

% barh places the first element at the bottom. This is the opposite of
% what we'd like to do.
h = barh(handles.eRightPanel,lambUBToPlot(end:-1:1));
hh = get(h,'children');
set(h, 'buttondownfcn', {@eBar_downfcn})
set(hh,'cdatamapping', 'direct', 'facevertexcdata', ones(nDisp,1)) %all the same color
handles.eBarUB = h;
handles.ePatchUB = hh;

h = barh(handles.eLeftPanel,lambLBToPlot(end:-1:1));
hh = get(h,'children');
set(h, 'buttondownfcn', {@eBar_downfcn})
set(hh,'cdatamapping', 'direct', 'facevertexcdata', ones(nDisp,1))
handles.eBarLB = h;
handles.ePatchLB = hh;

% Clear the y axis tick labels
set(handles.eRightPanel,'YTickLabel',{});
set(handles.eLeftPanel,'YTickLabel',{});

switch handles.exp.scale
    case 'linear'
        xlabel(handles.eRightPanel,'\lambda')
        xlabel(handles.eLeftPanel,'\lambda')
    case 'log'
        xlabel(handles.eRightPanel,'log10(\lambda)');
        xlabel(handles.eLeftPanel,'log10(\lambda)');
end


% The top is nDisp, the bottom is nDisp-maxShow (since the listbox
% is tall enough for max show and the bars to line up with it).
ylim = [nDisp-maxShow+0.5 nDisp+0.5];
%set([handles.eRightPanel handles.eLeftPanel],'Xscale',handles.exp.scale);
set([handles.eRightPanel handles.eLeftPanel],'XLim',xlim,'YLim',ylim);
if strcmp(handles.exp.scale,'log')
    %now relabel ticks to correct numbers
    xticks = get(handles.eRightPanel,'xtick');
    set([handles.eRightPanel handles.eLeftPanel],'xticklabel',xticks+x1);
end


% Both panels should have the same XTickLabels. We want to remove the first
% element (which is zero) to around crowding.
%tmp = get(handles.eRightPanel,'XTickLabel');
%tmp(1,:) = ' '; %set to a blank string
%set([handles.eRightPanel handles.eLeftPanel],'XTickLabel',tmp);

function handles = plotPBars(handles,top)
% This function draws the bar graphs. The second input is the number of the
% listbox element that is at the top of the currently displayed elements.
% I.e., if the top element is the list box is visible, top = 1.

maxShow = 11; %this is the maximum number of bars we can fit in the display
nShow = length(handles.param.lambUB); %this is the number of experiments
nDisp = min(maxShow,nShow); %this is the number of bars we will draw

% Order the multipliers in the order they are to be displayed then just
% grab the ones corresponding to the visible listbox elements.
order = handles.param.order;
switch get(handles.maxradio,'Value')
    case 1
        switch handles.exp.scale
            case 'linear'
                lambUB = handles.param.lambUBR(order);
                lambLB = handles.param.lambUB(order);
            case 'log'
                lambUB = log10(handles.param.lambUBR(order));
                lambLB = log10(handles.param.lambUB(order));
        end
    case 0
        switch handles.exp.scale
            case 'linear'
                lambUB = handles.param.lambLBR(order);
                lambLB = handles.param.lambLB(order);
            case 'log'
                lambUB = log10(handles.param.lambLBR(order));
                lambLB = log10(handles.param.lambLB(order));
        end
end

% Some funky code to determing the x-axis limits. First determine the
% largest mult (among ALL so the scale is constant)

mx = max(max([lambUB lambLB]));
mn = min(min([lambUB lambLB]));
lambUB1=lambUB; lambLB1=lambLB;
if ceil(mx)-floor(mn)>10
    lambUB1((mx-lambUB)>10) = mx-10;
    lambLB1((mx-lambLB)>10) = mx-10;
    mx = max(max([lambUB1 lambLB1]));
    mn = min(min([lambUB1 lambLB1]));
end
lambUBToPlot = lambUB(top:nDisp+top-1);
lambLBToPlot = lambLB(top:nDisp+top-1);


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
            lambUBToPlot = lambUBToPlot-x1;
            lambLBToPlot = lambLBToPlot-x1;
    end
end

% barh places the first element at the bottom. This is the opposite of
% what we'd like to do.
h = barh(handles.pRightPanel,lambUBToPlot(end:-1:1));
hh = get(h,'children');
set(h, 'buttondownfcn', {@pBar_downfcn})
set(hh,'cdatamapping', 'direct', 'facevertexcdata', ones(nDisp,1)) %all the same color
handles.pBarUB = h;
handles.pPatchUB = hh;

h = barh(handles.pLeftPanel,lambLBToPlot(end:-1:1));
hh = get(h,'children');
set(h, 'buttondownfcn', {@pBar_downfcn})
set(hh,'cdatamapping', 'direct', 'facevertexcdata', ones(nDisp,1))
handles.pBarLB = h;
handles.pPatchLB = hh;

% Clear the y axis tick labels
set(handles.pRightPanel,'YTickLabel',{});
set(handles.pLeftPanel,'YTickLabel',{});

switch handles.exp.scale
    case 'linear'
        xlabel(handles.pRightPanel,'\lambda')
        xlabel(handles.pLeftPanel,'\lambda')
    case 'log'
        xlabel(handles.pRightPanel,'log10(\lambda)');
        xlabel(handles.pLeftPanel,'log10(\lambda)');
end



% The top is nDisp, the bottom is nDisp-maxShow (since the listbox
% is tall enough for max show and the bars to line up with it).
ylim = [nDisp-maxShow+0.5 nDisp+0.5];
%set([handles.pRightPanel handles.pLeftPanel],'Xscale',handles.exp.scale);
try
set([handles.pRightPanel handles.pLeftPanel],'XLim',xlim,'YLim',ylim);
catch
    return
end
if strcmp(handles.exp.scale,'log')
    %now relabel ticks to correct numbers
    xticks = get(handles.pRightPanel,'xtick');
    set([handles.pRightPanel handles.pLeftPanel],'xticklabel',xticks+x1);
end

% Both panels should have the same XTickLabels. We want to remove the first
% element (which is zero) to around crowding.
%tmp = get(handles.pRightPanel,'XTickLabel');
%tmp(1,:) = ' '; %set to a blank string
%set([handles.pRightPanel handles.pLeftPanel],'XTickLabel',tmp);


function colorEBars(handles,barnumber)
% barnumber is the number of bars from the top of the visible listbox
% display, with barnumber = 1 meaning the first element is visible

maxShow = 11; %this is the maximum number of bars we can fit in the display
nShow = length(handles.exp.lambUB); %this is the number of experiments
nDisp = min(maxShow,nShow); %this is the number of bars we will draw

if barnumber < 1 || barnumber > nDisp
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
set([handles.ePatchUB handles.ePatchLB], 'facevertexcdata', faceCData);

function colorPBars(handles,barnumber)
% barnumber is the number of bars from the top of the visible listbox
% display, with barnumber = 1 meaning the first element is visible

maxShow = 11; %this is the maximum number of bars we can fit in the display
nShow = length(handles.param.lambUB); %this is the number of experiments
nDisp = min(maxShow,nShow); %this is the number of bars we will draw

if barnumber < 1 || barnumber > nDisp
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
set([handles.pPatchUB handles.pPatchLB], 'facevertexcdata', faceCData);


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
        %do nothing
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
    catch
        %do nothing
    end
end

%Call with the previous element
PredictionViewer(handles.predObjCell,handles.elem-1);

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
PredictionViewer(handles.predObjCell,handles.elem+1);



% --- Executes on button press in dsetorder_radio.
function dsetorder_radio_Callback(hObject, eventdata, handles)
% hObject    handle to dsetorder_radio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of dsetorder_radio


% We will make it impossible to unclick a checkbox by clicking on it. To
% unclick it, the user must click SortLeft or SortRight
if get(hObject,'Value')
    % Reset the order variable to 1:m and display the names in this order.
    % The timer callback will notice of the order was changed and update the
    % plots.
    handles.exp.order = (1:length(handles.exp.lambUB))';
    set(handles.ExperimentList,'String',handles.ExpList);

    handles.param.order = (1:length(handles.param.lambUB))';
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

    hParam = handles.param.highlighted;
    if ~isempty(hParam)
        LBval = find(handles.param.order == hParam);
        set(handles.ParameterList,'Value',LBval)
        uicontrol(handles.ParameterList)
    end

    guidata(hObject,handles);
else
    set(hObject,'Value',1);
end

% --- Executes on button press in sortright_radio.
function sortright_radio_Callback(hObject, eventdata, handles)
% hObject    handle to sortright_radio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of sortright_radio

% We will make it impossible to unclick a checkbox by clicking on it. To
% unclick it, the user must click SortByNum or SortLeft

if get(hObject,'Value')
    % Set the order variable and display the names in this order. The timer
    % callback will notice of the order was changed and update the plots.
    switch get(handles.maxradio,'Value')
        case 1
            [temp idx1] = sort(handles.exp.lambUB+handles.exp.lambUBR,1,'ascend'); %#ok
            [temp idx2] = sort(handles.param.lambUB+handles.param.lambUBR,1,'ascend'); %#ok
        case 0
            [temp idx1] = sort(handles.exp.lambLB+handles.exp.lambLBR,1,'ascend'); %#ok
            [temp idx2] = sort(handles.param.lambLB+handles.param.lambLBR,1,'ascend'); %#ok
    end
    handles.exp.order = idx1;
    set(handles.ExperimentList,'String',handles.ExpList(idx1,:));

    handles.param.order = idx2;
    set(handles.ParameterList,'String',handles.ParamList(idx2,:));


    % The timer callback will notice of the order was changed and update the
    % plots.

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
        uicontrol(handles.ParameterList)
    end

    guidata(hObject,handles);
else
    set(hObject,'Value',1);
end


% --- Executes on button press in sortleft_radio.
function sortleft_radio_Callback(hObject, eventdata, handles)
% hObject    handle to sortleft_radio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of sortleft_radio

% We will make it impossible to unclick a checkbox by clicking on it. To
% unclick it, the user must click SortByNum or SortRight
if get(hObject,'Value')
    % Set the order variable, display the names in this order, and
    % reorder the bar graphs. The timmer callback will notice of the order
    % was changed and update the plots.
    switch get(handles.maxradio,'Value')
        case 1
            [temp idx1] = sort(handles.exp.lambUB+handles.exp.lambUBR,1,'descend'); %#ok
            [temp idx2] = sort(handles.param.lambUB+handles.param.lambUBR,1,'descend'); %#ok
        case 0
            [temp idx1] = sort(handles.exp.lambLB+handles.exp.lambLBR,1,'descend'); %#ok
            [temp idx2] = sort(handles.param.lambLB+handles.param.lambLBR,1,'descend'); %#ok
    end
    handles.exp.order = idx1;
    set(handles.ExperimentList,'String',handles.ExpList(idx1,:));

    handles.param.order = idx2;
    set(handles.ParameterList,'String',handles.ParamList(idx2,:));

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
function Exp_menu_Callback(hObject, eventdata, handles)
% hObject    handle to Exp_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Param_Menu_Callback(hObject, eventdata, handles)
% hObject    handle to Param_Menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over ParameterList.
function ParameterList_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to ParameterList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over ExperimentList.
function ExperimentList_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to ExperimentList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function ShowXML_menu_Callback(hObject, eventdata, handles)
% hObject    handle to ShowXML_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
paramnum = handles.param.highlighted;
uDstruct = handles.paramAss(paramnum).userData;
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











% --- Executes when selected object is changed in uipanel13.
function uipanel13_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel13 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)




topE = get(handles.ExperimentList,'ListBoxTop');
topP = get(handles.ParameterList,'ListBoxTop');

handles = plotEBars(handles,topE);
colorEBars(handles,handles.exp.highlighted-(topE-1));
handles = plotPBars(handles,topP);
colorPBars(handles,handles.param.highlighted-(topP-1));

guidata(hObject,handles);

if get(handles.sortright_radio,'Value')
    sortright_radio_Callback(handles.sortright_radio,[],handles);
elseif get(handles.sortleft_radio,'Value')
    sortleft_radio_Callback(handles.sortleft_radio,[],handles);
else
    dsetorder_radio_Callback(handles.sortright_radio,[],handles);
end

