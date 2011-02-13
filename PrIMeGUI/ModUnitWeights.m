function varargout = ModUnitWeights(varargin)
% MODUNITWEIGHTS M-file for ModUnitWeights.fig
%      MODUNITWEIGHTS, by itself, creates a new MODUNITWEIGHTS or raises the existing
%      singleton*.
%
%      H = MODUNITWEIGHTS returns the handle to a new MODUNITWEIGHTS or the handle to
%      the existing singleton*.
%
%      MODUNITWEIGHTS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MODUNITWEIGHTS.M with the given input arguments.
%
%      MODUNITWEIGHTS('Property','Value',...) creates a new MODUNITWEIGHTS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ModUnitWeights_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ModUnitWeights_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
%
% Description of handles structure
%
%      ====graphics handles=====
%
%        ====additional data handles=====
%
%         existingNames: qx1 cell of chars
%                curIdx: scalar, index telling which element of existingNames/DUcell is currently viewable
%                DUcell: qx1 cell of DatasetUnit objects
%                DUinfo: qx1 struct with fields .def_d, .def_uncVect,
%                         .isdef, .curwt, .defwt. Structure generally
%                          tells us the default data/unc, and which are
%                          default for each of the q DatasetUnit objects.
%                          Here concerned only with the .curwt field
%               org_unc: {77x1 cell}


% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ModUnitWeights_OpeningFcn, ...
                   'gui_OutputFcn',  @ModUnitWeights_OutputFcn, ...
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

% --- Executes just before ModUnitWeights is made visible.
function ModUnitWeights_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ModUnitWeights (see VARARGIN)

%add fields to the handles structure
handles.existingNames = varargin{1};
handles.curIdx = varargin{2};
handles.DUcell = varargin{3};
handles.DUinfo = varargin{4};

handles.org_uncVect = cell2mat({handles.DUcell{handles.curIdx}.observationUncertaintyPlusMinus}');
handles.orgwt = handles.DUinfo(handles.curIdx).curwt;
handles.curwt = handles.orgwt;
handles.isdef = 2./diff(handles.org_uncVect,[],2)==handles.orgwt;
handles.goodwt = true;
handles.curSelection = 1;

% Update handles structure
guidata(hObject, handles);

%tweak ui
set(handles.figure1,'name',['Dataset: ' handles.existingNames{handles.curIdx}])
set(handles.DatasetUnitLB,'value',handles.curSelection) %initially select the first

%listbox string
str = {handles.DUcell{handles.curIdx}.name}';
if ~all(handles.isdef)
  pad = repmat({blanks(3)},size(str));
  pad(~handles.isdef) = {'*! '};
  str = cellstr([strvcat(pad{:}) strvcat(str{:})]); %#ok
  set(handles.DatasetUnitLB,'String',str)
else
  set(handles.DatasetUnitLB,'String',str)
end

%Update the edit displays
DatasetUnitLB_Callback(handles.DatasetUnitLB,eventdata,handles)

%position ui
Units = get(handles.figure1,'Units');
set(handles.figure1,'Units','pixels');
curpos = get(handles.figure1,'Position');
scrnsize = get(0,'ScreenSize');
newpos = curpos;
newpos(1) = (scrnsize(3)-curpos(3))/2;
newpos(2) = (scrnsize(4)-curpos(4))/2;
set(handles.figure1,'Position',newpos,'Units',Units);

%hide axis
hideAxis(handles)

% Make the GUI modal
%set(handles.figure1,'WindowStyle','modal')

% UIWAIT makes DatasetUnitImport wait for user response (see UIRESUME)
uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = ModUnitWeights_OutputFcn(hObject, eventdata, handles) %#ok
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

varargout{1} = handles.DUinfo;
%varargout{1} = handles.existingNames;
%varargout{2} = handles.curIdx;
%varargout{3} = handles.DUcell;
%varargout{4} = handles.DUinfo;

% The figure can be deleted now
delete(handles.figure1);

%===Create functions===
function UncLBEdit_CreateFcn(hObject, eventdata, handles) %#ok
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function UncUBEdit_CreateFcn(hObject, eventdata, handles) %#ok
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function UncEdit_CreateFcn(hObject, eventdata, handles) %#ok
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function WeightEdit_CreateFcn(hObject, eventdata, handles) %#ok
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function DatasetUnitLB_CreateFcn(hObject, eventdata, handles) %#ok
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%===Callback functions===
function WeightEdit_Callback(hObject, eventdata, handles) %#ok

idx = get(handles.DatasetUnitLB,'Value');
newwt = str2double(get(hObject,'String'));

if isnan(newwt) || newwt <= 0
  set(hObject,'BackgroundColor','red')
  
  %The listbox callback will disable it
  set([handles.SavePB handles.showDistPB],'enable','off')
  handles.goodwt = false;
else
  set(hObject,'BackgroundColor','white')
  handles.curwt(idx) = newwt;
  set([handles.DatasetUnitLB handles.SavePB handles.showDistPB],'enable','on')
  handles.goodwt = true;
  
  %Check if everything is default, and compare this to the old value. If it
  %changed, we need to update the display 
  isdef = newwt == 2/diff(handles.org_uncVect(idx,:));
  if isdef~=handles.isdef(idx)
    handles.isdef(idx) = isdef;
    get(handles.DatasetUnitLB,'Value');
    str = {handles.DUcell{handles.curIdx}.name}';
    if ~all(handles.isdef)
      pad = repmat({blanks(3)},size(str));
      pad(~handles.isdef) = {'*! '};
      str = cellstr([strvcat(pad{:}) strvcat(str{:})]); %#ok
      set(handles.DatasetUnitLB,'String',str)
    else
      set(handles.DatasetUnitLB,'String',str)
    end
  end
  if strcmp(get(handles.distAxis,'Visible'),'on')
    bar(handles.distAxis,handles.curwt)
  end

end
guidata(hObject,handles)

function RevertCurrentPB_Callback(hObject, eventdata, handles) %#ok
idx = get(handles.DatasetUnitLB,'Value');

% Is there anything to do?
if ~handles.isdef(idx) || ~handles.goodwt
  
  handles.curwt(idx) = 2/diff(handles.org_uncVect(idx,:));  

  if ~handles.goodwt
    handles.goodwt = true;
    set(handles.WeightEdit,'BackgroundColor','white')
    set([handles.DatasetUnitLB handles.SavePB handles.showDistPB],'enable','on')
  end
  if ~handles.isdef(idx)
    %if it wasn't default, make it so and fix the strings
    handles.isdef(idx) = true;
    str = {handles.DUcell{handles.curIdx}.name}';
    if ~all(handles.isdef)
      pad = repmat({blanks(3)},size(str));
      pad(~handles.isdef) = {'*! '};
      str = cellstr([strvcat(pad{:}) strvcat(str{:})]); %#ok
      set(handles.DatasetUnitLB,'String',str)
    else
      set(handles.DatasetUnitLB,'String',str)
    end
  end
  
  if strcmp(get(handles.distAxis,'Visible'),'on')
    bar(handles.distAxis,handles.curwt)
  end
  
  %save and call the listbox callback to refresh displays
  guidata(hObject,handles)
  DatasetUnitLB_Callback(handles.DatasetUnitLB, eventdata, handles)
end

function RevertAllPB_Callback(hObject, eventdata, handles) %#ok

if ~handles.goodwt
  handles.goodwt = true;
  set(handles.WeightEdit,'BackgroundColor','white')
  set([handles.DatasetUnitLB handles.SavePB handles.showDistPB],'enable','on')
end

handles.curwt = 2./diff(handles.org_uncVect,[],2);
handles.isdef(:) = true;
str = {handles.DUcell{handles.curIdx}.name}';
set(handles.DatasetUnitLB,'String',str)

if strcmp(get(handles.distAxis,'Visible'),'on')
  bar(handles.distAxis,handles.curwt)
end

%call the listbox callback to refresh displays
guidata(hObject,handles)
DatasetUnitLB_Callback(handles.DatasetUnitLB, eventdata, handles)

function DatasetUnitLB_Callback(hObject, eventdata, handles) %#ok

% if the current edit box entry is bad, listbox should reset to previous
% value and disable itself
if ~handles.goodwt
  set(hObject,'Value',handles.curSelection)
  set(hObject,'enable','off')
  return
end

val = get(hObject,'Value');
handles.curSelection = val;

% String in wt edit box
set(handles.WeightEdit,'String',num2str(handles.curwt(val)))

% Values and visibility of uncertainty edit boxes
set(handles.UncLBEdit,'String',num2str(handles.org_uncVect(val,1)))
set(handles.UncUBEdit,'String',num2str(handles.org_uncVect(val,2)))
if -handles.org_uncVect(val,1) == handles.org_uncVect(val,2)
  set([handles.UncLBEdit handles.UncUBEdit],'Visible','off')
  set(handles.UncEdit,'Visible','on','String',num2str(handles.org_uncVect(val,2)))
else
  set(handles.UncEdit,'Visible','off','String','inf')
  set([handles.UncLBEdit handles.UncUBEdit],'Visible','on')
end

%check for double click, which means we should display further details
if strcmp(get(handles.figure1,'SelectionType'),'open')
  LBval = get(hObject,'Value');
  
  %do nothing if the user double clicked and dragged at the same time
  if length(LBval) == 1
    hand = handles.DUcell{handles.curIdx}(LBval).dispCallback;
    %if a callback is available, invoke it
    if ~isempty(hand)
      hand(handles.DUcell{handles.curIdx}(LBval));
    end
  end
end

guidata(hObject,handles)

function SavePB_Callback(hObject, eventdata, handles) %#ok
%Save changes in DUinfo
handles.DUinfo(handles.curIdx).curwt = handles.curwt;
guidata(hObject,handles)
uiresume(handles.figure1)

function CancelPB_Callback(hObject, eventdata, handles) %#ok
%replace values with the originals, update the object, and return
handles.DUinfo(handles.curIdx).curwt = handles.orgwt;
guidata(hObject,handles)
uiresume(handles.figure1)

function hideAxis(handles)
%every thing we're playing with should be in characters, i won't check to
%save time...

set(handles.distAxis,'Visible','off')

%move a bunch of stuff down
offset = 9; %(in characters)

toMove = [handles.DatasetUnitLB handles.du_text handles.obsInfoPanel ...
    handles.showDistPB handles.RevertAllPB handles.SetAllOnePB];
for i1 = 1:length(toMove)
  pos = get(toMove(i1),'Position');
  pos(2) = pos(2)-offset;
  set(toMove(i1),'Position',pos)
end

%resize main panel
panelPos = get(handles.mainpanel,'Position');
panelPos(4) = panelPos(4)-offset;
set(handles.mainpanel,'Position',panelPos)

%now shorten the gui and move it up to be visually pleasing
curpos = get(handles.figure1,'Position');
curpos(2) = curpos(2)+offset;
curpos(4) = curpos(4)-offset;
set(handles.figure1,'Position',curpos)

drawnow

function showAxis(handles)
%every thing we're playing with should be in characters, i won't check to
%save time...

set(handles.distAxis,'Visible','on')

%move a bunch of stuff down
offset = 9; %(in characters)

toMove = [handles.DatasetUnitLB handles.du_text handles.obsInfoPanel ...
    handles.showDistPB handles.RevertAllPB handles.SetAllOnePB];
for i1 = 1:length(toMove)
  pos = get(toMove(i1),'Position');
  pos(2) = pos(2)+offset;
  set(toMove(i1),'Position',pos)
end

%resize main panel
panelPos = get(handles.mainpanel,'Position');
panelPos(4) = panelPos(4)+offset;
set(handles.mainpanel,'Position',panelPos)

%now shorten the gui and move it up to be visually pleasing
curpos = get(handles.figure1,'Position');
curpos(2) = curpos(2)-offset;
curpos(4) = curpos(4)+offset;
set(handles.figure1,'Position',curpos)

drawnow

function figure1_CloseRequestFcn(hObject, eventdata, handles) %#ok

curobj = get(handles.figure1,'CurrentObject');
if curobj==handles.WeightEdit
  %the user hitting exit does NOT trigger the edit callbacks. thus we call
  %small gui with a uicontrol in the center of the screen to flush all call
  %edit callbacks.
  scrnsz = get(0,'screensize');
  h = dummy('Position',[scrnsz(3) scrnsz(4) 1 1],'Units','pixels');
  delete(h)
  %update handles
  handles = guidata(hObject);
end

if isequal(get(handles.figure1, 'waitstatus'), 'waiting')
  if ~handles.goodwt
    %don't exit: there are problems
  elseif ~isequal(handles.DUinfo(handles.curIdx).curwt,handles.curwt)
    ButtonName = savedlg_nosaveas('string',['Do you want to save modified weights for Dataset ' handles.existingNames{handles.curIdx}],'title','User Input');
    switch ButtonName,
     case 'Save Changes and Return'
      SavePB_Callback(handles.SavePB,eventdata,handles)
     case 'Discard Changes and Return'
      CancelPB_Callback(handles.CancelPB,eventdata,handles)
     case 'Cancel'
      return
     otherwise
      %keyboard
    end % switch
  else
    % just close it
    CancelPB_Callback(handles.CancelPB,eventdata,handles)
  end

else
    % The GUI is no longer waiting, something bombed, just close it
    delete(handles.figure1);
end

% --- Executes on button press in showDistPB.
function showDistPB_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to showDistPB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

switch get(hObject,'String')
  case 'View Weight Distribution'
    set(hObject,'String','Hide Weight Distribution')
    showAxis(handles)
    bar(handles.distAxis,handles.curwt)
  case 'Hide Weight Distribution'
    set(hObject,'String','View Weight Distribution')
    hideAxis(handles)
  otherwise
    disp('show never get here, check code')
end

% --- Executes on button press in SetAllOnePB.
function SetAllOnePB_Callback(hObject, eventdata, handles)
% hObject    handle to SetAllOnePB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~handles.goodwt
  handles.goodwt = true;
  set(handles.WeightEdit,'BackgroundColor','white')
  set([handles.DatasetUnitLB handles.SavePB handles.showDistPB],'enable','on')
end

n = length(handles.org_uncVect(:,1));
handles.curwt = ones(n,1);
handles.isdef(:) = true;
str = {handles.DUcell{handles.curIdx}.name}';
set(handles.DatasetUnitLB,'String',str)

if strcmp(get(handles.distAxis,'Visible'),'on')
  bar(handles.distAxis,handles.curwt)
end

%call the listbox callback to refresh displays
guidata(hObject,handles)
DatasetUnitLB_Callback(handles.DatasetUnitLB, eventdata, handles)


% --------------------------------------------------------------------
function Exp_menu_Callback(hObject, eventdata, handles)
% hObject    handle to Exp_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function ShowXML_menu_Callback(hObject, eventdata, handles)
% hObject    handle to ShowXML_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function target_menu_Callback(hObject, eventdata, handles)
% hObject    handle to target_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
LBval = get(handles.DatasetUnitLB,'Value');
%do nothing if the user double clicked and dragged at the same time
if length(LBval) == 1
    uDstruct = handles.DUcell{handles.curIdx}(LBval).ResponseObservation.userData;
    pid = uDstruct.targetID;
    if iscell(pid) && isscalar(pid)
        uDstruct.callback2show(pid{1});
    elseif iscell(pid)
        uDstruct.callback2show(pid{1});
    elseif ischar(pid)
        uDstruct.callback2show(pid);
    end
end

% --------------------------------------------------------------------
function target_links_menu_Callback(hObject, eventdata, handles)
% hObject    handle to target_links_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
LBval = get(handles.DatasetUnitLB,'Value');

%do nothing if the user double clicked and dragged at the same time
if length(LBval) == 1
    uDstruct = handles.DUcell{handles.curIdx}(LBval).ResponseObservation.userData;
    pid = uDstruct.targetLinks;
    if iscell(pid) && isscalar(pid)
        uDstruct.callback2show(pid{1});
    elseif iscell(pid)
        uDstruct.callback2show(pid{1});
    elseif ischar(pid)
        uDstruct.callback2show(pid);
    end
end

% --------------------------------------------------------------------
function model_menu_Callback(hObject, eventdata, handles)
% hObject    handle to model_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
LBval = get(handles.DatasetUnitLB,'Value');

%do nothing if the user double clicked and dragged at the same time
if length(LBval) == 1
    uDstruct = handles.DUcell{handles.curIdx}(LBval).ResponseModel.userData;
    pid = uDstruct.modelID;
    if iscell(pid) && isscalar(pid)
        uDstruct.callback2show(pid{1});
    elseif iscell(pid)
        uDstruct.callback2show2(pid{1},pid{2});
    elseif ischar(pid)
        uDstruct.callback2show(pid);
    end
end

