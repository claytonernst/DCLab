function varargout = ModelImport(varargin)
% MODELIMPORT M-file for ModelImport.fig
%      MODELIMPORT, by itself, creates a new MODELIMPORT or raises the existing
%      singleton*.
%
%      H = MODELIMPORT returns the handle to a new MODELIMPORT or the handle to
%      the existing singleton*.
%
% ***  [idx name] = MODELIMPORT(existingNames,MACell) recieves a cell
%      array of existing dataset names that cannot be duplicated and a
%      cell array of ResponseModels from which the user will select the
%      one he wants Index is the index of the selected Model, while
%      name is a char array the user supplied as the name of the model.
%
%      MODELIMPORT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MODELIMPORT.M with the given input arguments.
%
%      MODELIMPORT('Property','Value',...) creates a new MODELIMPORT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ModelImport_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ModelImport_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ModelImport

% Last Modified by GUIDE v2.5 05-Feb-2009 14:00:12

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ModelImport_OpeningFcn, ...
                   'gui_OutputFcn',  @ModelImport_OutputFcn, ...
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

% --- Executes just before ModelImport is made visible.
function ModelImport_OpeningFcn(hObject, eventdata, handles, varargin) %#ok
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ModelImport (see VARARGIN)

%add fields to the handles structure
handles.existingNames = varargin{1};
handles.ModelAssertions = varargin{2};
handles.name = '';
handles.idx = [];

if isempty(handles.ModelAssertions)
  str = '';
else
  str = cell(length(handles.ModelAssertions),1);
  for i1 = 1:length(handles.ModelAssertions);
    if ~isempty(handles.ModelAssertions{i1}.name)
      str{i1} = handles.ModelAssertions{i1}.name;
    else
      str{i1} = num2str(i1);
    end
  end
end

%Tweak UI
set(handles.ModelLB,'String',str)
set(handles.figure1,'name','Add New Surrogate Model')

%Position UI
Units = get(handles.figure1,'Units');
set(handles.figure1,'Units','pixels');
curpos = get(handles.figure1,'Position');
scrnsize = get(0,'ScreenSize');
newpos = curpos;
newpos(1) = (scrnsize(3)-curpos(3))/2;
newpos(2) = (scrnsize(4)-curpos(4))/2;
set(handles.figure1,'Position',newpos,'Units',Units);

if length(handles.ModelAssertions) == 1 && ~isempty(handles.ModelAssertions)
  set(handles.ModelLB,'value',1); %select the first one automatically
  uicontrol(handles.ModelLB);
end

% Update handles structure
guidata(hObject, handles);

% Make the GUI modal
set(handles.figure1,'WindowStyle','modal')

% UIWAIT makes DatasetUnitImport wait for user response (see UIRESUME)
uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = ModelImport_OutputFcn(hObject, eventdata, handles) %#ok
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

varargout{1} = handles.idx;
varargout{2} = handles.name;
% The figure can be deleted now
delete(handles.figure1);

% --- Executes on selection change in ModelLB.
function ModelLB_Callback(hObject, eventdata, handles) %#ok

%check for double click
% if strcmp(get(handles.figure1,'SelectionType'),'open')
%   LBval = get(hObject,'Value');
%   %do nothing if the user double clicked and dragged at the same time
%   if length(LBval) == 1
%     hand = get(handles.ModelAssertions{LBval},'dispCallback');
%     %if a callback is available, invoke it
%     if ~isempty(hand)
%       hand(handles.unitArray(LBval));
%     end
%   end
% end

function ModelLB_CreateFcn(hObject, eventdata, handles) %#ok
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in Load.
function Load_Callback(hObject, eventdata, handles) %#ok
handles.idx = get(handles.ModelLB,'Value');

if isempty(handles.idx)
  errordlg('No Model was selected','Please retry','modal')
  return
end
%handles.name = GetNewName({'Please Supply Model Name'},'User Input',handles.existingNames);
handles.name = handles.ModelAssertions{handles.idx}.name;
if isempty(handles.name)
  return
else
  guidata(hObject,handles)
  uiresume(handles.figure1)
end

% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles) %#ok
if isequal(get(handles.figure1, 'waitstatus'), 'waiting')
    Cancel_Callback(handles.Cancel,eventdata,handles)
else
    % The GUI is no longer waiting, something bombed, just close it
    delete(handles.figure1);
end

% --- Executes on button press in Cancel.
function Cancel_Callback(hObject, eventdata, handles) %#ok
handles.idx = [];
handles.name = '';
guidata(hObject,handles)
uiresume(handles.figure1)



% --------------------------------------------------------------------
function model_menu_Callback(hObject, eventdata, handles)
% hObject    handle to ShowXML_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
LBval = get(handles.ModelLB,'Value');

%do nothing if the user double clicked and dragged at the same time
if length(LBval) == 1
    uDstruct = handles.ModelAssertions{LBval}.userData;
    pid = uDstruct.modelID;
    if iscell(pid) && isscalar(pid)
        uDstruct.callback2show(pid{1});
    elseif iscell(pid)
        uDstruct.callback2show2(pid{1},pid{2});
    elseif ischar(pid)
        uDstruct.callback2show(pid);
    end
end

% --------------------------------------------------------------------
function target_menu_Callback(hObject, eventdata, handles)
% hObject    handle to ShowXML_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
LBval = get(handles.ModelLB,'Value');

uDstruct = handles.ModelAssertions{LBval}.userData;
pid = uDstruct.targetID;
if iscell(pid) && isscalar(pid)
    uDstruct.callback2show(pid{1});
elseif iscell(pid)
    uDstruct.callback2show(pid{1});
elseif ischar(pid)
    uDstruct.callback2show(pid);
end


% --------------------------------------------------------------------
function target_links_menu_Callback(hObject, eventdata, handles)
% hObject    handle to target_links_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
LBval = get(handles.ModelLB,'Value');

uDstruct = handles.ModelAssertions{LBval}.userData;
pid = uDstruct.targetLinks;
if iscell(pid) && isscalar(pid)
    uDstruct.callback2show(pid{1});
elseif iscell(pid)
    uDstruct.callback2show(pid{1});
elseif ischar(pid)
    uDstruct.callback2show(pid);
end


% --------------------------------------------------------------------
function Exp_Menu_Callback(hObject, eventdata, handles)
% hObject    handle to Exp_Menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function ShowXML_menu_Callback(hObject, eventdata, handles)
% hObject    handle to ShowXML_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
