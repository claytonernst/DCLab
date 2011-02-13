function varargout = DatasetUnitImport(varargin)
% DATASETUNITIMPORT M-file for DatasetUnitImport.fig
%      DATASETUNITIMPORT, by itself, creates a new DATASETUNITIMPORT or raises the existing
%      singleton*.
%
%      H = DATASETUNITIMPORT returns the handle to a new DATASETUNITIMPORT or the handle to
%      the existing singleton*.
%
%      DATASETUNITIMPORT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DATASETUNITIMPORT.M with the given input arguments.
%
% *** [idx name] = DATASETUNITIMPORT(existingNames,dsetUnits) recieves a cell array of
%      existing dataset names that cannot be duplicated and a
%      multidimensional array of dataset units from which the user will
%      select those to include in his dataset. Index is a vector of indices
%      of the selected strings, while name is a char array the user
%      supplied as the name of the dataset.
%
%      DATASETUNITIMPORT('Property','Value',...) creates a new DATASETUNITIMPORT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before DatasetUnitImport_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to DatasetUnitImport_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help DatasetUnitImport

% Last Modified by GUIDE v2.5 05-Feb-2009 12:51:28

%pragmas
%#function Cancel_Callback DatasetUnitImport_OpeningFcn
%#function DatasetUnitImport_OutputFcn figure1_CloseRequestFcn 
%#function LoadSelected_Callback SelectAll_Callback TargetListing_Callback
%#function TargetListing_CreateFcn

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @DatasetUnitImport_OpeningFcn, ...
                   'gui_OutputFcn',  @DatasetUnitImport_OutputFcn, ...
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

% --- Executes just before DatasetUnitImport is made visible.
function DatasetUnitImport_OpeningFcn(hObject, eventdata, handles, varargin) %#ok
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to DatasetUnitImport (see VARARGIN)

%add fields to the handles structure
handles.existingNames = varargin{1};
handles.unitArray = varargin{2};
handles.name = '';
handles.idx = [];

%tweak ui
str = {handles.unitArray.name}';
set(handles.TargetListing,'String',str)
set(handles.TargetListing,'value',[]) %initially nothing should be selected
set(handles.figure1,'name','Add New Dataset')

%position ui
Units = get(handles.figure1,'Units');
set(handles.figure1,'Units','pixels');
curpos = get(handles.figure1,'Position');
scrnsize = get(0,'ScreenSize');
newpos = curpos;
newpos(1) = (scrnsize(3)-curpos(3))/2;
newpos(2) = (scrnsize(4)-curpos(4))/2;
set(handles.figure1,'Position',newpos,'Units',Units);

% Update handles structure
guidata(hObject, handles);

% Make the GUI modal
set(handles.figure1,'WindowStyle','modal')

% UIWAIT makes DatasetUnitImport wait for user response (see UIRESUME)
uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = DatasetUnitImport_OutputFcn(hObject, eventdata, handles) %#ok
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

varargout{1} = handles.idx;
varargout{2} = handles.name;
% The figure can be deleted now
delete(handles.figure1);

function TargetListing_Callback(hObject, eventdata, handles) %#ok

%check for double click
% if strcmp(get(handles.figure1,'SelectionType'),'open')
%   LBval = get(hObject,'Value');
%   
%   %do nothing if the user double clicked and dragged at the same time
%   if length(LBval) == 1
%     hand = handles.unitArray(LBval).dispCallback;
%     %if a callback is available, invoke it
%     if ~isempty(hand)
%       hand(handles.unitArray(LBval));
%     end
%   end
% end

% --- Executes during object creation, after setting all properties.
function TargetListing_CreateFcn(hObject, eventdata, handles) %#ok
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function SelectAll_Callback(hObject, eventdata, handles) %#ok
set(handles.TargetListing,'Value',1:length(handles.unitArray));
uicontrol(handles.TargetListing)

function LoadSelected_Callback(hObject, eventdata, handles) %#ok
handles.idx = get(handles.TargetListing,'Value');
if isempty(handles.idx)
  errordlg('No DatasetUnits have been selected','Please retry','modal')
  return
end
handles.name = GetNewName({'Please Supply Dataset Name'},'User Input',handles.existingNames);

if isempty(handles.name)
  return
else
  guidata(hObject,handles)
  uiresume(handles.figure1)
end

function figure1_CloseRequestFcn(hObject, eventdata, handles) %#ok
if isequal(get(handles.figure1, 'waitstatus'), 'waiting')
    Cancel_Callback(handles.Cancel,eventdata,handles)
else
    % The GUI is no longer waiting, something bombed, just close it
    delete(handles.figure1);
end

function Cancel_Callback(hObject, eventdata, handles) %#ok
handles.idx = [];
handles.name = '';
guidata(hObject,handles)
uiresume(handles.figure1)


% --------------------------------------------------------------------
function target_menu_Callback(hObject, eventdata, handles)
% hObject    handle to ShowXML_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
LBval = get(handles.TargetListing,'Value');

%do nothing if the user double clicked and dragged at the same time
if length(LBval) == 1
    uDstruct = handles.unitArray(LBval).ResponseObservation.userData;
    pid = uDstruct.targetID;
    if iscell(pid) && isscalar(pid)
        uDstruct.callback2show(pid{1});
    elseif iscell(pid)
        uDstruct.callback2show(pid{1});
    elseif ischar(pid)
        uDstruct.callback2show(pid);
    end
elseif length(LBval)>1
    errordlg('You may highlight many to build a Dataset, but please only select one before viewing Target/Model information.','Too many Highlighted','Modal')
else
    errordlg('Please select exactly one target before viewing Target/Model information.  You may select many before creating a Dataset.','None highlighted','Modal')
end


% --------------------------------------------------------------------
function target_links_menu_Callback(hObject, eventdata, handles)
% hObject    handle to target_links_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
LBval = get(handles.TargetListing,'Value');

%do nothing if the user double clicked and dragged at the same time
if length(LBval) == 1
    uDstruct = handles.unitArray(LBval).ResponseObservation.userData;
    pid = uDstruct.targetLinks;
    if iscell(pid) && isscalar(pid)
        uDstruct.callback2show(pid{1});
    elseif iscell(pid)
        uDstruct.callback2show(pid{1});
    elseif ischar(pid)
        uDstruct.callback2show(pid);
    end
elseif length(LBval)>1
    errordlg('You may highlight many to build a Dataset, but please only select one before viewing Target/Model information.','Too many Highlighted','Modal')
else
    errordlg('Please select exactly one target before viewing Target/Model information.  You may select many before creating a Dataset.','None highlighted','Modal')
end


% --------------------------------------------------------------------
function model_menu_Callback(hObject, eventdata, handles)
% hObject    handle to model_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
LBval = get(handles.TargetListing,'Value');

%do nothing if the user double clicked and dragged at the same time
if length(LBval) == 1
    uDstruct = handles.unitArray(LBval).ResponseModel.userData;
    pid = uDstruct.modelID;
    if iscell(pid) && isscalar(pid)
        uDstruct.callback2show(pid{1});
    elseif iscell(pid)
        uDstruct.callback2show2(pid{1},pid{2});
    elseif ischar(pid)
        uDstruct.callback2show(pid);
    end
elseif length(LBval)>1
    errordlg('You may highlight many to build a Dataset, but please only select one before viewing Target/Model information.','Too many Highlighted','Modal')
else
    errordlg('Please select exactly one target before viewing Target/Model information.  You may select many before creating a Dataset.','None highlighted','Modal')
end


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over TargetListing.
function TargetListing_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to TargetListing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on key press with focus on TargetListing and no controls selected.
function TargetListing_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to TargetListing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Exp_menu_Callback(hObject, eventdata, handles)
% hObject    handle to ShowXML_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function ShowXML_menu_Callback(hObject, eventdata, handles)
% hObject    handle to ShowXML_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


