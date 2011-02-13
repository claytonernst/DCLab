function varargout = PrIMeDatasetImport(varargin)

% PRIMEDATASETIMPORT M-file for PrIMeDatasetImport.fig

%      PRIMEDATASETIMPORT, by itself, creates a new PRIMEDATASETIMPORT or raises the existing

%      singleton*.

%

%      H = PRIMEDATASETIMPORT returns the handle to a new PRIMEDATASETIMPORT or the handle to

%      the existing singleton*.

%

%      PRIMEDATASETIMPORT('CALLBACK',hObject,eventData,handles,...) calls the local

%      function named CALLBACK in PRIMEDATASETIMPORT.M with the given input arguments.

%

%      PRIMEDATASETIMPORT('Property','Value',...) creates a new PRIMEDATASETIMPORT or raises the

%      existing singleton*.  Starting from the left, property value pairs are

%      applied to the GUI before PrIMeDatasetImport_OpeningFunction gets called.  An

%      unrecognized property name or invalid value makes property application

%      stop.  All inputs are passed to PrIMeDatasetImport_OpeningFcn via varargin.

%

%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one

%      instance to run (singleton)".

%

% See also: GUIDE, GUIDATA, GUIHANDLES



% Edit the above text to modify the response to help PrIMeDatasetImport



% Last Modified by GUIDE v2.5 22-Feb-2006 01:44:32



% Begin initialization code - DO NOT EDIT

gui_Singleton = 1;

gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @PrIMeDatasetImport_OpeningFcn, ...
                   'gui_OutputFcn',  @PrIMeDatasetImport_OutputFcn, ...
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





% --- Executes just before PrIMeDatasetImport is made visible.

function PrIMeDatasetImport_OpeningFcn(hObject, eventdata, handles, varargin)

% This function has no output args, see OutputFcn.

% hObject    handle to figure

% eventdata  reserved - to be defined in a future version of MATLAB

% handles    structure with handles and user data (see GUIDATA)

% varargin   command line arguments to PrIMeDatasetImport (see VARARGIN)



handles.Dataset = [];

handles.name = '';

handles.success = 0; 

handles.existingNames = varargin{1};

handles.plugin = varargin{2};



set(handles.TargetListing,'String',list_available_targets(handles.plugin))

set(handles.TargetListing,'value',[])



% Update handles structure

guidata(hObject, handles);



% UIWAIT makes PrIMeDatasetImport wait for user response (see UIRESUME)

% uiwait(handles.figure1);



% --- Outputs from this function are returned to the command line.

function varargout = PrIMeDatasetImport_OutputFcn(hObject, eventdata, handles) 

% varargout  cell array for returning output args (see VARARGOUT);

% hObject    handle to figure

% eventdata  reserved - to be defined in a future version of MATLAB

% handles    structure with handles and user data (see GUIDATA)



uiwait(handles.figure1);

handles = guidata(hObject); %load most recent handles structure

varargout{1} = handles.Dataset;

varargout{2} = handles.name;

varargout{3} = handles.success;

delete(handles.figure1);



% --- Executes on selection change in TargetListing.

function TargetListing_Callback(hObject, eventdata, handles)

% hObject    handle to TargetListing (see GCBO)

% eventdata  reserved - to be defined in a future version of MATLAB

% handles    structure with handles and user data (see GUIDATA)



% Hints: contents = get(hObject,'String') returns TargetListing contents as cell array

%        contents{get(hObject,'Value')} returns selected item from TargetListing



% --- Executes during object creation, after setting all properties.

function TargetListing_CreateFcn(hObject, eventdata, handles)

% hObject    handle to TargetListing (see GCBO)

% eventdata  reserved - to be defined in a future version of MATLAB

% handles    empty - handles not created until after all CreateFcns called



% Hint: listbox controls usually have a white background on Windows.

%       See ISPC and COMPUTER.

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))

    set(hObject,'BackgroundColor','white');

end





% --- Executes on button press in LoadAll.

function LoadAll_Callback(hObject, eventdata, handles)

% hObject    handle to LoadAll (see GCBO)

% eventdata  reserved - to be defined in a future version of MATLAB

% handles    structure with handles and user data (see GUIDATA)



handles.Dataset = upload_all_targets(handles.plugin);

handles.name = DatasetName(handles.existingNames);

if isempty(handles.name)

  return

else

  handles.success = 1;

  guidata(hObject,handles)

  figure1_CloseRequestFcn(handles.figure1,eventdata,handles)

  guidata(hObject,handles)

end



% --- Executes on button press in LoadSelected.



function LoadSelected_Callback(hObject, eventdata, handles)



if isempty(get(handles.TargetListing,'value'))

  errordlg('No Targets have been selected')

  return

else

  TrgList = get(handles.TargetListing,'String');

  handles.Dataset = upload_requested_targets(handles.plugin,TrgList(get(handles.TargetListing,'Value')));

  handles.name = DatasetName(handles.existingNames);

  if isempty(handles.name)

    return

  else

    handles.success = 1;

    guidata(hObject,handles)

    figure1_CloseRequestFcn(handles.figure1,eventdata,handles)

    guidata(hObject,handles)

  end

end



% --- Executes when user attempts to close figure1.

function figure1_CloseRequestFcn(hObject, eventdata, handles)

uiresume(handles.figure1)







