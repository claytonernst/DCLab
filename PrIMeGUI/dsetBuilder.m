function varargout = dsetBuilder(varargin)

% DSETBUILDER M-file for dsetBuilder.fig

%      DSETBUILDER, by itself, creates a new DSETBUILDER or raises the existing

%      singleton*.

%

%      H = DSETBUILDER returns the handle to a new DSETBUILDER or the handle to

%      the existing singleton*.

%

%      DSETBUILDER('CALLBACK',hObject,eventData,handles,...) calls the local

%      function named CALLBACK in DSETBUILDER.M with the given input arguments.

%

%      DSETBUILDER('Property','Value',...) creates a new DSETBUILDER or raises the

%      existing singleton*.  Starting from the left, property value pairs are

%      applied to the GUI before dsetBuilder_OpeningFunction gets called.  An

%      unrecognized property name or invalid value makes property application

%      stop.  All inputs are passed to dsetBuilder_OpeningFcn via varargin.

%

%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one

%      instance to run (singleton)".

%

% See also: GUIDE, GUIDATA, GUIHANDLES



% Edit the above text to modify the response to help dsetBuilder



% Last Modified by GUIDE v2.5 06-Apr-2006 02:08:05



% Begin initialization code - DO NOT EDIT

gui_Singleton = 1;

gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @dsetBuilder_OpeningFcn, ...
                   'gui_OutputFcn',  @dsetBuilder_OutputFcn, ...
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





% --- Executes just before dsetBuilder is made visible.

function dsetBuilder_OpeningFcn(hObject, eventdata, handles, varargin)

% This function has no output args, see OutputFcn.

% hObject    handle to figure

% eventdata  reserved - to be defined in a future version of MATLAB

% handles    structure with handles and user data (see GUIDATA)

% varargin   command line arguments to dsetBuilder (see VARARGIN)



% Choose default command line output for dsetBuilder

if nargin > 3

  handles.existingNames = varargin{1};

else

  handles.existingNames = {''};

end

handles.var.name = '';

handles.var.value = [];

handles.str = struct([]);

handles.success = 0;

% Update handles structure

guidata(hObject, handles);



% UIWAIT makes dsetBuilder wait for user response (see UIRESUME)

% uiwait(handles.figure1);





% --- Outputs from this function are returned to the command line.

function varargout = dsetBuilder_OutputFcn(hObject, eventdata, handles) 

% varargout  cell array for returning output args (see VARARGOUT);

% hObject    handle to figure

% eventdata  reserved - to be defined in a future version of MATLAB

% handles    structure with handles and user data (see GUIDATA)



% Get default command line output from handles structure

uiwait(handles.figure1);

handles = guidata(hObject); %load most recent handles structure

varargout{1} = handles.var.value;

varargout{2} = handles.var.name;

varargout{3} = handles.success;

delete(handles.figure1);



% --- Executes during object creation, after setting all properties.

function disp_CreateFcn(hObject, eventdata, handles)

% hObject    handle to disp (see GCBO)

% eventdata  reserved - to be defined in a future version of MATLAB

% handles    empty - handles not created until after all CreateFcns called





% --- Executes on button press in Cancel.

function Cancel_Callback(hObject, eventdata, handles)

% hObject    handle to Cancel (see GCBO)

% eventdata  reserved - to be defined in a future version of MATLAB

% handles    structure with handles and user data (see GUIDATA)

handles.var.name = '';

handles.var.value = [];

guidata(hObject,handles)

uiresume(handles.figure1);



% --- Executes on button press in addUnit.

function addUnit_Callback(hObject, eventdata, handles)

% hObject    handle to addUnit (see GCBO)

% eventdata  reserved - to be defined in a future version of MATLAB

% handles    structure with handles and user data (see GUIDATA)

[Unit str] = unitBuilder(handles.str);

handles.str = str;

handles.var.value = vertcat(handles.var.value,Unit);

set(handles.disp,'String',evalc('guidisplay(handles.var.value)'));

guidata(hObject,handles)



% --- Executes on button press in OK.

function OK_Callback(hObject, eventdata, handles)

% hObject    handle to OK (see GCBO)

% eventdata  reserved - to be defined in a future version of MATLAB

% handles    structure with handles and user data (see GUIDATA)

if isa(handles.var.value,'DCDataset') & length(handles.var.value) > 0

  handles.var.name = DatasetName(handles.existingNames);

  handles.success = 1;

  guidata(hObject,handles)

  uiresume(handles.figure1)

else

  errordlg('Unsucessful Object Creation, Please Retry!')

end





