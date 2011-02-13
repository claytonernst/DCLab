function varargout = unitBuilder(varargin)
% UNITBUILDER M-file for unitBuilder.fig
%      UNITBUILDER, by itself, creates a new UNITBUILDER or raises the existing
%      singleton*.
%
%      H = UNITBUILDER returns the handle to a new UNITBUILDER or the handle to
%      the existing singleton*.
%
%      UNITBUILDER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in UNITBUILDER.M with the given input arguments.
%
%      UNITBUILDER('Property','Value',...) creates a new UNITBUILDER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before unitBuilder_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to unitBuilder_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help unitBuilder

% Last Modified by GUIDE v2.5 06-Apr-2006 01:10:15

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @unitBuilder_OpeningFcn, ...
                   'gui_OutputFcn',  @unitBuilder_OutputFcn, ...
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


% --- Executes just before unitBuilder is made visible.
function unitBuilder_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to unitBuilder (see VARARGIN)

% Choose default command line output for unitBuilder
handles.output = hObject;

if nargin > 3
  handles.str = varargin{1};
else
  handles.str = struct([]);
end
handles.Dset = [];

% Update handles structure
guidata(hObject, handles);
str = handles.str;
if isfield(str,'code')
  set(handles.codeinput,'String',str.code);
end
if isfield(str,'pList')
  set(handles.pList,'String',str.pList);
end
if isfield(str,'subfunctions')
  set(handles.subfunctions,'String',str.subfunctions);
end

% UIWAIT makes unitBuilder wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = unitBuilder_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure

uiwait(handles.figure1);
handles = guidata(hObject); %load most recent handles structure
varargout{1} = handles.Dset;
varargout{2} = handles.str;
delete(handles.figure1);

function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function uinput_Callback(hObject, eventdata, handles)
% hObject    handle to uinput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of uinput as text
%        str2double(get(hObject,'String')) returns contents of uinput as a double


% --- Executes during object creation, after setting all properties.
function uinput_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uinput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function dinput_Callback(hObject, eventdata, handles)
% hObject    handle to dinput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dinput as text
%        str2double(get(hObject,'String')) returns contents of dinput as a double


% --- Executes during object creation, after setting all properties.
function dinput_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dinput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in CancelPB.
function CancelPB_Callback(hObject, eventdata, handles)
% hObject    handle to CancelPB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.Dataset = [];
guidata(hObject,handles)
uiresume(handles.figure1);

% --- Executes on button press in OKPB.
function OKPB_Callback(hObject, eventdata, handles)
str.description = get(handles.description,'String');
str.pList = get(handles.pList,'String');
str.code = get(handles.codeinput,'String');
str.subfunctions = get(handles.subfunctions,'String');
str.data = get(handles.dinput,'String');
str.unc = get(handles.uinput,'String');
str.mName = get(handles.mName,'String');
guidata(hObject,handles)

if isempty(str.pList) | isempty(str.code) | isempty(str.data) | isempty(str.unc) | isempty(str.mName)
  errordlg('Please Supply All Required Inputs!')
  return
else  
  [Dset,bool] = makeit(str);
  if bool
    handles.Dset = Dset;
    handles.str = str;
    guidata(hObject,handles)
    uiresume(handles.figure1)
  else
    errordlg('Unsucessful Object Creation, Please Check your Inputs!')
  end
end

function pList_Callback(hObject, eventdata, handles)
% hObject    handle to pList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pList as text
%        str2double(get(hObject,'String')) returns contents of pList as a double


% --- Executes during object creation, after setting all properties.
function pList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function codeinput_Callback(hObject, eventdata, handles)
% hObject    handle to codeinput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of codeinput as text
%        str2double(get(hObject,'String')) returns contents of codeinput as a double


% --- Executes during object creation, after setting all properties.
function codeinput_CreateFcn(hObject, eventdata, handles)
% hObject    handle to codeinput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function subfunctions_Callback(hObject, eventdata, handles)
% hObject    handle to subfunctions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of subfunctions as text
%        str2double(get(hObject,'String')) returns contents of subfunctions as a double


% --- Executes during object creation, after setting all properties.
function subfunctions_CreateFcn(hObject, eventdata, handles)
% hObject    handle to subfunctions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function description_Callback(hObject, eventdata, handles)
% hObject    handle to description (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of description as text
%        str2double(get(hObject,'String')) returns contents of description as a double


% --- Executes during object creation, after setting all properties.
function description_CreateFcn(hObject, eventdata, handles)
% hObject    handle to description (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function figure1_CloseRequestFcn(hObject, eventdata, handles)
uiresume(handles.figure1)

function [Dset,bool] = makeit(str);

Dset = [];
bool = 0;

EA = ResponseObservation(str2num(str.data),str2num(str.unc),{'NA'},str.description);
varName = str.mName;

%write m-file
fid = fopen([pwd filesep varName '.m'],'w');
fprintf(fid,['function OUT = %s(executionType,' ...
             'pValuesMatrix,flag)\n'],varName);
fprintf(fid,'ni = nargin;\nno = nargout;\n');
fprintf(fid,['error(nargchk(1,3,ni));\nerror(nargoutchk(0,1,no));\' ...
             'n']);
fprintf(fid,'if ~ischar(executionType)\n');
fprintf(fid,['  error([''First input to model simulation file must ' ...
             'be of type char.''])\n']);
fprintf(fid,'end\n\n');

fprintf(fid,'switch executionType\n');
fprintf(fid,' case ''description''\n');
fprintf(fid,'  OUT = [');
desc = str.description;
for i=1:size(desc,1)
  fprintf(fid,'''%s '', ...\n',desc(i,:));
end
fprintf(fid,'             ];\n');
fprintf(fid,' case ''featureList''\n');
fprintf(fid,'  OUT = {''NA''};\n');
fprintf(fid,' case ''featurePrecedence''\n');
fprintf(fid,'  OUT = 1;\n');
fprintf(fid,' case ''islinear''\n');
fprintf(fid,'  OUT = 0;\n');
fprintf(fid,' case ''isquadratic''\n');
fprintf(fid,'  OUT = 0;\n');
fprintf(fid,' case ''paramList''\n');
fprintf(fid,'  OUT = %s;\n',str.pList);
fprintf(fid,' case ''simulate''  %%assign parameter values\n');
fprintf(fid,'  N = size(pValuesMatrix,1);\n')
fprintf(fid,'  OUT = zeros(N,1);\n');
fprintf(fid,'  %%assign parameter values\n');
fprintf(fid,'  for i1 = 1:N\n');
fprintf(fid,'    pValues = pValuesMatrix(i1,:);\n');
simcode = str.code;
for i=1:size(simcode,1);
  fprintf(fid,'    %s\n',simcode(i,:));
end

fprintf(fid,'    OUT(i1,1) = out;\n');
fprintf(fid,['  end\n otherwise\n  error(''Incorrect input type to ' ...
             'model simulation file'')\nend\n\n']);

subfunc = str.subfunctions;
for i=1:size(subfunc,1)
  fprintf(fid,'%s\n',subfunc(i,:));
end

fclose(fid);

MA = ResponseModel(varName,'NA');
Dset = DClab.DCDataset(MA,EA);
bool = 1;

function mName_Callback(hObject, eventdata, handles)
% hObject    handle to mName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mName as text
%        str2double(get(hObject,'String')) returns contents of mName as a double


% --- Executes during object creation, after setting all properties.
function mName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


