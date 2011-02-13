function varargout = variableImport(varargin)

% [outputvalue, outputname, fighandle] = variableImport(varargin)

%

%

% VARIABLEIMPORT M-file for variableImport.fig

%      VARIABLEIMPORT, by itself, creates a new VARIABLEIMPORT or raises the existing

%      singleton*.

%

%      H = VARIABLEIMPORT returns the handle to a new VARIABLEIMPORT or the handle to

%      the existing singleton*.

%

%      VARIABLEIMPORT('CALLBACK',hObject,eventData,handles,...) calls the local

%      function named CALLBACK in VARIABLEIMPORT.M with the given input arguments.

%

%      VARIABLEIMPORT('Property','Value',...) creates a new VARIABLEIMPORT or raises the

%      existing singleton*.  Starting from the left, property value pairs are

%      applied to the GUI before variableImport_OpeningFunction gets called.  An

%      unrecognized property name or invalid value makes property application

%      stop.  All inputs are passed to variableImport_OpeningFcn via varargin.

%

%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one

%      instance to run (singleton)".

%

% See also: GUIDE, GUIDATA, GUIHANDLES



% Edit the above text to modify the response to help variableImport



% Last Modified by GUIDE v2.5 21-Feb-2006 01:44:36



% Begin initialization code - DO NOT EDIT

gui_Singleton = 1;

gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @variableImport_OpeningFcn, ...
                   'gui_OutputFcn',  @variableImport_OutputFcn, ...
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





% --- Executes just before variableImport is made visible.

function variableImport_OpeningFcn(hObject, eventdata, handles, varargin)

% This function has no output args, see OutputFcn.

% hObject    handle to figure

% eventdata  reserved - to be defined in a future version of MATLAB

% handles    structure with handles and user data (see GUIDATA)

% varargin   command line arguments to variableImport (see VARARGIN)



% %employ input variables

% if nargin > 3 && isstruct(varargin{1})

%   variables = varargin{1};

%   if isfield(variables,'name') && isfield(variables,'value')

%     handles.variables = variables(:);

%     set(handles.Display,'String',unique(cellstr(strvcat(handles.variables.name))));

%   end

% end

if isstruct(varargin{1})

  variables = varargin{1};

  if isfield(variables,'name') && isfield(variables,'value')

    handles.variables = variables(:);

    set(handles.Display,'String',unique(cellstr(strvcat(handles.variables.name))));

  end

else

  handles.variables = struct([]);

end



% Default outputs

handles.output = [];

handles.ivarName = '';

handles.success = 0; 





handles.lastruncommand = 0;

set(handles.input,'String',''); %make this an empty string so size = [0,0]



% Update handles structure

guidata(hObject, handles);



% --- Outputs from this function are returned to the command line.

function varargout = variableImport_OutputFcn(hObject, eventdata, handles) 

% varargout  cell array for returning output args (see VARARGOUT);

% hObject    handle to figure

% eventdata  reserved - to be defined in a future version of MATLAB

% handles    structure with handles and user data (see GUIDATA)



uiwait(handles.figure1);

handles = guidata(hObject); %load most recent handles structure

varargout{1} = handles.output;

varargout{2} = handles.ivarName;

varargout{3} = handles.success;

varargout{4} = handles.variables;



delete(handles.figure1);



function input_Callback(hObject, eventdata, handles)

% hObject    handle to input (see GCBO)

% eventdata  reserved - to be defined in a future version of MATLAB

% handles    structure with handles and user data (see GUIDATA)



% Hints: get(hObject,'String') returns contents of input as text

%        str2double(get(hObject,'String')) returns contents of input as a double



handles.commands = get(hObject,'String');

guidata(hObject,handles);



% --- Executes during object creation, after setting all properties.

function input_CreateFcn(hObject, eventdata, handles)

% hObject    handle to input (see GCBO)

% eventdata  reserved - to be defined in a future version of MATLAB

% handles    empty - handles not created until after all CreateFcns called



% Hint: edit controls usually have a white background on Windows.

%       See ISPC and COMPUTER.

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))

    set(hObject,'BackgroundColor','white');

end



% --- Executes on button press in Execute.

function Execute_Callback(hObject, eventdata, handles)

% hObject    handle to Execute (see GCBO)

% eventdata  reserved - to be defined in a future version of MATLAB

% handles    structure with handles and user data (see GUIDATA)



%Execute textedit callbacks incase the user clicked, tabbed, or otherwise

%didn't hit return

input_Callback(handles.input, eventdata, handles)

handles = guidata(hObject);



try

 

  %Define all 'existing' variables (use i1asdf cuz the user probably won't)

  for i1asdf = 1:length(handles.variables)

    eval([handles.variables(i1asdf).name ' = handles.variables(i1asdf).value;']);

  end



  %Run new commands

  for i1asdf = handles.lastruncommand+1:size(handles.commands,1)

    try

      eval([handles.commands(i1asdf,:) ';'])

    catch

      MessageBox(handles,lasterr)

    end

  end

  handles.lastruncommand = size(handles.commands,1);

  vars = who;

  idx = strmatch('handles',vars,'exact');

  vars(idx) = [];

  idx = strmatch('i1asdf',vars,'exact');

  vars(idx) = [];

  idx = strmatch('hObject',vars,'exact');

  vars(idx) = [];

  idx = strmatch('eventdata',vars,'exact');

  vars(idx) = [];



  %Add new variables to handles.variables structure

  if ~isempty(handles.variables)

    for i1 = 1:length(vars)

      handles.variables(end+1,1).name = vars{i1};

      eval(['handles.variables(end).value = ' vars{i1} ';']);

    end

  else

    for i1 = 1:length(vars)

      handles.variables(i1,1).name = vars{i1};

      eval(['handles.variables(i1).value = ' vars{i1} ';']);

    end

  end



  MessageBox(handles,'=== Execution Complete ===')



  %Update variable display

  if isfield(handles,'variables')

    set(handles.Display,'String',unique(cellstr(strvcat(handles.variables.name))));

  end



catch

  %keyboard

end



guidata(hObject,handles)



function VarName_Callback(hObject, eventdata, handles)

handles.varName = get(hObject,'String');

guidata(hObject,handles);





% --- Executes during object creation, after setting all properties.

function VarName_CreateFcn(hObject, eventdata, handles)

% hObject    handle to VarName (see GCBO)

% eventdata  reserved - to be defined in a future version of MATLAB

% handles    empty - handles not created until after all CreateFcns called



% Hint: edit controls usually have a white background on Windows.

%       See ISPC and COMPUTER.

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))

    set(hObject,'BackgroundColor','white');

end





function ImportName_Callback(hObject, eventdata, handles)

handles.ivarName = get(hObject,'String');

guidata(hObject,handles);





% --- Executes during object creation, after setting all properties.

function ImportName_CreateFcn(hObject, eventdata, handles)

% hObject    handle to ImportName (see GCBO)

% eventdata  reserved - to be defined in a future version of MATLAB

% handles    empty - handles not created until after all CreateFcns called



% Hint: edit controls usually have a white background on Windows.

%       See ISPC and COMPUTER.

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))

    set(hObject,'BackgroundColor','white');

end





% --- Executes on button press in Import.

function Import_Callback(hObject, eventdata, handles)

% hObject    handle to Import (see GCBO)

% eventdata  reserved - to be defined in a future version of MATLAB

% handles    structure with handles and user data (see GUIDATA)



ImportName_Callback(handles.ImportName, eventdata, handles);

handles = guidata(hObject);

VarName_Callback(handles.VarName,eventdata,handles);

handles = guidata(hObject);



if isempty(handles.varName)

  MessageBox(handles,'Variable Name must be present for import.')

  return

end

if isempty(handles.ivarName)

  MessageBox(handles,'Imported Name must be present for import.')

  return

end



if ~isfield(handles,'variables')

  MessageBox(handles,['You have not supplied the variable ' handles.varName]);

  return

end



idx = strmatch(handles.varName,strvcat(handles.variables.name),'exact');

if isempty(idx)

  MessageBox(handles,['You have not supplied the variable ' handles.varName]);

  return

else

  handles.output = handles.variables(idx(end)).value;

end



handles.success = 1;

guidata(hObject,handles)

figure1_CloseRequestFcn(handles.figure1,eventdata,handles)



% --- Executes during object creation, after setting all properties.

function Messages_CreateFcn(hObject, eventdata, handles)

% hObject    handle to Messages (see GCBO)

% eventdata  reserved - to be defined in a future version of MATLAB

% handles    empty - handles not created until after all CreateFcns called



%create an empty 3x1 cell array so lat

%set(hObject,'String',{'';'';''})a



function MessageBox(handles,str)



%number of lines available in display

availLines = 5;



message = get(handles.Messages,'String');

creturns = find(str == 10);



newLines = 1+length(creturns);



if size(message,1) > availLines-newLines

  message(1:(size(message,1)+newLines)-availLines) = [];

  message{end+1} = str;

else

  message = vertcat(message,{str});

end

set(handles.Messages,'String',message);





% --- Executes when user attempts to close figure1.

function figure1_CloseRequestFcn(hObject, eventdata, handles)

uiresume(handles.figure1)





