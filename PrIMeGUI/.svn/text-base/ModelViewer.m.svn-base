function varargout = ModelViewer(varargin)

% MODELVIEWER M-file for ModelViewer.fig

%      MODELVIEWER, by itself, creates a new MODELVIEWER or raises the existing

%      singleton*.

%

%      H = MODELVIEWER returns the handle to a new MODELVIEWER or the handle to

%      the existing singleton*.

%

%      MODELVIEWER('CALLBACK',hObject,eventData,handles,...) calls the local

%      function named CALLBACK in MODELVIEWER.M with the given input arguments.

%

%      MODELVIEWER('Property','Value',...) creates a new MODELVIEWER or raises the

%      existing singleton*.  Starting from the left, property value pairs are

%      applied to the GUI before ModelViewer_OpeningFunction gets called.  An

%      unrecognized property name or invalid value makes property application

%      stop.  All inputs are passed to ModelViewer_OpeningFcn via varargin.

%

%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one

%      instance to run (singleton)".

%

% See also: GUIDE, GUIDATA, GUIHANDLES



% Edit the above text to modify the response to help ModelViewer



% Last Modified by GUIDE v2.5 06-Feb-2009 14:20:35



% Begin initialization code - DO NOT EDIT

gui_Singleton = 1;

gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ModelViewer_OpeningFcn, ...
                   'gui_OutputFcn',  @ModelViewer_OutputFcn, ...
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



% --- Executes just before ModelViewer is made visible.

function ModelViewer_OpeningFcn(hObject, eventdata, handles, varargin) %#ok

% This function has no output args, see OutputFcn.

% hObject    handle to figure

% eventdata  reserved - to be defined in a future version of MATLAB

% handles    structure with handles and user data (see GUIDATA)

% varargin   command line arguments to ModelViewer (see VARARGIN)



% Choose default command line output for ModelViewer

handles.output = hObject;



if length(varargin) >= 1 && isa(varargin{1},'DClab.ResponseModel')

  handles.Model = varargin{1};

  if length(varargin) > 1

    handles.Name = varargin{2};

  else

    handles.Name = '';

  end

else

  handles.Model = ResponseModel;

  handles.Name = '';

end



%tweak ui



%active parameter list box

set(handles.activeParamLB,'value',1) %initially select the first

params = handles.Model.parameterList;

domain = handles.Model.domain;

param = {};

for i = 1:length(domain)

names = gate2primeData('getPreferredKey',{'primeID',domain(i).name});

param = [param,names];

end

params = param';



set(handles.activeParamLB,'String',params);



%key display

key = handles.Model.name;

set(handles.modelPrimeKey,'String',key(1:min(length(key),35)));



%gui title

set(handles.figure1,'name',['View Model -- ' handles.Name])



%position ui

guiunits = get(handles.figure1,'Units');

set(handles.figure1,'Units','pixels');

curpos = get(handles.figure1,'Position');

scrnsize = get(0,'ScreenSize');

newpos = curpos;

newpos(1) = (scrnsize(3)-curpos(3))/2;

newpos(2) = (scrnsize(4)-curpos(4))/2;

set(handles.figure1,'Position',newpos,'Units',guiunits);



% Update handles structure

guidata(hObject, handles);



% UIWAIT makes ModelViewer wait for user response (see UIRESUME)

% uiwait(handles.figure1);



% --- Outputs from this function are returned to the command line.

function varargout = ModelViewer_OutputFcn(hObject, eventdata, handles) %#ok

% varargout  cell array for returning output args (see VARARGOUT);

% hObject    handle to figure

% eventdata  reserved - to be defined in a future version of MATLAB

% handles    structure with handles and user data (see GUIDATA)



% Get default command line output from handles structure

varargout{1} = handles.output;



% --- CreateFcns execute during object creation, after setting all properties.

% Hint: listbox and edit controls usually have a white background on Windows.

%       See ISPC and COMPUTER.



function activeParamLB_CreateFcn(hObject, eventdata, handles) %#ok

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))

    set(hObject,'BackgroundColor','white');

end



function modelPrimeKey_CreateFcn(hObject, eventdata, handles) %#ok

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))

    set(hObject,'BackgroundColor','white');

end



% --- Callback functions



function activeParamLB_Callback(hObject, eventdata, handles) %#ok



function modelPrimeKey_Callback(hObject, eventdata, handles) %#ok

%keyboard





% --------------------------------------------------------------------

function parameterID_menu_Callback(hObject, eventdata, handles)

% hObject    handle to parameterID_menu (see GCBO)

% eventdata  reserved - to be defined in a future version of MATLAB

% handles    structure with handles and user data (see GUIDATA)



LBval = get(handles.activeParamLB,'Value');

uDstruct = handles.Model.userData;

pid = uDstruct.modelParamIDs(LBval);

if iscell(pid) && isscalar(pid)

    uDstruct.callback2show(pid{1});

elseif iscell(pid)

    uDstruct.callback2show(pid{1});

elseif ischar(pid)

    uDstruct.callback2show(pid);

end



% --------------------------------------------------------------------

function modelID_menu_Callback(hObject, eventdata, handles)

% hObject    handle to modelID_menu (see GCBO)

% eventdata  reserved - to be defined in a future version of MATLAB

% handles    structure with handles and user data (see GUIDATA)



uDstruct = handles.Model.userData;

pid = uDstruct.modelID;

if iscell(pid) && isscalar(pid)

    uDstruct.callback2show(pid{1});

elseif iscell(pid)

    uDstruct.callback2show2(pid{1},pid{2});

elseif ischar(pid)

    uDstruct.callback2show(pid);

end



% --------------------------------------------------------------------

function target_menu_Callback(hObject, eventdata, handles)

% hObject    handle to ShowXML_menu (see GCBO)

% eventdata  reserved - to be defined in a future version of MATLAB

% handles    structure with handles and user data (see GUIDATA)



uDstruct = handles.Model.userData;

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



uDstruct = handles.Model.userData;

pid = uDstruct.targetLinks;

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



LBval = get(handles.activeParamLB,'Value');

uDstruct = handles.Model.domain;

%pid = uDstruct(LBval).parameterLinks;



%uDstruct = handles.Model.userData;

%if iscell(pid) && isscalar(pid)

    %uDstruct.callback2show(pid{1});

%elseif iscell(pid)

    %uDstruct.callback2show2(pid{1},pid{2});

%elseif ischar(pid)

    %uDstruct.callback2show(pid);

%end



% --------------------------------------------------------------------

function trial_model_menu_Callback(hObject, eventdata, handles)

% hObject    handle to trial_model_menu (see GCBO)

% eventdata  reserved - to be defined in a future version of MATLAB

% handles    structure with handles and user data (see GUIDATA)



uDstruct = handles.Model.userData;

pid = uDstruct.trialModelID;

if iscell(pid) && isscalar(pid)

    uDstruct.callback2show(pid{1});

elseif ischar(pid)

    uDstruct.callback2show(pid);

end



% --------------------------------------------------------------------

function Model_Menu_Callback(hObject, eventdata, handles)

% hObject    handle to Param_Menu (see GCBO)

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



uDstruct = handles.Model.userData;

mid = uDstruct.trialModelID;



%y = gate2primeData('get',{mid});

%mname = getPreferredKey(y);

mname = gate2primeData('getPreferredKey',{'primeID',mid});



set(handles.trial_model_menu, 'Label', char(mname));



LBval = get(handles.activeParamLB,'Value');

uDstruct = handles.Model.domain;

%pid = uDstruct(LBval).parameterLinks;

%recMark = pid{2}(1:length(pid{2})-8);

%      switch recMark

%         case 'rk'

%            y = 'Rate Coefficient';

%         case 'thp'

%            y = 'Thermodynamic Polynomial';

%      end

%set(handles. Parameter_Links_Menu, 'Label', y);



