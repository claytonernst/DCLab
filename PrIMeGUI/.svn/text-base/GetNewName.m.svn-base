function varargout = GetNewName(varargin)

% GETNEWNAME M-file for GetNewName.fig

%      GETNEWNAME by itself, creates a new GETNEWNAME or raises the

%      existing singleton*.

%

%      H = GETNEWNAME returns the handle to a new GETNEWNAME or the handle to

%      the existing singleton*.

%

% ***  name = GETNEWNAME(prompt,figname,existingNames,defname) calls the GUI to

%      get a name from the user that is not a member of the cell array of

%      strings existingNames. The input figname is the title of the gui,

%      prompt is a cell array containing the string placed above the

%      edit box, and defname is the string initially in the box.

%

%      GETNEWNAME('CALLBACK',hObject,eventData,handles,...) calls the local

%      function named CALLBACK in GETNEWNAME.M with the given input arguments.

%

%      GETNEWNAME('Property','Value',...) creates a new GETNEWNAME or raises the

%      existing singleton*.  Starting from the left, property value pairs are

%      applied to the GUI before GetNewName_OpeningFunction gets called.  An

%      unrecognized property name or invalid value makes property application

%      stop.  All inputs are passed to GetNewName_OpeningFcn via varargin.

%

%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one

%      instance to run (singleton)".

%

% See also: GUIDE, GUIDATA, GUIHANDLES



% Edit the above text to modify the response to help GetNewName



% Last Modified by GUIDE v2.5 15-Mar-2007 17:58:27



% Begin initialization code - DO NOT EDIT

gui_Singleton = 1;

gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GetNewName_OpeningFcn, ...
                   'gui_OutputFcn',  @GetNewName_OutputFcn, ...
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



% --- Executes just before GetNewName is made visible.

function GetNewName_OpeningFcn(hObject, eventdata, handles, varargin) %#ok

% This function has no output args, see OutputFcn.

% hObject    handle to figure

% eventdata  reserved - to be defined in a future version of MATLAB

% handles    structure with handles and user data (see GUIDATA)

% varargin   command line arguments to GetNewName (see VARARGIN)



ni = length(varargin);

switch ni

  case 1

    prompt = varargin{1};

    tstr = '';

    handles.names = {};

    set(handles.edit1,'String','')

  case 2

    prompt = varargin{1};

    tstr = varargin{2};

    handles.names = {};

    set(handles.edit1,'String','')

  case 3

    prompt = varargin{1};

    tstr = varargin{2};

    handles.names = varargin{3};

    set(handles.edit1,'String','')

  case 4

    prompt = varargin{1};

    tstr = varargin{2};

    handles.names = varargin{3};

    set(handles.edit1,'String',varargin{4})

  otherwise

    error(nargchk(1,4,ni))

end



% Choose default command line output for GetNewName

handles.newName = '';



%Tweak the UI

set(handles.figure1,'name',tstr)

set(handles.prompt,'String',prompt)



% Update handles structure

guidata(hObject, handles);



% Determine the position of the dialog - centered on the callback figure

% if available, else, centered on the screen

FigPos=get(0,'DefaultFigurePosition');

OldUnits = get(hObject, 'Units');

set(hObject, 'Units', 'pixels');

OldPos = get(hObject,'Position');

FigWidth = OldPos(3);

FigHeight = OldPos(4);

if isempty(gcbf)

    ScreenUnits=get(0,'Units');

    set(0,'Units','pixels');

    ScreenSize=get(0,'ScreenSize');

    set(0,'Units',ScreenUnits);



    FigPos(1)=1/2*(ScreenSize(3)-FigWidth);

    FigPos(2)=2/3*(ScreenSize(4)-FigHeight);

else

    GCBFOldUnits = get(gcbf,'Units');

    set(gcbf,'Units','pixels');

    GCBFPos = get(gcbf,'Position');

    set(gcbf,'Units',GCBFOldUnits);

    FigPos(1:2) = [(GCBFPos(1) + GCBFPos(3) / 2) - FigWidth / 2, ...

                   (GCBFPos(2) + GCBFPos(4) / 2) - FigHeight / 2];

end

FigPos(3:4)=[FigWidth FigHeight];

set(hObject, 'Position', FigPos);

set(hObject, 'Units', OldUnits);



% Make the GUI modal

set(handles.figure1,'WindowStyle','modal')



% Give focus to the edit box

uicontrol(handles.edit1)

refresh %this probably does nothing since it's not visible...

drawnow

% UIWAIT makes DatasetUnitImport wait for user response (see UIRESUME)

uiwait(handles.figure1);



% --- Outputs from this function are returned to the command line.

function varargout = GetNewName_OutputFcn(hObject, eventdata, handles) %#ok

varargout{1} = handles.newName;



% The figure can be deleted now

delete(handles.figure1);



function figure1_CloseRequestFcn(hObject, eventdata, handles) %#ok



if isequal(get(handles.figure1, 'waitstatus'), 'waiting')

  Cancel_Callback(handles.Cancel,eventdata,handles)

else

  % The GUI is no longer waiting, something bombed, just close it

  delete(handles.figure1);

end



% --- Executes on key press over figure1 with no controls selected.

function figure1_KeyPressFcn(hObject, eventdata, handles) %#ok



% Check for "enter" or "escape"

if isequal(get(hObject,'CurrentKey'),'escape')

  % User said cancel by hitting escape

  Cancel_Callback(handles.Cancel,eventdata,handles)

end    

if isequal(get(hObject,'CurrentKey'),'return')

  % User said OK by hitting return

  OK_Callback(handles.Cancel,eventdata,handles) 

end    



function edit1_Callback(hObject, eventdata, handles) %#ok

%

%Maybe someday we will call OK if the user types enter



function edit1_CreateFcn(hObject, eventdata, handles) %#ok

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))

    set(hObject,'BackgroundColor','white');

end



function OK_Callback(hObject, eventdata, handles) %#ok



name = get(handles.edit1,'String');

if isempty(name)

  errordlg('No Name Provided','Please retry','modal')

elseif ismember(name,handles.names)

  errordlg('Specified Name already in use','Please retry','modal')

else

  handles.newName = name;

  guidata(hObject,handles)

  uiresume(handles.figure1)

end



function Cancel_Callback(hObject, eventdata, handles) %#ok

handles.newName = '';

guidata(hObject,handles)

uiresume(handles.figure1)



