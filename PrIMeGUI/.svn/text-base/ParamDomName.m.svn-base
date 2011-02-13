function varargout = ParamDomName(varargin)

% PARAMDOMNAME M-file for ParamDomName.fig

%      PARAMDOMNAME by itself, creates a new PARAMDOMNAME or raises the

%      existing singleton*.

%

%      H = PARAMDOMNAME returns the handle to a new PARAMDOMNAME or the handle to

%      the existing singleton*.

%

%      PARAMDOMNAME('CALLBACK',hObject,eventData,handles,...) calls the local

%      function named CALLBACK in PARAMDOMNAME.M with the given input arguments.

%

%      PARAMDOMNAME('Property','Value',...) creates a new PARAMDOMNAME or raises the

%      existing singleton*.  Starting from the left, property value pairs are

%      applied to the GUI before ParamDomName_OpeningFunction gets called.  An

%      unrecognized property name or invalid value makes property application

%      stop.  All inputs are passed to ParamDomName_OpeningFcn via varargin.

%

%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one

%      instance to run (singleton)".

%

% See also: GUIDE, GUIDATA, GUIHANDLES



% Edit the above text to modify the response to help ParamDomName



% Last Modified by GUIDE v2.5 01-Mar-2006 05:44:00



% Begin initialization code - DO NOT EDIT

gui_Singleton = 1;

gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ParamDomName_OpeningFcn, ...
                   'gui_OutputFcn',  @ParamDomName_OutputFcn, ...
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



% --- Executes just before ParamDomName is made visible.

function ParamDomName_OpeningFcn(hObject, eventdata, handles, varargin)

% This function has no output args, see OutputFcn.

% hObject    handle to figure

% eventdata  reserved - to be defined in a future version of MATLAB

% handles    structure with handles and user data (see GUIDATA)

% varargin   command line arguments to ParamDomName (see VARARGIN)



handles.output = '';



% Update handles structure

guidata(hObject, handles);



% Insert custom Title and Text if specified by the user

% Hint: when choosing keywords, be sure they are not easily confused 

% with existing figure properties.  See the output of set(figure) for

% a list of figure properties.

if(nargin > 3)

    for index = 1:2:(nargin-3),

        if nargin-3==index, break, end

        switch lower(varargin{index})

         case 'title'

          set(hObject, 'Name', varargin{index+1});

         case 'string'

          set(handles.text1, 'String', varargin{index+1});

        end

    end

end



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



% Show a question icon from dialogicons.mat - variables questIconData

% and questIconMap

% load dialogicons.mat



% IconData=questIconData;

% questIconMap(256,:) = get(handles.figure1, 'Color');

% IconCMap=questIconMap;

% 

% Img=image(IconData, 'Parent', handles.axes1);

% set(handles.figure1, 'Colormap', IconCMap);

% 

% set(handles.axes1, ...

%     'Visible', 'off', ...

%     'YDir'   , 'reverse'       , ...

%     'XLim'   , get(Img,'XData'), ...

%     'YLim'   , get(Img,'YData')  ...

%     );



% Make the GUI modal

set(handles.figure1,'WindowStyle','modal')



% UIWAIT makes ParamDomName wait for user response (see UIRESUME)

uiwait(handles.figure1);



% --- Outputs from this function are returned to the command line.

function varargout = ParamDomName_OutputFcn(hObject, eventdata, handles)

% varargout  cell array for returning output args (see VARARGOUT);

% hObject    handle to figure

% eventdata  reserved - to be defined in a future version of MATLAB

% handles    structure with handles and user data (see GUIDATA)



% Get default command line output from handles structure

varargout{1} = handles.output;



% The figure can be deleted now

delete(handles.figure1);



% --- Executes when user attempts to close figure1.

function figure1_CloseRequestFcn(hObject, eventdata, handles)

% hObject    handle to figure1 (see GCBO)

% eventdata  reserved - to be defined in a future version of MATLAB

% handles    structure with handles and user data (see GUIDATA)



if isequal(get(handles.figure1, 'waitstatus'), 'waiting')

  Cancel_Callback(handles.Cancel,eventdata,handles)

  % The GUI is still in UIWAIT, us UIRESUME

      

%    uiresume(handles.figure1);

else

    % The GUI is no longer waiting, just close it

    delete(handles.figure1);

end





% --- Executes on key press over figure1 with no controls selected.

function figure1_KeyPressFcn(hObject, eventdata, handles)

% hObject    handle to figure1 (see GCBO)

% eventdata  reserved - to be defined in a future version of MATLAB

% handles    structure with handles and user data (see GUIDATA)



% Check for "enter" or "escape"

if isequal(get(hObject,'CurrentKey'),'escape')

    % User said cancel by hitting escape

    Cancel_Callback(handles.Cancel,eventdata,handles)

end    

    

%if isequal(get(hObject,'CurrentKey'),'return')

%  OK_Callback(handles.Cancel,eventdata,handles) 

%  uiresume(handles.figure1);

%end    



function edit1_Callback(hObject, eventdata, handles)

% hObject    handle to edit1 (see GCBO)

% eventdata  reserved - to be defined in a future version of MATLAB

% handles    structure with handles and user data (see GUIDATA)



% Hints: get(hObject,'String') returns contents of edit1 as text

%        str2double(get(hObject,'String')) returns contents of edit1 as a double



handles.output = get(handles.edit1,'String');



if ~isempty(handles.output)

  % Use UIRESUME instead of delete because the OutputFcn needs

  % to get the updated handles structure.

  uiresume(handles.figure1);

%else

%  errordlg('Please Supply Dataset Name!')

end

guidata(hObject,handles)



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





% --- Executes on button press in OK.

function OK_Callback(hObject, eventdata, handles)

% hObject    handle to OK (see GCBO)

% eventdata  reserved - to be defined in a future version of MATLAB

% handles    structure with handles and user data (see GUIDATA)



handles.output = get(handles.edit1,'String');

guidata(hObject,handles)



if ~isempty(handles.output)

  % Use UIRESUME instead of delete because the OutputFcn needs

  % to get the updated handles structure.

  uiresume(handles.figure1);

else

  errordlg('Please Supply Parameter Domain Name!')

end



% --- Executes on button press in Cancel.

function Cancel_Callback(hObject, eventdata, handles)

% hObject    handle to Cancel (see GCBO)

% eventdata  reserved - to be defined in a future version of MATLAB

% handles    structure with handles and user data (see GUIDATA)



handles.output = '';

guidata(hObject,handles)

uiresume(handles.figure1);



