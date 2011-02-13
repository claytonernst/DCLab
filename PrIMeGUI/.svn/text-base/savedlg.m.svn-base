function varargout = savedlg(varargin)

% SAVEDLG M-file for savedlg.fig

%      SAVEDLG by itself, creates a new SAVEDLG or raises the

%      existing singleton*.

%

%      H = SAVEDLG returns the handle to a new SAVEDLG or the handle to

%      the existing singleton*.

%

%      SAVEDLG('CALLBACK',hObject,eventData,handles,...) calls the local

%      function named CALLBACK in SAVEDLG.M with the given input arguments.

%

%      SAVEDLG('Property','Value',...) creates a new SAVEDLG or raises the

%      existing singleton*.  Starting from the left, property value pairs are

%      applied to the GUI before savedlg_OpeningFunction gets called.  An

%      unrecognized property name or invalid value makes property application

%      stop.  All inputs are passed to savedlg_OpeningFcn via varargin.

%

%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one

%      instance to run (singleton)".

%

% See also: GUIDE, GUIDATA, GUIHANDLES



% Edit the above text to modify the response to help savedlg



% Last Modified by GUIDE v2.5 20-Mar-2007 00:39:07



% Begin initialization code - DO NOT EDIT

gui_Singleton = 1;

gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @savedlg_OpeningFcn, ...
                   'gui_OutputFcn',  @savedlg_OutputFcn, ...
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



% --- Executes just before savedlg is made visible.

function savedlg_OpeningFcn(hObject, eventdata, handles, varargin) %#ok



% Choose default command line output for savedlg

handles.output = 'Cancel';



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

load dialogicons.mat



IconData=questIconData;

questIconMap(256,:) = get(handles.figure1, 'Color');

IconCMap=questIconMap;



Img=image(IconData, 'Parent', handles.axes1);

set(handles.figure1, 'Colormap', IconCMap);



set(handles.axes1, ...
    'Visible', 'off', ...
    'YDir'   , 'reverse'       , ...
    'XLim'   , get(Img,'XData'), ...
    'YLim'   , get(Img,'YData')  ...
    );



% Make the GUI modal

set(handles.figure1,'WindowStyle','modal')



% UIWAIT makes savedlg wait for user response (see UIRESUME)

uiwait(handles.figure1);



% --- Outputs from this function are returned to the command line.

function varargout = savedlg_OutputFcn(hObject, eventdata, handles) %#ok



% Get default command line output from handles structure

varargout{1} = handles.output;



% The figure can be deleted now

delete(handles.figure1);



function PB_Callback(hObject, eventdata, handles) %#ok

handles.output = get(hObject,'String');



% Update handles structure

guidata(hObject, handles);



% Use UIRESUME instead of delete because the OutputFcn needs

% to get the updated handles structure.

uiresume(handles.figure1);



function figure1_CloseRequestFcn(hObject, eventdata, handles) %#ok

if isequal(get(handles.figure1, 'waitstatus'), 'waiting')

    % The GUI is still in UIWAIT, us UIRESUME

    uiresume(handles.figure1);

else

    % The GUI is no longer waiting, just close it

    delete(handles.figure1);

end



function figure1_KeyPressFcn(hObject, eventdata, handles) %#ok

% 

% % Check for "enter" or "escape"

% if isequal(get(hObject,'CurrentKey'),'escape')

%     % User said no by hitting escape

%     handles.output = 'Cancel';

%     

%     % Update handles structure

%     guidata(hObject, handles);

%     

%     uiresume(handles.figure1);

% end    

%     

% if isequal(get(hObject,'CurrentKey'),'return')

%     uiresume(handles.figure1);

% end    

