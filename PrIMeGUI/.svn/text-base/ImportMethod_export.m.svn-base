function varargout = ImportMethod_export(varargin)

% IMPORTMETHOD_EXPORT M-file for ImportMethod_export.fig

%      IMPORTMETHOD_EXPORT by itself, creates a new IMPORTMETHOD_EXPORT or raises the

%      existing singleton*.

%

%      H = IMPORTMETHOD_EXPORT returns the handle to a new IMPORTMETHOD_EXPORT or the handle to

%      the existing singleton*.

%

%      IMPORTMETHOD_EXPORT('CALLBACK',hObject,eventData,handles,...) calls the local

%      function named CALLBACK in IMPORTMETHOD_EXPORT.M with the given input arguments.

%

%      IMPORTMETHOD_EXPORT('Property','Value',...) creates a new IMPORTMETHOD_EXPORT or raises the

%      existing singleton*.  Starting from the left, property value pairs are

%      applied to the GUI before ImportMethod_export_OpeningFunction gets called.  An

%      unrecognized property name or invalid value makes property application

%      stop.  All inputs are passed to ImportMethod_export_OpeningFcn via varargin.

%

%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one

%      instance to run (singleton)".

%

% See also: GUIDE, GUIDATA, GUIHANDLES



% Edit the above text to modify the response to help ImportMethod_export



% Last Modified by GUIDE v2.5 06-Apr-2006 03:11:50



% Begin initialization code - DO NOT EDIT

gui_Singleton = 1;

gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ImportMethod_export_OpeningFcn, ...
                   'gui_OutputFcn',  @ImportMethod_export_OutputFcn, ...
                   'gui_LayoutFcn',  @ImportMethod_export_LayoutFcn, ...
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



% --- Executes just before ImportMethod_export is made visible.

function ImportMethod_export_OpeningFcn(hObject, eventdata, handles, varargin)

% This function has no output args, see OutputFcn.

% hObject    handle to figure

% eventdata  reserved - to be defined in a future version of MATLAB

% handles    structure with handles and user data (see GUIDATA)

% varargin   command line arguments to ImportMethod_export (see VARARGIN)



% Choose default command line output for ImportMethod_export

handles.output = 'Yes';



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

% Img=image(IconData);

% set(handles.figure1, 'Colormap', IconCMap);



% set(gca, ...

%     'Visible', 'off', ...

%     'YDir'   , 'reverse'       , ...

%     'XLim'   , get(Img,'XData'), ...

%     'YLim'   , get(Img,'YData')  ...

%     );



% Make the GUI modal

set(handles.figure1,'WindowStyle','modal')



% UIWAIT makes ImportMethod_export wait for user response (see UIRESUME)

uiwait(handles.figure1);



% --- Outputs from this function are returned to the command line.

function varargout = ImportMethod_export_OutputFcn(hObject, eventdata, handles)

% varargout  cell array for returning output args (see VARARGOUT);

% hObject    handle to figure

% eventdata  reserved - to be defined in a future version of MATLAB

% handles    structure with handles and user data (see GUIDATA)



% Get default command line output from handles structure

varargout{1} = handles.output;



% The figure can be deleted now

delete(handles.figure1);



% --- Executes on button press in pushbutton1.

function pushbutton1_Callback(hObject, eventdata, handles)

% hObject    handle to pushbutton1 (see GCBO)

% eventdata  reserved - to be defined in a future version of MATLAB

% handles    structure with handles and user data (see GUIDATA)



handles.output = get(hObject,'String');



% Update handles structure

guidata(hObject, handles);



% Use UIRESUME instead of delete because the OutputFcn needs

% to get the updated handles structure.

uiresume(handles.figure1);



% --- Executes on button press in pushbutton2.

function pushbutton2_Callback(hObject, eventdata, handles)

% hObject    handle to pushbutton2 (see GCBO)

% eventdata  reserved - to be defined in a future version of MATLAB

% handles    structure with handles and user data (see GUIDATA)



handles.output = get(hObject,'String');



% Update handles structure

guidata(hObject, handles);



% Use UIRESUME instead of delete because the OutputFcn needs

% to get the updated handles structure.

uiresume(handles.figure1);





% --- Executes when user attempts to close figure1.

function figure1_CloseRequestFcn(hObject, eventdata, handles)

% hObject    handle to figure1 (see GCBO)

% eventdata  reserved - to be defined in a future version of MATLAB

% handles    structure with handles and user data (see GUIDATA)



if isequal(get(handles.figure1, 'waitstatus'), 'waiting')

    % The GUI is still in UIWAIT, us UIRESUME

    uiresume(handles.figure1);

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

    % User said no by hitting escape

    handles.output = 'No';

    

    % Update handles structure

    guidata(hObject, handles);

    

    uiresume(handles.figure1);

end    

    

if isequal(get(hObject,'CurrentKey'),'return')

    uiresume(handles.figure1);

end    





% --- Creates and returns a handle to the GUI figure. 

function h1 = ImportMethod_export_LayoutFcn(policy)

% policy - create a new figure or use a singleton. 'new' or 'reuse'.



persistent hsingleton;

if strcmpi(policy, 'reuse') & ishandle(hsingleton)

    h1 = hsingleton;

    return;

end



appdata = [];

appdata.GUIDEOptions = struct(...
    'active_h', [], ...
    'taginfo', struct(...
    'figure', 2, ...
    'pushbutton', 3, ...
    'axes', 2, ...
    'text', 2), ...
    'override', 1, ...
    'release', 13, ...
    'resize', 'none', ...
    'accessibility', 'callback', ...
    'mfile', 1, ...
    'callbacks', 1, ...
    'singleton', 1, ...
    'syscolorfig', 1, ...
    'blocking', 0, ...
    'lastSavedFile', 'F:\Documents and Settings\feeley\My Documents\tbx10\gui\ImportMethod_export.m');

appdata.UsedByGUIData_m = struct(...
    'figure1', [], ...
    'text1', 2.00439453125, ...
    'axes1', 3.00439453125, ...
    'pushbutton2', [], ...
    'pushbutton1', [], ...
    'output', []);

appdata.lastValidTag = 'figure1';

appdata.GUIDELayoutEditor = [];



h1 = figure(...
'Units','characters',...
'CloseRequestFcn','ImportMethod_export(''figure1_CloseRequestFcn'',gcbf,[],guidata(gcbf))',...
'Color',[0.564705882352941 0.690196078431373 0.658823529411765],...
'Colormap',[0 0 0;1 1 1;0.984313725490196 0.956862745098039 0.6;0.984313725490196 0.952941176470588 0.6;0 0 0.6;0.988235294117647 0.956862745098039 0.603921568627451;0.988235294117647 0.956862745098039 0.6;0.690196078431373 0.662745098039216 0.666666666666667;0.0549019607843137 0.0509803921568627 0.0549019607843137;0.0627450980392157 0.0588235294117647 0.0627450980392157;0.0705882352941176 0.0666666666666667 0.0705882352941176;0.0156862745098039 0.0117647058823529 0.0196078431372549;0.0352941176470588 0.0313725490196078 0.0392156862745098;0.623529411764706 0.596078431372549 0.658823529411765;0.0196078431372549 0.0156862745098039 0.0274509803921569;0.501960784313725 0.482352941176471 0.647058823529412;0.447058823529412 0.427450980392157 0.643137254901961;0.388235294117647 0.372549019607843 0.63921568627451;0.270588235294118 0.258823529411765 0.627450980392157;0.294117647058824 0.282352941176471 0.627450980392157;0.309803921568627 0.298039215686275 0.631372549019608;0.352941176470588 0.341176470588235 0.635294117647059;0.125490196078431 0.12156862745098 0.611764705882353;0.149019607843137 0.145098039215686 0.615686274509804;0.192156862745098 0.184313725490196 0.619607843137255;0.223529411764706 0.215686274509804 0.619607843137255;0 0 0.00784313725490196;0 0 0.00392156862745098;0.0196078431372549 0.0196078431372549 0.603921568627451;0.0470588235294118 0.0431372549019608 0.603921568627451;0.0784313725490196 0.0745098039215686 0.607843137254902;0.0823529411764706 0.0784313725490196 0.607843137254902;0.105882352941176 0.101960784313725 0.611764705882353;0.00392156862745098 0.00392156862745098 0.0196078431372549;0.00784313725490196 0.00784313725490196 0.0196078431372549;0.0117647058823529 0.0117647058823529 0.0274509803921569;0.0235294117647059 0.0235294117647059 0.0352941176470588;0.0274509803921569 0.0274509803921569 0.0392156862745098;1 0.996078431372549 0.623529411764706;1 1 0.627450980392157;1 0.996078431372549 0.631372549019608;1 1 0.635294117647059;1 1 0.643137254901961;1 1 0.650980392156863;0.0705882352941176 0.0705882352941176 0.0509803921568627;0.305882352941176 0.305882352941176 0.227450980392157;0.16078431372549 0.16078431372549 0.12156862745098;0.0392156862745098 0.0392156862745098 0.0352941176470588;0.0705882352941176 0.0705882352941176 0.0666666666666667;0.0862745098039216 0.0862745098039216 0.0823529411764706;0.184313725490196 0.184313725490196 0.176470588235294;0.0941176470588235 0.0941176470588235 0.0901960784313725;0.101960784313725 0.101960784313725 0.0980392156862745;0.145098039215686 0.145098039215686 0.141176470588235;1 0.988235294117647 0.615686274509804;1 0.992156862745098 0.619607843137255;0.925490196078431 0.913725490196078 0.6;0.423529411764706 0.419607843137255 0.298039215686275;1 0.976470588235294 0.611764705882353;0.996078431372549 0.972549019607843 0.607843137254902;1 0.980392156862745 0.615686274509804;1 0.984313725490196 0.619607843137255;1 0.976470588235294 0.619607843137255;0.988235294117647 0.972549019607843 0.615686274509804;1 0.980392156862745 0.627450980392157;0.988235294117647 0.972549019607843 0.619607843137255;0.984313725490196 0.964705882352941 0.615686274509804;0.219607843137255 0.215686274509804 0.145098039215686;0.4 0.392156862745098 0.270588235294118;0.258823529411765 0.254901960784314 0.192156862745098;0.145098039215686 0.141176470588235 0.0862745098039216;0.992156862745098 0.96078431372549 0.603921568627451;0.988235294117647 0.96078431372549 0.6;0.96078431372549 0.929411764705882 0.584313725490196;0.996078431372549 0.968627450980392 0.607843137254902;0.988235294117647 0.96078431372549 0.603921568627451;0.96078431372549 0.933333333333333 0.588235294117647;0.945098039215686 0.913725490196078 0.576470588235294;0.996078431372549 0.964705882352941 0.611764705882353;0.984313725490196 0.952941176470588 0.603921568627451;0.964705882352941 0.941176470588235 0.592156862745098;0.964705882352941 0.937254901960784 0.592156862745098;0.956862745098039 0.925490196078431 0.588235294117647;0.949019607843137 0.92156862745098 0.584313725490196;0.984313725490196 0.96078431372549 0.607843137254902;0.952941176470588 0.925490196078431 0.588235294117647;0.972549019607843 0.949019607843137 0.607843137254902;0.956862745098039 0.929411764705882 0.6;0.937254901960784 0.909803921568627 0.588235294117647;0.929411764705882 0.901960784313726 0.584313725490196;0.92156862745098 0.898039215686275 0.584313725490196;0.909803921568627 0.882352941176471 0.576470588235294;0.850980392156863 0.827450980392157 0.541176470588235;0.611764705882353 0.596078431372549 0.4;0.407843137254902 0.396078431372549 0.270588235294118;0.458823529411765 0.447058823529412 0.309803921568627;0.368627450980392 0.36078431372549 0.258823529411765;0.329411764705882 0.32156862745098 0.235294117647059;0.231372549019608 0.227450980392157 0.176470588235294;0.988235294117647 0.952941176470588 0.6;0.988235294117647 0.952941176470588 0.603921568627451;0.984313725490196 0.949019607843137 0.6;0.92156862745098 0.890196078431373 0.580392156862745;0.819607843137255 0.792156862745098 0.52156862745098;0.83921568627451 0.811764705882353 0.537254901960784;0.8 0.772549019607843 0.509803921568627;0.764705882352941 0.737254901960784 0.494117647058824;0.713725490196078 0.690196078431373 0.462745098039216;0.741176470588235 0.713725490196078 0.482352941176471;0.580392156862745 0.56078431372549 0.380392156862745;0.215686274509804 0.207843137254902 0.141176470588235;0.698039215686274 0.674509803921569 0.458823529411765;0.619607843137255 0.6 0.407843137254902;0.682352941176471 0.658823529411765 0.450980392156863;0.450980392156863 0.435294117647059 0.301960784313725;0.262745098039216 0.254901960784314 0.176470588235294;0.584313725490196 0.564705882352941 0.396078431372549;0.486274509803922 0.470588235294118 0.329411764705882;0.6 0.580392156862745 0.407843137254902;0.470588235294118 0.454901960784314 0.32156862745098;0.505882352941176 0.490196078431373 0.349019607843137;0.388235294117647 0.376470588235294 0.274509803921569;0.403921568627451 0.392156862745098 0.290196078431373;0.266666666666667 0.258823529411765 0.192156862745098;0.180392156862745 0.176470588235294 0.137254901960784;0.72156862745098 0.694117647058824 0.470588235294118;0.6 0.576470588235294 0.392156862745098;0.101960784313725 0.0980392156862745 0.0705882352941176;0.309803921568627 0.298039215686275 0.215686274509804;0.313725490196078 0.301960784313725 0.219607843137255;0.250980392156863 0.243137254901961 0.180392156862745;0.141176470588235 0.137254901960784 0.105882352941176;0.156862745098039 0.152941176470588 0.12156862745098;0.0862745098039216 0.0823529411764706 0.0588235294117647;0.494117647058824 0.474509803921569 0.349019607843137;0.286274509803922 0.274509803921569 0.203921568627451;0.219607843137255 0.211764705882353 0.164705882352941;0.243137254901961 0.235294117647059 0.184313725490196;0.0627450980392157 0.0588235294117647 0.0392156862745098;0.192156862745098 0.184313725490196 0.145098039215686;0.443137254901961 0.43921568627451 0.419607843137255;0.0784313725490196 0.0745098039215686 0.0588235294117647;0.164705882352941 0.156862745098039 0.125490196078431;0.117647058823529 0.113725490196078 0.0980392156862745;0.152941176470588 0.145098039215686 0.117647058823529;0.850980392156863 0.815686274509804 0.682352941176471;0.835294117647059 0.8 0.67843137254902;0.0470588235294118 0.0431372549019608 0.0313725490196078;0.0862745098039216 0.0823529411764706 0.0705882352941176;0.803921568627451 0.772549019607843 0.67843137254902;0.23921568627451 0.235294117647059 0.223529411764706;0.513725490196078 0.505882352941176 0.482352941176471;0.568627450980392 0.56078431372549 0.537254901960784;0.56078431372549 0.552941176470588 0.529411764705882;0.556862745098039 0.549019607843137 0.525490196078431;0.552941176470588 0.545098039215686 0.52156862745098;0.270588235294118 0.266666666666667 0.254901960784314;0.607843137254902 0.6 0.576470588235294;0.576470588235294 0.568627450980392 0.545098039215686;0.290196078431373 0.286274509803922 0.274509803921569;0.498039215686275 0.490196078431373 0.470588235294118;0.482352941176471 0.474509803921569 0.454901960784314;0.47843137254902 0.470588235294118 0.450980392156863;0.533333333333333 0.525490196078431 0.505882352941176;0.529411764705882 0.52156862745098 0.501960784313725;0.513725490196078 0.505882352941176 0.486274509803922;0.505882352941176 0.498039215686275 0.47843137254902;0.501960784313725 0.494117647058824 0.474509803921569;0.552941176470588 0.545098039215686 0.525490196078431;0.772549019607843 0.741176470588235 0.674509803921569;0.662745098039216 0.650980392156863 0.623529411764706;0.647058823529412 0.635294117647059 0.607843137254902;0.701960784313725 0.690196078431373 0.662745098039216;0.686274509803922 0.674509803921569 0.647058823529412;0.670588235294118 0.658823529411765 0.631372549019608;0.0352941176470588 0.0313725490196078 0.0235294117647059;0.129411764705882 0.12156862745098 0.105882352941176;0.6 0.588235294117647 0.564705882352941;0.588235294117647 0.576470588235294 0.552941176470588;0.580392156862745 0.568627450980392 0.545098039215686;0.63921568627451 0.627450980392157 0.603921568627451;0.627450980392157 0.615686274509804 0.592156862745098;0.623529411764706 0.611764705882353 0.588235294117647;0.619607843137255 0.607843137254902 0.584313725490196;0.611764705882353 0.6 0.576470588235294;0.423529411764706 0.415686274509804 0.4;0.686274509803922 0.674509803921569 0.650980392156863;0.682352941176471 0.670588235294118 0.647058823529412;0.67843137254902 0.666666666666667 0.643137254901961;0.674509803921569 0.662745098039216 0.63921568627451;0.666666666666667 0.654901960784314 0.631372549019608;0.662745098039216 0.650980392156863 0.627450980392157;0.654901960784314 0.643137254901961 0.619607843137255;0.650980392156863 0.63921568627451 0.615686274509804;0.454901960784314 0.447058823529412 0.431372549019608;0.450980392156863 0.443137254901961 0.427450980392157;0.43921568627451 0.431372549019608 0.415686274509804;0.435294117647059 0.427450980392157 0.411764705882353;0.227450980392157 0.223529411764706 0.215686274509804;0.466666666666667 0.458823529411765 0.443137254901961;0.247058823529412 0.243137254901961 0.235294117647059;0.243137254901961 0.23921568627451 0.231372549019608;0.23921568627451 0.235294117647059 0.227450980392157;0.235294117647059 0.231372549019608 0.223529411764706;0.258823529411765 0.254901960784314 0.247058823529412;0.294117647058824 0.290196078431373 0.282352941176471;0.32156862745098 0.317647058823529 0.309803921568627;0.745098039215686 0.717647058823529 0.670588235294118;0.529411764705882 0.517647058823529 0.498039215686275;0.52156862745098 0.509803921568627 0.490196078431373;0.572549019607843 0.56078431372549 0.541176470588235;0.564705882352941 0.552941176470588 0.533333333333333;0.545098039215686 0.533333333333333 0.513725490196078;0.592156862745098 0.580392156862745 0.56078431372549;0.694117647058824 0.67843137254902 0.654901960784314;0.352941176470588 0.345098039215686 0.333333333333333;0.345098039215686 0.337254901960784 0.325490196078431;0.368627450980392 0.36078431372549 0.349019607843137;0.407843137254902 0.4 0.388235294117647;0.4 0.392156862745098 0.380392156862745;0.388235294117647 0.380392156862745 0.368627450980392;0.490196078431373 0.47843137254902 0.462745098039216;0.474509803921569 0.462745098039216 0.447058823529412;0.733333333333333 0.705882352941177 0.670588235294118;0.0588235294117647 0.0549019607843137 0.0509803921568627;0.101960784313725 0.0980392156862745 0.0941176470588235;0.133333333333333 0.129411764705882 0.125490196078431;0.701960784313725 0.682352941176471 0.662745098039216;0.27843137254902 0.270588235294118 0.262745098039216;0.145098039215686 0.141176470588235 0.137254901960784;0.333333333333333 0.325490196078431 0.317647058823529;0.317647058823529 0.309803921568627 0.301960784313725;0.309803921568627 0.301960784313725 0.294117647058824;0.164705882352941 0.16078431372549 0.156862745098039;0.203921568627451 0.2 0.196078431372549;0.0196078431372549 0.0156862745098039 0.0156862745098039;0.0470588235294118 0.0431372549019608 0.0431372549019608;0.0549019607843137 0.0509803921568627 0.0509803921568627;0.0705882352941176 0.0666666666666667 0.0666666666666667;0.0745098039215686 0.0705882352941176 0.0705882352941176;0.0784313725490196 0.0745098039215686 0.0745098039215686;0.12156862745098 0.117647058823529 0.117647058823529;0.113725490196078 0.109803921568627 0.109803921568627;0.172549019607843 0.168627450980392 0.168627450980392;0.109803921568627 0.109803921568627 0.109803921568627;0.105882352941176 0.105882352941176 0.105882352941176;0.0941176470588235 0.0941176470588235 0.0941176470588235;0.0823529411764706 0.0823529411764706 0.0823529411764706;0.0627450980392157 0.0627450980392157 0.0627450980392157;0.0588235294117647 0.0588235294117647 0.0588235294117647;0.0509803921568627 0.0509803921568627 0.0509803921568627;0.0392156862745098 0.0392156862745098 0.0392156862745098;0.0313725490196078 0.0313725490196078 0.0313725490196078;0.0274509803921569 0.0274509803921569 0.0274509803921569;0.00784313725490196 0.00784313725490196 0.00784313725490196;0.752941176470588 0.752941176470588 0.752941176470588],...
'IntegerHandle','off',...
'InvertHardcopy',get(0,'defaultfigureInvertHardcopy'),...
'KeyPressFcn','ImportMethod_export(''figure1_KeyPressFcn'',gcbo,[],guidata(gcbo))',...
'MenuBar','none',...
'Name','ImportMethod',...
'NumberTitle','off',...
'PaperPosition',get(0,'defaultfigurePaperPosition'),...
'Position',[131.2 40.8461538461538 56 7.92307692307692],...
'Resize','off',...
'HandleVisibility','callback',...
'Tag','figure1',...
'UserData',[],...
'Behavior',get(0,'defaultfigureBehavior'),...
'Visible',get(0,'defaultfigureVisible'),...
'CreateFcn', {@local_CreateFcn, '', appdata} );



appdata = [];

appdata.lastValidTag = 'pushbutton1';



h2 = uicontrol(...
'Parent',h1,...
'Units','normalized',...
'Callback','ImportMethod_export(''pushbutton1_Callback'',gcbo,[],guidata(gcbo))',...
'FontSize',10,...
'FontWeight','bold',...
'ListboxTop',0,...
'Position',[0.132142857142857 0.524271844660194 0.703571428571428 0.378640776699029],...
'String','Command Line Construction',...
'Tag','pushbutton1',...
'Behavior',get(0,'defaultuicontrolBehavior'),...
'CreateFcn', {@local_CreateFcn, '', appdata} );



appdata = [];

appdata.lastValidTag = 'pushbutton2';



h3 = uicontrol(...
'Parent',h1,...
'Units','normalized',...
'Callback','ImportMethod_export(''pushbutton2_Callback'',gcbo,[],guidata(gcbo))',...
'CData',[],...
'FontSize',10,...
'ListboxTop',0,...
'Position',[0.135714285714286 0.0679611650485437 0.692857142857143 0.368932038834951],...
'String','GUI Construction',...
'Tag','pushbutton2',...
'UserData',[],...
'Behavior',get(0,'defaultuicontrolBehavior'),...
'CreateFcn', {@local_CreateFcn, '', appdata} );





hsingleton = h1;





% --- Set application data first then calling the CreateFcn. 

function local_CreateFcn(hObject, eventdata, createfcn, appdata)



if ~isempty(appdata)

   names = fieldnames(appdata);

   for i=1:length(names)

       name = char(names(i));

       setappdata(hObject, name, getfield(appdata,name));

   end

end



if ~isempty(createfcn)

   eval(createfcn);

end





% --- Handles default GUIDE GUI creation and callback dispatch

function varargout = gui_mainfcn(gui_State, varargin)





%   GUI_MAINFCN provides these command line APIs for dealing with GUIs

%

%      IMPORTMETHOD_EXPORT, by itself, creates a new IMPORTMETHOD_EXPORT or raises the existing

%      singleton*.

%

%      H = IMPORTMETHOD_EXPORT returns the handle to a new IMPORTMETHOD_EXPORT or the handle to

%      the existing singleton*.

%

%      IMPORTMETHOD_EXPORT('CALLBACK',hObject,eventData,handles,...) calls the local

%      function named CALLBACK in IMPORTMETHOD_EXPORT.M with the given input arguments.

%

%      IMPORTMETHOD_EXPORT('Property','Value',...) creates a new IMPORTMETHOD_EXPORT or raises the

%      existing singleton*.  Starting from the left, property value pairs are

%      applied to the GUI before untitled_OpeningFunction gets called.  An

%      unrecognized property name or invalid value makes property application

%      stop.  All inputs are passed to untitled_OpeningFcn via varargin.

%

%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one

%      instance to run (singleton)".



%   Copyright 1984-2004 The MathWorks, Inc.

%   $Revision: 1.4 $ $Date: 2010/02/18 22:25:11 $



gui_StateFields =  {'gui_Name'

                    'gui_Singleton'

                    'gui_OpeningFcn'

                    'gui_OutputFcn'

                    'gui_LayoutFcn'

                    'gui_Callback'};

gui_Mfile = '';

for i=1:length(gui_StateFields)

    if ~isfield(gui_State, gui_StateFields{i})

        error('Could not find field %s in the gui_State struct in GUI M-file %s', gui_StateFields{i}, gui_Mfile);        

    elseif isequal(gui_StateFields{i}, 'gui_Name')

        gui_Mfile = [gui_State.(gui_StateFields{i}), '.m'];

    end

end



numargin = length(varargin);



if numargin == 0

    % IMPORTMETHOD_EXPORT

    % create the GUI

    gui_Create = 1;

elseif isequal(ishandle(varargin{1}), 1) && ispc && iscom(varargin{1}) && isequal(varargin{1},gcbo)

    % IMPORTMETHOD_EXPORT(ACTIVEX,...)    

    vin{1} = gui_State.gui_Name;

    vin{2} = [get(varargin{1}.Peer, 'Tag'), '_', varargin{end}];

    vin{3} = varargin{1};

    vin{4} = varargin{end-1};

    vin{5} = guidata(varargin{1}.Peer);

    feval(vin{:});

    return;

elseif ischar(varargin{1}) && numargin>1 && isequal(ishandle(varargin{2}), 1)

    % IMPORTMETHOD_EXPORT('CALLBACK',hObject,eventData,handles,...)

    gui_Create = 0;

else

    % IMPORTMETHOD_EXPORT(...)

    % create the GUI and hand varargin to the openingfcn

    gui_Create = 1;

end



if gui_Create == 0

    varargin{1} = gui_State.gui_Callback;

    if nargout

        [varargout{1:nargout}] = feval(varargin{:});

    else

        feval(varargin{:});

    end

else

    if gui_State.gui_Singleton

        gui_SingletonOpt = 'reuse';

    else

        gui_SingletonOpt = 'new';

    end

    

    % Open fig file with stored settings.  Note: This executes all component

    % specific CreateFunctions with an empty HANDLES structure.

    

    % Do feval on layout code in m-file if it exists

    if ~isempty(gui_State.gui_LayoutFcn)

        gui_hFigure = feval(gui_State.gui_LayoutFcn, gui_SingletonOpt);

        % openfig (called by local_openfig below) does this for guis without

        % the LayoutFcn. Be sure to do it here so guis show up on screen.

        movegui(gui_hFigure,'onscreen')

    else

        gui_hFigure = local_openfig(gui_State.gui_Name, gui_SingletonOpt);            

        % If the figure has InGUIInitialization it was not completely created

        % on the last pass.  Delete this handle and try again.

        if isappdata(gui_hFigure, 'InGUIInitialization')

            delete(gui_hFigure);

            gui_hFigure = local_openfig(gui_State.gui_Name, gui_SingletonOpt);            

        end

    end

    

    % Set flag to indicate starting GUI initialization

    setappdata(gui_hFigure,'InGUIInitialization',1);



    % Fetch GUIDE Application options

    gui_Options = getappdata(gui_hFigure,'GUIDEOptions');

    

    if ~isappdata(gui_hFigure,'GUIOnScreen')

        % Adjust background color

        if gui_Options.syscolorfig 

            set(gui_hFigure,'Color', get(0,'DefaultUicontrolBackgroundColor'));

        end



        % Generate HANDLES structure and store with GUIDATA. If there is

        % user set GUI data already, keep that also.

        data = guidata(gui_hFigure);

        handles = guihandles(gui_hFigure);

        if ~isempty(handles)

            if isempty(data)

                data = handles;

            else

                names = fieldnames(handles);

                for k=1:length(names)

                    data.(char(names(k)))=handles.(char(names(k)));

                end

            end

        end

        guidata(gui_hFigure, data);

    end

    

    % If user specified 'Visible','off' in p/v pairs, don't make the figure

    % visible.

    gui_MakeVisible = 1;

    for ind=1:2:length(varargin)

        if length(varargin) == ind

            break;

        end

        len1 = min(length('visible'),length(varargin{ind}));

        len2 = min(length('off'),length(varargin{ind+1}));

        if ischar(varargin{ind}) && ischar(varargin{ind+1}) && ...

                strncmpi(varargin{ind},'visible',len1) && len2 > 1

            if strncmpi(varargin{ind+1},'off',len2)

                gui_MakeVisible = 0;

            elseif strncmpi(varargin{ind+1},'on',len2)

                gui_MakeVisible = 1;

            end

        end

    end

    

    % Check for figure param value pairs

    for index=1:2:length(varargin)

        if length(varargin) == index || ~ischar(varargin{index})

            break;

        end

        try set(gui_hFigure, varargin{index}, varargin{index+1}), catch break, end

    end



    % If handle visibility is set to 'callback', turn it on until finished

    % with OpeningFcn

    gui_HandleVisibility = get(gui_hFigure,'HandleVisibility');

    if strcmp(gui_HandleVisibility, 'callback')

        set(gui_hFigure,'HandleVisibility', 'on');

    end

    

    feval(gui_State.gui_OpeningFcn, gui_hFigure, [], guidata(gui_hFigure), varargin{:});

    

    if ishandle(gui_hFigure)

        % Update handle visibility

        set(gui_hFigure,'HandleVisibility', gui_HandleVisibility);

        

        % Make figure visible

        if gui_MakeVisible

            set(gui_hFigure, 'Visible', 'on')

            if gui_Options.singleton 

                setappdata(gui_hFigure,'GUIOnScreen', 1);

            end

        end



        % Done with GUI initialization

        rmappdata(gui_hFigure,'InGUIInitialization');

    end

    

    % If handle visibility is set to 'callback', turn it on until finished with

    % OutputFcn

    if ishandle(gui_hFigure)

        gui_HandleVisibility = get(gui_hFigure,'HandleVisibility');

        if strcmp(gui_HandleVisibility, 'callback')

            set(gui_hFigure,'HandleVisibility', 'on');

        end

        gui_Handles = guidata(gui_hFigure);

    else

        gui_Handles = [];

    end

    

    if nargout

        [varargout{1:nargout}] = feval(gui_State.gui_OutputFcn, gui_hFigure, [], gui_Handles);

    else

        feval(gui_State.gui_OutputFcn, gui_hFigure, [], gui_Handles);

    end

    

    if ishandle(gui_hFigure)

        set(gui_hFigure,'HandleVisibility', gui_HandleVisibility);

    end

end    



function gui_hFigure = local_openfig(name, singleton)



% this application data is used to indicate the running mode of a GUIDE

% GUI to distinguish it from the design mode of the GUI in GUIDE.

setappdata(0,'OpenGuiWhenRunning',1);



% openfig with three arguments was new from R13. Try to call that first, if

% failed, try the old openfig.

try 

    gui_hFigure = openfig(name, singleton, 'auto');

catch

    % OPENFIG did not accept 3rd input argument until R13,

    % toggle default figure visible to prevent the figure

    % from showing up too soon.

    gui_OldDefaultVisible = get(0,'defaultFigureVisible');

    set(0,'defaultFigureVisible','off');

    gui_hFigure = openfig(name, singleton);

    set(0,'defaultFigureVisible',gui_OldDefaultVisible);

end

rmappdata(0,'OpenGuiWhenRunning');



