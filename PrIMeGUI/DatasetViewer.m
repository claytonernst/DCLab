function varargout = DatasetViewer(varargin)

% DATASETVIEWER M-file for DatasetViewer.fig

%      DATASETVIEWER, by itself, creates a new DATASETVIEWER or raises the existing

%      singleton*.

%

%      H = DATASETVIEWER returns the handle to a new DATASETVIEWER or the handle to

%      the existing singleton*.

%

%      DATASETVIEWER('CALLBACK',hObject,eventData,handles,...) calls the local

%      function named CALLBACK in DATASETVIEWER.M with the given input arguments.

%

%      DATASETVIEWER('Property','Value',...) creates a new DATASETVIEWER or raises the

%      existing singleton*.  Starting from the left, property value pairs are

%      applied to the GUI before DatasetViewer_OpeningFunction gets called.  An

%      unrecognized property name or invalid value makes property application

%      stop.  All inputs are passed to DatasetViewer_OpeningFcn via varargin.

%

%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one

%      instance to run (singleton)".

%

% See also: GUIDE, GUIDATA, GUIHANDLES



% Edit the above text to modify the response to help DatasetViewer



% Last Modified by GUIDE v2.5 21-Apr-2007 01:28:15



% Begin initialization code - DO NOT EDIT

gui_Singleton = 1;

gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @DatasetViewer_OpeningFcn, ...
                   'gui_OutputFcn',  @DatasetViewer_OutputFcn, ...
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





% --- Executes just before DatasetViewer is made visible.

function DatasetViewer_OpeningFcn(hObject, eventdata, handles, varargin)

% This function has no output args, see OutputFcn.

% hObject    handle to figure

% eventdata  reserved - to be defined in a future version of MATLAB

% handles    structure with handles and user data (see GUIDATA)

% varargin   command line arguments to DatasetViewer (see VARARGIN)



% Choose default command line output for DatasetViewer

handles.output = hObject;



if ~isempty(varargin)

  Dset = varargin{1};

else

  Dset = DClab.DCDataset;

end



handles.Units = Dset.ModelAndObservationPair;

handles.d = vertcat(handles.Units.observedValue);

handles.unc = {handles.Units.observationUncertainty}'; 

handles.uncVect = vertcat(handles.Units.uncVect);



handles.ParamAssns = Dset.FreeParameter;

handles.xNom = {handles.ParamAssns.nominal}';

handles.xBnds = vertcat(handles.ParamAssns.range);



%initially, prefer symmetric unc display where possible

handles.prefSym = true;



% Update handles structure

guidata(hObject, handles);



%tweak ui

set(handles.modeledObsLB,'value',1) %initially select the first

set(handles.uncParamLB,'value',1)



str = {Dset.ModelAndObservationPair.name}';

if ~isempty(str)

  prefix = repmat({blanks(4)},size(str));

  for i1 = 1:length(str)

    num = num2str(i1);

    prefix{i1}(1,1:length(num)) = num;

  end

  str = cellstr([strvcat(prefix{:}) strvcat(str{:})]); %#ok

end

set(handles.modeledObsLB,'String',str)



str = {Dset.FreeParameter.name}';

if ~isempty(str)

  prefix = repmat({blanks(4)},size(str));

  for i1 = 1:length(str)

    num = num2str(i1);

    prefix{i1}(1,1:length(num)) = num;

  end

  str = cellstr([strvcat(prefix{:}) strvcat(str{:})]); %#ok

end

set(handles.uncParamLB,'String',str)



% Visibility of uncertainty edit boxes, listbox will set strings

if isempty(handles.unc{1}) || numel(handles.unc{1}) == 1 || -handles.unc{1}(1) == handles.unc{1}(2)

  set(handles.uncModePU,'Value',1)

  set([handles.uncLBEdit handles.uncUBEdit],'Visible','off')

  set(handles.uncEdit,'Visible','on')

else

  set(handles.uncModePU,'Value',2,'enable','inactive');

  set(handles.uncEdit,'Visible','off')

  set([handles.uncLBEdit handles.uncUBEdit],'Visible','on')

end



%Update the edit displays

modeledObsLB_Callback(handles.modeledObsLB,eventdata,handles)

uncParamLB_Callback(handles.uncParamLB,eventdata,handles)



%position ui

guiunits = get(handles.figure1,'Units');

set(handles.figure1,'Units','pixels');

curpos = get(handles.figure1,'Position');

scrnsize = get(0,'ScreenSize');

newpos = curpos;

newpos(1) = (scrnsize(3)-curpos(3))/2;

newpos(2) = (scrnsize(4)-curpos(4))/2;

set(handles.figure1,'Position',newpos,'Units',guiunits);





% UIWAIT makes DatasetViewer wait for user response (see UIRESUME)

% uiwait(handles.figure1);





% --- Outputs from this function are returned to the command line.

function varargout = DatasetViewer_OutputFcn(hObject, eventdata, handles) %#ok

% varargout  cell array for returning output args (see VARARGOUT);

% hObject    handle to figure

% eventdata  reserved - to be defined in a future version of MATLAB

% handles    structure with handles and user data (see GUIDATA)



% Get default command line output from handles structure

varargout{1} = handles.output;





% --- Create functions. Execute during object creation, after setting all properties.

function Edit_CreateFcn(hObject, eventdata, handles) %#ok

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))

    set(hObject,'BackgroundColor','white');

end



function uncModePU_CreateFcn(hObject, eventdata, handles) %#ok

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))

    set(hObject,'BackgroundColor','white');

end



function modeledObsLB_CreateFcn(hObject, eventdata, handles) %#ok

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))

    set(hObject,'BackgroundColor','white');

end



function uncParamLB_CreateFcn(hObject, eventdata, handles) %#ok

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))

    set(hObject,'BackgroundColor','white');

end



% --- Callback functions.



% All the edit callbacks are unused since they are inactive. I kept them

% here incase we later want to use this gui to edit values.

function dataEdit_Callback(hObject, eventdata, handles) %#ok

function uncLBEdit_Callback(hObject, eventdata, handles) %#ok

function uncUBEdit_Callback(hObject, eventdata, handles) %#ok

function uncEdit_Callback(hObject, eventdata, handles) %#ok

function xNomEdit_Callback(hObject, eventdata, handles) %#ok

function xLBEdit_Callback(hObject, eventdata, handles) %#ok

function xUBEdit_Callback(hObject, eventdata, handles) %#ok



function uncModePU_Callback(hObject, eventdata, handles) %#ok

val = get(hObject,'Value');

if val == 1

  set([handles.uncLBEdit handles.uncUBEdit],'Visible','off')

  set(handles.uncEdit,'Visible','on')

  handles.prefSym = true;

else

  set(handles.uncEdit,'Visible','off')

  set([handles.uncLBEdit handles.uncUBEdit],'Visible','on')

  handles.prefSym = false;

end

guidata(hObject,handles)



function modeledObsLB_Callback(hObject, eventdata, handles) %#ok



%here we update the edit boxes

val = get(hObject,'Value');



if isempty(handles.uncVect)

  %assume everything is empty

  set(handles.uncEdit,'String','')

  set(handles.uncLBEdit,'String','')

  set(handles.uncUBEdit,'String','')

  set(handles.dataEdit,'String','')

else

  set(handles.uncLBEdit,'String',num2str(handles.uncVect(val,1)))

  set(handles.uncUBEdit,'String',num2str(handles.uncVect(val,2)))

  set(handles.dataEdit,'String',num2str(handles.d(val)))



  if -handles.uncVect(val,1) == handles.uncVect(val,2)

    symPoss = true;

    set(handles.uncEdit,'String',num2str(handles.uncVect(val,2)))

  else

    symPoss = false;

    set(handles.uncEdit,'String','inf') %this will be hidden from view

  end



  PUval = get(handles.uncModePU,'value');



  %if we'd rather use sym, we can, and we're not now, make the change

  if PUval == 2 && symPoss && handles.prefSym

    set(handles.uncModePU,'Value',1)

    set([handles.uncLBEdit handles.uncUBEdit],'Visible','off')

    set(handles.uncEdit,'Visible','on')

  end



  %if we're using sym, but we can't, make the change

  if PUval == 1 && ~symPoss

    set(handles.uncModePU,'Value',2)

    set([handles.uncLBEdit handles.uncUBEdit],'Visible','on')

    set(handles.uncEdit,'Visible','off')

  end



  %if we can't use sym, disable popup, but do this after we mess with it

  if symPoss

    set(handles.uncModePU,'enable','on')

  else

    set(handles.uncModePU,'enable','inactive')

  end

end



%check for double click, which means we should display further details

if strcmp(get(handles.figure1,'SelectionType'),'open')

  LBval = get(hObject,'Value');

  

  %do nothing if the user double clicked and dragged at the same time

  if length(LBval) == 1

    hand = handles.Units(LBval).dispCallback;

    %if a callback is available, invoke it

    if ~isempty(hand)

      hand(handles.Units(LBval));

    end

  end

end

guidata(hObject,handles)



function uncParamLB_Callback(hObject, eventdata, handles) %#ok

%here we update the edit boxes

val = get(hObject,'Value');



if isempty(handles.xNom)

  %assume everything is empty

  set(handles.xNomEdit,'String','')

  set(handles.xLBEdit,'String','')

  set(handles.xUBEdit,'String','')

else

  set(handles.xNomEdit,'String',num2str(handles.xNom(val)))

  set(handles.xLBEdit,'String',num2str(handles.xBnds(val,1)))

  set(handles.xUBEdit,'String',num2str(handles.xBnds(val,2)))

end



%check for double click, which means we should display further details

if strcmp(get(handles.figure1,'SelectionType'),'open')

  LBval = get(hObject,'Value');

  

  %do nothing if the user double clicked and dragged at the same time

  if length(LBval) == 1

    hand = handles.ParamAssns(LBval).dispCallback;

    %if a callback is available, invoke it

    if ~isempty(hand)

      hand(handles.ParamAssns(LBval));

    end

  end

end

guidata(hObject,handles)





