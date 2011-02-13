function varargout = ModDatasetUnits(varargin)
% MODDATASETUNITS M-file for ModDatasetUnits.fig
%      MODDATASETUNITS, by itself, creates a new MODDATASETUNITS or raises the existing
%      singleton*.
%
%      H = MODDATASETUNITS returns the handle to a new MODDATASETUNITS or the handle to
%      the existing singleton*.
%
%      MODDATASETUNITS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MODDATASETUNITS.M with the given input
%      arguments.
%
%      MODDATASETUNITS('Property','Value',...) creates a new MODDATASETUNITS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ModDatasetUnits_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ModDatasetUnits_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
%
% Description of handles structure
%
%      ====graphics handles=====
%
%
%               figure1: handle to main fig
%             mainpanel: 
%              CancelPB: 
%              SaveAsPB:
%                SavePB:
%          obsInfoPanel:
%           RevertAllPB:
%               du_text:
%         DatasetUnitLB:
%       RevertCurrentPB:
%             UncModePU: 
%               UncEdit: edit box for sym unc, value must be pos
%             UncUBEdit: edit box for asym unc, value must be pos
%             UncLBEdit: edit box for asym unc, value must be neg
%              DataEdit: edit box for data
%              unc_text: handle to static string
%             data_text: handle to static string
%
%        ====additional data handles=====
%
%         existingNames: qx1 cell of chars
%                curIdx: scalar, index telling which element of existingNames/DUcell is currently viewable
%                DUcell: qx1 cell of DatasetUnit objects
%                DUinfo: qx1 struct with fields .def_d, .def_uncVect, .isdef
%                          tells us the default data/unc, and which are
%                          default for each of the q DatasetUnit objects
%                 goodd: 1
%                goodlb: 1
%                 goodu: 1
%                goodub: 1
%                 org_d: [77x1 double]
%               org_unc: {77x1 cell}
%           org_uncVect: [77x2 double]
%                 cur_d: [77x1 double]
%           cur_uncVect: [77x2 double]
%                 isdef: [77x1] bool, tracks which of cur are nondef
%               prefSym: 1
%        curLBSelection: Bookkeep current listbox selction, in case we need
%                        to revert to it if entered data is bad.
%        curPUSelection: Bookkeep current PU selction, in case we need
%                        to revert to it if entered data is bad.


% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ModDatasetUnits_OpeningFcn, ...
                   'gui_OutputFcn',  @ModDatasetUnits_OutputFcn, ...
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

% --- Executes just before ModDatasetUnits is made visible.
function ModDatasetUnits_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ModDatasetUnits (see VARARGIN)

%blow away any existing timers saved in the handles structure
if isfield(handles,'timer')
  try
    stop(handles.timer);
    %delete(handles.timer);
    handles = rmfield(handles,'timer');
  catch
    %do nothing
  end
end 

%add fields to the handles structure
handles.existingNames = varargin{1};
handles.curIdx = varargin{2};
handles.DUcell = varargin{3};
handles.DUinfo = varargin{4};

%bools to keep track of illegitimate data entry
handles.goodd = true;
handles.goodlb = true;
handles.goodu = true;
handles.goodub = true;

%save the original data and uncertainty, so we know if we need to save and
%can revert upon cancel. unc and uncVect are not redundant, because we look
%at the dimensions of org_unc to see if we should save the uncertainty with
%a single bound or upper and lower bounds in the final object. We use the
%single bound method when possible if that was originally used. Generally
%this function works with the vector version, and decides right at the end
%how to store it in the object.
handles.org_d = vertcat(handles.DUcell{handles.curIdx}.observedValue);
handles.org_unc = {handles.DUcell{handles.curIdx}.observationUncertainty}';
handles.org_uncVect = vertcat(handles.DUcell{handles.curIdx}.observationUncertaintyPlusMinus);

%save the current data and uncertainty here, rather than always modify the
%objects
handles.cur_d = handles.org_d;
handles.cur_uncVect = handles.org_uncVect;

% keep track of which exhibit nondefault values
handles.isdef = handles.DUinfo(handles.curIdx).isdef;

%initially, prefer symmetric unc display where possible
handles.prefSym = true;
handles.curLBSelection = 1;
if numel(handles.org_unc{1}) == 1 || -handles.org_unc{1}(1) == handles.org_unc{1}(2)
  handles.curPUSelection = 1;
else
  handles.curPUSelection = 2;
end

% Update handles structure
guidata(hObject, handles);

%tweak ui
str = {handles.DUcell{handles.curIdx}.name}';
if ~all(handles.isdef)
  pad = repmat({blanks(3)},size(str));
  pad(~handles.isdef) = {'*! '};
  str = cellstr([strvcat(pad{:}) strvcat(str{:})]); %#ok
  set(handles.DatasetUnitLB,'String',str)
else
  set(handles.DatasetUnitLB,'String',str)
end

set(handles.DatasetUnitLB,'value',handles.curLBSelection) %initially select the first
set(handles.figure1,'name',['Modify Dataset -- ' handles.existingNames{handles.curIdx}])

% Visibility of uncertainty edit boxes, listbox will set strings
if numel(handles.org_unc{1}) == 1 || -handles.org_unc{1}(1) == handles.org_unc{1}(2)
  set(handles.UncModePU,'Value',1)
  set([handles.UncLBEdit handles.UncUBEdit],'Visible','off')
  set(handles.UncEdit,'Visible','on')
else
  set(handles.UncModePU,'Value',2,'enable','inactive');
  set(handles.UncEdit,'Visible','off')
  set([handles.UncLBEdit handles.UncUBEdit],'Visible','on')
end

%Update the edit displays
DatasetUnitLB_Callback(handles.DatasetUnitLB,eventdata,handles)

%position ui
Units = get(handles.figure1,'Units');
set(handles.figure1,'Units','pixels');
curpos = get(handles.figure1,'Position');
scrnsize = get(0,'ScreenSize');
newpos = curpos;
newpos(1) = (scrnsize(3)-curpos(3))/2;
newpos(2) = (scrnsize(4)-curpos(4))/2;
set(handles.figure1,'Position',newpos,'Units',Units);

handles.timer = timer('executionmode','fixedrate',...
  'TimerFcn',@timer_Callback,'Period',0.1,'startdelay',0.25,'UserData',hObject);
start(handles.timer);

% Make the GUI modal
%set(handles.figure1,'WindowStyle','modal')

% UIWAIT makes DatasetUnitImport wait for user response (see UIRESUME)
uiwait(handles.figure1);

function timer_Callback(varargin)
% This is the callback of the timer. Everytime it executes it checks if the
% listbox has been moved since the previous execution (by comparing the
% current listboxtop to the saved prev. value).

timer = varargin{1};
try
  handles = guidata(timer.UserData); %.UserData is a graphics handle
catch
  return
end

LBval = get(handles.DatasetUnitLB,'Value');
uDstruct = handles.DUcell{handles.curIdx}(LBval).ResponseModel.userData;
targetUnits = uDstruct.targetUnits;
if strcmp(targetUnits,'dimensionless')
    targetUnits = '';
end
set(handles.DataEdit1,'String',[num2str(10^str2double(get(handles.DataEdit,'String')),'%0.4g') ' ' targetUnits])
set(handles.lbEdit,'String',num2str(10^str2double(get(handles.UncLBEdit,'String')),'%0.5g'))
set(handles.ubEdit,'String',num2str(10^str2double(get(handles.UncUBEdit,'String')),'%0.5g'))


guidata(timer.UserData,handles);

% --- Outputs from this function are returned to the command line.
function varargout = ModDatasetUnits_OutputFcn(hObject, eventdata, handles) %#ok
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

varargout{1} = handles.existingNames;
varargout{2} = handles.curIdx;
varargout{3} = handles.DUcell;
varargout{4} = handles.DUinfo;

% The figure can be deleted now
delete(handles.figure1);

%===Create functions===
function DataEdit_CreateFcn(hObject, eventdata, handles) %#ok
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function UncLBEdit_CreateFcn(hObject, eventdata, handles) %#ok
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function UncUBEdit_CreateFcn(hObject, eventdata, handles) %#ok
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function UncEdit_CreateFcn(hObject, eventdata, handles) %#ok
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function UncModePU_CreateFcn(hObject, eventdata, handles) %#ok
% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function DatasetUnitLB_CreateFcn(hObject, eventdata, handles) %#ok
% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%===Callback functions===
function DataEdit_Callback(hObject, eventdata, handles) %#ok
idx = get(handles.DatasetUnitLB,'Value');
newd = str2double(get(hObject,'String'));

if isnan(newd)
  set(hObject,'BackgroundColor','red')
  handles.goodd = false;
  %disable stuff; listbox will take care of itself
  set([handles.SavePB handles.SaveAsPB],'enable','off')
else
  set(hObject,'BackgroundColor','white')
  handles.cur_d(idx) = newd;
  handles.goodd = true;
  if handles.goodd && handles.goodlb && handles.goodu && handles.goodub
    set([handles.DatasetUnitLB handles.SavePB handles.SaveAsPB],'enable','on')
  end

  %Check if everything is default, and compare this to the old value. If it
  %changed, we need to update the display 
  isdef = handles.cur_d(idx) == handles.DUinfo(handles.curIdx).def_d(idx) && isequal(handles.cur_uncVect(idx,:),handles.DUinfo(handles.curIdx).def_uncVect(idx,:));
  if isdef~=handles.isdef(idx)
    handles.isdef(idx) = isdef;
    str = {handles.DUcell{handles.curIdx}.name}';
    if ~all(handles.isdef)
      pad = repmat({blanks(3)},size(str));
      pad(~handles.isdef) = {'*! '};
      str = cellstr([strvcat(pad{:}) strvcat(str{:})]); %#ok
      set(handles.DatasetUnitLB,'String',str)
    else
      set(handles.DatasetUnitLB,'String',str)
    end
  end
    
end
guidata(hObject,handles)

function UncLBEdit_Callback(hObject, eventdata, handles) %#ok

idx = get(handles.DatasetUnitLB,'Value');
newlb = str2double(get(hObject,'String'));

if isnan(newlb) || newlb >= handles.cur_uncVect(idx,2)
  set(hObject,'BackgroundColor','red')
  handles.goodlb = false;
  %disable stuff; listbox and popup will take care of themselves
  set([handles.SavePB handles.SaveAsPB],'enable','off')
else
  set(hObject,'BackgroundColor','white')
  handles.cur_uncVect(idx,1) = newlb;
  handles.goodlb = true;
  if handles.goodd && handles.goodlb && handles.goodu && handles.goodub
    set([handles.DatasetUnitLB handles.SavePB handles.SaveAsPB],'enable','on')
  end
  
  %update symmetric unc string and popup enable
  if handles.goodlb && handles.goodub
    if -newlb == handles.cur_uncVect(idx,2)
      set(handles.UncEdit,'String',num2str(-newlb))
      set(handles.UncModePU,'enable','on')
    else
      set(handles.UncEdit,'String','inf')
      set(handles.UncModePU,'enable','inactive')
    end
  end
  
  %Check if everything is default, and compare this to the old value. If it
  %changed, we need to update the display 
  isdef = handles.cur_d(idx) == handles.DUinfo(handles.curIdx).def_d(idx) && isequal(handles.cur_uncVect(idx,:),handles.DUinfo(handles.curIdx).def_uncVect(idx,:));
  if isdef~=handles.isdef(idx)
    handles.isdef(idx) = isdef;
    str = {handles.DUcell{handles.curIdx}.name}';
    if ~all(handles.isdef)
      pad = repmat({blanks(3)},size(str));
      pad(~handles.isdef) = {'*! '};
      str = cellstr([strvcat(pad{:}) strvcat(str{:})]); %#ok
      set(handles.DatasetUnitLB,'String',str)
    else
      set(handles.DatasetUnitLB,'String',str)
    end
  end
 
end
guidata(hObject,handles)

function UncUBEdit_Callback(hObject, eventdata, handles) %#ok
idx = get(handles.DatasetUnitLB,'Value');
newub = str2double(get(hObject,'String'));

if isnan(newub) || newub <= handles.cur_uncVect(idx,1)
  set(hObject,'BackgroundColor','red')
  handles.goodub = false;
  %disable stuff; listbox and popup will take care of themselves
  set([handles.SavePB handles.SaveAsPB],'enable','off')
else
  set(hObject,'BackgroundColor','white')
  handles.cur_uncVect(idx,2) = newub;
  handles.goodub = true;
  if handles.goodd && handles.goodlb && handles.goodu && handles.goodub
    set([handles.DatasetUnitLB handles.SavePB handles.SaveAsPB],'enable','on')
  end

  %update symmetric unc string and popup enable
  if handles.goodlb && handles.goodub
    if -handles.cur_uncVect(idx,1) == newub
      set(handles.UncEdit,'String',num2str(newub))
      set(handles.UncModePU,'enable','on')
    else
      set(handles.UncEdit,'String','inf')
      set(handles.UncModePU,'enable','inactive')
    end
  end
  
  %Check if everything is default, and compare this to the old value. If it
  %changed, we need to update the display 
  isdef = handles.cur_d(idx) == handles.DUinfo(handles.curIdx).def_d(idx) && isequal(handles.cur_uncVect(idx,:),handles.DUinfo(handles.curIdx).def_uncVect(idx,:));
  if isdef~=handles.isdef(idx)
    handles.isdef(idx) = isdef;
    str = {handles.DUcell{handles.curIdx}.name}';
    if ~all(handles.isdef)
      pad = repmat({blanks(3)},size(str));
      pad(~handles.isdef) = {'*! '};
      str = cellstr([strvcat(pad{:}) strvcat(str{:})]); %#ok
      set(handles.DatasetUnitLB,'String',str)
    else
      set(handles.DatasetUnitLB,'String',str)
    end
  end
end
guidata(hObject,handles)

function UncEdit_Callback(hObject, eventdata, handles) %#ok
idx = get(handles.DatasetUnitLB,'Value');
newu = str2double(get(hObject,'String'));

if isnan(newu) || newu <= 0
  set(hObject,'BackgroundColor','red')
  handles.goodu = false;
  %disable stuff; listbox and popup will take care of themselves
  set([handles.SavePB handles.SaveAsPB],'enable','off')
else
  set(hObject,'BackgroundColor','white')
  handles.cur_uncVect(idx,:) = [-newu,newu];
  handles.goodu = true;
  if handles.goodd && handles.goodlb && handles.goodu && handles.goodub
    set([handles.DatasetUnitLB handles.UncModePU, handles.SavePB handles.SaveAsPB],'enable','on')
  end
  set(handles.UncLBEdit,'String',num2str(-newu))
  set(handles.UncUBEdit,'String',num2str(newu))
  
  %Check if everything is default, and compare this to the old value. If it
  %changed, we need to update the display 
  isdef = handles.cur_d(idx) == handles.DUinfo(handles.curIdx).def_d(idx) && isequal(handles.cur_uncVect(idx,:),handles.DUinfo(handles.curIdx).def_uncVect(idx,:));
  if isdef~=handles.isdef(idx)
    handles.isdef(idx) = isdef;
    str = {handles.DUcell{handles.curIdx}.name}';
    if ~all(handles.isdef)
      pad = repmat({blanks(3)},size(str));
      pad(~handles.isdef) = {'*! '};
      str = cellstr([strvcat(pad{:}) strvcat(str{:})]); %#ok
      set(handles.DatasetUnitLB,'String',str)
    else
      set(handles.DatasetUnitLB,'String',str)
    end
  end
end
guidata(hObject,handles)

function UncModePU_Callback(hObject, eventdata, handles) %#ok
%we should only be able to call this if both symmetric and asymmetric
%displaying is possible

if ~all([handles.goodlb handles.goodu handles.goodub])
  set(hObject,'Value',handles.curPUSelection)
  set(hObject,'enable','off')
  return
end

val = get(hObject,'Value');
handles.curPUSelection = val;
if val == 1
  set([handles.UncLBEdit handles.UncUBEdit],'Visible','off')
  set(handles.UncEdit,'Visible','on')
  handles.prefSym = true;
else
  set(handles.UncEdit,'Visible','off')
  set([handles.UncLBEdit handles.UncUBEdit],'Visible','on')
  handles.prefSym = false;
end
guidata(hObject,handles)

function RevertCurrentPB_Callback(hObject, eventdata, handles) %#ok
idx = get(handles.DatasetUnitLB,'Value');
handles.cur_d(idx) = handles.DUinfo(handles.curIdx).def_d(idx);
handles.cur_uncVect(idx,:) = handles.DUinfo(handles.curIdx).def_uncVect(idx,:);
handles.isdef(idx) = true;
handles.goodd = true;
handles.goodu = true;
handles.goodlb = true;
handles.goodub = true;

set([handles.DataEdit handles.UncLBEdit handles.UncEdit handles.UncUBEdit],'BackgroundColor','white')

% The listbox call will take care of the UncModePU enable property
set([handles.DatasetUnitLB handles.SavePB handles.SaveAsPB],'enable','on')

str = {handles.DUcell{handles.curIdx}.name}';
if ~all(handles.isdef)
  pad = repmat({blanks(3)},size(str));
  pad(~handles.isdef) = {'*! '};
  str = cellstr([strvcat(pad{:}) strvcat(str{:})]); %#ok
  set(handles.DatasetUnitLB,'String',str)
else
  set(handles.DatasetUnitLB,'String',str)
end

%save and call the listbox callback to refresh displays
guidata(hObject,handles)
DatasetUnitLB_Callback(handles.DatasetUnitLB, eventdata, handles)

function RevertAllPB_Callback(hObject, eventdata, handles) %#ok
handles.cur_d = handles.DUinfo(handles.curIdx).def_d;
handles.cur_uncVect = handles.DUinfo(handles.curIdx).def_uncVect;
handles.isdef(:) = true;

handles.goodd = true;
handles.goodu = true;
handles.goodlb = true;
handles.goodub = true;

set([handles.DataEdit handles.UncLBEdit handles.UncEdit handles.UncUBEdit],'BackgroundColor','white')

% The listbox call will take care of the UncModePU enable property
set([handles.DatasetUnitLB handles.SavePB handles.SaveAsPB],'enable','on')

str = {handles.DUcell{handles.curIdx}.name}';
set(handles.DatasetUnitLB,'String',str)

%call the listbox callback to refresh displays
guidata(hObject,handles)
DatasetUnitLB_Callback(handles.DatasetUnitLB, eventdata, handles)

function DatasetUnitLB_Callback(hObject, eventdata, handles) %#ok

if ~all([handles.goodd handles.goodlb handles.goodu handles.goodub])
  set(hObject,'Value',handles.curLBSelection)
  set(hObject,'enable','off')
  return
end

%here we update the edit boxes
val = get(hObject,'Value');
handles.curLBSelection = val;

set(handles.UncLBEdit,'String',num2str(handles.cur_uncVect(val,1)))
set(handles.UncUBEdit,'String',num2str(handles.cur_uncVect(val,2)))
set(handles.DataEdit,'String',num2str(handles.cur_d(val)))

LBval = get(handles.DatasetUnitLB,'Value');
uDstruct = handles.DUcell{handles.curIdx}(LBval).ResponseObservation.userData;
targetLabel = uDstruct.targetLabel;
targetUnits = uDstruct.targetUnits;
if strcmp(targetUnits,'dimensionless')
    targetUnits = '';
end
set(handles.DataEdit1,'String',[num2str(10^handles.cur_d(val),'%0.4g') ' ' targetUnits])
set(handles.lbEdit,'String',num2str(10^handles.cur_uncVect(val,1)))
set(handles.ubEdit,'String',num2str(10^handles.cur_uncVect(val,2)))
set(handles.labelEdit,'String',targetLabel)

if -handles.cur_uncVect(val,1) == handles.cur_uncVect(val,2)
  symPoss = true;
  set(handles.UncEdit,'String',num2str(handles.cur_uncVect(val,2)))
else
  symPoss = false;
  set(handles.UncEdit,'String','inf')
end

PUval = get(handles.UncModePU,'value');

%if we'd rather use sym, we can, and we're not now, make the change
if PUval == 2 && symPoss && handles.prefSym
  set(handles.UncModePU,'Value',1)
  set([handles.UncLBEdit handles.UncUBEdit],'Visible','off')
  set(handles.UncEdit,'Visible','on')
end

%if we're using sym, but we can't, make the change
if PUval == 1 && ~symPoss
  set(handles.UncModePU,'Value',2)
  set([handles.UncLBEdit handles.UncUBEdit],'Visible','on')
  set(handles.UncEdit,'Visible','off')
end

%if we can't use sym, disable popup, but do this after we mess with it
if symPoss
  set(handles.UncModePU,'enable','on')
else
  set(handles.UncModePU,'enable','inactive')
end

%check for double click, which means we should display further details
if strcmp(get(handles.figure1,'SelectionType'),'open')
  LBval = get(hObject,'Value');
  
  %do nothing if the user double clicked and dragged at the same time
  if length(LBval) == 1
    %TODO: This doesn't seem to work
    %hand = handles.DUcell{handles.curIdx}(LBval).dispCallback;
    %if a callback is available, invoke it
    if ~isempty(hand)
      hand(handles.DUcell{handles.curIdx}(LBval));
    end
  end
end
guidata(hObject,handles)

function SavePB_Callback(hObject, eventdata, handles) %#ok
%update the object, return
handles = updateValues(handles);
guidata(hObject,handles)
uiresume(handles.figure1)

function SaveAsPB_Callback(hObject, eventdata, handles) %#ok
%get a new Dataset name
newName = GetNewName({'Please Supply New Dataset Name'},'User Input',handles.existingNames,handles.existingNames{handles.curIdx});

if isempty(newName)
  return
else
  %create new DatasetUnit object
  handles.existingNames{end+1,1} = newName;
  handles.DUcell{end+1,1} = handles.DUcell{handles.curIdx};
  handles.DUinfo(end+1,1) = handles.DUinfo(handles.curIdx);
  handles.curIdx = length(handles.existingNames);
  
  %place the current data and uncertainty in the DatasetUnit object, and
  %update its isdef flags.
  handles = updateValues(handles);
  
  %cache the "original" data and uncertainty of this new dataset
  handles.org_d = handles.cur_d;
  handles.org_unc = {handles.DUcell{handles.curIdx}.observationUncertainty}';
  handles.org_uncVect = handles.cur_uncVect;
  
  set(handles.figure1,'name',['Dataset: ' handles.existingNames{handles.curIdx}])
  guidata(hObject,handles)
end

function CancelPB_Callback(hObject, eventdata, handles) %#ok
%replace values with the originals, update the object, and return
handles.cur_d = handles.org_d;
handles.cur_uncVect = handles.org_uncVect;
tmp1 = handles.DUinfo(handles.curIdx).def_d == handles.cur_d;
tmp2 = handles.DUinfo(handles.curIdx).def_uncVect == handles.cur_uncVect;
handles.isdef = tmp1 & tmp2(:,1) & tmp2(:,2);

handles = updateValues(handles);
guidata(hObject,handles)
uiresume(handles.figure1)

function handles = updateValues(handles)
%function to take all values in handles.cur_d and handles.cur_uncVect, and
%use them to modify the object in handles.DUcell{handles.curIdx}
DU = handles.DUcell{handles.curIdx};
m = nPairs(DU);
for i1 = 1:m
  EA = DU(i1).ResponseObservation;

  %if symmetric uncertainty was originally used and that representation is
  %still ok, use it.
  if numel(handles.org_unc{i1}) == 1 && -handles.cur_uncVect(i1,1) == handles.cur_uncVect(i1,2)
    EA.uncertainty = handles.cur_uncVect(i1,2);
  else
    EA.uncertainty = handles.cur_uncVect(i1,:);
  end
  EA.observedValue = handles.cur_d(i1);
  DU(i1).ResponseObservation = EA;
end
handles.DUcell{handles.curIdx} = DU;
handles.DUinfo(handles.curIdx).isdef = handles.isdef;

function figure1_CloseRequestFcn(hObject, eventdata, handles) %#ok

curobj = get(handles.figure1,'CurrentObject');
if ismember(curobj,[handles.UncLBEdit handles.UncUBEdit handles.UncEdit handles.DataEdit]) 
  %the user hitting exit does NOT trigger the edit callbacks. thus we call
  %small gui with a uicontrol in the center of the screen to flush all call
  %edit callbacks.
  scrnsz = get(0,'screensize');
  h = dummy('Position',[scrnsz(3) scrnsz(4) 1 1],'Units','pixels');
  delete(h)
  %update handles
  handles = guidata(hObject);
end

if isequal(get(handles.figure1, 'waitstatus'), 'waiting')
  if ~(handles.goodd && handles.goodlb && handles.goodu && handles.goodub)
    %don't exit: there are problems
  elseif ~isequal(handles.cur_d,handles.org_d) || ~isequal(handles.cur_uncVect,handles.org_uncVect)
    ButtonName = savedlg('string',['Do you want to save changes to Dataset ' handles.existingNames{handles.curIdx}],'title','User Input');
    switch ButtonName,
     case 'Save Changes and Return'
      SavePB_Callback(handles.SavePB,eventdata,handles)
     case 'Save As New'
      SaveAsPB_Callback(handles.SaveAsPB,eventdata,handles)
     case 'Discard Changes and Return'
      CancelPB_Callback(handles.CancelPB,eventdata,handles)
     case 'Cancel'
      return
     otherwise
      %keyboard
    end % switch
  else
    % just close it
    CancelPB_Callback(handles.CancelPB,eventdata,handles)
  end

else
    % The GUI is no longer waiting, something bombed, just close it
    delete(handles.figure1);
end


% --------------------------------------------------------------------
function target_menu_Callback(hObject, eventdata, handles)
% hObject    handle to ShowXML_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
LBval = get(handles.DatasetUnitLB,'Value');
%do nothing if the user double clicked and dragged at the same time
if length(LBval) == 1
    uDstruct = handles.DUcell{handles.curIdx}(LBval).ResponseObservation.userData;
    pid = uDstruct.targetID;
    if iscell(pid) && isscalar(pid)
        uDstruct.callback2show(pid{1});
    elseif iscell(pid)
        uDstruct.callback2show(pid{1});
    elseif ischar(pid)
        uDstruct.callback2show(pid);
    end
end


% --------------------------------------------------------------------
function target_links_menu_Callback(hObject, eventdata, handles)
% hObject    handle to target_links_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
LBval = get(handles.DatasetUnitLB,'Value');

%do nothing if the user double clicked and dragged at the same time
if length(LBval) == 1
    uDstruct = handles.DUcell{handles.curIdx}(LBval).ResponseObservation.userData;
    pid = uDstruct.targetLinks;
    if iscell(pid) && isscalar(pid)
        uDstruct.callback2show(pid{1});
    elseif iscell(pid)
        uDstruct.callback2show(pid{1});
    elseif ischar(pid)
        uDstruct.callback2show(pid);
    end
    
end


% --------------------------------------------------------------------
function model_menu_Callback(hObject, eventdata, handles)
% hObject    handle to model_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
LBval = get(handles.DatasetUnitLB,'Value');

%do nothing if the user double clicked and dragged at the same time
if length(LBval) == 1
    uDstruct = handles.DUcell{handles.curIdx}(LBval).ResponseModel.userData;
    pid = uDstruct.modelID;
    if iscell(pid) && isscalar(pid)
        uDstruct.callback2show(pid{1});
    elseif iscell(pid)
        uDstruct.callback2show2(pid{1},pid{2});
    elseif ischar(pid)
        uDstruct.callback2show(pid);
    end
end


% --------------------------------------------------------------------
function Exp_menu_Callback(hObject, eventdata, handles)
% hObject    handle to Exp_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function ShowXML_menu_Callback(hObject, eventdata, handles)
% hObject    handle to ShowXML_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function DataEdit1_Callback(hObject, eventdata, handles)
% hObject    handle to DataEdit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DataEdit1 as text
%        str2double(get(hObject,'String')) returns contents of DataEdit1 as a double


% --- Executes during object creation, after setting all properties.
function DataEdit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DataEdit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function lbEdit_Callback(hObject, eventdata, handles)
% hObject    handle to lbEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lbEdit as text
%        str2double(get(hObject,'String')) returns contents of lbEdit as a double


% --- Executes during object creation, after setting all properties.
function lbEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lbEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ubEdit_Callback(hObject, eventdata, handles)
% hObject    handle to ubEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ubEdit as text
%        str2double(get(hObject,'String')) returns contents of ubEdit as a double


% --- Executes during object creation, after setting all properties.
function ubEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ubEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function labelEdit_Callback(hObject, eventdata, handles)
% hObject    handle to labelEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of labelEdit as text
%        str2double(get(hObject,'String')) returns contents of labelEdit as a double


% --- Executes during object creation, after setting all properties.
function labelEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to labelEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


