function varargout = FitBuilder1(varargin)
%
% FitBuilder1(handles) takes an existing handles structure and fills out
% the gui with any available data/runs
%
% FitBuilder1(ResponseModel,FreeParameter) launches the gui.
%
% FITBUILDER1 M-file for FitBuilder1.fig
%      FITBUILDER1, by itself, creates a new FITBUILDER1 or raises the existing
%      singleton*.
%
%      H = FITBUILDER1 returns the handle to a new FITBUILDER1 or the handle to
%      the existing singleton*.
%
%      FITBUILDER1('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FITBUILDER1.M with the given input arguments.
%
%      FITBUILDER1('Property','Value',...) creates a new FITBUILDER1 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before FitBuilder1_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to FitBuilder1_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run2d (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Description of handles structure
%
%             .MA: the model assertion
%             .PA: the parameter assertions
%           .range: nx2 matrix of bounds on the model parameters
%            .nom: nx1 vector containing the mean value of the bounds
%           .twoD: structure containing results from all 2-d fits
%         .global: structure containing results from all global fits
%           .ymin: scalar, current minimum value of the model over its domain
%           .ymax: scalar, current maximum value of the model over its domain
%           .logY: boolean, equals 1 if logY transformations are possible
%                  (min of model positive)
%    .avgEvalTime: scalar, current average model evaluation time
%              .N: scalar, number of sample points per dimension for the 2-D fits
%
%
% Descriptions of handles.twoD:
%   This is a column structure array with an element for each 2-D run. It
%   has the following fields:
%         .settings: scalar, the number of sample points per dimension
%   .surrogateCoeff: structure with fields .inf.t1, ... , .two.t4. Each field
%                    contains a square cell array. Only the upper triangular
%                    portion is filled out, and it contains the coefficients of the
%                    corresponding fit. For example. .C.inf.t1{2,3} contains the
%                    coefficients of the fit in x2 and x3.
%          .storedY: structure with fields .t1, .t2, .t3, and .t4. Each field
%                    contains a square cell array. Only the upper triangular
%                    portion is filled out, and it contains y values used to
%                    develop the corresponding fit. These are in a vector
%                    for a 1-D fit (diagonal portion) or a matrix suitable
%                    for plotting with surf for a 2-D fit.

% Last Modified by GUIDE v2.5 27-Jun-2007 16:25:18

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @FitBuilder1_OpeningFcn, ...
                   'gui_OutputFcn',  @FitBuilder1_OutputFcn, ...
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


% --- Executes just before FitBuilder1 is made visible.
function FitBuilder1_OpeningFcn(hObject, eventdata, handles, varargin) %#ok
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to FitBuilder1 (see VARARGIN)

% Choose default command line output for FitBuilder1
handles.output = hObject;

switch length(varargin)
  case 1
    %received a saved handles structure
    tmp = varargin{1};
    
    %establish problem data
    handles.MA = tmp.MA;
    handles.PA = tmp.PA;
    handles.range = tmp.range;
    handles.nom = tmp.nom;
    handles.twoD = tmp.twoD;
    handles.global = tmp.global;
    handles.ymin = tmp.ymin;
    handles.ymax = tmp.ymax;
    handles.ynom = tmp.ynom;
    handles.avgEvalTime = tmp.avgEvalTime;
    handles.logY = tmp.logY;
    handles.variesIdx = tmp.variesIdx;
    handles.activeIdx = tmp.activeIdx;

  case 2
    %received a model assertion and parameter assertions
    
    %establish problem data
    handles.MA = varargin{1};
    handles.PA = varargin{2};
    tempbnds = vertcat(handles.PA.range);
    tempnom = mean(tempbnds,2);

    [trash idx1 idx2] = intersect({handles.PA.name}',handles.MA.parameterList); %#ok
    if length(trash) < length(handles.MA.parameterList)
      error('Insufficient parameters in FreeParameters')
    end
    
    handles.range(idx2,:) = tempbnds(idx1,:);
    handles.nom(idx2,1) = tempnom(idx1,1);
    handles.twoD = [];
    handles.global = [];
    handles.ymin = inf;
    handles.ymax = -inf;

    %establish baseline model evaluation time
    tym = cputime;
    handles.ynom = eval(handles.MA,handles.nom);
    handles.avgEvalTime =cputime-tym;
    
    %do we think y is alway positive?
    if handles.ynom > 0
      handles.logY = true;
    else
      handles.logY = false;
    end
    tmp = (1:size(handles.range,1))';
    handles.variesIdx = tmp(diff(handles.range,[],2) > 0);
    handles.activeIdx = handles.variesIdx;
    
  otherwise
    %no inputs
    handles.MA = ResponseModel;
    handles.PA = DClab.FreeParameter;
end

% --- tweak ui ---

%set slider dimensions so there are n available positions, where n = number
%of model parameters with nonsingleton range
n = length(handles.variesIdx);
set(handles.x1Slide,'value',1,'max',n,'min',1,'SliderStep',inv(n-1)*[1 1]);
set(handles.x2Slide,'value',1,'max',n,'min',1,'SliderStep',inv(n-1)*[1 1]);

%set the pop-ups relevent to the 2D to default values
set(handles.twoDNormSelectPU,'value',1)
set(handles.twoDPntsPerDimPU,'value',1)
set(handles.twoDTransPU,'value',1)
set(handles.Fig1DispOptPU,'value',1)

%refresh the 2D display
refreshTwoDDisplay(handles.twoDNormSelectPU, eventdata, handles)

%set pop-ups relevent to the global to default values
set(handles.availGlobalRunsPU,'Value',1)
set(handles.methodSelectPU,'Value',1)

% %call a uicontrol to refresh the display since we've messed with them in software
% %uicontrol(handles.availGlobalRunsPU)
% %refresh(handles.figure1)

if length(handles.global) > 0
  set(handles.availGlobalRunsPU,'String',{'Current work' , handles.global.settings}')
  %call the popup to fill out the edit boxes with strings corresponding to
  %these results and pull off the current active parameters.
  availGlobalRunsPU_Callback(handles.availGlobalRunsPU,eventdata,handles)
else
  set(handles.availGlobalRunsPU,'String','Current work')
  set(handles.Npnts4FitM,'string','5')
  set(handles.Npnts4ValidM,'string','3')
  set(handles.Nrestart,'string','0')
end

%refresh the global display
refreshGlobalDisplay(handles.Npnts4FitM, eventdata, handles)

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes FitBuilder1 wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = FitBuilder1_OutputFcn(hObject, eventdata, handles) %#ok
% varargout  cell array for returning output args (see VARARGOUT);
varargout{1} = [];

% --- Create Functions: execute during object creation, after setting all
% properties
function Pulldown_CreateFcn(hObject, eventdata, handles) %#ok
% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function Slide_CreateFcn(hObject, eventdata, handles) %#ok
% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
function Edit_CreateFcn(hObject, eventdata, handles) %#ok
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on slider movement.
function x1Slide_Callback(hObject, eventdata, handles) %#ok
v1 = get(handles.x1Slide,'value');
if get(handles.SyncSlidersCB,'value')
  v2 = v1;
  set(handles.x2Slide,'value',v2);
end
refreshTwoDDisplay(hObject, eventdata, handles)

function x2Slide_Callback(hObject, eventdata, handles)
v2 = get(handles.x2Slide,'value');
if get(handles.SyncSlidersCB,'value')
  v1 = v2;
  set(handles.x1Slide,'value',v1);
end
refreshTwoDDisplay(hObject, eventdata, handles)

function refreshTwoDDisplay(hObject, eventdata, handles) %#ok
% This function is called after a change to a slider value, a change to a
% pop-up in the 2D settings, or a new 2D run. 
% The function
% a) determines if results are available for the current settings
% b) indicates they are not available or displays, as appropriate

% If results for the current selected options are not available, make the
% run button red to alert the user to press it.

[fitNorm trans resultsIdx] = findTwoDFields(handles);

if isempty(resultsIdx)
  set(handles.Run2D,'BackgroundColor','red')
  set(handles.oneDactivePB,'enable','off')
  set(handles.twoDactivePB,'enable','off')
  set([handles.ToFigPB1 handles.ToFigPB2 handles.ToFigPB3],'enable','off')
  set([handles.twoDAvgErrSpecPB handles.twoDPeakErrSpecPB],'enable','off')
  set([handles.twoDNormSelectPU handles.twoDTransPU],'enable','off')
  
  set(handles.firstVar,'String',' ')
  set(handles.secondVar,'String',' ')
  
  set(handles.twoDPeakErr_text,'String','N.A.')
  set(handles.twoDAvgErr_text,'String','N.A.')  
  set(handles.twoDVars_text,'String','N.A.')
  set(handles.twoDQ_text,'String','N.A.')
  
  %clear plots
  axes(handles.axes1);
  plot(0,0)
  axes(handles.axes2);
  plot(0,0);
  axes(handles.axes3);
  plot(0,0);
else
  
  % If we got this far, results must be available. Use them
  set(handles.Run2D,'BackgroundColor',[0.831373 0.815686 0.784314]);
  set([handles.oneDactivePB handles.twoDactivePB],'enable','on')
  set([handles.ToFigPB1 handles.ToFigPB2 handles.ToFigPB3],'enable','on')
  set([handles.twoDAvgErrSpecPB handles.twoDPeakErrSpecPB],'enable','on')
  set([handles.twoDNormSelectPU handles.twoDTransPU],'enable','on')
  
  % Bunch of monkey business to determine the pop-up list of available
  % transformations and keep the previous setting when possible.
  cord1 = handles.variesIdx(round(get(handles.x1Slide,'value')));
  cord2 = handles.variesIdx(round(get(handles.x2Slide,'value')));
  
  if handles.range(cord1,1) > 0 && handles.range(cord2,1) > 0
    logX = true;
  else
    logX = false;
  end
  logY = handles.logY;
  
  currAvailTrans = get(handles.twoDTransPU,'String');
  t = {'linXlinY';'logXlinY';'linXlogY';'logXlogY'}; %considered transformations
  newAvailTrans = t([true logX logY logX&&logY]);
  
  if ~isequal(currAvailTrans,newAvailTrans)
    %we need to do something.
    currentTrans = currAvailTrans{get(handles.twoDTransPU,'Value')};
    newVal = strmatch(currentTrans,newAvailTrans,'exact');
    if isempty(newVal)
      set(handles.twoDTransPU,'Value',1) %linXlinY is alway an acceptible options
    else
      set(handles.twoDTransPU,'Value',newVal)
    end
    set(handles.twoDTransPU,'String',newAvailTrans)
  end
  % End monkey business
  
  % Get the x's, yM and yS for the present case. Since some results are
  % available, none of the returned values should be empty since. These
  % should be vectors for the 1d case and matrices for the 2d case.
  [cord1 cord2 yM yS x1m x2m] = getPlotInfo(handles);
  
  set(handles.firstVar,'string',int2str(cord1));
  set(handles.secondVar,'string',int2str(cord2));
  
  if cord1==cord2
    % 1D case. Plot y^M and y^S on the same plot
    axes(handles.axes1);
    plot(x1m,yM,'b',x1m,yS,'r')
    axes(handles.axes2);
    plot(0,0);
    axes(handles.axes3);
    err = yS-yM;
    plot(x1m,err);

  else
    % Must be 2D case
    axes(handles.axes1);
    surf(x1m,x2m,yM);
    axes(handles.axes2);
    surf(x1m,x2m,yS);
    axes(handles.axes3);
    err = yS-yM;
    surf(x1m,x2m,err);
  end
  high = max(max(err));
  low = min(min(err));

  set(handles.twoDPeakErr_text,'string',[ '[' num2str(low) ', ' num2str(high) ']']);
  set(handles.twoDAvgErr_text,'string',num2str(mean(mean(abs(err)))));
  set(handles.twoDVars_text,'string',[num2str(cord1) ' and ' num2str(cord2)]);
  ymin = min(min(yM));
  ymax = max(max(yM));

  %compute the fit quality. It is a function of the peak error and the total
  %swing of this slice of the function
  if ymax-ymin <= 0
    set(handles.twoDQ_text,'String',num2str(1))
  else
    set(handles.twoDQ_text,'string',num2str(1-(high-low)/(ymax-ymin)))
  end
end

function ToFigPB1_Callback(hObject, eventdata, handles) %#ok
% This push button should be disabled unless there are results to view

%get surfaces and x vectors/matrices
[cord1 cord2 yM yS x1m x2m] = getPlotInfo(handles);

figure
if cord1==cord2
  plot(x1m,yM,'b',x1m,yS,'r')
else
  surf(x1m,x2m,yM);
end

function ToFigPB2_Callback(hObject, eventdata, handles) %#ok
% This push button should be disabled unless there are results to view

[cord1 cord2 yM yS x1m x2m] = getPlotInfo(handles);

if cord1~=cord2
  figure
  surf(x1m,x2m,yS);
end

function ToFigPB3_Callback(hObject, eventdata, handles) %#ok
% This push button should be disabled unless there are results to view

[cord1 cord2 yM yS x1m x2m] = getPlotInfo(handles);

if cord1~=cord2
  figure
  surf(x1m,x2m,yS-yM);
end

% --- Executes on button press in Run2D.
function Run2D_Callback(hObject, eventdata, handles) %#ok

[fitNorm trans resultsIdx] = findTwoDFields(handles);

if ~isempty(resultsIdx)
  %if these results already exist, we don't need to recompute them
  return
end

% Determine the current setting of points per dimensions
switch get(handles.twoDPntsPerDimPU,'value')
  case 1
    N = 5;
  case 2
    N = 11;
  case 3
    N = 21;
  case 4
    N = 51;
  case 5
    N = 101;
  otherwise
    errordlg('Case not expected, check code')
    return
end

%pull some useful stuff off the handles structure
ymin = handles.ymin;
ymax = handles.ymax;
bnds = handles.range;
nom = handles.nom;
variesIdx = handles.variesIdx;

h = waitbar(0,'Performing single variable calcs');

%First do single variable calcs. Loop thru all n variables
for cord1 = variesIdx'

  if bnds(cord1,1) > 0
    logX = true;
  else
    logX = false;
  end
  
  %Evaluate y = M(x). x will be linspaced over its range or linspaced in
  %log cordinates
  x.lin = repmat(nom,1,N);
  x.lin(cord1,:) = linspace(bnds(cord1,1),bnds(cord1,2),N);
  y.t1 = eval(handles.MA,x.lin);
  if logX
    x.log = repmat(nom,1,N);
    x.log(cord1,:) = 10.^linspace(log10(bnds(cord1,1)),log10(bnds(cord1,2)),N);
    y.t2 = eval(handles.MA,x.log);
  end
  
  %If y is ever negative, we can't do logY transformations. Store this occurance in the handles structure  
  if min(y.t1) <= 0 || handles.logY == 0 || (isfield(y,'t2') && min(y.t2) <= 0)
    logY = false;
    handles.logY = 0;
  else
    logY = true;
  end
  
  if logY
    y.t3 = log10(y.t1);
    if logX
      y.t4 = log10(y.t2);
    end
  end
  %Finished evaluating y = M(x)
  
  
  t = {'t1';'t2';'t3';'t4'}; %considered transformations
  t = t([true logX logY logX&&logY]); %reduce to possible transformations
  norm = {'inf';'two'}; %considered norms
  
  %form 1-d regression matrix X in normalized variables for all possible
  %transformations
  normalizedx1.lin = 2*(linspace(bnds(cord1,1),bnds(cord1,2),N)-bnds(cord1,1))./diff(bnds(cord1,:),[],2)-1;
  normalizedX.t1 =  [ones(1,N); normalizedx1.lin; normalizedx1.lin.^2];
  normalizedX.t3 = normalizedX.t1;
  if logX
      normalizedx1.log = 2*(linspace(log10(bnds(cord1,1)),log10(bnds(cord1,2)),N)-log10(bnds(cord1,1)))./diff(log10(bnds(cord1,:)),[],2)-1;
      normalizedX.t2 = [ones(1,N); normalizedx1.log; normalizedx1.log.^2];
      normalizedX.t4 = normalizedX.t2;
  end

  for j1 = 1:length(t)
    currenty = y.(t{j1});
    storedY.(t{j1}){cord1,cord1} = currenty;
    ymin = min(ymin,min(currenty));
    ymax = max(ymax,max(currenty));
    for j2 = 1:length(norm)
      fitFctnHand = str2func([norm{j2} 'NormRegression']);
      coeff = fitFctnHand(normalizedX.(t{j1})',currenty','off');
      surrogateCoeff.(norm{j2}).(t{j1}){cord1,cord1} = coeff;
    end
  end
end

runIdx = 1;
nRuns = (length(variesIdx)-1)*length(variesIdx)/2;
for cord1 = variesIdx'
  for cord2 = variesIdx(variesIdx > cord1)'
  
    % Waitbar update
    waitbar(runIdx/nRuns,h,'Performing two variable calcs')
    runIdx = runIdx+1;
    
    if bnds(cord1,1) > 0 && bnds(cord2,1) > 0
      logX = true;
    else
      logX = false;
    end
    
    %Evaluate y = M(x). x will be linspaced over its range or linspaced in
    %log cordinates
    x.lin = repmat(nom,1,N^2);
    [x1m.lin x2m.lin] = meshgrid(linspace(bnds(cord1,1),bnds(cord1,2),N),linspace(bnds(cord2,1),bnds(cord2,2),N));
    x.lin(cord1,:) = x1m.lin(:)';
    x.lin(cord2,:) = x2m.lin(:)';
    y.t1 = eval(handles.MA,x.lin);

    normalizedx1.lin = 2*(x1m.lin(:)'-bnds(cord1,1))./diff(bnds(cord1,:),[],2)-1;
    normalizedx2.lin = 2*(x2m.lin(:)'-bnds(cord2,1))./diff(bnds(cord2,:),[],2)-1;
    normalizedX.t1 = [ones(1,N^2); normalizedx1.lin; normalizedx2.lin; normalizedx1.lin.^2; normalizedx1.lin.*normalizedx2.lin; normalizedx2.lin.^2];
    normalizedX.t3 = normalizedX.t1;
    if logX
      x.log = repmat(nom,1,N^2);
      [x1m.log x2m.log] = meshgrid(10.^linspace(log10(bnds(cord1,1)),log10(bnds(cord1,2)),N),10.^linspace(log10(bnds(cord2,1)),log10(bnds(cord2,2)),N));
      x.log(cord1,:) = x1m.log(:)';
      x.log(cord2,:) = x2m.log(:)';
      y.t2 = eval(handles.MA,x.log);

      normalizedx1.log = 2*(log10(x1m.log(:)')-log10(bnds(cord1,1)))./diff(log10(bnds(cord1,:)),[],2)-1;
      normalizedx2.log = 2*(log10(x2m.log(:)')-log10(bnds(cord2,1)))./diff(log10(bnds(cord2,:)),[],2)-1;
      normalizedX.t2 = [ones(1,N^2); normalizedx1.log; normalizedx2.log; normalizedx1.log.^2; normalizedx1.log.*normalizedx2.log; normalizedx2.log.^2];
      normalizedX.t4 = normalizedX.t2;
    end
    
    %If y is ever negative, we can't do logY transformations. Store this occurance in the handles structure  
    if handles.logY == 0 || min(y.t1) <= 0 || (isfield(y,'t2') && min(y.t2) <= 0)
      logY = false;
      handles.logY = 0;
    else
      logY = true;
    end
    
    if logY
      y.t3 = log10(y.t1);
      if logX
        y.t4 = log10(y.t2);
      end
    end
    %Finished evaluating y = M(x)
    
    t = {'t1';'t2';'t3';'t4'}; %considered transformations
    t = t([true logX logY logX&&logY]); %reduce to possible transformations
    
    for j1 = 1:length(t)
      currenty = y.(t{j1});
      ymin = min(ymin,min(currenty));
      ymax = max(ymax,max(currenty));
      storedY.(t{j1}){cord1,cord2} = reshape(currenty,N,N);
      for j2 = 1:length(norm)
        fitFctnHand = str2func([norm{j2} 'NormRegression']);
        coeff = fitFctnHand(normalizedX.(t{j1})',currenty','off');
        surrogateCoeff.(norm{j2}).(t{j1}){cord1,cord2} = coeff;
      end
    end
  end
end

delete(h) %remove waitbar

m = size(handles.twoD,1)+1;
handles.twoD(m,1).settings = int2str(N);
handles.twoD(m).surrogateCoeff = surrogateCoeff;
handles.twoD(m).storedY = storedY; %TODO:, why store both y and log10(y)

% Update handles structure
guidata(hObject, handles);

refreshTwoDDisplay(handles.x1Slide, eventdata, handles)

% --- Executes on button press in oneDactivePB.
function oneDactivePB_Callback(hObject, eventdata, handles) %#ok

n = size(handles.range,1);
ht = zeros(n,1);

[fitNorm trans resultsIdx] = findTwoDFields(handles);

for i1 = handles.variesIdx'
  ht(i1) = max(handles.twoD(resultsIdx).storedY.(trans){i1,i1}) - min(handles.twoD(resultsIdx).storedY.(trans){i1,i1});
end
figure
subplot(2,1,1)
bar(ht);
title('Swing of slice by number')
subplot(2,1,2)
bar(sort(ht,'descend'))
title('Swing of slice, sorted')

% --- Executes on button press in twoDactivePB.
function twoDactivePB_Callback(hObject, eventdata, handles) %#ok

n = size(handles.range,1);
ht = zeros(n);

[fitNorm trans resultsIdx] = findTwoDFields(handles);
variesIdx = handles.variesIdx;

for i1 = variesIdx'
  ht(i1,i1) = max(handles.twoD(resultsIdx).storedY.(trans){i1,i1}) - min(handles.twoD(resultsIdx).storedY.(trans){i1,i1});
end

for i1 = variesIdx'
  for i2 = variesIdx(variesIdx > i1)'
    ym = handles.twoD(resultsIdx).storedY.(trans){i1,i2};
    ht(i1,i2) = max(ym(:)) - min(ym(:));
  end
end

ht = ht + ht' - diag(diag(ht));
figure
subplot(2,1,1)
bar3(ht);
title('Swing of slice')
subplot(2,1,2)
bar(sort(max(ht),'descend'));
title('Column max')

function RunAll_Callback(hObject, eventdata, handles) %#ok

% Find what the current settings are
[fitNorm trans resultsIdx] = findGlobalFields(handles);

% If results exist for these settings and we are supposed to use existing
% results when possible, make a quick exit.
if ~isempty(resultsIdx) && get(handles.useCachedCB,'Value')
  return
end

% Grab current settings
Npnts4FitM = str2double(get(handles.Npnts4FitM,'string'));
Npnts4ValidM = str2double(get(handles.Npnts4ValidM,'string'));
activeIdx = handles.activeIdx;
bnds = handles.range;


% Form x vectors at which to evaluate the model. Notice we use the same
% design for the linearly and logrithmically spaced points.
if all(bnds(activeIdx,1) > 0)
  logX = true;
else
  logX = false;
end

n = length(activeIdx);
if n == 0
  Ncoeff = 1;
elseif n == 1
  Ncoeff = 3;
else
  Ncoeff = 1 + 2*n + nchoosek(n,2);
end

NFit = Ncoeff*Npnts4FitM;
normxFit = lhsdesign(NFit,n,'criterion','none');
xFit.lin = repmat(handles.nom,1,NFit);
xFit.lin(activeIdx,:) = (normxFit*diag(diff(bnds(activeIdx,:),[],2)) + repmat(bnds(activeIdx,1)',NFit,1))';
NVal = Ncoeff*Npnts4ValidM;
normxVal = lhsdesign(NVal,n,'criterion','none');
xVal.lin = repmat(handles.nom,1,NVal);
xVal.lin(activeIdx,:) = (normxVal*diag(diff(bnds(activeIdx,:),[],2)) + repmat(bnds(activeIdx,1)',NVal,1))';

if logX
  xFit.log = repmat(handles.nom,1,NFit);
  xFit.log(activeIdx,:) = (10.^(normxFit*diag(diff(log10(bnds(activeIdx,:)),[],2)) + repmat(log10(bnds(activeIdx,1))',NFit,1)))';
  xVal.log = repmat(handles.nom,1,NVal);
  xVal.log(activeIdx,:) = (10.^(normxVal*diag(diff(log10(bnds(activeIdx,:)),[],2)) + repmat(log10(bnds(activeIdx,1))',NVal,1)))';
end

% Evaluate the model y = M(x)

yFit.t1 = eval(handles.MA,xFit.lin);
yVal.t1 = eval(handles.MA,xVal.lin);

if logX
  yFit.t2 = eval(handles.MA,xFit.log);
  yVal.t2 = eval(handles.MA,xVal.log);
end

% Use optimization to try to determine the extrema of M(x).
% TODO: why not use some of the sample pnts as seeds?
NRE = str2double(get(handles.Nrestart,'string'));
errLowT = zeros(NRE+1,1);
errHighT = zeros(NRE+1,1);
for i3 = 1:NRE+1
  x = bnds(:,1) + rand(size(bnds,1),1).*diff(bnds,[],2);
  [errLowT(i3) trash errHighT(i3)] = optimError(handles.MA, bnds, zeros(n+1), 't1', activeIdx, x, x); %#ok
end
funLow = -errHighT;
funHigh = -errLowT;
handles.ymin = min([handles.ymin; funLow]);
handles.ymax = max([handles.ymax; funHigh]);

if handles.ymin <= 0
  handles.logY = false;
end

% Have we seen a negative y value?
if handles.logY == 0 || min([yFit.t1 yVal.t1]) <= 0 || (isfield(yFit,'t2') && min([yFit.t2 yVal.t2]) <= 0)
  logY = false;
  handles.logY = false;
else
  logY = true;
end

if logY
  yFit.t3 = log10(yFit.t1);
  yVal.t3 = log10(yVal.t1);
  if logX
    yFit.t4 = log10(yFit.t2);
    yVal.t4 = log10(yVal.t2);
  end
end

% Scale parameters to be betw +/-1 and create regression matrices
XFit = buildQuadX(2*normxFit-1);
XVal = buildQuadX(2*normxVal-1);

% Determine the legal transformations and specify the fitting norms
t = {'t1';'t2';'t3';'t4'};
t = t([true logX logY logX&&logY]);
norm = {'inf';'two'};

% Build approximations
switch get(handles.methodSelectPU,'value')
 case {1;2}
  for i1 = 1:length(t)
    for i2 = 1:length(norm)
      hand = str2func([norm{i2} 'NormRegression']);
      [C.(norm{i2}).(t{i1}) errFit.(norm{i2}).(t{i1})] = hand(XFit,yFit.(t{i1})','off');
      errVal.(norm{i2}).(t{i1}) = XVal*C.(norm{i2}).(t{i1}) - yVal.(t{i1})';

      % S is used for the optimization routines, but we don't save it
      S.(norm{i2}).(t{i1}) = coeff2quadform(C.(norm{i2}).(t{i1}),n);
    end
  end
% case 4
%  y = y - sum([ones(Nsamples,1) normx]*handles.inf.twoD.*[ones(Nsamples,1) normx],2)';
%  [coeff err] = twoNormRegression(XX,y');
%  poly = handles.inf.twoD + coeff2quadform(coeff,n);
 otherwise
  errordlg('code not complete')
end


% Use optimization to improve estimate of the peak fit error
NfromSamp = min([floor((NRE+1)/2) 50]); %Half of the restarts (up to a max of 50) will be from the worst performing sample points

for i1 = 1:length(t)
  for i2 = 1:length(norm)
    [trash idX] = sort([errFit.(norm{i2}).(t{i1}); errVal.(norm{i2}).(t{i1})],'descend'); %#ok
    [trash idx] = sort([errFit.(norm{i2}).(t{i1}); errVal.(norm{i2}).(t{i1})]); %#ok
    NN = min([50 length(idX)]); %We'll do up to 50 tries based on the sample points
    switch t{i1}
      case {'t1';'t3'}
        x = [xFit.lin xVal.lin];
      otherwise
        x = [xFit.log xVal.log];
    end
    xmin = x(:,idx(1:NN)); %these x's give small (neg.) y^S-y^M
    xmax = x(:,idX(1:NN)); %these x's give big y^S-y^M
    
    for i3 = 1:NRE+1
      x = bnds(:,1) + rand(size(bnds,1),1).*diff(bnds,[],2); %all parameters, not just active
      if i3 <= NfromSamp
        [errLowT(i3) trash errHighT(i3)] = optimError(handles.MA, bnds, S.(norm{i2}).(t{i1}), t{i1}, activeIdx, xmin(:,i3), xmax(:,i3)); %#ok
      else
        [errLowT(i3) trash errHighT(i3)] = optimError(handles.MA, bnds, S.(norm{i2}).(t{i1}), t{i1}, activeIdx, x, x); %#ok
      end
    end
    errLow.(norm{i2}).(t{i1}) = errLowT;
    errHigh.(norm{i2}).(t{i1}) = errHighT;

    %TODO: only save NRE+1 seeds, not all 50
    errXLow.(norm{i2}).(t{i1}) = xmin; %Note these are the seed points, not the optimizers 
    errXHigh.(norm{i2}).(t{i1}) = xmax; %Note these are the seed points, not the optimizers
  end
end

m = size(handles.global,1)+1;
activeIdxStr = num2str(activeIdx');
activeIdxStr = strrep(activeIdxStr,'  ',',');
activeIdxStr = strrep(activeIdxStr,',,',',');
handles.global(m,1).settings = ['Settings: ' int2str(Npnts4FitM) ',' int2str(Npnts4ValidM) ',' int2str(NRE) '  Active: ' activeIdxStr];
handles.global(m).activeIdx = activeIdx;
handles.global(m).surrogateCoeff = C;
handles.global(m).storedY_Train = yFit; %TODO: why save y if you don't save the x's?
handles.global(m).storedY_Validate = yVal;
handles.global(m).errOnTraining = errFit;
handles.global(m).errOnValidation = errVal;
handles.global(m).errFromOptimL = errLow;
handles.global(m).errFromOptimH = errHigh;
handles.global(m).errFromOptimLSeedPnt = errXLow;
handles.global(m).errFromOptimHSeedPnt = errXHigh;
handles.global(m).fctnMinFromOptim = funLow;
handles.global(m).fctnMaxFromOptim = funHigh;
    
str = cellstr(get(handles.availGlobalRunsPU,'String'));
str{m+1,1} = handles.global(m).settings;
set(handles.availGlobalRunsPU,'String',str,'Value',m+1);

% Update handles structure
guidata(hObject, handles);

refreshGlobalDisplay(hObject, eventdata, handles)

function [errLow xLow errHigh xHigh] = optimError(MA, bnds, poly, trans, activeIdx, xmin, xmax) 
% errLow = min(y^S - y^M)
% errHigh = max(y^S - y^M)

% We solve the optimization in normalized variables so that everything is
% nicely conditioned. Always normalize linearly, since some nonactive
% variables may have a negative lower bound.

options = optimset('fmincon');
options = optimset(options,'GradConstr','off','GradObj','off', ...
      'LargeScale','off','TolFun',5e-5,'Display','iter','TolCon',1e-6,'DiffMinChange',1e-3);

doesVaryBool = diff(bnds,[],2) ~= 0; 

%Eliminate constant components and normalize to keep fmincon sane
xmin = 2*(xmin(doesVaryBool) - bnds(doesVaryBool,1))./diff(bnds(doesVaryBool,:),[],2) - 1;
xmax = 2*(xmax(doesVaryBool) - bnds(doesVaryBool,1))./diff(bnds(doesVaryBool,:),[],2) - 1;

switch trans
 case 't1'
  hand1 = @(x) linXlinYMIN(x,MA,bnds,poly,activeIdx,doesVaryBool);
 case 't2'
  hand1 = @(x) logXlinYMIN(x,MA,bnds,poly,activeIdx,doesVaryBool);
 case 't3'
  hand1 = @(x) linXlogYMIN(x,MA,bnds,poly,activeIdx,doesVaryBool);
 case 't4'
  hand1 = @(x) logXlogYMIN(x,MA,bnds,poly,activeIdx,doesVaryBool);
 otherwise
  error('Internal inconsistency: condition should never occur')
end
hand2 = @(x) -hand1(x);

n = length(xmin);

[x,fval,exitflg]=fmincon(hand1,xmin,[],[],[],[],-ones(n,1),ones(n,1),[],options);

% Convert back to the nonnormalized variables
xLow = mean(bnds,2);
xLow(doesVaryBool) = 0.5*(x+1).*diff(bnds(doesVaryBool,:),[],2) + bnds(doesVaryBool,1);

if exitflg < 0
  disp('we had a bad exit flag')
  disp(exitflg)
  errLow = NaN;
  %keyboard
elseif exitflg == 0
  disp('Max iterations or function evaluations occured')
  errLow = fval;
else
  errLow = fval;
end

[x,fval,exitflg]=fmincon(hand2,xmax,[],[],[],[],-ones(n,1),ones(n,1),[],options);

% Convert back to the nonnormalized variables
xHigh = mean(bnds,2);
xHigh(doesVaryBool) = 0.5*(x+1).*diff(bnds(doesVaryBool,:),[],2) + bnds(doesVaryBool,1);

if exitflg < 0
  disp('we had a bad exit flag')
  disp(exitflg)
  errHigh = NaN;
  %keyboard
elseif exitflg == 0
  disp('Max iterations or function evaluations occured')
  errHigh = -fval;
else
  errHigh = -fval;
end

function err = linXlinYMIN(x,MA,bnds,poly,activeIdx,doesVary)
%doesVary is a logical array
realx = mean(bnds,2);
realx(doesVary) = 0.5*(x+1).*diff(bnds(doesVary,:),[],2) + bnds(doesVary,1);
normx = 2*(realx(activeIdx) - bnds(activeIdx,1))./diff(bnds(activeIdx,:),[],2) - 1;
err = [1 normx']*poly*[1;normx] - eval(MA,realx);

function err = logXlinYMIN(x,MA,bnds,poly,activeIdx,doesVary)
realx = mean(bnds,2);
realx(doesVary) = 0.5*(x+1).*diff(bnds(doesVary,:),[],2) + bnds(doesVary,1);
normx = 2*log10(realx(activeIdx)./bnds(activeIdx,1))./log10(bnds(activeIdx,2)./bnds(activeIdx,1)) - 1;
err = [1 normx']*poly*[1;normx] - eval(MA,realx);

function err = linXlogYMIN(x,MA,bnds,poly,activeIdx,doesVary)
realx = mean(bnds,2);
realx(doesVary) = 0.5*(x+1).*diff(bnds(doesVary,:),[],2) + bnds(doesVary,1);
normx = 2*(realx(activeIdx) - bnds(activeIdx,1))./diff(bnds(activeIdx,:),[],2) - 1;
err = [1 normx']*poly*[1;normx] - log10(eval(MA,realx));

function err = logXlogYMIN(x,MA,bnds,poly,activeIdx,doesVary)
realx = mean(bnds,2);
realx(doesVary) = 0.5*(x+1).*diff(bnds(doesVary,:),[],2) + bnds(doesVary,1);
normx = 2*log10(realx(activeIdx)./bnds(activeIdx,1))./log10(bnds(activeIdx,2)./bnds(activeIdx,1)) - 1;
err = [1 normx']*poly*[1;normx] - log10(eval(MA,realx));

function quadform = coeff2quadform(coeff,Nactive)

%force coeff to row vector
coeff = coeff(:)';

if Nactive >= 2
  Ncoeff = 1 + 2*Nactive + nchoosek(Nactive,2);
else
  Ncoeff = 1 + 2*Nactive;
end

if length(coeff) ~= Ncoeff
  error('Number of coefficients mismatches number of parameters')
end

matSize = Nactive+1;

quadform = zeros(matSize);
for i1 = 1:Nactive+1
  startidx = (i1-1)*matSize - sum(0:i1-1) + i1;
  endidx = i1*matSize - sum(0:i1-1);
  quadform(i1,i1:matSize) = coeff(startidx:endidx);
end
quadform = 0.5*(quadform + quadform');

% --- Executes on button press in SyncSlidersCB.
function SyncSlidersCB_Callback(hObject, eventdata, handles) %#ok
if get(hObject,'value')
  x2Slide_Callback(handles.x2Slide, eventdata, handles);
end

% --- Executes on button press in useCachedCB.
function useCachedCB_Callback(hObject, eventdata, handles) %#ok

% Callback does nothing. If this uicontrol is checked with we make a global run, we
% use existing results when possible.

function [cord1 cord2 yM yS x1m x2m] = getPlotInfo(handles)

%outputs
% x1m: matrix of x1 values from meshgrid
% x2m: matrix of x2 values from meshgrid

%find what the current settings are
[n t resultsIdx] = findTwoDFields(handles);

if isempty(resultsIdx)
  cord1 = [];
  cord2 = [];
  yM = [];
  yS = [];
  x1m = [];
  x2m = [];
else
  cord1 = handles.variesIdx(round(get(handles.x1Slide,'value')));
  cord2 = handles.variesIdx(round(get(handles.x2Slide,'value')));
  bnds = handles.range;
  N = str2double(handles.twoD(resultsIdx).settings);
  
  switch t  
   case {'t1';'t3'}
    x1 = linspace(bnds(cord1,1),bnds(cord1,2),N);
    x2 = linspace(bnds(cord2,1),bnds(cord2,2),N);
    [x1m x2m] = meshgrid(x1,x2);

    normx1 = 2*(x1-bnds(cord1,1))/diff(bnds(cord1,:),[],2)-1;
    normx1v = 2*(x1m(:)-bnds(cord1,1))/diff(bnds(cord1,:),[],2)-1;
    normx2v = 2*(x2m(:)-bnds(cord2,1))/diff(bnds(cord2,:),[],2)-1;

   case {'t2';'t4'}
    x1 = 10.^linspace(log10(bnds(cord1,1)),log10(bnds(cord1,2)),N);
    x2 = 10.^linspace(log10(bnds(cord2,1)),log10(bnds(cord2,2)),N);
    [x1m x2m] = meshgrid(x1,x2);

    normx1 = 2*(log10(x1)-log10(bnds(cord1,1)))/diff(log10(bnds(cord1,:)),[],2)-1;
    normx1v = 2*(log10(x1m(:))-log10(bnds(cord1,1)))/diff(log10(bnds(cord1,:)),[],2)-1;
    normx2v = 2*(log10(x2m(:))-log10(bnds(cord2,1)))/diff(log10(bnds(cord2,:)),[],2)-1;

   otherwise
      error('Condition should never occur')
  end
  
  if cord1==cord2
    yM = handles.twoD(resultsIdx).storedY.(t){cord1,cord1};
    coeff = handles.twoD(resultsIdx).surrogateCoeff.(n).(t){cord1,cord1};
    yS = coeff(1) + coeff(2)*normx1 + coeff(3)*normx1.^2;
    x1m = x1;
    x2m = [];
  elseif cord2 > cord1 %only the upper triangular part is saved.
    yM = handles.twoD(resultsIdx).storedY.(t){cord1,cord2};
    coeff = handles.twoD(resultsIdx).surrogateCoeff.(n).(t){cord1,cord2};
    yS = reshape(coeff(1) + [normx1v normx2v normx1v.^2 normx1v.*normx2v normx2v.^2]*coeff(2:end),N,N);
  else
    yM = handles.twoD(resultsIdx).storedY.(t){cord2,cord1}';
    coeff = handles.twoD(resultsIdx).surrogateCoeff.(n).(t){cord2,cord1};
    yS = reshape(coeff(1) + [normx2v normx1v normx2v.^2 normx1v.*normx2v normx1v.^2]*coeff(2:end),N,N);
  end
end

function refreshGlobalDisplay(hObject, eventdata, handles) %#ok

[fitNorm trans resultsIdx] = findGlobalFields(handles);

activeIdxStr = strrep(num2str(handles.activeIdx'),'  ',',');
activeIdxStr = strrep(activeIdxStr,',,',',');
set(handles.curActive_text,'String',activeIdxStr)

if isempty(resultsIdx)

  set(handles.RunAll,'BackgroundColor','red')
  set([handles.optimResultsPB],'enable','off')
  set([handles.avgTrainErr_text handles.avgValidErr_text handles.peakTrainErr_text handles.peakValidErr_text handles.peakOptimizedErr_text],'String','N.A.')
  set([handles.rangeOnSamples_text handles.rangeFromOptim_text handles.globalQ_text],'String','N.A.')
  set(handles.globalTransPU,'enable','off')

else  
  
  % If we got this far, results must be available. Use them.
  % The availGlobalRunsPU should take care of reseting Npnts4FitM, etc when
  % its value is changed.
  set(handles.RunAll,'BackgroundColor',[0.831373 0.815686 0.784314]);
  set([handles.optimResultsPB handles.globalTransPU],'enable','on')
  
  % Bunch of monkey business to determine the pop-up list of available
  % transformations and keep the previous setting when possible.
  activeIdx = handles.global(resultsIdx).activeIdx;
  if all(handles.range(activeIdx,1) > 0)
    logX = true;
  else
    logX = false;
  end
  logY = handles.logY;
  currAvailTrans = get(handles.globalTransPU,'String');
  t = {'linXlinY';'logXlinY';'linXlogY';'logXlogY'}; %considered transformations
  newAvailTrans = t([true logX logY logX&&logY]);
  
  if ~isequal(currAvailTrans,newAvailTrans)
    %we need to do something.
    currentTrans = currAvailTrans{get(handles.globalTransPU,'Value')};
    newVal = strmatch(currentTrans,newAvailTrans,'exact');
    if isempty(newVal)
      set(handles.globalTransPU,'Value',1) %linXlinY is alway an acceptible options
    else
      set(handles.globalTransPU,'Value',newVal)
    end
    set(handles.globalTransPU,'String',newAvailTrans)
  end
  % End monkey business
  
  err = handles.global(resultsIdx).errOnTraining.(fitNorm).(trans);
  set(handles.avgTrainErr_text,'string',num2str(mean(abs(err))));
  set(handles.peakTrainErr_text,'string',['[' num2str(min(err)) ', ' num2str(max(err)) ']']);

  err = handles.global(resultsIdx).errOnValidation.(fitNorm).(trans);
  set(handles.avgValidErr_text,'string',num2str(mean(abs(err))));
  set(handles.peakValidErr_text,'string',['[' num2str(min(err)) ', ' num2str(max(err)) ']']);

  low = min(handles.global(resultsIdx).errFromOptimL.(fitNorm).(trans));
  high = max(handles.global(resultsIdx).errFromOptimH.(fitNorm).(trans));
  set(handles.peakOptimizedErr_text,'string',['[' num2str(low) ', ' num2str(high) ']']);

  ally = [handles.global(resultsIdx).storedY_Train.(trans) handles.global(resultsIdx).storedY_Validate.(trans)];
  ymin = min(ally);
  ymax = max(ally);
  set(handles.rangeOnSamples_text,'string',['[' num2str(ymin) ', ' num2str(ymax) ']']);

  switch trans
    case {'t1';'t2'}
      range = max([handles.ymax ymax]) - min([handles.ymin ymin]);
      set(handles.rangeFromOptim_text,'string',['[' num2str(handles.ymin) ', ' num2str(handles.ymax) ']']);
    case {'t3';'t4'}
      range = max([log10(handles.ymax) ymax]) - min([log10(handles.ymin) ymin]);
      set(handles.rangeFromOptim_text,'string',['[' num2str(log10(handles.ymin)) ', ' num2str(log10(handles.ymax)) ']']);
  end

  set(handles.globalQ_text,'string',num2str(1-(high-low)/range));
  
end

% --- Executes on button press in twoDAvgErrSpecPB.
function twoDAvgErrSpecPB_Callback(hObject, eventdata, handles) %#ok
hand = @mean;
err = twoDErrSpec(handles,hand);
figure
bar3(err)

% --- Executes on button press in twoDPeakErrSpecPB.
function twoDPeakErrSpecPB_Callback(hObject, eventdata, handles) %#ok
hand = @max;
err = twoDErrSpec(handles,hand);
figure
bar3(err)

function err = twoDErrSpec(handles,hand)

[norm t m] = findTwoDFields(handles);
if isempty(m)
  return
end
bnds = handles.range;
n = size(bnds,1);
N = str2double(handles.twoD(m).settings);

err = zeros(n);
for i1 = 1:n
  for i2 = i1:n
    switch t  
     case {'t1';'t3'}
      x1 = linspace(bnds(i1,1),bnds(i1,2),N);
      x2 = linspace(bnds(i2,1),bnds(i2,2),N);
      [x1m x2m] = meshgrid(x1,x2);

      normx1 = 2*(x1-bnds(i1,1))/diff(bnds(i1,:),[],2)-1;
      normx1v = 2*(x1m(:)-bnds(i1,1))/diff(bnds(i1,:),[],2)-1;
      normx2v = 2*(x2m(:)-bnds(i2,1))/diff(bnds(i2,:),[],2)-1;

     case {'t2';'t4'}
      x1 = 10.^linspace(log10(bnds(i1,1)),log10(bnds(i1,2)),N);
      x2 = 10.^linspace(log10(bnds(i2,1)),log10(bnds(i2,2)),N);
      [x1m x2m] = meshgrid(x1,x2);

      normx1 = 2*(log10(x1)-log10(bnds(i1,1)))/diff(log10(bnds(i1,:)),[],2)-1;
      normx1v = 2*(log10(x1m(:))-log10(bnds(i1,1)))/diff(log10(bnds(i1,:)),[],2)-1;
      normx2v = 2*(log10(x2m(:))-log10(bnds(i2,1)))/diff(log10(bnds(i2,:)),[],2)-1;

     otherwise
        error('Condition should never occur')
    end

    if i1==i2
      yM = handles.twoD(m).ym.(norm).(t){i1,i1};
      coeff = handles.twoD(m).C.(norm).(t){i1,i1};
      yS = coeff(1) + coeff(2)*normx1 + coeff(3)*normx1.^2;
    else
      yM = handles.twoD(m).ym.(norm).(t){i1,i2};
      coeff = handles.twoD(m).C.(norm).(t){i1,i2};
      yS = reshape(coeff(1) + [normx1v normx2v normx1v.^2 normx1v.*normx2v normx2v.^2]*coeff(2:end),N,N);
    end
    err(i1,i2) = hand(hand(abs(yS-yM)));
  end
end
    
err = err+err'-diag(diag(err));

function [fitNorm trans resultsIdx] = findTwoDFields(handles)

% Determine which index of handles.twoD corresponds to the current setting
% for "pnts per dim"
switch get(handles.twoDPntsPerDimPU,'value')
  case 1
    N = 5;
  case 2
    N = 11;
  case 3
    N = 21;
  case 4
    N = 51;
  case 5
    N = 101;
  otherwise
    errordlg('Case not expected, check code')
    return
end
resultsIdx = [];
if length(handles.twoD) > 0
  resultsIdx = strmatch(int2str(N),{handles.twoD.settings});
end
%resultsIdx should be scalar since we don't allow multiple runs with the
%same settings to be performed, but just in case...
if ~isempty(resultsIdx)
  resultsIdx = resultsIdx(end);
end

%determine the currently selected surrogate fit norm
if get(handles.twoDNormSelectPU,'value') == 1
  fitNorm = 'inf';
else
  fitNorm = 'two';
end

% Determine the currently selected variable transformation
str = get(handles.twoDTransPU,'String');
currTrans = str{get(handles.twoDTransPU,'Value')};
switch currTrans
  case 'linXlinY'
    trans = 't1';
  case 'logXlinY'
    trans = 't2';
  case 'linXlogY'
    trans = 't3';
  case 'logXlogY';
    trans = 't4';
  otherwise
    error('unexpected value')
end

function [fitNorm trans resultsIdx] = findGlobalFields(handles)

% Determine which (if any) index of handles.global corresponds to the current setting
% for "pnts4Fit, pnts4Val, Nrestarts, and active

%find what the current settings are
Npnts4FitM = get(handles.Npnts4FitM,'string');
Npnts4ValidM = get(handles.Npnts4ValidM,'string');
NRE = get(handles.Nrestart,'string');
tmp = num2str(handles.activeIdx');
activeIdx = strrep(tmp,'  ',',');

%check for exact match
if length(handles.global) > 0
  resultsIdx = strmatch(['Settings: ' Npnts4FitM ',' Npnts4ValidM ',' NRE , '  Active: ' activeIdx],{handles.global.settings},'exact');

  if isempty(resultsIdx)
    %Check for partial match whose NRE setting is greater than or equal to
    %the current setting.
    possibleMatches = strmatch(['Settings: ' Npnts4FitM ',' Npnts4ValidM ','],{handles.global.settings});
    for i1 = 1:length(possibleMatches)
      %These match in NpntsPer. Now compare active parameters
      idx1 = strfind(handles.global(possibleMatches(i1)).settings,'Active: ');
      storedActive = handles.global(possibleMatches(i1)).settings(idx1+8:end);

      idx2 = length(['Settings: ' Npnts4FitM ',' Npnts4ValidM ',']) + 1;
      storedNRE = str2double(handles.global(possibleMatches(i1)).settings(idx2:(idx1-2)));
      if strcmp(storedActive,activeIdx) && NRE <= storedNRE
        resultsIdx = possibleMatches(i1);
      end
    end
  end
  %make sure its scalar
  if ~isempty(resultsIdx)
    resultsIdx = resultsIdx(end);
  end

else
  resultsIdx = [];
end

switch get(handles.methodSelectPU,'value')
  case 1
    fitNorm = 'inf';
  case 2
    fitNorm = 'two';
  otherwise
    error('code note complete')
end

% Determine the currently selected variable transformation
str = get(handles.globalTransPU,'String');
currTrans = str{get(handles.globalTransPU,'Value')};
switch currTrans
  case 'linXlinY'
    trans = 't1';
  case 'logXlinY'
    trans = 't2';
  case 'linXlogY'
    trans = 't3';
  case 'logXlogY';
    trans = 't4';
  otherwise
    error('unexpected value')
end

% --- Executes on selection change in Fig1DispOptPU.
function Fig1DispOptPU_Callback(hObject, eventdata, handles) %#ok

function fileName_Callback(hObject, eventdata, handles) %#ok


function savePB_Callback(hObject, eventdata, handles) %#ok

out.MA = handles.MA;
out.PA = handles.PA;
out.nom = handles.nom;
out.range = handles.range;
out.twoD = handles.twoD;
out.global = handles.global;
out.ymin = handles.ymin;
out.ymax = handles.ymax;
out.ynom = handles.ynom;
out.logY = handles.logY;
out.avgEvalTime = handles.avgEvalTime;
out.variesIdx = handles.variesIdx;
out.activeIdx = handles.activeIdx;

fName = get(handles.fileName,'string');
save(fName,'out')

%TODO:
%make optimized range use smart initial points
%make fig1 options functional
%add comp times, segregates betw optim and samp
%add buttons to view the optim+restarts results
%document handles structure
%estimate comp times
%list existing runs

function availGlobalRunsPU_Callback(hObject, eventdata, handles) %#ok
% This uicontrol should be disabled until results are available

val = get(handles.availGlobalRunsPU,'value');
if val ~= 1
  %we didn't select the current working setup. If val==1, this callback
  %leaves the gui unchanged.
  val = val-1;
  settings = handles.global(val).settings;

  % Change the text boxes as if the user had done it, and then refresh the
  % display.

  % Settings looks like "Settings: NptsPerFit,NptsPerVal,NRestarts  Active:
  % p1,p2,...,pn

  % toss out the first 10 characters, i.e., "Settings: "
  settings(1:10) = [];

  idx1 = strfind(settings,',');
  idx2 = strfind(settings,' ');
  set(handles.Npnts4FitM,'string',settings(1:idx1(1)-1));
  set(handles.Npnts4ValidM,'string',settings(idx1(1)+1:idx1(2)-1));
  set(handles.Nrestart,'string',settings(idx1(2)+1:idx2(1)-1));

  % toss out up to the comma separated active parameters
  settings(1:idx2(2)+8) = [];

  activeIdx = strrep(settings,',','  ');
  activeIdx = str2num(activeIdx)'; %#ok
  
  %save as current active
  handles.activeIdx = activeIdx;
  guidata(hObject,handles)
  refreshGlobalDisplay(hObject, eventdata, handles);
end



% --- Executes on button press in optimResultsPB.
function optimResultsPB_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to optimResultsPB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[norm t M] = findGlobalFields(handles);

y{1} = handles.global(M).fctnMaxFromOptim;
y{2} = handles.global(M).fctnMinFromOptim;
y{3} = handles.global(M).errFromOptimH.(norm).(t);
y{4} = handles.global(M).errFromOptimL.(norm).(t);
settings = handles.global(M).settings;
idx1 = strfind(settings,',');
idx2 = strfind(settings,'  ');
NRE = str2double(settings(idx1(2)+1:idx2(1)-1));
NfromSamp = min([floor((NRE+1)/2) 50]);

n = length(y{1});
h = zeros(4,1);
figure
for i1 = 1:4
  h(i1) = subplot(4,1,i1);
  xl = [0.5 0.5+n];
  if i1==1 || i1==3
    [yy idx] = sort(y{i1},'descend');
  else
    [yy idx] = sort(y{i1});
  end
  if n > 1
    yFromSampSeeds = yy(idx(2:NfromSamp));
    xFromSampSeeds = idx(2:NfromSamp);
  else
    yFromSampSeeds = [];
    xFromSampSeeds = [];
  end
  if n > NfromSamp
    yFromRandSeeds = yy(idx(NfromSamp+1:end));
    xFromRandSeeds = idx(NfromSamp+1:n);
  else
    yFromRandSeeds = [];
    xFromRandSeeds = [];
  end
  if i1==1 || i1==2
    plot(yy,'ko');
  else
    plot(idx(1),yy(idx(1)),'ro',xFromSampSeeds,yFromSampSeeds,'bo',xFromRandSeeds,yFromRandSeeds,'ko');
  end
  if n == 1
    yl = y{i1} + [-0.1 0.1]*y{i1}*sign(y{i1});
  else
    yl = [min(y{i1}) max(y{i1})]+ (max(y{i1})-min(y{i1}))*[-0.1 0.1];
  end
  set(h(i1),'YLim',yl,'XLim',xl)
end
set(h,'XTick',(1:n))

legend('Best Sample Seed','Sample Seed','Random Seed','Location','Best')

function SelActivePB_Callback(hObject, eventdata, handles) %#ok

n = size(handles.nom,1);
Units = DClab.ModelAndObservationPair(n);
for i1 = 1:n
  Units(i1).name = num2str(i1);
end
idx = DatasetUnitImport({},Units);

if isempty(idx)
  %user hit cancel
else
  handles.activeIdx = idx';
  guidata(hObject,handles)
end

% Refresh the display
refreshGlobalDisplay(hObject, eventdata, handles) %#ok


%h = uibuttongroup('Visible','off','Units','Characters','Position',[2 2 20 20]);
%u0 = uicontrol('Style','Radio','String','asdf','Units','Characters','pos',[1 3 10 1.45],'parent','h','HandleVisibility','off');
%set(h,'Visible','on')



