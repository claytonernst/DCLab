function out =  PDBuild(flag,varargin)
%PDBuild   ParameterDomain Builder
%
%  PDBuild opens the graphical user interface
%  to help construct ParameterDomain objects for the 
%  Data Collaboration tools produced by the Berkeley Center
%  for Control and Identification (BCCI) at UC Berkeley.
%  
%  PDBuild -ABOUT opens GUI's about window.
%
%  PDBuild -HELP opens GUI's HTML help in a web browser
%
%  PDBuild(ParameterAssObj,name) loads the gui with a given
%  parameterassertion object.

if nargin==1
  switch flag
    case {'about','-about','-ABOUT','-About'}
      menuAbout;
    case {'help','-help','-HELP','-Help'}
      menuHelp;
    otherwise
      error('Invalid input')
  end
  return
end

%create figure
figH = figure('position',[300 100 700 500], ...
              'menubar','none', ...
              'numbertitle','off', ...
              'name','ParameterDomain Builder - untitled.mat*', ...
              'renderer','openGL',...
              'resize','off',...
              'closeRequestFcn',@quitIt);

%new ParameterDomain
newParamDom([],[],1);

%menus
file = uimenu(figH,'label','File');
uimenu(file,'label','New ParameterDomain',...
            'Accelerator','n',...
            'callback',@newParamDom)
uimenu(file,'label','Load',...
            'Accelerator','l',...
            'callback',@loadIt)
sMenu = uimenu(file,'label','Save',...
            'Accelerator','s', ...
            'callback',@saveIt);
uimenu(file,'label','Save As',...
            'callback',{@saveIt,0})
export=uimenu(file,'label','Export');
uimenu(export,'label','To workspace',...
            'callback',@export2WS);
uimenu(export,'label','To m-file',...
            'callback',@export2mFile);
uimenu(file,'label','Quit',...
              'Accelerator','Q',...
              'callback',@quitIt);
Edit = uimenu(figH,'label','Edit');
uimenu(Edit,'label','Add Parameter',...
                'Accelerator','a', ...
                'callback',@newParameter)
uimenu(Edit,'label','Delete Parameter',...
                'Accelerator','d', ...
                'callback',@deleteParam)
rmenu=uimenu(Edit,'label','Revert', ...
                'Accelerator','z',...
                'callback',@revertIt,...
                'enable','off');
remenu=uimenu(Edit,'label','Redo', ...
                'Accelerator','r', ...
                'callback',@redoIt,...
                'enable','off');
Help = uimenu(figH,'label','Help');
uimenu(Help,'label','GUI Help',...
                'callback',@menuHelp);
uimenu(Help,'label','About',...
            'callback',@menuAbout);


%setappdata
setappdata(figH,'highlight',0);
setappdata(figH,'saved',0);
setappdata(figH,'revertMenu',rmenu)
setappdata(figH,'redoMenu',remenu)
setappdata(figH,'smenu',sMenu)
setappdata(figH,'nout',nargout)

%rpf add
if nargin==2
  P = flag; %get parameter obj
  name = varargin{1}; %get name
  newParamDom([],[],1);
  setappdata(gcf,'varName',name);
  varName = getappdata(gcf,'varNameHdl');
  set(varName,'string',name)
  setappdata(gcf,'paramObj',P);

  %setappdata(gcf,'saved',dirName);

  buildFromAppData;
  %set(gcf,'name',['ParameterDomain Builder - ' fileName])

  setappdata(gcf,'revertList',{P});
  setappdata(gcf,'redo',[]);

  rmenu=getappdata(gcf,'revertMenu');
  remenu=getappdata(gcf,'redoMenu');
  set(remenu,'enable','off');
  set(rmenu,'enable','off');
end
%end rpf add
  
if nargout
  uiwait(figH)
  out = getappdata(figH,'paramObj');
  delete(figH)
end



%-----------------------------------------------------
function newParamDom(obj,eventD,flag)

if nargin<3
    a = areYouSure;
    switch a
        case 'Yes'
            saveIt;
        case 'Cancel'
            return
    end
end
set(gcf,'name','ParameterDomain Builder - untitled.mat*');
            

for i1=get(gcf,'children')';
    if ~strcmp(get(i1,'type'),'uimenu')
        delete(i1)
    end
end

setappdata(gcf,'paramObj',[]);

%Variable Name
uicontrol(gcf, ...
          'style','text', ...
          'string','Variable Name: ', ...
          'position',[50 450 200 20], ...
          'HorizontalAlignment','left', ...
          'FontSize',16, ...
          'backgroundcolor',[0.8 0.8 0.8]);
varNameC=uicontrol(gcf, ...
          'style','edit', ...
          'string','P', ...
          'position',[220 450 400 20], ...
          'callback',@varName, ...
          'HorizontalAlignment','left', ...
          'FontSize',16, ...
          'backgroundcolor',[0.8 0.8 0.8]);
setappdata(gcf,'varName','P');
setappdata(gcf,'varNameHdl',varNameC);
      
%panel
ph = uipanel(gcf, ...
             'units','pixels', ...
             'position',[20 20 630 400]);

%scroll bar
slide=uicontrol(gcf, ...
                'style','slider', ...
                'position', [650 20 20 400], ...
                'sliderstep',[1 1], ...
                'max',1,...
                'min',0,...
                'value',1, ...
                'enable','off', ...
                'callback',@scrollIt);
     
%setappdata
setappdata(gcf,'panel',ph);
setappdata(gcf,'orderSize',zeros(0,2));
setappdata(gcf,'slide',slide);
setappdata(gcf,'lastvalue',1);
if nargin<3
    setappdata(gcf,'revertList',{});
    setappdata(gcf,'redo',[]);

    rmenu=getappdata(gcf,'revertMenu');
    remenu=getappdata(gcf,'redoMenu');
    set(remenu,'enable','off');
    set(rmenu,'enable','off');
end


%-----------------------------------------------------
function [ph1,button] = newParameter(obj,eventData,flag)

P = getappdata(gcf,'paramObj');
ph = getappdata(gcf,'panel');
indx = length(get(ph,'children'));
orderSize = getappdata(gcf,'orderSize');

if isempty(orderSize)
    pos = 400;
else
    pos = get(orderSize(end,2),'position');
    pos = pos(2);
end

if mod(indx+1,2)
    bgcolor=[0.85 0.85 0.85];
else
    bgcolor=[0.9 0.9 0.9];
end

ph1=uipanel('parent',ph, ...
            'units','pixels', ...
            'position',[0,pos-25,630,25], ...
            'backgroundcolor',bgcolor, ...
            'borderType','none', ...
            'borderwidth',0,...
            'buttonDownFcn',@highlight);
        
%button
button = uicontrol('parent',ph1, ...
          'style','pushbutton', ...
          'position',[5,2,21,21], ...
          'string','+',...
          'callback',@expandCollapse,...
          'buttonDownFcn',@highlight);
      
%paramName
uicontrol('parent',ph1, ...
          'style','text', ...
          'string','Name:', ...
          'HorizontalAlignment','left', ...
          'position',[52,2,45,21], ...
          'backgroundColor',bgcolor,...
          'buttonDownFcn',@highlight,...
          'enable','inactive');
name = uicontrol('parent',ph1, ...
          'style','edit', ...
          'string','', ...
          'HorizontalAlignment','left', ...
          'position',[102,2,140,21], ...
          'backgroundColor',bgcolor,...
          'buttonDownFcn',@highlight, ...
          'callback',{@updateField,'parameter_name'});
      
%nominal
uicontrol('parent',ph1, ...
          'style','text', ...
          'string','Nominal:', ...
          'HorizontalAlignment','left', ...
          'position',[258,2,60,21], ...
          'backgroundColor',bgcolor,...
          'buttonDownFcn',@highlight,...
          'enable','inactive');
nom = uicontrol('parent',ph1, ...
          'style','edit', ...
          'string','', ...
          'HorizontalAlignment','left', ...
          'position',[323,2,60,21], ...
          'backgroundColor',bgcolor,...
          'buttonDownFcn',@highlight, ...
          'callback',{@updateField,'nominal'});
       
%range
uicontrol('parent',ph1, ...
          'style','text', ...
          'string','Range:', ...
          'HorizontalAlignment','left', ...
          'position',[415,2,45,21], ...
          'backgroundColor',bgcolor,...
          'buttonDownFcn',@highlight,...
          'enable','inactive');
range = uicontrol('parent',ph1, ...
          'style','edit', ...
          'string','', ...
          'HorizontalAlignment','left', ...
          'position',[465,2,160,21], ...
          'backgroundColor',bgcolor,...
          'buttonDownFcn',@highlight, ...
          'callback',{@updateField,'range'});
       
%description
uicontrol('parent',ph1, ...
          'style','text', ...
          'string','Description:', ...
          'HorizontalAlignment','left', ...
          'position',[52,-25,80,21], ...
          'backgroundColor',bgcolor,...
          'visible','off',...
          'buttonDownFcn',@highlight,...
          'enable','inactive');
desc = uicontrol('parent',ph1, ...
          'style','edit', ...
          'max',4, ...
          'min',1, ...
          'string','', ...
          'HorizontalAlignment','left', ...
          'position',[137,-49,160,45], ...
          'backgroundColor',bgcolor,...
          'visible','off',...
          'buttonDownFcn',@highlight, ...
          'callback',{@updateField,'description'});
       
%uncertainty value
uicontrol('parent',ph1, ...
          'style','text', ...
          'string','Uncertainty:', ...
          'HorizontalAlignment','left', ...
          'position',[307,-49,80,45], ...
          'backgroundColor',bgcolor,...
          'visible','off',...
          'buttonDownFcn',@highlight,...
          'enable','inactive');
uvalue = uicontrol('parent',ph1, ...
          'style','edit', ...
          'string','', ...
          'HorizontalAlignment','left', ...
          'position',[392,-25,60,21], ...
          'backgroundColor',bgcolor,...
          'visible','off',...
          'buttonDownFcn',@highlight, ...
          'callback',{@updateField,'uncertainty_percent'}, ...
          'enable','off');
       
%uncertainty type
uicontrol('parent',ph1, ...
          'style','text', ...
          'string','Uncertainty Type: ', ...
          'HorizontalAlignment','right', ...
          'position',[475,-49,80,45], ...
          'backgroundColor',bgcolor,...
          'visible','off',...
          'buttonDownFcn',@highlight,...
          'enable','inactive');
umethod = uicontrol('parent',ph1, ...
          'style','popupmenu', ...
          'string',['bnds ';'log10';'  %  '], ...
          'HorizontalAlignment','left', ...
          'position',[560,-25,65,21], ...
          'backgroundColor',bgcolor,...
          'visible','off',...
          'buttonDownFcn',@highlight, ...
          'callback',{@updateField,'uncertainty_method'});

%setappdata
setappdata(ph1,'name',name);
setappdata(ph1,'nom',nom);
setappdata(ph1,'range',range);
setappdata(ph1,'desc',desc);
setappdata(ph1,'uvalue',uvalue);
setappdata(ph1,'umethod',umethod);

   
%add ParameterDomain
if nargin<3
    len = length(P);
    if (isempty(P) && ~isa(P,'FreeParameter')) || len==0
        P = DClab.FreeParameter;
    else
        P(len+1,1) = DClab.FreeParameter;
    end
    setappdata(gcf,'paramObj',P);
end
if nargin<3
    revertList(P);
end

%add to orderSize
orderSize = getappdata(gcf,'orderSize');
orderSize(end+1,:) = [1 ph1];
setappdata(gcf,'orderSize',orderSize);

%highlight
highlight(ph1,[],0);
setappdata(gcf,'lastvalue',1);
adjustScroll([],[],1);

name = get(gcf,'name');
if name(end)~='*'
    set(gcf,'name',[name '*']);
    saveMenu = getappdata(gcf,'smenu');
    set(saveMenu,'enable','on')
end

expandCollapse(button,[]);

%-----------------------------------------------------
function highlight(obj,event,flag) %#ok

if nargin<3
    flag=0;
end

orderSize = getappdata(gcf,'orderSize');
if nargin==3
    cbo = obj;
else
    cbo = gcbo;
end
if ~flag
    for i1=orderSize(:,2)'
        highlight(i1,[],1);
    end
end

type = get(cbo,'type');
switch type
    case 'uipanel'
        panel = cbo;
    otherwise
        panel = get(cbo,'parent');
end
toChange = [];
toChange2 = [];
for i1=get(panel,'children')'
    if ~strcmp(get(i1,'style'),'pushbutton') && ~strcmp(get(i1,'enable'),'off')
        toChange = [toChange; i1];
    elseif strcmp(get(i1,'enable'),'off')
        toChange2 = i1;
    end
end

indx = find(orderSize(:,2)==panel);

if get(panel,'backgroundcolor')==[0.7 0.7 0.9] | (nargin==3 & flag)
    if mod(indx,2)
        bgcolor=[0.85 0.85 0.85];
    else
        bgcolor=[0.9 0.9 0.9];
    end
    bgcolor2 = [0.7 0.7 0.7];
else
    bgcolor=[0.7 0.7 0.9];
    bgcolor2=[0.5 0.5 0.7];
    setappdata(gcf,'highlight',indx);
end

set(panel,'backgroundcolor',bgcolor);
set(toChange,'backgroundcolor',bgcolor);
set(toChange2,'backgroundcolor',bgcolor2);

%------------------------------------------------------
function expandCollapse(obj,eventD)

highlight(obj,[],0);

orderSize = getappdata(gcf,'orderSize');
indx = getappdata(gcf,'highlight');
slide = getappdata(gcf,'slide');

switch orderSize(indx,1)
    case 1
        set(obj,'string','-');
        for i1=orderSize(indx:end,2)'
            pos = get(i1,'position');
            pos = pos + [0 -50 0 0];
            if orderSize(indx,2)==i1
                pos(4)=75;
                for i2 = get(i1,'children')'
                    set(i2,'position',get(i2,'position')+[0 50 0 0]);
                    set(i2,'visible','on');
                end
            end
            set(i1,'position',pos);
        end
        orderSize(indx,1)=3;
    case 3
        set(obj,'string','+');
        for i1=orderSize(indx:end,2)'
            pos = get(i1,'position');
            pos = pos +[0 50 0 0];
            if orderSize(indx,2)==i1
                pos(4)=25;
                for i2 = get(i1,'children')'
                    pos2 = get(i2,'position')+[0 -50 0 0];
                    set(i2,'position',pos2);
                    if pos2(2)<0
                        set(i2,'visible','off');
                    end
                end
            end
            set(i1,'position',pos);
        end
        orderSize(indx,1)=1;
end

len = sum(orderSize(:,1));
if len>16
    set(slide,'max',len-16);
    count=0;
    for i1=orderSize(indx:end,2)'
        pos = get(i1,'position');
        count = count+(pos(2)<0);
    end
    set(slide,'value',count);
end
setappdata(gcf,'lastvalue',get(slide,'value'));

setappdata(gcf,'orderSize',orderSize);
adjustScroll([],[],1);

%---------------------------------------------------------
function varName(obj,eventD)

value = get(obj,'string');
if ~isempty(regexp(value,'\W')) | ~isempty(regexp(value(1),'[0-9]'))
    set(obj,'foregroundcolor',[1 0 0])
    errordlg('Invalid variable name.','Error');
else
    set(obj,'foregroundcolor',[0 0 0]);
    setappdata(gcf,'varName',value);
end

%---------------------------------------------------------
function updateField(obj,eventD,property,flag)

if nargin<4
  flag = 1;
end

P = getappdata(gcf,'paramObj');
orderSize = getappdata(gcf,'orderSize');
indx = find(orderSize(:,2)==get(gcbo,'parent'));

switch property
    case 'parameter_name'
        value = get(gcbo,'string');
        
        indx2 = find(strcmp(P.name,value));
        if ~isempty(indx2) && any(~strcmp(P(indx2).name,''))
            errordlg('Parameter name must be unique','Error');
            set(obj,'foregroundcolor',[1 0 0]);
            return
        end

        P(indx).name = value;
        set(obj,'foregroundcolor',[0 0 0]);

    case 'nominal_value'
        tmp = get(gcbo,'string');
        if isempty(tmp)
            value = [];
        else
            try
                value = eval(tmp);
            catch
                errordlg('Invalid nominal value','Error');
                set(obj,'foregroundcolor',[1 0 0]);
                return
            end
        end
        if ~isnumeric(value) | (~isscalar(value) & ~isempty(value))
            errordlg('Nominal value must be scalar numeric','Error');
            set(obj,'foregroundcolor',[1 0 0]);
            return
        end
        tmp = getappdata(get(gcbo,'parent'),'umethod');
        if get(tmp,'value')==2 & value<=0
            errordlg('Nominal must be positive for log10 uncertainty method','Error');
            set(obj,'foregroundcolor',[1 0 0]);
            return
        end
        minn = P(indx).lower_bound;
        maxx = P(indx).upper_bound;
        if value < minn
            minn = value;
        elseif value > maxx
            maxx = value;
        end
        umethod = P(indx).uncertainty_method;
        P(indx).nominal = value;
        P(indx).lower_bound = minn;
        P(indx).upper_bound = maxx;
        P(indx).uncertainty_method = umethod;
        tmp = getappdata(get(gcbo,'parent'),'uvalue');
        setappdata(gcf,'paramObj',P);
        if ~isempty(get(tmp,'string'))
            updateField(tmp,[],'uncertainty_percent',0);
        end
        set(gcbo,'foregroundcolor',[0 0 0]);

    case 'range'
        tmp = get(gcbo,'string');
        if isempty(tmp)
            value = zeros(0,2);
        else
            try
                value = eval(tmp);
            catch
                errordlg('Invalid range value','Error');
                set(obj,'foregroundcolor',[1 0 0]);
                return
            end
        end
        if ~isnumeric(value) | size(value,2)~=2
            errordlg('Range value must be 1x2 numeric','Error');
            set(obj,'foregroundcolor',[1 0 0]);
            return
        end
        nom = P(indx).nominal_value;
        if isempty(value)
            val1=[]; val2=[];
        else
            val1=value(1);val2=value(2);
        end
        if val1>nom
            nom = val1;
        elseif nom>val2
            nom = val2;
        end
        P(indx).nominal_value = nom;
        P(indx).lower_bound = val1;
        P(indx).upper_bound = val2;

        tmp = getappdata(get(gcbo,'parent'),'umethod');
        set(tmp,'value',1);
        tmp = getappdata(get(gcbo,'parent'),'uvalue');
        set(tmp,'enable','off');
        set(gcbo,'foregroundcolor',[0 0 0]);

    case 'description'
        value = get(gcbo,'string');
        P(indx).description = value;

    case 'uncertainty_percent'
        tmp = get(obj,'string');
        if isempty(tmp)
            value = [];
        else
            try
                value = eval(tmp);
            catch
                errordlg('Invalid uncertainty value','Error');
                set(obj,'foregroundcolor',[1 0 0]);
                return
            end
        end
        if ~isnumeric(value) | (~isscalar(value) & ~isempty(value))
            errordlg('Uncertainty Value must be scalar numeric','Error');
            set(obj,'foregroundcolor',[1 0 0]);
            return
        end
        P(indx).uncertainty_percent = value;
        
        if isempty(value)
            return
        end

        nom = P(indx).nominal_value;
        switch P(indx).uncertainty_method
            case 'bnds'
                return
            case '%'
                minn = (1-value/100)*nom;
                maxx = (1+value/100)*nom;
            case 'log10'
                minn = 10^(log10(nom)-value/2);
                maxx = 10^(log10(nom)+value/2);
        end

        tmp = getappdata(get(gcbo,'parent'),'range');
        if ~isempty(minn) & ~isempty(maxx)
            set(tmp,'string',sprintf('[%g, %g]',minn,maxx));
        end
        if flag
          P(indx).lower_bound = minn;
          P(indx).upper_bound = maxx;
        end
        set(obj,'foregroundcolor',[0 0 0]);

    case 'uncertainty_method'
        value = get(gcbo,'value');
        tmp = getappdata(get(gcbo,'parent'),'uvalue');
        switch value
            case 1
                set(tmp,'enable','off')
                P(indx).uncertainty_method = 'bnds';
            case 2
                set(tmp,'enable','on')
                P(indx).uncertainty_method = 'log10';
                setappdata(gcf,'paramObj',P);
                if ~isempty(get(tmp,'string'))
                    updateField(tmp,[],'uncertainty_percent',0);
                end
            case 3
                set(tmp,'enable','on')
                P(indx).uncertainty_method = '%';
                setappdata(gcf,'paramObj',P);
                if ~isempty(get(tmp,'string'))
                    updateField(tmp,[],'uncertainty_percent',0);
                end
        end

end

setappdata(gcf,'paramObj',P);
revertList(P);
highlight;

name = get(gcf,'name');
if name(end)~='*'
    set(gcf,'name',[name '*']);
    saveMenu = getappdata(gcf,'smenu');
    set(saveMenu,'enable','on')
end

%---------------------------------------------------------
function deleteParam(obj,eventD)

indx = getappdata(gcf,'highlight');
if indx==0
    return
end
orderSize = getappdata(gcf,'orderSize');
paramObj = getappdata(gcf,'paramObj');

toMove = get(orderSize(indx,2),'position');
toMove = toMove(4);
delete(orderSize(indx,2));

if indx<size(orderSize,1)
    for i=orderSize(indx+1:end,2)'
        pos = get(i,'position');
        pos = pos + [0 toMove 0 0];
        set(i,'position',pos);
        if pos(2)>0
            set(i,'visible','on');
        end
    end
end

orderSize(indx,:) = [];
paramObj(indx) = [];

setappdata(gcf,'orderSize',orderSize);
setappdata(gcf,'paramObj',paramObj);
revertList(paramObj);
if indx>size(orderSize,1)
    indx=indx-1;
end
highlight(orderSize(indx,2),[],0);
adjustScroll;

name = get(gcf,'name');
if name(end)~='*'
    set(gcf,'name',[name '*']);
    saveMenu = getappdata(gcf,'smenu');
    set(saveMenu,'enable','on')
end

%-------------------------------------------------------
function adjustScroll(obj,eventD,snap);

if nargin<3
    snap=0;
end

slide = getappdata(gcf,'slide');
orderSize = getappdata(gcf,'orderSize');
indx = getappdata(gcf,'highlight');
len = sum(orderSize(:,1));

if len<=16
    pos = get(orderSize(1,2),'position');
    if pos(2)+pos(4)<=400
        set(slide,'enable','off');
        set(slide,'value',1);
        setappdata(gcf,'lastvalue',1);
    else
        set(slide,'enable','on');
        set(slide,'value',0);
        maxx = (pos(2)+pos(4)-400)/25;
        set(slide,'max',maxx);
        set(slide,'sliderstep',[1/maxx 1/maxx]);
        setappdata(gcf,'lastvalue',0);
    end
else
    pos = get(orderSize(end,2),'position');
    if pos(2)<=0
        pos(2)=0;
    end
    maxx=len-16+pos(2)/25;
    set(slide,'max',maxx);
    set(slide,'sliderstep',[1/maxx 1/maxx]);
    if strcmp(get(slide,'enable'),'off')
        setappdata(gcf,'lastvalue',maxx);
    end
    set(slide,'enable','on');
    
    pos = get(orderSize(indx,2),'position');
    if (pos(2)>=0 && pos(2)+pos(4)<=400) | ~snap
        len=size(orderSize,1);
        moveTo=0;
        for i1=orderSize(:,2)'
            pos = get(i1,'position');
            moveTo = moveTo+(pos(2)<0);
        end
        setappdata(gcf,'lastvalue',moveTo);
    elseif sum(orderSize(1:indx,1))>16
        len = size(orderSize,1);
        if indx==len
            moveTo = 0;
        else
            moveTo = sum(orderSize(indx+1:end,1));
        end
    end
    set(slide,'value',moveTo);
end
scrollIt(slide);

%-------------------------------------------------------
function scrollIt(obj,eventD)

value = get(obj,'value');
orderSize = getappdata(gcf,'orderSize');
lastvalue = getappdata(gcf,'lastvalue');

value = round(value);
set(obj,'value',value);
setappdata(gcf,'lastvalue',value);

toMove = 25*(lastvalue-value);

for i1=orderSize(:,2)'
    pos = get(i1,'position');
    pos = pos + [0 toMove 0 0];
    set(i1,'position',pos);
    if pos(2)<0 | pos(2)+pos(4)>400
        set(i1,'visible','off');
    else
        set(i1,'visible','on');
    end
end


if value==get(obj,'max') & sum(orderSize(:,1))<16
    set(obj,'enable','off');
    set(obj,'value',1);
end

len = sum(orderSize(:,1));
if len>16
    pos = get(orderSize(end,2),'position');
    if pos(2)<0
        pos(2)=0;
    else
        set(obj,'value',0);
        setappdata(gcf,'lastvalue',0);
    end
    maxx=len-16+pos(2)/25;
    set(obj,'max',maxx);
    set(obj,'sliderstep',[1/maxx 1/maxx]);
else
    pos = get(orderSize(end,2),'position');
    if pos(2)>=0
        set(obj,'value',0);
        setappdata(gcf,'lastvalue',0);
    end
    maxx=pos(2)/25-16+len;
    if maxx==0
        set(obj,'max',1);
        set(obj,'value',1);
        set(obj,'enable','off');
    else
        set(obj,'max',maxx);
        set(obj,'sliderstep',[1/maxx 1/maxx]);
    end
end
    

%-------------------------------------------------------
function saveIt(obj,eventD,flag)

if nargin<3
    flag = getappdata(gcf,'saved');
end
P = getappdata(gcf,'paramObj');
varName = getappdata(gcf,'varName');

name = get(gcf,'name');
if (name(end)~='*' & flag) | isempty(P)
    return
end
P = prepareObj(P);
if ~P
    return
end
if name(end)=='*'
    name=name(1:end-1);
    saveMenu = getappdata(gcf,'smenu');
    set(saveMenu,'enable','off')
end

eval([varName '=P;']);
if ~flag
    [dirName,fileName]=saveItAs(name(27:end-4));
    if dirName==0 && fileName==0
        return
    end
    setappdata(gcf,'saved',dirName);
else
    dirName=flag;
    fileName=name(27:end-4);
end

save(fullfile(dirName,fileName),varName);
set(gcf,'name',['ParameterDomain Builder - ' fileName ]);

%-------------------------------------------------------
function [dirName,fileName]=saveItAs(name);

pos = get(gcf,'position');

[fileName,dirName]=uiputfile('*.mat','Save',[name '.mat'],'location',pos(1:2)+[50 50]);

if all(fileName==0) || all(dirName==0)
    return
end

while ~exist(dirName) || regexp(fileName,'[^a-z_A-Z0-9.]')
    
    if regexp(fileName,'[^a-z_A-Z0-9.]')
        errordlg('Invalid file name','Error')
    end
    if ~exist(dirName)
        errordlg('Not a valid directory','Error')
    end
    
    [fileName,dirName]=uiputfile('*.mat','Save',[name '.mat'],'location',pos(1:2)+[50 50]);
    
    if fileName==0 | dirName==0
         return
    end
end

%-------------------------------------------------------
function a=areYouSure

name = get(gcf,'name');
if name(end)=='*' && ~isempty(getappdata(gcf,'paramObj'))
    a=questdlg('Would you like to save first?','Save?');
else
    a='No';
end

%-------------------------------------------------------
function loadIt(obj,eventD)

a = areYouSure;
switch a
    case 'Yes'
        saveIt;
    case 'Cancel'
        return
end

pos = get(gcf,'position');

[fileName,dirName]=uigetfile('.mat','Load','location',pos(1:2)+[50 50]);
if all(fileName==0) || all(dirName==0)
    return
end

while ~exist(fullfile(dirName,fileName))
    errordlg('Not a valid filename','Error')
    [fileName,dirName]=uigetfile('.mat','Load','location',pos(1:2)+[50 50]);
    if all(fileName==0) || all(dirName==0)
        return
    end
end
load(fullfile(dirName,fileName))

tmp=who;
for i1=1:length(tmp)
    if eval(['isa(' tmp{i1} ',''FreeParameter'')'])
        break
    end
    if i1==length(tmp)
        errordlg('No FreeParameter in mat-file','Error')
        return
    end
end

eval(['P=' tmp{i1} ';']);

newParamDom([],[],1);
setappdata(gcf,'varName',tmp{i1});
varName = getappdata(gcf,'varNameHdl');
set(varName,'string',tmp{i1})

setappdata(gcf,'paramObj',P);

setappdata(gcf,'saved',dirName);

buildFromAppData;
set(gcf,'name',['ParameterDomain Builder - ' fileName])

setappdata(gcf,'revertList',{P});
setappdata(gcf,'redo',[]);

rmenu=getappdata(gcf,'revertMenu');
remenu=getappdata(gcf,'redoMenu');
set(remenu,'enable','off');
set(rmenu,'enable','off');


%------------------------------------------------
function buildFromAppData

P = getappdata(gcf,'paramObj');

for i1=1:length(P)
    [ph1(i1),button] = newParameter([],[],1);
    expandCollapse(button,[]);
    name = getappdata(ph1(i1),'name');
    desc = getappdata(ph1(i1),'desc');
    uvalue = getappdata(ph1(i1),'uvalue');
    umethod = getappdata(ph1(i1),'umethod');
    nom = getappdata(ph1(i1),'nom');
    range = getappdata(ph1(i1),'range');
    set(name,'string',P(i1).name);
    set(desc,'string',P(i1).description);
    set(uvalue,'string',num2str(P(i1).uncertainty_percent));
    set(nom,'string',num2str(P(i1).nominal_value));
    if ~isempty(P(i1).lower_bound) && ~isempty(P(i1).upper_bound)
        set(range,'string',sprintf('[%g,%g]',P(i1).lower_bound,P(i1).upper_bound));
    end

    switch P(i1).uncertainty_method
        case 'bnds'
            value=1;
        case 'log10'
            value=2;
            set(uvalue,'enable','on')
        case '%'
            value=3;
            set(uvalue,'enable','on')
    end
    set(umethod,'value',value);
end
highlight(ph1(1),[],0);

%----------------------------------------------------
function revertList(P)

PList = getappdata(gcf,'revertList');
rmenu = getappdata(gcf,'revertMenu');
remenu = getappdata(gcf,'redoMenu');

set(rmenu,'enable','on');
set(remenu,'enable','off');

if length(PList)==40
    PList = [{P} PList(1:end-1)];
else
    PList = [{P} PList];
end

if length(PList)==1
    set(rmenu,'enable','off');
end

setappdata(gcf,'revertList',PList);
setappdata(gcf,'redo',[]);

%-----------------------------------------------------
function revertIt(obj,eventD)

PList = getappdata(gcf,'revertList');
newParamDom([],[],1);
setappdata(gcf,'paramObj',PList{2});
buildFromAppData;
name = get(gcf,'name');
if name(end)~='*'
    name(end+1)='*';
    saveMenu = getappdata(gcf,'smenu');
    set(saveMenu,'enable','on')
end
set(gcf,'name',name)

rmenu=getappdata(gcf,'revertMenu');
remenu=getappdata(gcf,'redoMenu');

redo = PList{1};
setappdata(gcf,'redo',redo);
set(remenu,'enable','on');
PList = PList(2:end);
if length(PList)==1
    set(rmenu,'enable','off');
end
setappdata(gcf,'revertList',PList);

%------------------------------------------------------
function redoIt(obj,eventD)

redo = getappdata(gcf,'redo');
newParamDom([],[],1);
setappdata(gcf,'paramObj',redo);
buildFromAppData;
name = get(gcf,'name');
if name(end)~='*'
    name(end+1)='*';
    saveMenu = getappdata(gcf,'smenu');
    set(saveMenu,'enable','on')
end
set(gcf,'name',name)
revertList(redo);

%-------------------------------------------------------
function obj=prepareObj(obj)

names = obj.name;

if any(strcmp(obj.name,''))
    errordlg('All parameters must have a name first.','Error');
    obj=0;
    return
end
[trash,idx]=sort(names);
obj = obj(idx);

%-------------------------------------------------------
function export2WS(obj,eventD)

varName = getappdata(gcf,'varName');
obj = prepareObj(getappdata(gcf,'paramObj'));
if ~obj
    return
end
assignin('base',varName,obj);

%-------------------------------------------------------
function export2mFile(obj,eventD)

varName = getappdata(gcf,'varName');
P = getappdata(gcf,'paramObj');
P = prepareObj(P);
if ~P
    return
end

name = get(gcf,'name');
if name(end)=='*'
    name=name(1:end-1);
    saveMenu = getappdata(gcf,'smenu');
    set(saveMenu,'enable','off')
end
[dirName,fileName]=saveItAs(name(27:end-4));

fid = fopen(fullfile(dirName,[fileName '.m']),'w');

fprintf(fid,'%%%s\n',fullfile(dirName,[fileName '.m']));

fprintf(fid,'\n');


for i1=1:length(P)
    fprintf(fid,'%s%d = ParameterDomain;\n',varName,i1);
    name = P(i1).parameter_name;
    nom = P(i1).nominal_value;
    lb = P(i1).lower_bound;
    ub = P(i1).upper_bound;
    desc = P(i1).description;
    uvalue = P(i1).uncertainty_percent;
    umethod = P(i1).uncertainty_method;
    
    if isempty(name)
        fprintf(fid,'%s%d.parameter_name = '''';\n',varName,i1);
    else
        fprintf(fid,'%s%d.parameter_name = ''%s'';\n',varName,i1,name);
    end
    if isempty(nom)
        fprintf(fid,'%s%d.nominal_value = [];\n',varName,i1);
    else
        fprintf(fid,'%s%d.nominal_value = %0.5f;\n',varName,i1,nom);
    end
    if isempty(lb)
        fprintf(fid,'%s%d.lower_bound = [];\n',varName,i1);
    else
        fprintf(fid,'%s%d.lower_bound = %0.5f;\n',varName,i1,lb);
    end
    if isempty(ub)
        fprintf(fid,'%s%d.upper_bound = [];\n',varName,i1);
    else
        fprintf(fid,'%s%d.upper_bound = %0.5f;\n',varName,i1,ub);
    end
    if isempty(desc)
        fprintf(fid,'%s%d.description = '''';\n',varName,i1);
    else
        fprintf(fid,'%s%d.description = ''%s'';\n',varName,i1,desc);
    end
    if isempty(uvalue)
        fprintf(fid,'%s%d.uncertainty_percent = [];\n',varName,i1);
    else
        fprintf(fid,'%s%d.uncertainty_percent = %0.5f;\n',varName,i1,uvalue);
    end
    if isempty(umethod)
        fprintf(fid,'%s%d.uncertainty_method = '''';\n',varName,i1);
    else
        fprintf(fid,'%s%d.uncertainty_method = ''%s'';\n',varName,i1,umethod);
    end
    fprintf(fid,'\n');
end

fprintf(fid,'%s = [',varName);
for i1=1:length(P)
    fprintf(fid,'%s%d;',varName,i1);
end
fprintf(fid,'];\n');
fprintf(fid,'clear');
for i1=1:length(P)
    fprintf(fid,' %s%d',varName,i1);
end

fclose(fid);

%-------------------------------------------------------
function quitIt(obj,eventD)

a = areYouSure;
switch a
    case 'Yes'
        saveIt;
    case 'Cancel'
        return
end

if getappdata(gcf,'nout')
    uiresume(gcf);
else
    closereq;
end

%------------------------------------------------------
function menuAbout(obj,eventD)

s = ['ParameterDomain Builder and associated software were created at the Berkeley Center for Control and Identification (BCCI) at UC Berkeley.' 10 10 ...
     'Supported by the NSF under grant no. 0113985' 10 10 ...
     'Data Collaboration: http://jagger.me.berkeley.edu/~pack/nsfuncertainty' 10 ...
     'BCCI: http://jagger.me.berkeley.edu'];
msgbox(s,'About','none','modal');

%-------------------------------------------------------
function menuHelp(obj,eventD)

if str2num(version('-release'))>=14
  web PDBuildHelp.html -helpbrowser
else
  web PDBuildHelp.html
end
