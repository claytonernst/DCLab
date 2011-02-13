function objOut = OptimOptions(obj)

ss = get(0,'screensize');
posx = max([0,ss(3)-670])/2;
posy = max([0,ss(4)-450])/2;

ni = nargin;
no = nargout;
error(nargoutchk(0,1,no));
error(nargchk(0,1,ni));

OptimOptns = {};
OptimOptns.method = 1;
OptimOptns.tolFun = 1.0e-6;
OptimOptns.objectiveWeightRatio = 1;

if (ni==0) || isempty(obj)
    obj = OptimOptns;
end

%figure
figH = figure('MenuBar','none',...
    'NumberTitle','off',...
    'name','Optimization Options',...
    'ToolBar','none',...
    'resize','off',...
    'dockcontrols','off',...
    'closerequestfcn','uiresume(gcf)',...
    'position',[posx posy 600 260]);

%Buttons
uicontrol(figH,'style','pushbutton',...
    'string','Resort to Defaults',...
    'callback',@local_defaults,...
    'position',[15 10 125 20],...
    'selected','off');
uicontrol(figH,'style','pushbutton',...
    'string','OK',...
    'callback','uiresume(gcf)',...
    'position',[470 10 40 20 ],...
    'selected','off');
uicontrol(figH,'style','pushbutton',...
    'string','Cancel',...
    'callback','setappdata(gcf,''opt'',getappdata(gcf,''original''));uiresume(gcf)',...
    'position',[525 10 60 20],...
    'selected','off');

%title
uicontrol(figH,'style','text',...
    'string','Options',...
    'fontsize',16,...
    'horizontalalignment','left',...
    'backgroundcolor',[0.8 0.8 0.8],...
    'position',[20 225 80 25]);

%panel
mainpan = uipanel(figH,'title','',...
    'units','pixels',...
    'position',[10 40 580 175]);

uicontrol(mainpan,'style','text',...
    'string','Optimization method:',...
    'Horizontalalignment','left',...
    'position',[20 130 110 20]);
method = uicontrol(mainpan,'style','popupmenu',...
    'position',[130 135 170 20],'Enable','off',...
    'string',{'LS-H','LS-F','1N-F','MO-F'},...
    'callback',{@local_pulldowns,'method'});

owrtext = uicontrol(mainpan,'style','text',...
    'horizontalalignment','left',...
    'position',[20,90,220,20],'Visible','off',...
    'string','Weight ratio of objective functions 2 to 1:');
objectiveWeightRatio = uicontrol(mainpan,'style','edit',...
    'position',[240 95 60 20],'Visible','off',...
    'callback',{@local_generic_set,'objectiveWeightRatio'});

uicontrol(mainpan,'style','text',...
    'string','Objective Function Convergence Tolerence:',...
    'horizontalalignment','left',...
    'position',[20 50 220 20]);
tolFun = uicontrol(mainpan,'style','edit',...
    'position',[240 55 60 20],...
    'callback',{@local_generic_set,'tolFun'});


%--------save handles-------------
handles.fig = figH;
handles.tolFun = tolFun;
handles.method = method;
handles.objectiveWeightRatio = objectiveWeightRatio;

setappdata(figH,'handles',handles);
setappdata(figH,'isSaved',0);
setappdata(figH,'original',obj);

local_setForm([],[],obj);

uiwait(figH);
objOut = getappdata(figH,'opt');
delete(figH);

%-----------import data---------------------
    function local_setForm(obj,eventD,opt)
        
        handles = getappdata(gcf,'handles');
        
        switch opt.method
            case 1
                set(handles.method,'value',1);
            case 2
                set(handles.method,'value',2);
            case 3
                set(handles.method,'value',3);
            case 4
                set(handles.method,'value',4);
        end
        if (opt.method == 4)
            set(objectiveWeightRatio,'Visible','on')
            set(owrtext,'Visible','on')
        else
            set(objectiveWeightRatio,'Visible','off')
            set(owrtext,'Visible','off')
        end

        set(handles.tolFun,'string',num2str(opt.tolFun));
        set(handles.objectiveWeightRatio,'string',num2str(opt.objectiveWeightRatio));
        setappdata(gcf,'opt',opt);
    end

%----generic----
    function local_generic_set(obj,eventD,prop)
        opt = getappdata(gcf,'opt');
        try
            try
                a = eval(get(obj,'string'));
            catch
                a = get(obj,'string');
            end
            opt.(prop)=a;
            setappdata(gcf,'isSaved',1);
        catch
            msgbox('Invalid Entry: Option was not modified','Error','Error','modal');
            set(obj,'ForegroundColor',[1 0 0]);
            return
        end
        %change color to black:
        set(obj,'ForegroundColor',[0 0 0]);
        setappdata(gcf,'opt',opt);
    end

%----pull-downs----
    function local_pulldowns(obj,eventD,prop)
        opt = getappdata(gcf,'opt');
        switch prop
            case 'method'
                types = {1,2,3,4};
        end
        opt.(prop)=types{get(obj,'value')};
        if (opt.(prop) == 4)
            set(objectiveWeightRatio,'Visible','on')
            set(owrtext,'Visible','on')
        else
            set(objectiveWeightRatio,'Visible','off')
            set(owrtext,'Visible','off')
        end
        setappdata(gcf,'opt',opt);
        setappdata(gcf,'isSaved',1);
    end

%----reset options----
    function local_defaults(obj,eventD)
        optDef = OptimOptns;
        local_setForm([],[],optDef);
    end
end
