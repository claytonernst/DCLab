function EABuild



ssize = get(0,'screensize');



figHan = figure('name','ExperimentAssertion Builder', ...
                'numbertitle','off', ...
                'menubar', 'none');



x=1/7.5;

p1 = uipanel(figHan,'title','Description', ...
                    'position',x*[0.5 5 6.5 2]);

p2 = uipanel(figHan,'title','Measured Data Point', ...
                    'position',x*[0.5 4 6.5 1]);

p3 = uipanel(figHan,'title','Data Uncertainty', ...
                    'position',x*[0.5 3 6.5 1]);

p4 = uipanel(figHan,'title','Variable Name', ...
                    'position',x*[0.5 2 6.5 1]);

bt = uicontrol(figHan,'style','pushbutton', ...
                      'string','Create Object',...
                      'units','normalized',...
                      'position',x*[2 0.5 3.5 1],...
                      'callback',@makeIt);



desc = uicontrol(p1,'style','edit',...
                    'units','normalized',...
                    'position',[0.1 0.1 0.8 0.8],...
                    'HorizontalAlignment','left',...
                    'tooltipstring','Description of experiment',...
                    'max',3,'min',1,...
                    'backgroundcolor',[0.9 0.9 0.9]);



d = uicontrol(p2,'style','edit',...
                 'units','normalized',...
                 'position',[0.1 0.1 0.8 0.8],...
                 'HorizontalAlignment','left',...
                 'tooltipstring','d',...
                 'backgroundcolor',[0.9 0.9 0.9]);



u = uicontrol(p3,'style','edit',...
                 'units','normalized',...
                 'position',[0.1 0.1 0.8 0.8],...
                 'HorizontalAlignment','left',...
                 'tooltipstring','u',...
                 'backgroundcolor',[0.9 0.9 0.9]);



var = uicontrol(p4,'style','edit',...
                   'units','normalized',...
                   'position',[0.1 0.1 0.8 0.8],...
                   'HorizontalAlignment','left',...
                   'tooltipstring','Name of variable to be save in workspace',...
                   'backgroundcolor',[0.9 0.9 0.9]);



setappdata(figHan,'hans',{desc,d,u,var});



%-------------------------------------------------

function makeIt(obj,eventData)



tmp = getappdata(gcf,'hans');

[deschan,dhan,uhan,varhan] = deal(tmp{:});



varName = get(varhan,'string');

d = eval(get(dhan,'string'));

u = eval(get(uhan,'string'));

desc = get(deschan,'string');



Obj = ExperimentAssertion(d,u,desc);



assignin('base',varName,Obj);

