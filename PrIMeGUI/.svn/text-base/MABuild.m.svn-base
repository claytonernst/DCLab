function MABuild



ssize = get(0,'screensize');



figHan = figure('name', 'ModelAssertion Builder', ...
                'numbertitle', 'off', ...
                'menubar', 'none',...
                'position',[100 100 0.6*ssize(3) 0.8*ssize(4)]);



x=1/19;

p1 = uipanel(figHan,'title','Description',...
                    'position',x*[1 16 17 2]);

p2 = uipanel(figHan,'title','Parameter List', ...
                    'position',x*[1 15 17 1]);

p3 = uipanel(figHan,'title','Simulation Code', ...
                    'position',x*[1 9 17 6]);

p4 = uipanel(figHan,'title','Sub-Functions', ...
                    'position',x*[1 4 17 5]);

p5 = uipanel(figHan,'title','Variable Name', ...
                    'position',x*[1 3 17 1]);

bt = uicontrol(figHan,'style','pushbutton',...
                    'string','Create Object and m-file',...
                    'units','normalized',...
                    'position',x*[6 1 7 1],...
                    'callback',@makeIt);



desc = uicontrol(p1,'style','edit',...
                    'units','normalized',...
                    'position',[0.1 0.1 0.8 0.8],...
                    'HorizontalAlignment','left',...
                    'tooltipstring','Description of model',...
                    'max',3,'min',1,...
                    'backgroundcolor',[0.9 0.9 0.9]);



param = uicontrol(p2,'style','edit',...
                     'units','normalized',...
                     'position',[0.1 0.1 0.8 0.8],...
                     'HorizontalAlignment','left',...
                     'tooltipstring','{''name1'',''name2'',...}',...
                     'max',3,'min',1,...
                     'backgroundcolor',[0.9 0.9 0.9]);


simtext = uicontrol(p3,'style','text',...
                       'units','normalized',...
                       'position',[0.05 0.6 0.9 0.4],...
                       'HorizontalAlignment','left',...
                       'string',['Assume the variable "paramValues", ', ...
                                 'a vector, is the input and is ', ...
                                 'already defined.  Your output should ',...
                                 'be assigned to the variable ' ...
                                 '"out".']);

simc = uicontrol(p3,'style','edit',...
                    'units','normalized',...
                    'position',[0.05 0.05 0.9 0.8],...
                    'HorizontalAlignment','left',...
                    'tooltipstring','Simulation Code',...
                    'max',4,'min',2,...
                    'backgroundcolor',[0.9 0.9 0.9]);



substext = uicontrol(p4,'style','text',...
                       'units','normalized',...
                       'position',[0.05 0.6 0.9 0.4],...
                       'HorizontalAlignment','left',...
                       'string',['Type in any needed subfunctions ' ...
                                 'here (including declarations).']);

subs = uicontrol(p4,'style','edit',...
                    'units','normalized',...
                    'position',[0.05 0.05 0.9 0.8],...
                    'HorizontalAlignment','left',...
                    'tooltipstring','Sub Functions',...
                    'max',4,'min',2 ,...
                    'backgroundcolor',[0.9 0.9 0.9]);



var = uicontrol(p5,'style','edit',...
                   'units','normalized',...
                   'position',[0.1 0.1 0.8 0.8],...
                   'HorizontalAlignment','left',...
                   'tooltipstring','Name of variable to save in workspace',...
                   'backgroundcolor',[0.9 0.9 0.9]);



setappdata(figHan,'hans',{desc,param,simc,subs,var});



%--------------------------------------------------



function makeIt(obj,eventD)



tmp = getappdata(gcf,'hans');

[desc,param,simc,subs,var] = deal(tmp{:});



varName = get(var,'string');



%write m-file

fid = fopen([pwd filesep varName '.m'],'w');

fprintf(fid,['function response = %s(executionType,' ...

             'paramValuesMatrix,flag)\n'],varName);

fprintf(fid,'ni = nargin;\nno = nargout;\n');

fprintf(fid,['error(nargchk(1,3,ni));\nerror(nargoutchk(0,1,no));\' ...

             'n']);

fprintf(fid,'if ~ischar(executionType)\n');

fprintf(fid,['  error([''First input to model simulation file must ' ...

             'be of type char.''])\n']);

fprintf(fid,'end\n\n');



fprintf(fid,'switch executionType\n');

fprintf(fid,' case ''description''\n');

descrip = get(desc,'string');

fprintf(fid,'  response = [');

for i=1:size(descrip,1)

  fprintf(fid,'''%s '', 10, ...\n',descrip(i,:));

end

fprintf(fid,'             ];\n');

fprintf(fid,' case ''observable''\n');

fprintf(fid,'  response = '''';\n');

fprintf(fid,' case ''featureList''\n');

fprintf(fid,'  response = {''NA''};');

fprintf(fid,' case ''featurePrecedence''\n');

fprintf(fid,'  response = 1;\n');

fprintf(fid,' case ''paramList''\n');

fprintf(fid,'  response = %s;\n',get(param,'string'));

fprintf(fid,' case ''simulate''  %%assign parameter values\n');

fprintf(fid,['  N = size(paramValuesMatrix,1);\n  response = ' ...

             'zeros(N,1);\n\n']);

fprintf(fid,'  %assign parameter values\n');

fprintf(fid,'  for i1 = 1:N\n');

fprintf(fid,'    paramValues = paramValuesMatrix(i1,:);\n');

simcode = get(simc,'string');

for i=1:size(simcode,1);

  fprintf(fid,'%s\n',simcode(i,:));

end

fprintf(fid,'    response(i1) = out;\n');

fprintf(fid,['  end\n otherwise\n  error(''Incorrect input type to ' ...

             'model simulation file'')\nend\n\n']);



subfunc = get(subs,'string');

for i=1:size(subfunc,1)

  fprintf(fid,'%s\n',subfunc(i,:));

end



fclose(fid);



Obj = ModelAssertion(varName,'NA');



assignin('base',varName,Obj);

