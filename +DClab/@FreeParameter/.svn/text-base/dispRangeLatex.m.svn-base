function out = dispRangeLatex(obj)
%DISPLAYRANGELATEX creates a char array formatted for a LaTeX table
%
%   STR = DISPLAYRANGELATEX(OBJ) creates a character array STR displaying
%   the parameter ranges and nominal values of the FreeParameter object OBJ
%   with appropriate formatting and line breaks for a LaTeX table.


names = {obj.name};
ranges = {obj.range};
nom = [obj.nominal];

n = length(names);

formStr = '%0.5g';

%formating latex table
[myand{1:n,1}] = deal('  &  ');
[slash{1:n,1}] = deal('\\');
out = [char(names) char(myand)  num2str(ranges(:,1),formStr) char(myand) ...
       num2str(nom,formStr) char(myand) num2str(ranges(:,2),formStr) char(slash)];
out = strvcat('Parameter Name & Lower Bound & Nominal & Upper Bound',out); 
