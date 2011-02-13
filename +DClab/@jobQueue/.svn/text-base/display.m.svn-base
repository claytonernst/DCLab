function display(obj)
%DISPLAY  display function for jobQueue objects
%

tmp = obj.jobs;

disp(' ')

if strcmp(inputname(1),'')
  disp('ans =')
else
  disp([inputname(1) ' ='])
end

disp(' ')

disp('     jobQueue Object:')
disp(' ')

disp(['               Total number of jobs: ' ...
      num2str(length(tmp))])
if length(tmp)
  disp(['                 Total jobs running:  ' num2str(sum(~ ...
                                                  strcmp({tmp.status},'done')))])
end
disp(['              Most recent job index: ' obj.indx])
disp(['    Number of computers using queue: ' num2str(obj.numComps)])
disp(['          Queue creation time stamp: ' num2str(fix(obj.starttime))])

disp(' ')
