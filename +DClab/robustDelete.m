function [] = robustDelete(fileName,initialPause)
% function [] = robustDelete(fileName)
%
% a function to delete a mat file. File name should be on the path. 
% fileName should include full path to make things smoother. if not, EXIST
% can successfully see a file, but DELETE won't remove it.

ni = nargin;
no = nargout;

error(nargchk(1,2,ni));
error(nargoutchk(0,0,no));

%==check input types==
if ~ischar(fileName)
  error('fileName must be a char array')
end

if ni==1
  initialPause = 20;
end
%start off with a long pause, so if robustSave is still tinkering with a
%file, it finished before we try to delete it. This is the only crack in
%the methodology.
pause(initialPause)

pauseTime = 4;
maxTries = 10;

notDeleted = 1;
attempt = 1;
while notDeleted & attempt <= maxTries
  if ~exist(fileName)
    fprintf('Warning, attempting to delete noexistant file %s, retrying anyway \n', fileName);
    fileExists = 0;
  else
    fileExists = 1;
  end
  delete(fileName)
  pause(pauseTime);
  if fileExists & exist(fileName) %if we originally found the file, 
                                  %and still do after delete attempt
    fprintf('Unsuccessful delete of %s, attempt %d of %d \n', fileName, attempt, maxTries);
  else
    notDeleted = 0;
  end
  attempt = attempt+1;
end % while notDeleted & attempt <= maxTries

    
  
  
