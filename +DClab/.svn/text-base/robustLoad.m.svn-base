function out = robustLoad(fileName,initialPause)
% function out = robustLoad(fileName)
%
% a function to load a mat file. File name should have a .mat
% extension. The function will load the saved variables the check for
% proper loading. If no = 0, all loaded variables will be assigned to the
% caller workspace. If no=1, loaded variables will be returned in a
% structure array, mimicing the matlab load.

ni = nargin;
no = nargout;

error(nargchk(1,2,ni));
error(nargoutchk(0,1,no));

if ni == 1
  initialPause = 10;
end

pause(initialPause)

%==check input types==
if ~ischar(fileName)
  error('fileName must be a char array')
end

pauseTime = 0.05;
maxTries = 50;

notLoaded = 1;
attempt = 1;
while notLoaded & attempt <= maxTries
  try
    loadedVars = load(fileName);
    notLoaded = 0;
    pause(pauseTime)
  catch
    %eventually log error, for now just dump to screen
    fprintf('Error loading %s, attempt: %s of %s \n',fileName,...
      num2str(attempt), num2str(maxTries));
    pause(pauseTime)
  end
  attempt = attempt+1;
end %while notLoaded & attempt <= maxTries

if notLoaded
  fprintf('Unsuccessful load of %s \n', fileName);
  out = [];
else
  if no==1 
    out = loadedVars;
  elseif no==0
    varNames = fieldnames(loadedVars);
    for i1 = 1:length(varNames)
      assignin('caller',varNames{i1},loadedVars.(varNames{i1}));
    end
  else
    error('robustLoad requires 0 or 1 output argument')
  end
end
