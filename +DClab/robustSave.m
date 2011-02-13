function bool = robustSave(fileName, varargin)
% function [] = robustSave(fileName,varargin)
%
% a function to save the variables in the cell array 
% varargin to the char array fileName. File name should have a .mat
% extension. The function will load the saved variables the check for
% proper saving and will retry in an error is encountered. Additionally and
% file saving errors are logged in the file robustSave.log

ni = nargin;
no = nargout;

error(nargchk(2,27,ni));
error(nargoutchk(0,1,no));

%==check input types==
if ~ischar(fileName)
  error('fileName must be a char array')
end

if length(varargin) > 26
  error('im not smart enough to deal with more than 26 variables to save')
end

pauseTime = 0.75;
maxTries = 10;

%assign varargin to nice variables
saveStr = '';
for i1 = 1:length(varargin)
  eval([char(96+i1) ' = varargin{i1};'])
  saveStr = [saveStr ',''' char(96+i1) ''''];
end

%==try to save a maximum of maxTries times===
notSaved = 1;
attempt = 1;
while notSaved & attempt <= maxTries
  try
    eval(['save(fileName' saveStr ')']);
    notSaved = 0;
    pause(pauseTime)
  catch
    %eventually log error, for now just dump to screen
    fprintf('Error saving %s, attempt: %s of %s \n',fileName,...
      num2str(attempt), num2str(maxTries));
    pause(pauseTime)
  end
  attempt = attempt+1;
end

%now check if it loaded
if notSaved == 0 
  %matlab thinks it saved sucessfully. now verify.
  notLoaded = 1;
  attempt = 1;
  while notLoaded & attempt <= maxTries
    try
      saveVars = load(fileName);
      savedOk = 1; %tenatively we saved ok since the load on the prev line didn't bomb
      for i1 = 1:length(varargin)
        temp = [];
        eval(['temp = ~isequal(saveVars.(char(i1+96)),' char(i1+96) ');']); 
        if temp
          fprintf('In robustSave, mismatch between loaded and saved variables, attempt: %s of %s \n', ...
            num2str(attempt), num2str(maxTries));
        end
      end
        
      if savedOk == 1
        notLoaded = 0;
      else
        pause(pauseTime)
      end
    
    catch
      %log error
      savedOk = 0; %code in the try bombed, so we didn't sucessfully save
      fprintf('In robustSave, error loading saved file %s, attempt %s of %s \n', ...
            fileName, num2str(attempt), num2str(maxTries));
      pause(pauseTime)
    end
    attempt = attempt+1;
  end %while notLoaded & attempt <= maxTries
  
  if savedOk == 1
    bool = 1;
  else
    %log error
    fprintf('Unsuccessful save of %s \n', fileName);
    bool = 0;
  end

else
  %log error
  fprintf('Unsuccessful save of %s \n', fileName);
  bool = 0;
end
  
  
  
