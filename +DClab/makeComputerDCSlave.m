function makeComputerDCSlave(compNum)
%MAKECOMPUTERDCSLAVE runs an infinite loop and waits for tasks
%
%   Run this on each computer that will be doing work
%
%   MAKECOMPUTERDCSLAVE(COMPNUM) starts an infinite loop that will perform
%   tasks designated to the COMPNUM computer by a jobQueue object. Use this
%   when DCOptions.nComputer > 0. Each computer that you execute this file
%   on should be given a distinct integer COMPNUM. No integer should be
%   skipped (i.e. if there are 4 computers, they should be labeled 
%   1,2,3, & 4).

notDone=true;

while notDone

  file = which('multiComputerTempDir');
  instrDir = fileparts(file);
  
  filename = [instrDir filesep 'Instr' num2str(compNum) '.mat'];
  
  if exist(filename,'file')
    initialPause = 1; 
    DClab.robustLoad(filename, initialPause);
    dataforjob = a;

    switch dataforjob.distribute
     case 'modelEvals'
      LOCALmodeleval(filename,dataforjob,instrDir,compNum);
     case 'consistency'
      LOCALconsistency(filename,dataforjob,instrDir,compNum);
     otherwise
      error('Don''t know distribution type')
    end
  else
    pause(5)
  end
  
end

%---------------------------------------------------------------
function LOCALmodeleval(filename,dataforjob,instrDir,compNum)

RM = dataforjob.ResponseModel;
x = dataforjob.designPoints;

% Evaluate model
switch RM.type
 case {'linearModel','quadraticModel'}
  y = eval(RM,x);  
 case 'dcModel'
  try
    addl = RM.additionalInputs;
    y = feval(RM.model,'simulate',x,addl{:});
  catch
    fprintf('error evaluating model handle %s in distribute_evaluations\n',func2str(RM.model));
    keyboard
  end
  otherwise 
    error('Internal inconsistency: condition should never occur')
end

initialPause = 2;
DClab.robustDelete(filename,initialPause);

savename = [instrDir filesep 'Evals' num2str(compNum) 'job' ...
           dataforjob.indx '.mat'];

s.y = y;
s.x = x;
DClab.robustSave(savename,s);
    
%--------------------------------------------------------------
%TODO, this function is dated

function LOCALconsistency(filename,dataforjob,instrDir,compNum)

[uu,lu,ul,ll,l] = checkConsistency(    dataforjob.PDataset, ...
                                     dataforjob.startPt, ...
                                     dataforjob.nodeIndices, ...
                                       dataforjob.opt);

% since we're already loaded the file, save is hopefully not still
% dicking with it. opt for a shorter initial pause
initialPause = 2;
DClab.robustDelete(filename,initialPause);
indices = dataforjob.indices;
    
savename = [instrDir filesep 'Cmeas' num2str(compNum) 'job' ...
           dataforjob.indx '.mat'];
DClab.robustSave(savename,uu,lu,ul,ll,l,indices);
