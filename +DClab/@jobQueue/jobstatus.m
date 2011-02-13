function out = jobstatus(obj,indxstr,status)
%JOBSTATUS  Get/set job status
%
%   STATUS = JOBSTATUS(OBJ,JOBLABEL) get the status
%   of job JOBLABEL according to jobQueue OBJ.
%
%   OBJ = JOBSTATUS(OBJ,JOBLABEL,STATUS) sets the status
%   of job JOBLABEL and returns the updated jobQueue OBJ.
%   It then updates 'jobQueueObj.mat'.
%
%   See also JOBQUEUE, ISFREE,ISRUNNING.

indx = findjob(obj,indxstr);

if nargin==2
  
  out = obj.jobs(indx).status;
  
elseif nargin==3
  if ~any(strcmp({'running','done','old'},status))
    error(['Invalid status.  Must be ''running'', ''done'', or ' ...
           '''old'''])
    
  end
  tmp = obj.jobs(indx).status;
  obj.jobs(indx).status = status;

  if strcmp(status,'done') && ~strcmp(tmp,'done')
    obj.jobs(indx).endTime = clock;
  end
  
  out = obj;
  
  jobObj = obj;
  instrDir = fileparts(which('multiComputerTempDir'));
  DClab.robustSave([instrDir filesep 'jobQueueObj.mat'],jobObj);
  
else
  
  error('Incorrect number of inputs.')
  
end

