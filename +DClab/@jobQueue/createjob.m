function [obj,exitflag,indx] = createjob(obj, dataforjob)
%CREATEJOB  Adds job to jobQueue, and creates instruction file.
%
%   OBJ = CREATEJOB(OBJ,BNDS,MHANDLE,N)
%   creates the job corresponding to the given inputs.  OBJ is the
%   jobQueue object.  N is the number of points to compute.  The
%   function returns the updated jobQueue object as well as updates
%   the jobQueue mat file.
%
%   [OBJ,EXITFLAG] = CREATEJOB(...) returns an exitflag.  EXITFLAG
%   is 1 if the job successfully was created.  EXITFLAG is 0 if
%   there are no computers currently available to perform the job.
%
%   See also JOBQUEUE, FINDJOB.

indx = obj.indx;

indx = double(indx);

%update indx
notDone = 1;
where = length(indx);
while notDone
  indx(where) = indx(where)+1;
  switch indx(where)
   case 58
    indx(where) = 65;
    notDone=0;
   case 91
    indx(where) = 97;
    notDone=0;
   case 123
    indx(where) = 48;
    where = where-1;
   otherwise
    notDone=0;
  end
  if where==0
    indx = [49 indx];
    notDone=0;
  end
end

obj.indx = char(indx);
dataforjob.indx = char(indx);

compindx = min(isfree(obj));

if isempty(compindx)
  exitflag = 0;
else
  file = which('multiComputerTempDir');
  instrDir = fileparts(file);
  filename = [instrDir filesep 'Instr' num2str(compindx) '.mat'];

  DClab.robustSave(filename,dataforjob);

  obj.jobs(end+1).indx = char(indx);
  obj.jobs(end).compindx = compindx;
  obj.jobs(end).status = 'running';
  obj.jobs(end).startTime = clock;
  obj.jobs(end).endTime = [];

  jobObj = obj;
  DClab.robustSave([instrDir filesep 'jobQueueObj.mat'],jobObj);

  exitflag = 1;
  indx = char(indx);

end
