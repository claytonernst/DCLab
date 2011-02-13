function out = jobtime(obj,indxstr)
%JOBTIME   Elapsed time for a job
%
%   TIME = JOBTIME(OBJ,JOBLABEL) returns the elapsed time
%   for the job JOBLABEL according to the jobQueue OBJ.
%   If the job has finished, it returns the total time the 
%   job ran.
%
%   See also JOBQUEUE.

indx = findjob(obj,indxstr);

if isempty(obj.jobs(indx).endTime)
  
  out = etime(clock,obj.jobs(indx).startTime);
  
else
  
  out = etime(obj.jobs(indx).endTime,obj.jobs(indx).startTime);
  
end
