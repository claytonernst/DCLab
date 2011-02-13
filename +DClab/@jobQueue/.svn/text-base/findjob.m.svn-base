function indx = findjob(obj,indxstr)
%FINDJOB    Attempts to find job indx based on job label.
%
%   INDX = FINDJOB(OBJ,JOBLABEL) returns the index of where
%   job JOBLABEL is located in the jobQueue OBJ, such that 
%   OBJ.jobs(INDX).indx == JOBLABEL is true.
%
%   NOTE: this is mostly just used by other jobQueue methods
%
%   See also JOBQUEUE, CREATEJOB.

tmp = obj.jobs;
indx = find(strcmp({tmp.indx},indxstr));

if isempty(indx)
  error('No such job')
end

