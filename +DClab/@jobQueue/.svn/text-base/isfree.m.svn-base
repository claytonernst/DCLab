function out = isfree(obj,compindx)
%ISFREE  Looks for computers not running a job.
%
%   BOOLEAN = ISFREE(OBJ,COMPINDX) specifies whether
%   the computer defined by COMPINDX is available to 
%   process a job, according to the jobQueue OBJ
%
%   ARRAY = ISFREE(OBJ) returns an array of the computer
%   indices that are currently not processing a job.
%
%   See also JOBQUEUE, ISRUNNING, JOBSTATUS.

if isempty(obj.jobs)
  if nargin==2
    out = true;
  else
    out = 1:6;
  end
else
  if nargin==2
    tmp = obj.jobs;
    indx = find(~strcmp({tmp.status},'done'));
    comps = [tmp(indx).compindx];
    out = isempty(find(comps==compindx,1));
  elseif nargin==1
    tmp = obj.jobs;
    indx = find(~strcmp({tmp.status},'done'));
    comps = [tmp(indx).compindx];
    out = setdiff(1:obj.numComps,comps);
  else
    error('Incorrect number of inputs.')
  end
end
