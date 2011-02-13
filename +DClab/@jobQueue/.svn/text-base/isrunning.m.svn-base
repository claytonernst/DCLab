function out = isrunning(obj)
%ISRUNNING  Finds jobs with status 'running'
%
%   JOBLABELS = ISRUNNING(OBJ) looks in the jobQueue
%   OBJ for jobs with the status 'running' and returns
%   a cell array of the job labels such that,
%   jobstatus(OBJ,JOBLABEL{i})=='running' is true for
%   every JOBLABEL i.
%
%   See also JOBQUEUE, JOBSTATUS, ISFREE.

tmp = obj.jobs;
indx = find(strcmp({tmp.status},'running'));
out = {tmp(indx).indx};
