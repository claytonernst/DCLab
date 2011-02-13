function out = isold(obj)
%ISOLD  Finds jobs with status 'old'
%
%   JOBLABELS = ISOLD(OBJ) looks in the jobQueue
%   OBJ for jobs with the status 'old' and returns
%   a cell array of the job labels such that,
%   jobstatus(OBJ,JOBLABEL{i})=='old' is true for
%   every JOBLABEL i.
%
%   See also JOBQUEUE, JOBSTATUS, ISFREE, ISRUNNING.

tmp = obj.jobs;
indx = find(strcmp({tmp.status},'old'));
out = {tmp(indx).indx};
