function obj = clean(obj,howOld)
%CLEAN  Cleans out old jobs
%
%  OBJ = CLEAN(OBJ) cleans out all jobs with status 'old' that
%  were created over an hour ago.  It does so by creating a false
%  output mat-file so the computer can receive instruction files
%  again.  OBJ is the jobQueue object.  The function just returns
%  the updated jobQueue object OBJ.
%
%  OBJ = CLEAN(OBJ,HOWOLD) lets you specify how old a job has to be to be
%  cleaned out.  The function will clean out all jobs that have
%  status 'old' and were created over HOWOLD seconds ago.


if nargin==1
  howOld = 36000;
end

file = which('instructionFlag');
instrDir = fileparts(file);
filename = [instrDir filesep 'Out'];

jobLabels = isold(obj);
dummy = [];
tmp = obj.jobs;

try
  for i = jobLabels
    if jobtime(obj,i{1})>howOld  %if older than an hour, dump it.
      comp = tmp(findjob(obj,i{1})).compindx;
      save([filename num2str(comp) 'job' i{1}],'dummy');
      obj = jobstatus(obj,i{1},'done');
    end
  
  end
catch

end
