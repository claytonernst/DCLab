classdef jobQueue < DClab.DCObject
    %JOBQUEUE  Constructor function for jobQueue objects
    %
    %   OBJ = JOBQUEUE(NUMCOMPS) does one of two things.
    %   If 'jobQueueObj.mat' exists, it is loaded, the number of
    %   computers is set to NUMCOMPS, and it is returned.  If
    %   the mat file does not exist, it will create it, and
    %   return the object.
    %
    %   See also FINDJOB, JOBSTATUS, JOBTIME, CREATEJOB, ISFREE, ISRUNNING.

    properties
        jobs;
        indx;
        starttime;
        numComps;
    end

    methods
        
        function jobObj = jobQueue(numComps)

            if nargin>0
                %create new object
                jobObj.jobs = struct('indx',[],'compindx',[],'status','done','startTime',[],'endTime',[]);
                jobObj.indx = '0';
                jobObj.starttime = clock;
                jobObj.numComps = numComps;
                if exist('jobQueueObj.mat','file')
                    %load it
                    initialPause = 0; % we shouldn't be trying to load this
                    % unless we've already sucessfully saved it
                    DClab.robustLoad('jobQueueObj.mat', initialPause);
                    jobObj = a;
                    jobObj.numComps = numComps;
                end
            end

        end

        function [list,sz] = displayProps(obj)
            list = [];
            sz = '';
        end

    end %public methods

end %classdef
