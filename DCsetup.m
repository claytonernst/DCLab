%function DCsetup(dispOn)

% DCSETUP
%
%   DCsetup
%
%   adds needed directories to the path for the Data Collaboration
%   tools.

% make a script so imports work

%if nargin == 0
    dispOn = true;
%end

if ~isdeployed
    mainpath = fileparts(which('DCsetup'));

    DCpath = fullfile(mainpath,'source');

    P = {mainpath, ...
        fullfile(mainpath, 'examples', 'GRI'), ...
        fullfile(mainpath, 'examples', 'massSpringDamper'), ...
        fullfile(mainpath, 'savedEvaluations'), ...
        fullfile(mainpath, 'multiComputerTemp'), ...
        fullfile(mainpath, 'PrIMeGUI'), ...
        };

    if ~exist('sedumi','file')
        P = [P,{fullfile(mainpath, 'SeDuMi_1_21'),fullfile(mainpath, 'SeDuMi_1_21', 'conversion')}];
    else

        %Do nothing in this case, we detected a version of SeDuMi. If it
        %is outdated, the test script should detect this.
    end

    if dispOn

        disp('Appending directories to Matlab path:');

        for k = 1:length(P)
            path(P{k}, path);
            fprintf(1,'  %s\n', P{k});
        end

        disp('Path updated: run dctest.m to confirm successful installation')

    else
        for k = 1:length(P)
            path(P{k}, path);
        end
    end
    
    import DClab.*
    if dispOn
        disp('')
        disp('Imported the DClab package')
    end

    clear P k

end
