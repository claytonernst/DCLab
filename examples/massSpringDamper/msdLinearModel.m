function out = msdLinearModel(flag,paramMat,feature,inputForce)
% Function to simulate the mass-spring-damper system and compute
% the requested feature of the displacement of the mass. A simple
% schematic of the system is pictured below
%
%            \|       b
%            \|-------|]------|===|
%            \|               | m |---> F
%            \|----\/\/\/\/---|===|
%            \|       k
%                               |---> x
%
% Inputs:
%   flag: A char array indicating the task to be performed: return the
%         domain of the model; communicate if model evaluations
%         should be saved; or simulate the model at a particular
%         combination of parameter values to produce the desired feature.
%   paramMat: A 3x1 vector of values for the parameters b,k,and m. Can also
%         be a horizontal concatenation of such vectors.
%   feature: A char array
%     feature == 'peakDis' implies out is the peak displacement over the
%        150 second simulation window.
%     feature == 'rt95' implies out is the time when the displacement first
%        reaches 95% of its peak value.
%     feature == 'ssMean' implies out is the mean displacement over the
%        last 50 seconds of the simulation.
%   inputForce: a scalar indicating the value of the step input of force
%        applied to the mass.
%
% Outputs:
%   out: Depending on flag, either the model domain, a bool indicating if
%        model evaluations should be saved, or the model prediction for the
%        requested output. If paramMat is N horizontally concatenated 3x1
%        vectors (i.e., 3xN), out is 1xN.

if ~ischar(flag)
    error('Inputs: The first input to a DCModelFun must be a char array');
end
switch flag
    case 'getModelDomain'
        % Initialize the 3x1 model domain structure array with fields .name and .range
        modelDomain = cell2struct(cell(3,2),{'name','range'},2);
        modelDomain(1).name = 'b';
        modelDomain(1).range = [0.01 100];
        modelDomain(2).name = 'k';
        modelDomain(2).range = [0.01 100];
        modelDomain(3).name = 'm';
        modelDomain(3).range = [0.01 100];
        out = modelDomain;
        
    case 'isSaveEnabled'
        % Specify that model evaluations should be saved to speed analysis
        out = true;
        
        
    case 'isMultipleResponsesEnabled'
        % There are several features
        out = true;
        
    case 'getResponseList'
        out = {'peakDis';'rt95';'ssMean'};
        
        
    case 'simulate'
        % Determine how many horizontally concatenated parameter vectors
        % were supplied
        N = size(paramMat,2);
        
        tspan = [0 150];  % Define the time interval of simulation
        IC = [0;0];  % Define the initial conditions
        
        % For each of the N parameter vectors, integrate the equations of motion
        % with ode45.
        
        out = zeros(3,N);
        for i1 = 1:N
            % Define the parameter values
            paramVect = paramMat(:,i1);
            b = paramVect(1);
            k = paramVect(2);
            m = paramVect(3);
            
            odefun = @(t,x) localEOM(t,x,inputForce,b,k,m);
            [t, traj] = ode45(odefun,tspan,IC);
            displ = traj(:,1); % Displacement is the first state
            
            % Compute the maximum displacement during the 150 second simulation.
            out(1,i1) = max(displ);
            
            % Determine when the displacement first achieves 95% of its peak
            % value. Since displacement is only available from the simulation at
            % discrete points, linearly interpolate between points on either side
            % of peak95 to estimate the 95% rise time.
            [peak peakIdx] = max(displ);
            peak95 = 0.95*peak;
            idx1 = find(displ(1:peakIdx) < peak95);
            idx2 = find(displ(1:peakIdx) > peak95);
            
            t95 = interp1(displ(idx1(end):idx2(1)),t(idx1(end):idx2(1)),peak95);
            out(2,i1) = t95;
            
            % Average the displacement over the last 50 seconds.
            out(3,i1) = mean(displ(t >= 100));
        end
    otherwise
        error(['Inputs: Behavior for flag value ' flag ' is not defined']);
end
%==subfunction for ODE45==
function xdot = localEOM(t,x,F,b,k,m) %#ok
xdot = [0 1;-k/m -b/m]*x + [0; F/m];
