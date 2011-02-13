function [t displ] = simTrueSystem(F)
% This function plots the displacement (position) of the mass of a 
% one-degree of freedom mass-spring-damper system. The component 
% properties are described by
%
% (nonlinear) spring: F(x) = -k1x + k2x^3 
% damper: F(x') = -bx'
%
% We will simulate this system starting from zero initial
% conditions with a step-input of F units of force

IC = [0;0];
tspan = (0:0.1:150);
[t, traj] = ode45(@EOM,tspan,IC,[],F);
displ = traj(:,1); % Displacement is the first state variable and

meas = displ+0.1; % Add bias
% doctor the data
idx = find(meas <= 0.25);
scale = meas(idx(end))/displ(idx(end));
meas(idx) = displ(idx)*scale; % Correct at t=0

% Plot the displacement
h = plot(t,meas);
xl = xlabel('time (sec)');
yl = ylabel('displacement (m)');
set(h,'linewidth',2);
set([xl yl],'fontname','times','fontsize',14)

function xdotvect = EOM(t,xvect,F) %#ok
% Equations of motion for the nonlinear mass-spring-damper.
% The input t is time and the input xvect is the 2x1 vector
% [x;xdot]. For clarity, some operations that could easily be vectorized 
% have been explicitly written out.

% Model parameters
k1 = 1;
k2 = 0.025;
b = 0.75;
m = 8;

% System state
x = xvect(1);
xdot = xvect(2);

% Sum of forces = ma
xddot = inv(m)*(F - k1*x + k2*x^3 - b*xdot);

% Define the time derivative of xvect as the output
xdotvect = [xdot; xddot];
