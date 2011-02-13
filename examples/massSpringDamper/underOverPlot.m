function underOverPlot

F = 1;
b = 0.5;
k = 1;
m = 5;
IC = [0;0];
tspan = [0 100];

[t1, traj] = ode45(@EOM,tspan,IC,[],F,b,k,m);
pos1 = traj(:,1); %grab the state trajectory corresponding to position

F = 1;
b = 10;
k = 1;
m = 5;
IC = [0;0];
tspan = [0 100];

[t2, traj] = ode45(@EOM,tspan,IC,[],F,b,k,m);
pos2 = traj(:,1); %grab the state trajectory corresponding to position
h = plot(t1,pos1,'r',t2,pos2,'b');
legend(h,'underdamped','overdamped')

xl = xlabel('time (sec)');
yl = ylabel('position (m)');
set(h,'linewidth',2);
set([xl yl],'fontname','times','fontsize',14)

%==subfunction for ODE45==
function xdot = EOM(t,x,F,b,k,m)
xdot = [0 1;-k/m -b/m]*x + [0; F/m]; 
