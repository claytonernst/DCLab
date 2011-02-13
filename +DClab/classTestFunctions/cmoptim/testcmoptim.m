% This is NOT comprehensive. It is intended only to find big problems.

% Lets pose a very simple problem.

% min \gamma s.t. 
% -1 <= x <= 1;
% -1(1+\gamma) <= 5*x - 2 <= 1*(1+\gamma).

Zquad = {[0 2.5; 2.5 0]};
d = 2;
u = [-1 1];
LB = -1;
UB = 1;

obj = cmoptim(Zquad,d,u,LB,UB);

lb = lowerBnd(obj);
ub = upperBnd(obj);

% Previous works. Lets make the constraint quadratic.

profile on
Zquad = {[0 2.5; 2.5 3]};
obj = cmoptim(Zquad,d,u,LB,UB);

lb = lowerBnd(obj);
ub = upperBnd(obj);
profile viewer

% previous works, now lets do some funny business.

% let x2 = log10(x1);

Zquad = rand(3);
Zquad = {Zquad + Zquad'};
d = 0.5;
u = [-1 1];
LB = [1;0];
UB = [10;1];

% Let [i j c1 c2 c3] denote a row of LINXLOGX. This row imposes the
%     equality constraint 
%             UB(i)-LB(i)
%      X(i) = ---------- c2 ( c1^[(X(j)-LB(j))/(UB(j)-LB(j))] - 1 ) + LB(i)
%                c3

linXlogX = [1 2 10 1 9];
obj = cmoptim(Zquad,d,u,LB,UB,[0 0],1,[0 0],{'none'},linXlogX);







