function passed = testBuildRelaxedEqualities(dispBool,plotBool)
%TestBuildRelaxedEqualities returns true if it appears OK
%
%   PASSED = TestBuildRelaxedEqualities returns true or false.
%   Several messages will be displayed to the screen.
%
%   PASSED = TestBuildRelaxedEqualities(DISPBOOL) will suppress
%   all screen displays if DISPBOOL==false.
%
%   PASSED = TestBuildRelaxedEqualities(DISPBOOL,PLOTBOOL) will display a
%   plot showing the orginal equality constraint and the quadratic
%   relaxation when PLOTBOOL==TRUE. The default is FALSE.

if nargin == 0
  dispBool = true;
  plotBool = false;
end
if nargin == 1
  plotBool = false;
end

passed = true;

if dispBool
  disp('===Testing nqcqp/private/buildRelaxedEqualities function===')
end

% suppose \rho \in [a,b] and xl = \rho, xg = log10(\rho)
a = 1;
b = 100;

% Create an object that should call the private function

n = 5; % number of decision variables
Znot = eye(n+1);
ZiCell = [];
%LB = [a; log10(a)];
%UB = [b; log10(b)];
LB = -2*rand(n,1);
UB = 2*rand(n,1);

idxLin = 4; 
idxLog = 2;

c1 = b/a;
c2 = a;
c3 = b-a;
c4 = LB(idxLin);
c5 = UB(idxLin) - LB(idxLin);
c6 = LB(idxLog);
c7 = UB(idxLog) - LB(idxLog);

linXlogX = [idxLin idxLog c1 c2 c3];

obj = nqcqp(Znot,ZiCell,LB,UB,linXlogX);

% points used for ploting
N = 200; 
xl = linspace(LB(idxLin),UB(idxLin),N);
xg = linspace(LB(idxLog),UB(idxLog),N);

%fxg = c5*c2/c3 * ( c1.^((xg-c6)./c7) - 1 ) + c4;
finvxl = c7/log10(c1) * log10( c3/(c2*c5)*(xl-c4) + 1 ) + c6;
t = struct(obj);
Zcell = t.ZrelaxedEq;

% % plot fxg and the quadratic approximations that should over and under bound it.
% qUB = -Zcell{1};
% coeffUB0 = qUB(1,1);
% coeffUB1 = 2*qUB(idxLog+1,1);
% coeffUB2 = qUB(idxLog+1,idxLog+1);
% 
% qLB = Zcell{2};
% coeffLB0 = qLB(1,1);
% coeffLB1 = 2*qLB(idxLog+1,1);
% coeffLB2 = qLB(idxLog+1,idxLog+1);
% 
% figure
% plot(xg,fxg,'b',xg,coeffLB0 + xg*coeffLB1 + xg.^2*coeffLB2,'r',xg,coeffUB0 + xg*coeffUB1 + xg.^2*coeffUB2,'k')
% 
% % plot finvxl and the quadratic approximations that should over and under bound it.
% pUB = -Zcell{3};
% coeffUB0 = pUB(1,1);
% coeffUB1 = 2*pUB(idxLin+1,1);
% coeffUB2 = pUB(idxLin+1,idxLin+1);
% 
% pLB = Zcell{4};
% coeffLB0 = pLB(1,1);
% coeffLB1 = 2*pLB(idxLin+1,1);
% coeffLB2 = pLB(idxLin+1,idxLin+1);
% 
% figure
% plot(xl,finvxl,'b',xl,coeffLB0 + xl*coeffLB1 + xl.^2*coeffLB2,'r',xl,coeffUB0 + xl*coeffUB1 + xl.^2*coeffUB2,'k')

% Create images of the points (xl,xg) that satisfy the relaxed constraint
[XL XG] = meshgrid(xl,xg);

const1 = zeros(N);
const2 = zeros(N);
const3 = zeros(N);
const4 = zeros(N);

for i1 = 1:N
  for i2 = 1:N
    xx = [1; zeros(n,1)];
    xx(idxLin+1) = XL(i1,i2);
    xx(idxLog+1) = XG(i1,i2);
    const1(i1,i2) = xx'*Zcell{1}*xx;
    const2(i1,i2) = xx'*Zcell{2}*xx;
    const3(i1,i2) = xx'*Zcell{3}*xx;
    const4(i1,i2) = xx'*Zcell{4}*xx;
  end
end

% Feasible sets for quadratic constraints
feas1 = const1<=0 & const2<=0;
feas2 = const3<=0 & const4<=0;
feasAll = feas1 & feas2;

% Actual feasible set
FCXG = c5*c2/c3 * ( c1.^((XG-c6)./c7) - 1 ) + c4;
feasAct = abs(XL-FCXG) <= 1e-2/(UB(idxLin)-LB(idxLin));

if plotBool
  figure
  subplot(2,2,1)
  image(xl,xg,100*feas1)
  title('Constraint from fitting x_l-f_c(x_g)')
  xlabel('x_l')
  ylabel('x_g')
  set(gca,'YDir','normal')

  subplot(2,2,2)
  image(xl,xg,100*feas2)
  title('Constraint from fitting x_g - f^{-1}_c(x^l)')
  xlabel('x_l')
  ylabel('x_g')
  set(gca,'YDir','normal')

  subplot(2,2,3)
  image(xl,xg,100*feasAll)
  title('Both constraints')
  xlabel('x_l')
  ylabel('x_g')
  set(gca,'YDir','normal')

  subplot(2,2,4)
  plot(xl,finvxl)
  set(gca,'xlim',[xl(1) xl(end)])
  title('Actual contraint x^l = f_c(x_g)')
  xlabel('x_l')
  ylabel('x_g')
  set(gca,'xlim',[xl(1) xl(end)])
  set(gca,'ylim',[finvxl(1) finvxl(end)])
end

if ~isequal(feasAct & feasAll,feasAct)
  passed = false;
  if dispBool
    disp('    Error, relaxation does not appear to contain original constraint')
  end
end
