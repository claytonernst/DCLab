function response = G_protein(executionType,paramValuesMatrix,u)

ni = nargin; 
no = nargout; 
error( nargchk(1,3,ni) ); 
error( nargoutchk(0,1,no) ); 

if ~ischar(executionType) 
  error('First input to model simulation file must be of type char.') 
end 

switch executionType
    case 'getModelDomain'
      tmp = cell(2,9);
      tmp(1,:) = {'k1';'k2';'k3';'k4';'k5';'k6';'k7';'k8';'Gt'};
      [tmp{2,:}] = deal([-inf inf]);
      response = cell2struct(tmp,{'name';'range'});
  case 'simulate'
        options = odeset('RelTol',5e-6,'AbsTol',5e-9);
        Neval = size(paramValuesMatrix,2);
        response = zeros(1,Neval);
        for i1=1:Neval
            k3 = paramValuesMatrix(3,i1);
            k5 = paramValuesMatrix(5,i1);
            Gt = paramValuesMatrix(end,i1);
            [t,x] = ode15s(@LOCAL_ode,[0 60],[k5/k3;0;Gt;0],options,paramValuesMatrix(:,i1),u);
            Gt = paramValuesMatrix(end,i1);
            y = (Gt - x(end,3))./Gt;
            response(1,i1) = y;
            if response(1,i1) < eps
              options = odeset('RelTol',5e-5,'AbsTol',5e-10);
              [t,x] = ode15s(@LOCAL_ode,[0 60],[k5/k3;0;Gt;0],options,paramValuesMatrix(:,i1),u);
              Gt = paramValuesMatrix(end,i1);
              y = (Gt - x(end,3))./Gt;
              response(1,i1) = y;
              if response(1,i1) < eps        
                disp('the output may be negative, you need to tighten the ode15s tolerances')
                keyboard
              end
            end
        end
    case 'gimme'
        response = @LOCAL_ode;
    case 'isSaveEnabled'
        response = 1;
    
    otherwise
        error('Improper Execution Type')
end

function x_dot = LOCAL_ode(t,x,k,u) %#ok

x1d = -k(1)*x(1)*u + k(2)*x(2)-k(3)*x(1)+k(5);
x2d = k(1)*x(1)*u-k(2)*x(2)-k(4)*x(2);
x3d = -k(6)*x(2)*x(3)+k(8)*(k(9)-x(3)-x(4))*(k(9)-x(3));
x4d = k(6)*x(2)*x(3)-k(7)*x(4);

x_dot = [x1d;x2d;x3d;x4d];
