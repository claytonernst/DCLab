function passed = testResponseModelEval(dispBool)
%TESTRESPONSEMODELEVAL returns true if it appears OK
%
%   PASSED = TESTRESPONSEMODELEVAL returns true or false. Several
%   messages will be displayed to the screen.
%
%   PASSED = TESTRESPONSEMODELEVAL(DISPBOOL) will suppress all
%   screen displays if DISPBOOL==false.

if nargin == 0
  dispBool = true;
end

passed = true;

if dispBool
  disp('===Testing ResponseModel eval method===')
end

if dispBool
  disp('    Linear model, no outputUncertainty')
end

coeffMatrix = [3 4 5]';
domain = struct('name',{'p1';'p2'},'range',{[-10 10];[-inf inf]});
linRM = ResponseModel(coeffMatrix,domain);

% Evaluate model at p1 = 4, p2 = 7. Should get 3+4*4+5*7 = 54
[y yrange] = eval(linRM,[4; 7]);
if y ~= 54 || yrange(1) ~= yrange(2)
  passed = false;
  if dispBool
    disp('   ERROR: computation incorrect')
  end
end

if dispBool
  disp('    Linear model in transformed coordinates, no outputUncertainty')
end

structWithTrans.value = [-4 0.5 1]';
structWithTrans.variableTransformations = {'none';'log10'};
structWithTrans.responseTransformation = 'log10';
domain = struct('name',{'p1';'p2'},'range',{[1 10];[5 6]});
linRM = ResponseModel(structWithTrans,domain);

% Evaluate model at p1 = 5, p2 = 6. Should get 10^(-4+0.5*5+log10(6)) = 0.189
[y yrange] = eval(linRM,[5; 6]);
if abs(y-10^(-4+0.5*5+log10(6))) > 1e-10 || yrange(1) ~= yrange(2)
  passed = false;
  if dispBool
    disp('   ERROR: computation incorrect')
  end
end

if dispBool
  disp('    Linear model, outputUncertainty [-1 5]')
end
coeffMatrix = [3 4 5]';
domain = struct('name',{'p1';'p2'},'range',{[-10 10];[-inf inf]});
linRM = ResponseModel(coeffMatrix,domain,[-1 5]);
[y yrange] = eval(linRM,[4; 7]);
if y ~= 54 || yrange(1) ~= y-1 || yrange(2) ~= y+5
  passed = false;
  if dispBool
    disp('   ERROR: computation incorrect')
  end
end

if dispBool
  disp('    Linear model, log10 transformed relative uncertainty 0.1')
end
linRM = set(linRM,'outputUncertainty',0.1,'outputUncertaintyType','relative','outputUncertaintyTransformation','log10');
[y yrange] = eval(linRM,[4; 7]);
if abs(y-54) > 1e-10 || abs(yrange(1)-10^(log10(y)*0.9)) > 1e-10 || abs(yrange(2)-10^(log10(y)*1.1)) > 1e-10
  passed = false;
  if dispBool
    disp('   ERROR: computation incorrect')
  end
end

if dispBool
  disp('    Quadratic model, log10 transformed absolute uncertainty [-0.2 0.1]')
end
coeffMatrix = diag([3 4 5]);
domain = struct('name',{'p1';'p2'},'range',{[-10 10];[-inf inf]});
unc.value = [-0.2 0.1];
unc.type = 'absolute';
unc.transformation = 'log10';
quadRM = ResponseModel(coeffMatrix,domain,unc);

% Evaluate model at p1 = 4, p2 = 7. Should get 3+4*4*4+5*7*7 = 242
[y yrange] = eval(quadRM,[4; 7]);
if abs(y-312) > 1e-10 || abs(yrange(1)-10^(log10(312)-0.2)) > 1e-10 || abs(yrange(2)-10^(log10(312)+0.1)) > 1e-10
  passed = false;
  if dispBool
    disp('   ERROR: computation incorrect')
  end
end

if dispBool
  disp('    Simple dcmodel, no output uncertainty')
end
simpleDCObj = ResponseModel(@simpleDCModel); %#ok
% Evaluate model at p1 = 2, p2 = 4. Should get 10^(2+4) = 10^6
[y yrange] = eval(simpleDCObj,[0.5; 4]);
if y ~= 10^4.5 || yrange(1) ~= yrange(2)
  passed = false;
  if dispBool
    disp('   ERROR: computation incorrect')
  end
end

if dispBool
  disp('    Simple dcmodel, supplying insufficient parameters')
end
try %#ok
  eval(simpleDCObj,2);
  passed = false;
  if dispBool
    disp('   ERROR: previous call should have bombed')
  end
end

if dispBool
  disp('    Simple dcmodel, x exceeding the model domain')
end
try %#ok
  eval(simpleDCObj,[2;4]);
  passed = false;
  if dispBool
    disp('   ERROR: previous call should have bombed')
  end
end

if dispBool
  disp('    Simple dcmodel, supplying extra parameters and paramList')
end
y = eval(simpleDCObj,[5;4;0.5],{'ptrash';'p2';'p1'});
if y ~= 10^4.5
  passed = false;
  if dispBool
    disp('   ERROR: computation incorrect')
  end
end

if dispBool
  disp('    Complex dcmodel, resp2')
end
cplxDCObj = ResponseModel(@complexDCModel,'resp2'); %#ok
[y yrange] = eval(cplxDCObj,[5;4;0.5],{'ptrash';'p2';'p1'});
if y ~= exp(4.5) || abs(yrange(1) - 0.9*exp(4.5)) > 1e-8 || abs(yrange(2) - 1.1*exp(4.5)) > 1e-8
  passed = false;
  if dispBool
    disp('   ERROR: computation incorrect')
  end
end

if dispBool
  disp('    Complex dcmodel, resp3')
end
cplxDCObj = ResponseModel(@complexDCModel,'resp3'); %#ok
[y yrange] = eval(cplxDCObj,[5;4;0.5],{'ptrash';'p2';'p1'});
if y ~= 10^4.5 || yrange(1) ~= 0.9*10^4.5 || yrange(2) ~= 1.1*10^4.5
  passed = false;
  if dispBool
    disp('   ERROR: computation incorrect')
  end
end
