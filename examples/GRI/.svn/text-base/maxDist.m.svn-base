% This file obtains an outer bound on the maximum distance between two
% points in the feasible set of the GRI dataset. Note that the modified
% values are used for the two Stanford targets.

% Note, it would be rather easy to modify createGRI to realize the 'old'
% version, where H was [-1 1]. Trent knows how...

% Note, this file takes around a half hour to run. 

D = createGRI;
n = getPrivate(D,'nParameters');
m = getPrivate(D,'nPairs');

% Make an additional copy of each parameter to realize a 204 dim 
% parameter domain
H = D.FreeParameter;
H2 = H;
for i1 = 1:102
  H2(i1).name = [H2(i1).name 'copy2'];
end

Hall = vertcat(H,H2);

% Make an additional copy of each MOPair to realize a 154 total.
% Change the parameter names and pair names to reflect copy2 variables.
Pairs = D.ModelAndObservationPair;
Pairs2 = Pairs;
for i1 = 1:m
  oldDom = Pairs(i1).ResponseModel.domain;
  newDom = oldDom;
  for i2 = 1:length(oldDom)
    newDom(i2).name = [oldDom(i2).name 'copy2'];
  end
  Pairs2(i1).ResponseModel.domain = newDom;
  Pairs2(i1).name = [Pairs2(i1).name 'copy2'];
end

PairsAll = vertcat(Pairs,Pairs2);

% Create the joint dataset
DAll = Dataset(PairsAll,Hall);

%create euclidean norm objective

% First its domain.

modelDom = cell(2*n,2);
modelDom(:,1) = Hall.name;
modelDom(:,2) = mat2cell(Hall.range,ones(2*n,1),2);
modelDom = cell2struct(modelDom,{'name','range'},2);

%create objective matrix
%Nominally, the objective is sum_{i=1:102} (xi-xicopy2)^2.
%However, we what to make it consistent with the previous [-1 1]^n
%definition of H.
%
% Let rho_i be a dataset variable and let lrho_i = log10(lrho_i). To make
% the log10-transformed version between +/- 1, 
% x_i = 2*(lrho_i)-log10(lb))/(log10(ub)-log10(lb))-1;

% The objective function in terms of x_i:
objMat = speye(2*n);
for i1 = 1:n;
  objMat(i1,i1+n) = -1;
  objMat(i1+n,i1) = -1;
end
objMat = blkdiag(0,objMat);

% To make this a quadratic function of lrho_i, we need to compose it with
% the affine function
%x_i = [2/(log10(ub)-log10(lb))]*lrho_i + {-[2*log10(lb)/(log10(ub)-log10(lb))] - 1}
bnds = Hall.range;
A = zeros(2*n,1);
B = zeros(2*n,1);
for i1 = 1:2*n;
  A(i1) = 2/log10(bnds(i1,2)/bnds(i1,1));
  B(i1) = -2*log10(bnds(i1,1))/(log10(bnds(i1,2)/bnds(i1,1))) - 1;
end

% We know from the definition of the GRI parameter ranges, B should be 0.
% Let's enforce this to avoid any numerical troubles...
B = zeros(size(B));
A = diag(A);

% This creates a quadratic form such that Q(log10(rho)) = the euclidean
% distance between two points in the old H = [-1 1]^n
newObjMat = composeQuadWithAffine(objMat,A,B);

quadform.value = newObjMat;
quadform.responseTransformation = 'none';
quadform.variableTransformations = repmat({'log10'},2*n,1);

% Call the constructor
EuclObj = ResponseModel(quadform,modelDom);
EuclObj.name = 'Euclidean objective in old H vars';

opt = DCOptions('display','iter','nRestart',0,'maxBranchBoundIter',1);
obj = ResponsePrediction(EuclObj,DAll,opt);

distLowerBnd = sqrt(obj.UBi);
distUpperBnd = sqrt(obj.UBo);

xopt = obj.UBx;
x1 = xopt(1:n);
x2 = xopt(n+1:end);

d = D.ModelAndObservationPair.observedValue;

y = evalResponseModels(D,x1);
if max(abs(log10(y)-log10(d))>0.1001)
  x1feas = false;
else
  x1feas = true;
end

y = evalResponseModels(D,x2);
if max(abs(log10(y)-log10(d))>0.1001)
  x2feas = false;
else
  x2feas = true;
end

