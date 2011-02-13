function X = buildQuadX(normParamValues)
% given a matrix normParamValues = [x11 x21 ... xn1
%                                   x12 x22 ... xn2
%                                    :   :       :
%                                   x1N x2N ... xnN]
%
% this function creates a matrix suitable for creating a leastsquare
% quadratic approximation
%
% X = [1 x11 x21 ... xn1 | x11^2 x11*x21 ... x11*xn1 | x21^2 x21*x31 ... xn1^2
%      1 x12 x22 ... xn2 | x12^2 x12*x22 ... x12*xn2 | x22^2 x22*x32 ... xn2^2
%            :                     :                         :            :
%      1 x1N x2N ... xnN | x1N^2 x1N*x2N ... x1N*xnN | x2N^2 x2N*x3N ... xnN^2]
%
%
% normParamValues should have size [Nsamples, Nparams];
%
[N, n] = size(normParamValues);
Ncoeff = (n+1)*(n+2)/2;
X = zeros(N,Ncoeff);

%affine term
X(:,1) = ones(N,1);

%linear terms
X(:,[2:n+1]) = normParamValues;

%bilinear terms
placeHolder = 1;
for i1 = 1:n
  for i2 = i1:n
    X(:,n+1+placeHolder) = normParamValues(:,i1).*normParamValues(:,i2);
    placeHolder = placeHolder+1; 
  end
end
