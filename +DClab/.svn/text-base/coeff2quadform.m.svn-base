function quadform = coeff2quadform(coeff,Nactive)
%COEFF2QUADFORM takes the coefficients of a quadratic polynomial and returns
%a symmetric matrix.
%
%  QUADFORM = COEFF2QUADFORM(COEFF,NACTIVE) does XXX let q(x1,...,xn) be
%  defined by
%  q(x1,...,xn) = c00 + c10*x1 + ... + cn0*xn + c11*x1^2 + c12*x1*x2 + ... +
%  c1n*x1*xn + c22*x2^2 + c23*x2*x3 + ... + c2n*x2*xn + ... + cnn*xn^2.
%
%  Then COEFF should be [c00 c10 ... cn0 c11 c12 ... c1n c22 c23 ... c2n ... cnn]'
%
%  


%force coeff to row vector
coeff = coeff(:)';

if Nactive >= 2
  Ncoeff = 1 + 2*Nactive + nchoosek(Nactive,2);
else
  Ncoeff = 1 + 2*Nactive;
end

if length(coeff) ~= Ncoeff
  error('Number of coefficients mismatches number of parameters')
end

matSize = Nactive+1;

quadform = zeros(matSize);
for i1 = 1:Nactive+1
  startidx = (i1-1)*matSize - sum([0:i1-1]) + i1;
  endidx = i1*matSize - sum([0:i1-1]);
  quadform(i1,[i1:matSize]) = coeff(startidx:endidx);
end
quadform = 0.5*(quadform + quadform');
