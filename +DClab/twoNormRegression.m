function [coeff,fitError,EXITFLAG] = ...
    twoNormRegression(x,pred,displayLevel)
% function to fit a quadratic to function evaluations in two norm.
%
%  min_{theta,gamma} gamma
%  s.t. |A*theta - b| <= gamma
%
%  the inputs are pred = b, A = x
%
%  we use this to fit a response surface to input output data.
%  

[n, m] = size(x);
[o, p] = size(pred);

if n ~=o || p ~= 1
  error('you have dimension issues')
end

coeff = x\pred;
fitError = x*coeff - pred;
EXITFLAG = [];
