function MM = uminus(MM)

% Special case:
if isempty(MM.model)
  return
end

for i1 = 1:length(MM)
  if strcmp(MM(i1).model.type,'svmd')
    MM(i1).model.coef = - MM(i1).model.coef;
    MM(i1).model.bias = - MM(i1).model.bias;
    MM(i1).model.quadCoeff = - MM(i1).model.quadCoeff;
    MM(i1).model.yScale(1) = - MM(i1).model.yScale(1);
    tmp = MM(i1).fitInfo.peakErrorOnSamplePoints4Fit;
    MM(i1).fitInfo.peakErrorOnSamplePoints4Fit = [-tmp(2) -tmp(1)];
    tmp = MM(i1).fitInfo.peakErrorOnSamplePoints4Validation;
    MM(i1).fitInfo.peakErrorOnSamplePoints4Validation= [-tmp(2) -tmp(1)];
    tmp = MM(i1).fitInfo.peakErrorFromOptimization;
    MM(i1).fitInfo.peakErrorFromOptimization = [-tmp(2) -tmp(1)];
    tmp = MM(i1).fitInfo.optimFailed;
    MM(i1).fitInfo.optimFailed = [tmp(2) tmp(1)];

    tmp2 = MM(i1).fitInfo.xGivingErrMin;
    MM(i1).fitInfo.xGivingErrMin = MM(i1).fitInfo.xGivingErrMax;
    MM(i1).fitInfo.xGivingErrMax = tmp2;
  else
    error('Code not complete')
    %TODO
  end
end
