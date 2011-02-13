function surf = uminus(surf)
%UMINUS method takes a surface created for M(x) and makes it a surface created for -M(x)
%
%   SURF = uminus(SURF) takes the DCSurface SURF created for a model
%   M(x) and converts it into a surface created for -M(x). If the surface
%   is fit in log10(y), i.e., surf.responseTransformation = 'log10', then
%   this function takes a surface created for log10(M(x)) and converts it
%   into a surface for -log10(M(x)). 
%
%   SURF = -SURF
%
%   See DCSurface

%   This function is used so we can maximize functions (in
%   ResponsePrediction) using code that minimizes.
%


% Last modified 8/15/07


if isstruct(surf.surrogateModel)
  surf.surrogateModel.num = - surf.surrogateModel.num;
else
  surf.surrogateModel = - surf.surrogateModel;
end
surf.surrogateFitInfo.peakErrorOnSamplePoints4Fit = -surf.surrogateFitInfo.peakErrorOnSamplePoints4Fit([2 1]);
surf.surrogateFitInfo.peakErrorOnSamplePoints4Validation = -surf.surrogateFitInfo.peakErrorOnSamplePoints4Validation([2 1]);
surf.surrogateFitInfo.peakErrorFromOptimization = -surf.surrogateFitInfo.peakErrorFromOptimization([2 1]);
tmp = surf.surrogateFitInfo.xInModelOrderGivingErrMax;
surf.surrogateFitInfo.xInModelOrderGivingErrMax = surf.surrogateFitInfo.xInModelOrderGivingErrMin;
surf.surrogateFitInfo.xInModelOrderGivingErrMin = tmp;
surf.trueModelRangeEstimate = -surf.trueModelRangeEstimate([2 1]);

%Make sure we're not associating this flipped surface with anything
%concrete.
surf.MOPairIdx = [];
