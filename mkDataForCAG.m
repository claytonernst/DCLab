%% makedataforCAG class
% October 3, 2010

Ms = struct('Qform',[],'Pidx',[],'Interval',[]);

D = createGRI;
E = D;
E.ModelAndObservationPair(37) = [];

for i=1:76
   data = E.ModelAndObservationPair(i);
   Ms(i).Qform = data.ResponseModel.model;
   [a,b,c] = intersect(data.ResponseModel.parameterList,E.parameterList);
   [~,sbi] = sort(b);
   Ms(i).Pidx = c(sbi);
   Ms(i).Interval = data.observedValue + data.observationUncertaintyPlusMinus;
end
 
opts = DCOptions;
opts.maxBranchBoundIter = 1;
C = ConsistencyTest(E,opts);

X = C.iterStruct.lower.xfeas{1};

for i=1:76
   xi = [1;X(Ms(i).Pidx)];
   Val = xi'*Ms(i).Qform*xi;
   [i (Ms(i).Interval(1)<=Val&Val<=Ms(i).Interval(2)) Ms(i).Interval Val]
end

D = Ms';
Xfeas = X;

save CAGdata D Xfeas