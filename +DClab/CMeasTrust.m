function [iter, subiter] = CMeasTrust(Dset,opt,absBnds)
% TRUSTREGION
%
%   [iter subiter] = trustRegion(Dset,opt,stopCriteria,absBnds)
%
%   This is the trust region algorithm (or "magnifying glass"
%   algorithm).  By looking at a subset of the feasible set that is
%   small enough for the surrogate fits to be sufficiently "good",
%   the optimization is performed locally.  Based on the result of
%   that optimization, a new subset (trust region) is chosen, and
%   the optimization is repeated.
%
%   Inputs:
%     Dset: This is a DCDataset object.
%     opt: Options object
%     absBnds: an optional nx2 vector to represent hard bounds on the
%           feasible space (n is the number of parameters listed in
%           Dset.ParamAssn).  If these are given, the bounds in
%           .ParamAssn are used for scaling the problem, and the
%           initial griding size, but the algorithm will search
%           over the space defined by absBnds for the optimal
%           point.  This is designed for when the uncertainty bnds
%           are uncertain, or fuzzy.
%
%   Outputs:
%     iter: a structure containing information about each iteration
%           of the algorithm which successfully caused the trust
%           region to move.  It has the following fields:
%           bestPt: the current best point in the parameter space
%           CmeasReal: the consistency measure cost function
%               evaluated at bestPt
%           CmeasSurrLow: a lower bound on the consistency measure
%               as calculated using the surrogate fits
%           CmeasSurrUp: an upper bound on the consistency measure
%               as calculated using the surrogate fits
%           PDset: the PolyDataset object constructed at this
%               iteration
%     subiter: this is a structure containing information about
%           each iteration of the algorithm regardless of whether
%           it was succesful or not.  If an iteration was
%           unsuccessful the trust region is shrunk, and the
%           optimization repeated.  This structure has all the
%           fields that iter has (with xFeas in place of bestPt)
%           plus the following fields:
%           G: The current shrinking or growth factor of the trust
%               region
%           wasGood: this field is 1 if the bounds on the
%               consistency measure improved.  If so the
%               information of this index should be in iter

ni = nargin;

%extract the needed parameters from the parameter object
%[ParamObj,idx] = extract(ParamObj, getParamList(Dset));

switch ni
 case 1
  opt = DClab.DCOptions;
  absBnds = [];
 case 2 
  absBnds = [];
 otherwise
  error(nargchk(1,3,ni))
end

if isempty(opt)
  opt = DClab.DCOptions;
end

%Overwrite some of the fitting options
opt.nPnts4ValidationM=1;
opt.nLocalValidationSearches=1;
opt.maxFitIter=4;
opt.fitType=2;

%number of parameters
%n = length(Dset.ParamAssn);

%define the absolute bounds
if isempty(absBnds)
  absBnds = vertcat(Dset.FreeParameter.range);
end

LB = absBnds(:,1);
UB = absBnds(:,2);

if strcmp(opt.diagnostics,'on')
  if isempty(opt.diagName)
    dbSaveName = 'magnifyingGlassDiagnostics';
  else
    dbSaveName = opt.diagName;
  end
else
  dbSaveName = '';
end

%======================================================
% Section to determine the initial grid size. This exists because we make
% more dratic change to the cube size here that we do in the trust-region
% iterations, and also, because we'd like to initially fit reasonably well
% to preserve the "local-search" nature of the algorithm

if any(strmatch(opt.display,{'iter';'all';'ALL'}))
  disp([blanks(2) 'Finding initial grid size...']);
end

notDone = true;
i1 = 1;

perIterScaleFactor = 0.5;
Units = Dset.Unit;
tmpDset = Dset;

while notDone
  if i1 == 1
    iter(1).cubeSize = 1;
    if any(strmatch(opt.display,{'iter';'all';'ALL'}))
      disp([blanks(4) 'Fitting over initial domain'])
    end
  else
    iter(1).cubeSize = perIterScaleFactor*iter(1).cubeSize;    
    if any(strmatch(opt.display,{'iter';'all';'ALL'}))
      disp([blanks(4) 'Shrinking ' num2str(iter(1).cubeSize) ' of domain'])
    end
    Dset.ParamAssn = scaleUnc(Dset.ParamAssn,perIterScaleFactor);
  end
  
  PDset = DClab.PolyDataset;
  
  passed = zeros(nUnits(Dset),1);
  
  for i2=1:nUnits(Dset)
    tmpDset.Unit = Units(i2);
    
    PDset = merge(PDset, DClab.PolyDataset(tmpDset,opt));
    initializingPDset = PDset; %#ok this is here just for saving
    if dbSaveName
      save(dbSaveName);
    end
    tmp = struct(PDset.surrogateTree(1).surf(i2));

    if strmatch('linXlinY',tmp.trans)
      flag1 = tmp.fitInfo.linXlinY.maxPntsFlag;
    else
      flag1 = 0;      
    end
    if strmatch('logXlinY',tmp.trans)
      flag2 = tmp.fitInfo.logXlinY.maxPntsFlag;
    else
      flag2 = 0;
    end
    if strmatch('linXlogY',tmp.trans)
      flag3 = tmp.fitInfo.linXlogY.maxPntsFlag;
    else
      flag3 = 0;
    end
    if strmatch('logXlogY',tmp.trans)
      flag4 = tmp.fitInfo.logXlogY.maxPntsFlag;
    else
      flag4 = 0;
    end
    
    if any([flag1 flag2 flag3 flag4])
      if any(strmatch(opt.display,{'iter';'all';'ALL'}))
        disp([blanks(4) 'Fit failed on surface ' num2str(i2) ' of ' num2str(nUnits(Dset))])
      end
      break %break the for loop prematurely and decrease grid size
    else
      passed(i2) = true;
    end
  end %i2=1:nUnits(Dset)
  
  %if we didn't break the for loop, it must have worked
  if passed(end)
    notDone = false;
  else
    i1 = i1+1;
  end
  
end

opt.maxFitIter = 3;

if any(strmatch(opt.display,{'iter';'all';'ALL'}))
  disp([blanks(2) 'Initial grid size found: ' num2str(iter(1).cubeSize) ' of original'])
end

initialPDset = PDset; %#ok this is defined just for saving
if dbSaveName
  save(dbSaveName);
end

%=====================================================
%section to perform the trust-region iterations

%the structure array iter will track outer loop iterations.
%stuff from the inner trials will get packed in subiter

%xFeasTry, 
%CmeasSurrUpTry, 
%CmeasSurrLowTry,
%cubeSizeTry, and
%G

%set starting cube point and evaluate at point
iter(1).bestPnt = vertcat(Dset.FreeParameter.nominal);
%iter(1).CmeasReal = ... 
%  findRealCmeasLB(PDataset,get(PDset.ParamAssn,'nominal_value')); 
iter(1).CmeasSurrLow = -inf;
iter(1).CmeasSurrUp = inf;
iter(1).PDset = PDset;
% iter(1).cubeSize is already defined;

%begin consitency check loop
ii = 1;
subiter(ii).G = 1;
keepLooking = true;
numShrinks = 0;
numMoves = 0;

% if strcmp(opt.stopGui,'on')
%   figH = stopGui;
%   drawnow;
% end

while keepLooking
    
  if ii>1
    %update cubeSize
    %maybe could adjust cube so not beyond bounds... later

%    %Make fit
%    if strcmp(opt.display,'on')
%      disp('=======Making Fit....======');
%    end

    if length(iter) >= 2
      Dset.ParamAssn = shiftNominal(Dset.ParamAssn,mean([iter(end).bestPnt iter(end-1).bestPnt],2));
    else
      %we're still at the initial point
      Dset.ParamAssn = shiftNominal(Dset.ParamAssn,iter(1).bestPnt);
    end

    Dset.ParamAssn = scaleUnc(Dset.ParamAssn,subiter(ii).G);

    %fix any bounds that are out of the absolute bnd domain
    range = vertcat(Dset.FreeParameter.range);
    lb = range(1);
    ub = range(2);
    tidx = find(lb < LB);
    lb(tidx) = LB(tidx);
    tidx = find(ub < LB);
    ub(tidx) = LB(tidx);

    tidx = find(ub > UB);
    ub(tidx) = UB(tidx);
    tidx = find(lb > UB);
    lb(tidx) = UB(tidx);

    for i1=1:length(lb)
        Dset.FreeParameter(i1).range = [lb(i1) ub(i1)];
    end

    subiter(ii).PDset = DClab.PolyDataset(Dset,opt);
  else
    subiter(ii).PDset = iter(1).PDset;
  end

%  %check consistency
%  if strcmp(opt.display,'on')
%    disp('=====Checking Consistency=====');
%  end

  try
    s1 = ilbCons2(subiter(ii).PDset,1,iter(end).bestPnt);
  catch
    disp('failed call to ilbCons1 in TRUSTREGION')
    keyboard
  end

  subiter(ii).CmeasSurrUp = -(s1.pound.bnd+s1.gap);
  subiter(ii).CmeasSurrLow = -(s1.pound.bnd);
  subiter(ii).xFeas = s1.pound.xfeas{1};

  if any(strmatch(opt.display,{'iter';'all';'ALL'}))
    disp([blanks(6) 'Trial lower bound: ' num2str(subiter(ii).CmeasSurrLow)]);
  end
    
  % At this point three things may happen
  % a) the lower bound has improved, in which case we save the current
  %    subiteration at this iteration, display stuff, and grow the cube
  % b) we've stopped improving. this happens when 3 subiteration
  %    tries have not increased the lower bound by more then tolFun of the
  %    previous full iteration lower bound. the 3 is arbitrary.
  % c) the lower bound has not improved, in which case we shrink the cube
  %    and try again

  if subiter(ii).CmeasSurrLow > iter(end).CmeasSurrLow
    %we sucessfully improved the lower bound, save last trial results in iter
    subiter(ii).wasGood = 1;
    Niter = length(iter)+1;
    iter(Niter).bestPnt = subiter(ii).xFeas;
    iter(Niter).PDset = subiter(ii).PDset;
    iter(Niter).CmeasSurrLow = subiter(ii).CmeasSurrLow;
    iter(Niter).CmeasSurrUp = subiter(ii).CmeasSurrUp;

    if any(strmatch(opt.display,{'iter';'all';'ALL'}))
      disp([blanks(4) 'Successful iteration!']);
      disp([blanks(6) 'Current lower bound: ' num2str(iter(Niter).CmeasSurrLow)]);
    end

%   if any(strmatch(opt.display,{'iter';'all';'ALL'}))
%     disp('====Growing Cube size====');
%   end
    numMoves = numMoves+1;

    %increase cube volume by 25 percent
    subiter(ii+1).G = 1.05;
    %subiter(ii+1).G = 1.25^(1/n);
    %G(ii+1,1) = abs(1.05*G(ii));

  elseif length(subiter) >= 3 && min([subiter(ii-2:ii).CmeasSurrLow]) >= iter(end).CmeasSurrLow - opt.tolFun^2
    % We're stuck. In the last three subiteration, we never did much worse
    % than the previous best bound, so the fits are probably good, etc.
    if any(strmatch(opt.display,{'final';'iter';'all';'ALL'}))
       disp([blanks(2) 'Optimization converged'])
       disp([blanks(4) 'Lower bound: ' num2str(subiter(ii).CmeasSurrLow)]);
    end
    
    Niter = length(iter)+1;
    iter(Niter).bestPnt = subiter(ii).xFeas;
    iter(Niter).PDset = subiter(ii).PDset;
    iter(Niter).CmeasSurrLow = subiter(ii).CmeasSurrLow;
    iter(Niter).CmeasSurrUp = subiter(ii).CmeasSurrUp;
    
  else
    %We failed to improve the lower bound---shrink the cube size
    if any(strmatch(opt.display,{'iter';'all';'ALL'}))
      disp([blanks(6) 'Unsuccessful try: shrinking cube size']);
    end

    numShrinks = numShrinks+1;
    %shrink cube volume by 70 percent
    %subiter(ii+1).G = 0.50^(1/n);
    subiter(ii+1).G = 0.8;        

    if dbSaveName
      save(dbSaveName);
    end
  end

  if iter(end).CmeasSurrLow <= iter(end-1).CmeasSurrLow + opt.tolFun^2;
    if any(strmatch(opt.display,{'final';'iter';'all';'ALL'}))
       disp([blanks(2) 'Optimization converged'])
       disp([blanks(4) 'Lower bound: ' num2str(subiter(ii).CmeasSurrLow)]);
    end
    keepLooking = false;
  else
    ii = ii+1;
  end
end
