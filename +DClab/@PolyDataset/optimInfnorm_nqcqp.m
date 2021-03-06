function [linXlinYnqcqp,logXlinYnqcqp,activeIdxCell] = optimInfnorm_nqcqp(PD,opt,trgWeights,approxType,nodes)
% function [linXnqcqp,logXnqcqp,activeIdxCell] = optimInfnorm_nqcqp(PD,trgWeights,approxType,nodes)
%
% This function creates nqcqp optimization objects in x and
% log10(x), depending on the surrogate models that are
% available. Each object is a column, of length corresponding to
% length nodes, i.e., a different optimization is formulated for
% each subdivision over which we're posing the best fit
% optimization. More info can be found in PolyDataset/infNormFitX
%
%  Inputs:
%    PD: a PolyDataset object 
%    opt: DCOptions object
%    trgWeights A mx1 vector of weights w_e in the above
%      optimization. Default is usually 1/u_e (see below).
%    approxType: inner, crude, or outer: indicates how fitting
%      errors are treated
%    nodes: a column vector specifying node numbers of the PolyDataset
%      object to employ in the nqcqp objects.
%  Outputs:
%    linXnqcqp: a rx1 nqcqp object using surrogate models that 
%      are quadratic in the parameters. r = length(nodes).
%    logXnqcqp: a rx1 nqcqp object using surrogate models that
%      are quadratic in the log10 of parameters.
%    activeIdxCell: a rx1 cell array. Each cell contains the
%      indicies into PD.ParamAssns that are "active" in the
%      corresponding optimization. We pass this back because we
%      need it later, and it is conveniently available here.

% TODO currently logXnqcqp is empty, and nodes is set to 1.

% code assumes activeIdx is the same for linX and logX constraints

nargchk(nargin,5,5);

% Determine which surrogate models are available. At most it should
% be linXlinY and logXlinY. The logY case is handled elsewhere.

t = {'linXlinY'};
m = nUnits(PD);
r = length(nodes);
QlinXlinY = cell(r,1);
QlogXlinY = cell(r,1);
optimBnds = cell(r,1);
objFctn = cell(r,1);
activeIdxCell = cell(r,1);

for i0 = 1:r
  
  activeIdx = [];
  for i1 = 1:m
     activeIdx = union(activeIdx,PD.surrogateTree(nodes(i0)).surf(i1).activeParams);
  end
  activeIdx = activeIdx(:); %force to column
  n = length(activeIdx);
  const.linXlinY = {}; const.logXlinY = {}; %cell arrays to hold
                                            %the constraint matrices

  for i1 = 1:m
    surf = PD.surrogateTree(nodes(i0)).surf(i1);
   
    d = PD.ModelAndObservationPair(i1).observedValue;
    u = PD.ModelAndObservationPair(i1).observationUncertaintyPlusMinus;    
    %center data such that u is symmetric
    if -u(1) ~= u(2)
      tmp = mean(u);
      d = d-tmp;
    end

    mm = length(surf.activeParams);
    %build constraints

    for i2 = 1:length(t)

      [trash gIdx] = intersect(activeIdx,surf.activeParams); %#ok
      str1 = [t{i2} 'PeakErr'];
      str2 = [t{i2} 'Poly'];
      fitErr = surf.(str1);
      poly = surf.(str2);

      %if approxType = outer, we want to implement the constraints
      % -gammai + wi*fitErr(1) <= wi*(S(x) - d) <= wi*fitErr(2) + gammai
      %if approxType = crude we want to implement
      % -gammai <= wi*(S(x) - d) <= + gammai
      %if approxType = inner we want to implement
      % -gammai + wi*fitErr(2) <= wi*(S(x) - d) <= wi*fitErr(1) + gammai

      switch approxType
       case 'outer'
        fitu = fitErr;
       case 'crude'
        fitu = [0 0]; 
       case 'inner'
        fitu = [fitErr(2) fitErr(1)];
      end
      %allocate for (n+m+1)x(n+m+1) matrices to include gammas as the last variables      
      temp = spalloc(1+n+1,1+n+1,(mm+1)^2+2); 

      %create right constraint
      if isstruct(poly)
        error('code not complete')
        %temp([1;gIdx+1],[1;gIdx+1]) = poly.num - (d+fitu+u(1,2))*poly.den;
      else
        temp([1;gIdx+1],[1;gIdx+1]) = trgWeights(i1)*poly;
        temp(1,1) = temp(1,1) - trgWeights(i1)*(d + fitu(2));
      end
      temp(1,1+n+1) = -0.5; 
      temp(1+n+1,1) = -0.5;
      if isempty(const.(t{i2}))
        const.(t{i2}){1} = temp;
      else
        const.(t{i2}){end+1,1} = temp;
      end
      
      %create left constraint
      temp = spalloc(1+n+1,1+n+1,(mm+1)^2+2);
      if isstruct(poly)
        error('code not complete')
        %temp([1;gIdx+1],[1;gIdx+1]) = (d + fitu + u(1,1))*poly.den - poly.num;
      else
        temp([1;gIdx+1],[1;gIdx+1]) = -trgWeights(i1)*poly;
        temp(1,1) = temp(1,1) + trgWeights(i1)*(d + fitu(1));
      end
      temp(1,1+n+1) = -0.5; 
      temp(1+n+1,1) = -0.5;
      const.(t{i2}){end+1,1} = temp;
    end
  end

  QlinXlinY{i0,1} = const.linXlinY;
  QlogXlinY{i0,1} = const.logXlinY;
  
  %constrain bounds on x's and gammas. since everything is normalized, this
  %is trivial.
  optimBnds{i0,1} = [-ones(n,1) ones(n,1); [-inf inf]];
  
  %build objective function consisting of gamma (the last variables)
  Z0 = zeros(n+2);
  Z0(1,n+2) = 0.5;
  Z0(n+2,1) = 0.5;
  
  objFctn{i0,1} = Z0;
  activeIdxCell{i0,1} = activeIdx;  
  
end

linXlinYnqcqp = cell(r,1);
logXlinYnqcqp = cell(r,1);
for i1=1:r
    % Call the nqcqp constructor
    if opt.constraints(1)
        linXlinYnqcqp = DClab.nqcqp(objFctn{i1},QlinXlinY{i1},optimBnds{i1}(:,1),optimBnds{i1}(:,2));
    else
        linXlinYnqcqp = DClab.nqcqp;
    end

    if opt.constraints(2)
        logXlinYnqcqp = DClab.nqcqp(objFctn{i1},QlogXlinY{i1},optimBnds{i1}(:,1),optimBnds{i1}(:,2));
    else
        logXlinYnqcqp = DClab.nqcqp;
    end
end
