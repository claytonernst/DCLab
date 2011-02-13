function [fval xfeas] = maximizeLB_meta(MMCell,d,u,uncCase,FP,domRng)

if ~all(uncCase==1)
  error('Code not complete')
end

%Our objective is to solve
% - min\gamma s.t. el(1+\gamma) \leq S(x) - d \leq u(1+\gamma)
%
% Should be a piece of cake.

% Step one. Find all active variables.

m = length(MMCell);
n = nParameters(FP);

modParam = cell(m,1);
actModParam = cell(m,1);
varTrans = cell(m,1);
for i1 = 1:m
  modParam{i1} = MMCell{i1}.parameterList;
  varTrans{i1} = MMCell{i1}.variableTransformations;
  actModParam{i1} = modParam{i1}(varTrans{i1}~=0);
end

active = unique(vertcat(actModParam{:}));
activeIdx = findIdx(FP,active); % index to parameters active in the metamodels.
actBool = false(n,1);
actBool(activeIdx) = true;
doesVaryBool = diff(domRng,[],2) > 1e-10;
activeIdx = find(actBool & doesVaryBool);

activeFP = FP(activeIdx);

active2model = cell(m,1);
dset2model = cell(m,1);
for i1 = 1:m
  % this will place NaNs in locations where the corresponding parameter is
  % not active.
  active2model{i1} = findIdx(activeFP,actModParam{i1}); %#ok
  dset2model{i1} = findIdx(FP,modParam{i1});
end

actUncTrans = vertcat(activeFP.uncertaintyTransformation);
actBnds = domRng(activeIdx,:);

nom = vertcat(FP.nominal); %#ok

% create array of handles that we can call with the normalized optimization
% variable and they return the metamodel output.
MMhandles = cell(m,1);

for i1 = 1:m
  % assume x is only a column vector. i don't believe fmincon takes
  % advantage of vectorizing.
  
  %create string to build anon function on the fly.

  nP = length(varTrans{i1});
  tmp = cell(nP,1);
  
  act2mod = active2model{i1};
  dset2mod = dset2model{i1}; %#ok
  for i2 = 1:nP
    if varTrans{i1}(i2) == 0
      tmp{i2} = ['nom(dset2mod(' num2str(i2) '));'];
    else
      % need to undo transformation into normalized variables.
      if strcmp(actUncTrans{act2mod(i2)},'log10')
        tmp{i2} = ['actBnds(act2mod(' num2str(i2) '),1)*(actBnds(act2mod(' num2str(i2) '),2)/actBnds(act2mod(' num2str(i2) '),1))^(0.5*(x(act2mod(' num2str(i2) '))+1));'];
      else
        tmp{i2} = ['0.5*(x(act2mod(' num2str(i2) '))+1)*diff(actBnds(act2mod(' num2str(i2) '),:))+actBnds(act2mod(' num2str(i2) '),1);'];
      end
    end
  end
  
  eval(['hand = @(x) eval(MMCell{' num2str(i1) '},[' horzcat(tmp{:}) ']);']);
  MMhandles{i1} = hand;
end

confun = @(x) fmconfun(x,MMhandles,d,u);
objfun = @(x) x(end);

opt.tolFun = 1e-5;
opt.tolCon = 1e-5;

fmoptions = optimset('fmincon');
fmoptions = optimset(fmoptions,'GradConstr','off','GradObj','off','TolX',opt.tolCon,...
   'LargeScale','off','MaxIter',500,'MaxFunEval',12500, ...
   'TolFun',opt.tolFun,'MaxSQPIter',10000,'Display','off','TolCon',opt.tolFun);

nact = length(activeIdx);
xinit = 2*rand(nact+1,1)-1; 

[xopt fval eflg] = fmincon(objfun,xinit,[],[],[],[],[-ones(nact,1);-inf],[ones(nact,1);inf],confun,fmoptions);

if eflg <= 0
  disp('nonpositive eflg')
  %keyboard
end

% convert xopt back to normal cordinates.

xfeas = nom;
for i1 = 1:length(activeIdx)
  if strcmp(actUncTrans(i1),'log10')
    xfeas(activeIdx(i1)) = actBnds(i1,1)*(actBnds(i1,2)/actBnds(i1,1))^(0.5*(xopt(i1)+1));
  else
    xfeas(activeIdx(i1)) = 0.5*(xopt(i1)+1)*diff(actBnds(i1,:))+actBnds(i1,1);
  end
end

fval = -fval;

% local function to implement the constraints.

% TODO we don't provide gradients. however, we should be able to.
function [iconval,econval] = fmconfun(x,MMhandles,d,u)

% gamma is the last variable;
m = length(MMhandles);
iconval = zeros(2*m,1);
econval = 0;
y = zeros(m,1);
xx = x(1:end-1);
for i1 = 1:m
  y(i1) = MMhandles{i1}(xx);

  % constraints are
  % -M(x) + d + l <= -l*\gamma
  % M(x) - d - u <= u*gamma;
  
  iconval(2*i1-1) = -y(i1) + d(i1) + u(i1,1)*(1+x(end));
  iconval(2*i1) = y(i1) - d(i1) - u(i1,2)*(1+x(end));

end


















