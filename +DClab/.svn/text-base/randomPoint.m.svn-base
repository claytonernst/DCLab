function x = randomPoint(PD)

%only works on non partitioned cubes and for noLog surface
tempPD = PD;
root = PD.cube.root;
bnds = PD.cube.partitionStruct(root).bnds;
%opt = DClab.DCOptions('omitInnerBound',1);
opt = DClab.DCOptions('nRestart',0,'doNotBB',1);

n = size(bnds,1);
x = zeros(n,1);
pNames = PD.cube.paramList;

for i1 = 1:n
  %find feasible parameter bounds
  model = [0 0.5;0.5 0];
  makeQuadModel(model,pNames(i1));
  
  MA = ModelAssertion('tempModel');
  
  [lb ub] = prediction(MA,tempPD,[],opt)
  
  if lb(1,2) > -1+1e-4 | ub(1,1) < 1-1e-4
    ff = tempPD.cube.fitterData;
    ff.options.nRestart = 4;
    tempPD.cube.fitterData = ff;
    [lb2 ub2] = prediction(MA,tempPD,[],opt)
    LB(i1,[1 2]) = min(lb,lb2);
    UB(i1,[1 2]) = max(ub,ub2);
  else 
    LB(i1,[1 2]) = lb;
    UB(i1,[1 2]) = ub;
  end

  %rm(fileNameFull);
  clear tempModel
  
  if LB(i1,1) > UB(i1,2)
      keyboard
  end
  
%  bnds(i1,:) = [LB(i1,1) UB(i1,2)];
  bnds(i1,:) = [LB(i1,2) UB(i1,1)];

  x(i1) = bnds(i1,1) + rand*diff(bnds(i1,:));
  
  for i2 = 1:length(PD)
    if ismember(i1,tempPD.cube.partitionStruct(root).surf(i2).surfInToCubeIdx)
        matrix = tempPD.cube.partitionStruct(root).surf(i2).noLog.Surf;
        matrix(1,1) = matrix(1,1) + 2*matrix(1,2)*x(i1) + matrix(2,2)*x(i1)^2;
        matrix(1,3:end) = matrix(1,3:end) + matrix(2,3:end)*x(i1);
        matrix(3:end,1) = matrix(3:end,1) + matrix(3:end,2)*x(i1);
        
        %remove 2nd column and row
        matrix(:,2) = [];
        matrix(2,:) = [];
        
        %update surface
        ss = tempPD.cube.partitionStruct(root).surf(i2);
        ss.noLog.Surf = matrix;
        ss.surfInToCubeIdx(1) = [];
        ss.fitBnds(1,:) = [];
        
        tempPD.cube.partitionStruct(root).surf(i2) = ss;

    end
  end
  i1
end
 
keyboard

  
function makeQuadModel(model,pNames)
  
fileName = 'tempModel';
fileNameFull = [pwd filesep 'tempModel' '.m'];

%open file
fid = fopen(fileNameFull,'w+');

%begin writing...
fprintf(fid,['function response=%s(executionType,paramValues,' ...
             'flag)\n\n'],fileName);
fprintf(fid,'ni = nargin;\nno = nargout;\n');
fprintf(fid,['error(nargchk(1,3,ni));\nerror(nargoutchk(0,1,no));\n\' ...
             'n']);

fprintf(fid,'if ~ischar(executionType)\n');
fprintf(fid,['  error(''First input to model simulation file must be ' ...
             'of type char.'')\n']);
fprintf(fid,'end\n\n');

fprintf(fid,'if ni == 2\n  flag = 0; %%don''t plot\nend\n\n');

fprintf(fid,'Q = [');
for i2=1:size(model,1)
  for i3=1:size(model,2)
    fprintf(fid,'%d, ',model(i2,i3));
  end
  fprintf(fid,'; ...\n     ');
end
fprintf(fid,'     ];\n\n');

fprintf(fid,'switch executionType\n');
fprintf(fid,'  case ''description''\n');
fprintf(fid,'   response = ''%s'';\n','model for testing');
%fprintf(fid,'  case ''observable''\n');
%fprintf(fid,'   response = ''%s'';\n',dataStruct(i1).type);
fprintf(fid,'  case ''featureList''\n');
fprintf(fid,'   response = {''NA''};\n');
fprintf(fid,'  case ''featurePrecedence''\n');
fprintf(fid,'   response = 1;\n');
fprintf(fid,'  case ''isquadratic''\n');
fprintf(fid,'   response = 1;\n');
fprintf(fid,'  case ''islinear''\n');
fprintf(fid,'   response = 0;\n');
fprintf(fid,'  case ''paramList''\n');
fprintf(fid,'   response={''%s''',pNames{1});
for i2=2:length(pNames)
  fprintf(fid,'; ''%s''',pNames{i2});
end
fprintf(fid,'};\n\n');

fprintf(fid,'  case ''getPoly''\n');
fprintf(fid,'   response = Q;\n');

fprintf(fid,'  case ''simulate''\n\n');
fprintf(fid,'    N = size(paramValues,1); \n');
fprintf(fid,'    x = [ones(N,1) paramValues];\n\n');

fprintf(fid,'    response = sum(x*Q.*x,2);\n\n');

fprintf(fid,'  otherwise\n');
fprintf(fid,['    error(''Incorrect input type to model simulation ' ...
             'file'')\n']);
fprintf(fid,'end');

% close file
fclose(fid);
