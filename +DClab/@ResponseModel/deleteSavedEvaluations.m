function deleteSavedEvaluations(RM)
%DELETESAVEDEVALUATIONS method deletes the files that store saved evalutions
%
%   DELETESAVEDEVALUATIONS(RM) deletes the files (if any exist) that
%   store the saved evalutaions for the ResponseModel RM. These files
%   are stored in a subdirectory of DClabVxxx/savedModelEvalutaions/

if strcmp(RM.type,'dcModel')

  name = func2str(RM.model);
  for i1 = 1:length(RM.additionalInputs)
    if ischar(RM.additionalInputs{i1})
      name = [name '_' RM.additionalInputs{i1}]; %#ok
    elseif isnumeric(RM.additionalInputs{i1})
      name = [name '_' num2str(RM.additionalInputs{i1})]; %#ok
    else
      error('Internal inconsistency: condition should never occur')
    end
  end

  name = strrep(name,'.','p');
  
  file = which('savedEvaluationsDir');
  evalPath = fileparts(file);
  dirName = fullfile(evalPath,name);
  if exist(dirName,'dir')
    rmdir(dirName,'s');
  end
end
