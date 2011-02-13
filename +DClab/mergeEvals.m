function [] = mergeEvals(toLoad)
% function [] = mergeEvals(toLoad)
%
% this function is used for combining saved evaluations from
% different sources. toLoad should be the name of the temporary
% directory that contains the files from a different source. the
% function then increases file names an moves files into subfolders
% of the base EVALUATIONS directory.

file = which('evaluationsFlag');
evalPath = fileparts(file);

toLoadDir = fullfile(evalPath,toLoad);
dirStruct = dir(toLoadDir);
folders = find(vertcat(dirStruct.isdir));
for i1 = folders'
  if strcmp(dirStruct(i1).name,'.') | strcmp(dirStruct(i1).name, ...
                                             '..')
    %do nothing
  else
    name = dirStruct(i1).name;
    if exist(fullfile(evalPath,name),'dir')
      %need to update file names
      matFiles = what(fullfile(toLoadDir,name));
      matFiles = matFiles.mat;
      temp = strmatch('lean',matFiles);
      
      i2 = 1;
      for i3 = 1:length(temp)
    
        while exist(fullfile(evalPath,name,[num2str(i2) '.mat']))
          i2 = i2+1;
        end
                
        %new file will be called i2
        str1 = fullfile(toLoadDir,name,['lean' num2str(i3) '.mat']);
        str2 = fullfile(toLoadDir,name,[num2str(i3) '.mat']);
        fileName1 = fullfile(evalPath,name,['lean' num2str(i2)]);
        fileName2 = fullfile(evalPath,name,[num2str(i2)]);
        try
          load(str1);
          slean = s;
          load(str2);
         
          save(fileName1,'slean');
          save(fileName2,'s');
        catch
          %do nothing
        end
        delete(str1);
        delete(str2);
        
        %increment file name counter
        i2 = i2+1;   
        
      end
    
    else

      str1 = ['mv ' fullfile(toLoadDir,name) ' ' fullfile(evalPath,name)];
      system(str1);
    end
  end
end
