currentDir = dir;
tot = 0;
for i1 = 1:length(currentDir)
  if exist(currentDir(i1).name,'dir')
   [trash L] = countcode(currentDir(i1).name);
   tot = tot+L;
  end
end
