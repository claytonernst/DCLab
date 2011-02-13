function domain = listBnds2domain(list,bnds)


n = length(list);
tmp = cell(n,2);
tmp(:,1) = list;
tmp(:,2) = mat2cell(bnds,ones(n,1),2);
domain = cell2struct(tmp,{'name';'range'},2);
