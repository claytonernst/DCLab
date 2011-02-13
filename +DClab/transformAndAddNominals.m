function xout = transformAndAddNominals(xin,PD,activeIdx,node)
%Transform xfeas and fill in nominal values as necessary

%x should come in rxp where r = # nodes and p = 1+nRestarts


if isempty(xin)
  xout = xin;
  return
end

%tolerance!
xin = round(xin*10^5)/(10^5);

n = PD.nParameters;
linAct = activeIdx <= n;
logAct = activeIdx > n;

linIdx = activeIdx(linAct);
logIdx = activeIdx(logAct)-n;

%if both lin and log exist for a given parameter, lin takes precedence.
logIdx = setdiff(logIdx,linIdx);
[trash logAct] = intersect(activeIdx,logIdx+n);

domrange = PD.PiecewiseSurrogateModelTree(node).domainRange;   
xout = vertcat(PD.FreeParameter.nominal);      

% AKP edit, Jan 22, 2011, failing on Prajna Example
if ~isempty(linIdx)
   xout(linIdx) = 0.5*(diff(domrange(linIdx,:),[],2)).*(xin(linAct)+1) + domrange(linIdx,1);
end
if ~isempty(logIdx)
   xout(logIdx) = domrange(logIdx,1).*((domrange(logIdx,2)./domrange(logIdx,1)).^(0.5*(xin(logAct)+1)));
end