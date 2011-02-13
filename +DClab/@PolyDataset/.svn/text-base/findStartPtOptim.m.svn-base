function xinit = findStartPtOptim(PD,activeIdx,nodeIndices,startPt,trans,Q,Ntry)
% xinit = findStartPtOptim(PD,nodeIndices,startPt,trans,Q,Ntry)
%
% Inputs
%   PD: A PolyDataset object
%   nodeIndices: rx1 vector of indices to the nodes of the optimization
%   startPt: The start point in the problem's natural cordinates. This
%     input must be empty or a rx1 cell array. Each element of the cell
%     array should be a nx1 vector.
%   trans: char array, either 'linXlinY', 'logXlinY', 'linXlogY', or
%   'logXlogY' 
%   Q: a rx1 cell array containing the primal variable from the S-procedure.
%   Ntry: if Ntry = 1, xinit will be derived from startPt, Q, or
%   randomly, in that order of preference. If Ntry > 1 a
%     random xinit will be created.
%
% Ouputs
%   xinit: a rx1 cell array of n'x1 vectors, where n' is the number of
%     active variables.

r = length(nodeIndices);
xinit = cell(r,1);
m = PD.nPairs;

if Ntry == 1 %use perferably user supplied startPt, or Q derived point
  if ~isempty(startPt)

    n = nParams(PD);
    for i1 = 1:r
      if ~isequal(size(startPt{i1}),[n,1])
          error('dimensions of startPt do not match number of parameters in PolyDataset');
      end

      tempx = startPt{i1}(activeIdx{i1});
      bnds = PD.PiecewiseSurrogateModelTree(nodeIndices(i1)).domainRange(activeIdx{i1},:);
      if ismember(trans,{'m2';'m4'})
        tempx = log10(tempx);
        bnds = log10(bnds);
      end

      %now transform tempx appropriately.
      tempx = 2*(tempx - bnds(:,1))./diff(bnds,[],2) - 1;
      xinit{i1} = [tempx; 0.1*rand(m,1)]; %add last element for gammas
    end
  elseif ~isempty(Q)
      for i1 = 1:r
          xinit{i1} = Q{i1}(2:end,1);
      end	
  else
      %call this function recursively, but set Ntry = 2 so we use a random point.
      xinit = findStartPtOptim(PD,activeIdx,nodeIndices,startPt,trans,Q,2);
  end
  
else
  %in this case, we're just going to generate random points 
  %(although we could use sample from distribution idea...)
  for i1 = 1:r
    n = length(activeIdx{i1});
    %THIS ASSUMES +/- 1 CUBE BOUNDS
    xinit{i1} = [2*rand(n,1)-1 ; 0.1*rand(m,1)];
  end
end

