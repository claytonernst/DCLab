function x=lhsample(nPnts,n,xOld)
% Latin Hypercube Sample
% 
% x = lhsample(nPnts,n)
% Produces a latin hypercube sampling of nPnts on the unit hypercube of
% dimension n.  The resulting x is has size [nPnts n].
%
% A latin hypercube sample is such that for each dimension (each column of
% x) points are taken from a uniform sampling over each uniformly spaced
% bin of [0,1].  Specifically, there will be one point from a uniform
% distribution of [0, 1/nPnts], one from [1/nPnts, 2/nPnts], ... and one from
% [1-1/nPnts,1].  These are then randomly permuted.  
%
% xNew = lhsample(nPnts,n,xOld)
% Produces a latin hypercube sampling of nPnts on the unit hypercube of
% dimension n.  It does so such that the new sampling is distributed
% properly amongst the sampling of points xOld.  Specifically, along each
% dimension, if nPnts number of equally sized bins are produced, the new
% sampling only samples from the bins not filled by xOld.  The function
% returns just the new points xNew.  nPnts represents the number of new
% points to calculate.  Requirement: size(xOld,2)=n.


%error checking
ni = nargin;
error(nargchk(2,3,ni,'struct'))

%check integer values
if ~mod(nPnts,1)
    error('Input ''nPnts'' needs to be an integer')
end
if ~mod(n,1)
    error('Input ''n'' needs to be an integer')
end

%check xOld proper sizes
if ni==3
    if size(xOld,2)~=n
        error('Size mismatch.  size(xOld,2) needs to be n')
    end
end

%check nPnts
if nPnts<0
    error('Input nPnts must be a positive integer')
end
if nPnts==0
    warning('Inputs nPnts is zero')
    x = zeros(0,n);
    return
end


x = zeros(nPnts,n);

switch ni
    
    case 2
        %Basic case
        
        x = rand(nPnts,n)/nPnts;
        toAdd = 0:(1/n):(1-1/n);
        
        for i1=1:n
            idx = randPerm(nPnts);
            x(i1,:) = x(i1,:) + toAdd(idx)';
        end
        
    case 3
        %Addnew points
        
        nOld = size(xOld,1);
        x = rand(nPnts,n)/nPnts;
        
        for i1=1:n
            
            intervals = ones(nPnts+nOld,1);
            
            for i2=1:nOld
                
                intervals(ceil(xOld(i2,i1)*(nOld+nPnts))) = 0;
                
            end
            
            idx = randperm(nPnts);
            toAdd = 0:(1/(nPnts+nOld)):(1-1/(nPnts+nOld));
            toAdd = toAdd(intervals);
            x(i1,:) = x(i1,:) + toAdd(idx)';
            
        end
        
            
    otherwise
        %you should never reach this point
        error('Something''s wrong.  Shouldn''t get this error')
        
end

