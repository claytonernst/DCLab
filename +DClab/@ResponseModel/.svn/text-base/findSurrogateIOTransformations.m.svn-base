function trans = findSurrogateIOTransformations(RM,domain,opt)
%
% If trans.variable{i,j} will be 'none' if the parameter of domain(i) is
% not employed by RM.

% Returns a structure. trans.response is a 1-by-nsurfs cell; trans.variable
% is an n-by-nsurfs cell.

%TODO make this actually do something.

if strcmp(RM.type,'quadraticModel') || strcmp(RM.type,'linearModel') || ~opt.findIOTransforms
    % Use the transformations of the Response model.
    
    trans.response = {RM.responseTransformation};
    trans.variable = repmat({'none'},size(domain));
    
    % for parameters that are used by the model, overwrite with the proper transformations.
    mNames = RM.parameterList;
    idx = NaN(size(mNames));
    domNames = {domain.name}';
    [trash domNames2sortFoundNames mNames2sortFoundNames] = intersect(domNames,mNames); 
    idx(mNames2sortFoundNames) = domNames2sortFoundNames;
    
    trans.variable(idx,1) = RM.variableTransformations;
    
else
    
    n = length(domain);
    
    %defaults
    trans.response = {'none'};
    trans.variable = repmat({'none'},n,1);
    
    %some constants
    Ngrid = 30; %grid points per dimension
    n = size(domain,1); %number of parameters
    t = 0.9; %threshold
    
    %initialize
    X = cell(1,n);
    Y = cell(1,n);
    
    rng = vertcat(domain.range);
    
    for i=1:n %for each parameter
        if rng(i,1)>0 %don't want to do logs for negative numbers
            
            %points for gridding
            gridLin = linspace(rng(i,1),rng(i,2),round(2*Ngrid/3));
            gridLog = logspace(log10(rng(i,1)),log10(rng(i,2)),round(2*Ngrid/3));
            grid = sort([gridLin gridLog(1:(Ngrid-round(2*Ngrid/3)))]);
            N = length(grid);
            
            %evaluate
            x1 = repmat(rng(:,1),1,N);
            x1(i,:) = grid;
            yL = RM.rapideval(x1);
            x2 = repmat(mean(rng,2),1,N);
            x2(i,:) = grid;
            yN = RM.rapideval(x2);
            x3 = repmat(rng(:,2),1,N);
            x3(i,:) = grid;
            yH = RM.rapideval(x3);
            X{i} = [x1 x2 x3];
            Y{i} = [yL yN yH];
            
            %check monotonicity
            isMonoL = all(diff(yL)>=0) || all(diff(yL)<=0);
            isMonoN = all(diff(yN)>=0) || all(diff(yN)<=0);
            isMonoH = all(diff(yH)>=0) || all(diff(yH)<=0);
            
            if isMonoL+isMonoN+isMonoH>=2
                
                %do some estimates for linear
                yLest = [ones(N,1) 2*grid'/diff(rng(i,:)) - sum(rng(i,:))/diff(rng(i,:))]*([ones(N,1) 2*grid'/diff(rng(i,:)) - sum(rng(i,:))/diff(rng(i,:))]\yL');
                yNest = [ones(N,1) 2*grid'/diff(rng(i,:)) - sum(rng(i,:))/diff(rng(i,:))]*([ones(N,1) 2*grid'/diff(rng(i,:)) - sum(rng(i,:))/diff(rng(i,:))]\yN');
                yHest = [ones(N,1) 2*grid'/diff(rng(i,:)) - sum(rng(i,:))/diff(rng(i,:))]*([ones(N,1) 2*grid'/diff(rng(i,:)) - sum(rng(i,:))/diff(rng(i,:))]\yH');
                errLin = std(yL'-yLest)/std(yL) + std(yN'-yNest)/std(yN) + std(yH'-yHest)/std(yH);
                
                %...and for log
                yLest = [ones(N,1) log10(grid)']*([ones(N,1) log10(grid)']\yL');
                yNest = [ones(N,1) log10(grid)']*([ones(N,1) log10(grid')]\yN');
                yHest = [ones(N,1) log10(grid)']*([ones(N,1) log10(grid)']\yH');
                errLog = std(yL'-yLest)/std(yL) + std(yN'-yNest)/std(yN) + std(yH'-yHest)/std(yH);
                
                if errLog <= t*errLin
                    trans.variable{i} = 'log10';
                end
                
            end  % end monotonicity check
            
        end  %end positive check
        
    end  %end for-loop
        
    
    %now to check for the output transform
%     X = horzcat(X{:});
%     Y = horzcat(Y{:});
%     N = length(Y);
%     for i=1:n
%         if strcmp(trans.variable{i},'log10')
%             X(i,:) = log10(X(i,:));
%         end
%     end
%     YestLin = [ones(N,1) X']*([ones(N,1) X']\Y');
%     YestLog = [ones(N,1) X']*([ones(N,1) X']\log10(Y'));
%     errLin = std(YestLin'-Y)/std(Y);
%     errLog = std(YestLog'-log10(Y))/std(log10(Y));
%     if errLog<=t*errLin
%         trans.response = {'log10'};
%     end
    
    if any(RM.uncertaintyCase==[3 4])
        trans.response = {'log10'};
    end
    
end

%TODO, check sign of stuff. In certain cases, this will disallow
%log10(y). I guess we should take L and U as inputs as well in order to
%determine this.



