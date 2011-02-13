classdef cmoptim < DClab.DCObject
    %CMOPTIM sets up a problem defining the dataset consistency measure
    %
    %   OBJ = CMOPTIM(ZQUAD,D,U,LB,UB) sets up the problem
    %      min \gamma s.t.
    %          LB <= x <= UB
    %          D+U(1)*(1+\gamma) <= ZQUAD(x) <= D+U(2)*(1+\gamma)
    %   ZQUAD is a m-by-1 cell array of sparce symmetric (n+1)-by-(n+1)
    %   matrices. D is a m-by-1 vector, and U is a m-by-2 vector, with the
    %   elements of the 1st column all negative. LB and UB are vectors of
    %   finite upper and lower bounds on the x variables.
    %
    %   OBJ = CMOPTIM(ZQUAD,D,U,LB,UB,OU,UNCCASE,SE,RESPTRANS) sets up a much
    %   nastier problem due to transformations in the constraints.
    %
    %   OBJ = CMOPTIM(ZQUAD,D,U,LB,UB,OU,UNCCASE,SE,RESPTRANS,LINXLOGX)
    %   provides for additional equality constraints of a special form. These
    %   can only be used to constrain variables whose range (as specified by LB
    %   and UB) is a bounded interval. LINXLOGX should be p-by-5.
    %   Let [i j c1 c2 c3] denote a row of LINXLOGX. This row imposes the
    %   equality constraint
    %           UB(i)-LB(i)
    %    X(i) = ---------- c2 ( c1^[(X(j)-LB(j))/(UB(j)-LB(j))] - 1 ) + LB(i)
    %              c3
    %   This equality constraint is relaxed using quadratic inequalities by the
    %   lowerBnd method. The upperBnd method handles it directly.
    %
    %   Inputs:
    %      ZQUAD: is an m-by-1 cell array of sparce symmetric (n+1)-by-(n+1)
    %         matrices.
    %      D: is a m-by-1 vector
    %      U: is a m-by-2 vector, first column is negative
    %      OU: is a m-by-2 vector, first column is negative
    %      UNCASE: is a m-by-1 vector of the integers 1:16
    %      SE: is a m-by-2 vector, first column is negative
    %      RESPTRANS: is a m-by-1 cell array of 'log10' or 'none'.
    %      LINXLOGX is a p-by-5 matrix. See above.
    
    properties
        type = 'empty';
        ZquadFinal;
        tau;
        ZrelaxedEq;
        nVars;
        variableBnds;
        linXlogX;
        xfeas;
    end
    
    methods
        
        function obj = cmoptim(Zquad,d,u,LB,UB,ou,uncCase,se,respTrans,linXlogX)
            
            ni = nargin;
            
            switch ni
                case 0
                    obj = initdefaultobj;
                    return
                case 5
                    m = length(Zquad);
                    ou = zeros(m,2);
                    uncCase = ones(m,1);
                    se = zeros(m,2);
                    respTrans = repmat({'none'},m,1);
                    linXlogX = [];
                case 9
                    linXlogX = [];
                case 10
                    %do nothing
                otherwise
                    error('Incorrect number of input arguments.')
            end
            
            % ==Input checking==
            
            n = length(UB);
            % variable bounds
            if isempty(UB) || isempty(LB)
                error('Inputs: LB and UB must be nonempty n-by-1 vectors.')
            else
                if ~isequal(size(LB),[n 1]) || ~isequal(size(UB),[n 1])
                    error('Inputs: LB and UB must be n-by-1 vectors.')
                end
                
                % Make sure they are finite.
                if any(any(isinf([LB UB]))) || any(any(isnan([LB UB])))
                    error('Inputs: LB and UB cannot contain infinite or NaN elements.')
                end
                
                % TODO: this 100*eps is a hack. I'm not sure how close these two bounds
                % can be before we start having numerical difficulties.
                if any(UB - LB < 100*eps)
                    error('Inputs: box constraint implied by LB and UB appears to be empty or have empty interior.')
                end
            end
            
            % constraint functions:
            if isempty(Zquad)
                error('Inputs: ZQUAD cannot be empty.')
            elseif ~iscell(Zquad) || size(Zquad,2) ~= 1
                error('Inputs: ZQUAD must be a column cell array')
            else
                if ~all(cellfun('prodofsize',Zquad) == (n+1)^2)
                    error('Inputs: ZQUAD must contain matrices of identical dimension as ZNOT')
                end
                for i1 = 1:length(Zquad)
                    if ~DClab.isrealsymmatrix(Zquad{i1})
                        Zquad{i1} = 0.5*(Zquad{i1}+Zquad{i1}');
                    end
                end
            end
            
            % d and u
            m = length(Zquad);
            if ~isequal(size(d),[m 1])
                error('Inputs: D of improper dimension')
            end
            if ~isequal(size(u),[m 2])
                error('Inputs: U of improper dimension')
            end
            if any(u(:,1) >= 0) || any(u(:,2) <= 0)
                error('Inputs: Components of U have improper sign')
            end
            
            % ou
            if ~isequal(size(ou),[m 2])
                error('Inputs: OU of improper dimension')
            end
            % 1st column should be all nonpositive, 2nd all nonnegative (outer approx
            % of F), OR visa versa.
            if all(sign(ou(:,1)) <= 0) && all(sign(ou(:,2)) >= 0)
                % must be outer approx.
            elseif all(sign(ou(:,1)) >= 0) && all(sign(ou(:,2)) <= 0)
                % must be inner approx.
            else
                error('Inputs: Components of OU have improper sign')
            end
            
            % uncCase
            if ~isequal(size(uncCase),[m 1])
                error('Inputs: UNCCASE of improper dimension')
            end
            if ~isequal(uncCase,round(uncCase)) || min(uncCase) < 1 || max(uncCase) > 16
                error('Inputs: UNCCASE must be a vector of the integers 1:16.')
            end
            
            % se
            if ~isequal(size(se),[m 2])
                error('Inputs: SE of improper dimension')
            end
            % 1st column should be all nonpositive, 2nd all nonnegative (outer approx
            % of F), OR visa versa.
            if all(sign(se(:,1)) <= 0) && all(sign(se(:,2)) >= 0)
                % must be outer approx.
            elseif all(sign(se(:,1)) >= 0) && all(sign(se(:,2)) <= 0)
                % must be inner approx.
            else
                %TODO: legally this can occur if on critical range, quad is totally below
                %or above model.
                %error('Inputs: Components of SE have improper sign')
            end
            
            % respTrans
            if ~isequal(size(respTrans),[m 1])
                error('Inputs: RESPTRANS of improper dimension')
            end
            if ~iscell(respTrans)
                error('Inputs: RESPTRANS must be a column cell array')
            end
            idx1 = strmatch('none',respTrans,'exact');
            idx2 = strmatch('log10',respTrans,'exact');
            
            if ~isequal(union(idx1,idx2),(1:m)')
                error('Inputs: each element of RESPTRANS must be ''log10'' or ''none''.')
            end
            
            % We will have a "simple" case, and a "messy case"
            
            %Simple cases:
            % no response transformation, uncCase =
            %   1,2,5,6,9,10,
            %   13 and 14 with no output uncertainty
            % log10 response transformation, uncCase =
            %   3 and 4 with no output uncertainty
            %   7,8,11,12,15,16.
            
            % Determine whether the present case is simple or messy. Simple cases are
            % those where u enters linearly in the right-hand side of the quadratic
            % constraints. This prevent us from having to bisect, etc, do to
            % nonlinearities.
            simple = true;
            for i1 = 1:m
                if strcmp(respTrans{i1},'none')
                    if ~ ( ismember(uncCase(i1),[1 2 5 6 9 10]) || ( ismember(uncCase(i1),[13 14]) && all(ou(i1,:)==0)) )
                        simple = false;
                    end
                else
                    if ~ ( ismember(uncCase(i1),[7 8 11 12 15 16]) || ( ismember(uncCase(i1),[3 4]) && all(ou(i1,:)==0)) )
                        simple = false;
                    end
                end
            end
            
            % Create quadratic relaxations of any equality constraints. We employ four
            % quadratic constraints per equality constraint.
            if isempty(linXlogX)
                Z4relaxedEqualities = [];
            else
                % Use the nqcqp code to create the matrices to avoid having multiple
                % copies of buildRelaxedEqualities floating around
                tmpobj = DClab.nqcqp(zeros(n+1),[],LB,UB,linXlogX);
                
                Z4relaxedEqualities = tmpobj.ZrelaxedEq;
                %  Z4relaxedEqualities = [];
                %  linXlogX = [];
            end
            
            %disp('deleting linlogconstraints===================')
            %Z4relaxedEqualities = [];
            %linXlogX = [];
            
            % Create object
            
            %If the object is simple, the resulting problem looks like
            %  min \gamma s.t.
            %    LB <= x <= UB
            %    Quad(x) <= c\gamma
            
            if simple
                obj.type = 'simple';
                tau = zeros(2*m,1);
                ZquadFinal = cell(2*m,1);
                for i1 = 1:m
                    switch uncCase(i1)
                        case {1,13} % recall, by the conditions for 'simple', ou will be zero if 13
                            tempL = -Zquad{i1};
                            tempL(1,1) = tempL(1,1) + d(i1) - ou(i1,2) + u(i1,1) + se(i1,1);
                            tau(i1*2-1) = -u(i1,1);
                            
                            tempR = Zquad{i1};
                            tempR(1,1) = tempR(1,1) - d(i1) + ou(i1,1) - u(i1,2) - se(i1,2);
                            tau(i1*2) = u(i1,2);
                            
                            ZquadFinal{i1*2-1} = tempL;
                            ZquadFinal{i1*2} = tempR;
                            
                        case {2,14} % recall, by the conditions for 'simple', ou will be zero if 14
                            tempL = -Zquad{i1};
                            tempL(1,1) = tempL(1,1) + d(i1)*(1+u(i1,1)) - ou(i1,2) + se(i1,1);
                            tau(i1*2-1) = -d(i1)*u(i1,1);
                            
                            tempR = Zquad{i1};
                            tempR(1,1) = tempR(1,1) - d(i1)*(1+u(i1,2)) + ou(i1,1) - se(i1,2);
                            tau(i1*2) = d(i1)*u(i1,2);
                            
                            ZquadFinal{i1*2-1} = tempL;
                            ZquadFinal{i1*2} = tempR;
                        case 3
                            % will be logY with no output uncertainty
                            tempL = -Zquad{i1};
                            tempL(1,1) = tempL(1,1) + log10(d(i1)) + u(i1,1) + se(i1,1);
                            tau(i1*2-1) = -u(i1,1);
                            
                            tempR = Zquad{i1};
                            tempR(1,1) = tempR(1,1) - log10(d(i1)) - u(i1,2) - se(i1,2);
                            tau(i1*2) = u(i1,2);
                            
                            ZquadFinal{i1*2-1} = tempL;
                            ZquadFinal{i1*2} = tempR;
                        case 4
                            % will be logY with no output uncertainty
                            tempL = -Zquad{i1};
                            tempL(1,1) = tempL(1,1) + log10(d(i1))*(1+u(i1,1)) + se(i1,1);
                            tau(i1*2-1) = -log10(d(i1))*u(i1,1);
                            
                            tempR = Zquad{i1};
                            tempR(1,1) = tempR(1,1) - log10(d(i1))*(1+u(i1,2)) - se(i1,2);
                            tau(i1*2) = log10(d(i1))*u(i1,2);
                            
                            ZquadFinal{i1*2-1} = tempL;
                            ZquadFinal{i1*2} = tempR;
                        case 5
                            tempL = -(1+ou(i1,2))*Zquad{i1};
                            tempL(1,1) = tempL(1,1) + d(i1) + u(i1,1) + se(i1,1)*(1+ou(i1,2));
                            tau(i1*2-1) = -u(i1,1);
                            
                            tempR = (1+ou(i1,1))*Zquad{i1};
                            tempR(1,1) = tempR(1,1) - d(i1) - u(i1,2) - se(i1,2)*(1+ou(i1,1));
                            tau(i1*2) = u(i1,2);
                            
                            ZquadFinal{i1*2-1} = tempL;
                            ZquadFinal{i1*2} = tempR;
                        case 6
                            tempL = -(1+ou(i1,2))*Zquad{i1};
                            tempL(1,1) = tempL(1,1) + d(i1)*(1 + u(i1,1)) + se(i1,1)*(1+ou(i1,2));
                            tau(i1*2-1) = -d(i1)*u(i1,1);
                            
                            tempR = (1+ou(i1,1))*Zquad{i1};
                            tempR(1,1) = tempR(1,1) - d(i1)*(1 + u(i1,2)) - se(i1,2)*(1+ou(i1,1));
                            tau(i1*2) = d(i1)*u(i1,2);
                            
                            ZquadFinal{i1*2-1} = tempL;
                            ZquadFinal{i1*2} = tempR;
                        case 7
                            % will be logY
                            tempL = -Zquad{i1};
                            tempL(1,1) = tempL(1,1) + log10(d(i1)) + u(i1,1) - log10(1 + ou(i1,2)) + se(i1,1);
                            tau(i1*2-1) = -u(i1,1);
                            
                            tempR = Zquad{i1};
                            tempR(1,1) = tempR(1,1) - log10(d(i1)) - u(i1,2) + log10(1 + ou(i1,1)) - se(i1,2);
                            tau(i1*2) = u(i1,2);
                            
                            ZquadFinal{i1*2-1} = tempL;
                            ZquadFinal{i1*2} = tempR;
                        case 8
                            % will be logY
                            tempL = -Zquad{i1};
                            tempL(1,1) = tempL(1,1) + log10(d(i1))*(1+u(i1,1)) - log10(1 + ou(i1,2)) + se(i1,1);
                            tau(i1*2-1) = -log10(d(i1))*u(i1,1);
                            
                            tempR = Zquad{i1};
                            tempR(1,1) = tempR(1,1) - log10(d(i1))*(1+u(i1,2)) + log10(1 + ou(i1,1)) - se(i1,2);
                            tau(i1*2) = log10(d(i1))*u(i1,2);
                            
                            ZquadFinal{i1*2-1} = tempL;
                            ZquadFinal{i1*2} = tempR;
                            
                        case 9
                            tempL = -(10^ou(i1,2))*Zquad{i1};
                            tempL(1,1) = tempL(1,1) + d(i1) + u(i1,1) + se(i1,1)*10^ou(i1,2);
                            tau(i1*2-1) = -u(i1,1);
                            
                            tempR = (10^ou(i1,1))*Zquad{i1};
                            tempR(1,1) = tempR(1,1) - d(i1) - u(i1,2) - se(i1,2)*10^ou(i1,1);
                            tau(i1*2) = u(i1,2);
                            
                            ZquadFinal{i1*2-1} = tempL;
                            ZquadFinal{i1*2} = tempR;
                        case 10
                            tempL = -(10^ou(i1,2))*Zquad{i1};
                            tempL(1,1) = tempL(1,1) + d(i1)*(1 + u(i1,1)) + se(i1,1)*10^ou(i1,2);
                            tau(i1*2-1) = -d(i1)*u(i1,1);
                            
                            tempR = (10^ou(i1,1))*Zquad{i1};
                            tempR(1,1) = tempR(1,1) - d(i1)*(1 + u(i1,2)) - se(i1,2)*10^ou(i1,1);
                            tau(i1*2) = d(i1)*u(i1,2);
                            
                            ZquadFinal{i1*2-1} = tempL;
                            ZquadFinal{i1*2} = tempR;
                            
                        case 11
                            % will be logY
                            tempL = -Zquad{i1};
                            tempL(1,1) = tempL(1,1) + log10(d(i1)) + u(i1,1) - ou(i1,2) + se(i1,1);
                            tau(i1*2-1) = -u(i1,1);
                            
                            tempR = Zquad{i1};
                            tempR(1,1) = tempR(1,1) - log10(d(i1)) - u(i1,2) + ou(i1,1) - se(i1,2);
                            tau(i1*2) = u(i1,2);
                            
                            ZquadFinal{i1*2-1} = tempL;
                            ZquadFinal{i1*2} = tempR;
                        case 12
                            % will be logY
                            tempL = -Zquad{i1};
                            tempL(1,1) = tempL(1,1) + log10(d(i1))*(1+u(i1,1)) - ou(i1,2) + se(i1,1) ;
                            tau(i1*2-1) = -log10(d(i1))*u(i1,1);
                            
                            tempR = Zquad{i1};
                            tempR(1,1) = tempR(1,1) - log10(d(i1))*(1+u(i1,2)) + ou(i1,1) - se(i1,2);
                            tau(i1*2) = log10(d(i1))*u(i1,2);
                            
                            ZquadFinal{i1*2-1} = tempL;
                            ZquadFinal{i1*2} = tempR;
                            
                            %cases 13 and 14 are combined with 1 and 2
                        case 15
                            % will be logY
                            tempL = -(1+ou(i1,2))*Zquad{i1};
                            tempL(1,1) = tempL(1,1) + log10(d(i1)) + u(i1,1) + se(i1,1)*(1+ou(i1,2));
                            tau(i1*2-1) = -u(i1,1);
                            
                            tempR = (1+ou(i1,1))*Zquad{i1};
                            tempR(1,1) = tempR(1,1) - log10(d(i1)) - u(i1,2) - se(i1,2)*(1+ou(i1,1));
                            tau(i1*2) = u(i1,2);
                            
                            ZquadFinal{i1*2-1} = tempL;
                            ZquadFinal{i1*2} = tempR;
                        case 16
                            % will be logY
                            tempL = -(1+ou(i1,2))*Zquad{i1};
                            tempL(1,1) = tempL(1,1) + log10(d(i1))*(1 + u(i1,1)) + se(i1,1)*(1+ou(i1,2));
                            tau(i1*2-1) = -log10(d(i1))*u(i1,1);
                            
                            tempR = (1+ou(i1,1))*Zquad{i1};
                            tempR(1,1) = tempR(1,1) - log10(d(i1))*(1 + u(i1,2)) - se(i1,2)*(1+ou(i1,1));
                            tau(i1*2) = log10(d(i1))*u(i1,2);
                            
                            ZquadFinal{i1*2-1} = tempL;
                            ZquadFinal{i1*2} = tempR;
                            
                        otherwise
                            error('Internal code inconsistency: condition should not occur')
                    end
                end
                
                obj.ZquadFinal = ZquadFinal;
                obj.tau = tau;
            else
                error('Code not complete')
            end
            
            obj.ZrelaxedEq = Z4relaxedEqualities;
            obj.nVars = n;
            obj.variableBnds = [LB UB];
            obj.linXlogX = linXlogX;
            obj.xfeas = [];
            
        end
        
        function [list,sz] = displayProps(obj)
            list = [];
            sz = '';
        end
        
        
        function [lb,lambda,Q] = lowerBnd(cmObj,opt)
            %LOWERBND method of cmoptim produces a lower bound on the optimal value
            %
            %   An cmoptim object describes an optimization with box constraints,
            %   nonlinear constraints that are quadratic in x, but may have some nasty
            %   dependence on \gamma. If .type = 'simple', the dependence on \gamma is
            %   linear. Additionally there may be some special purpose linXlogX
            %   equality constraints. This function determines a lower bound on the
            %   optimal value by calling SeDuMi. If .type = 'messy', bisection is
            %   required to solve the problem.
            %
            %   REMAINDER of help is not necessarily accurate
            %
            %   First each dimension of
            %   the box constraints is either a bounded interval x_i \in [A,B] or
            %   x_i \in (-Inf,Inf). Those that are a bounded interval are converted
            %   into the equivalent quadratic constraint (B-x_i)(A-x_i) <= 0. Let Zbox
            %   denote the collection of such constraints, and let
            %   Zconst = [Zbox; Zquad; ZrelaxedEq]
            %
            %   This function finds a lower bound on the following nonconvex,
            %   quadratic program (x is an n-dimensional vector):
            %      p=min [1 x']*Znot*[1; x]
            %         s.t.  [1 x']*Zconst_k*[1; x] <= 0  (k=1,...,N)
            %   The S-procedure gives a semidefinite program (SDP) that lower
            %   bounds the optimal cost (lb<=p):
            %      lb = max gamma
            %         s.t.  Znot - [gamma 0; 0 0] + sum_k lambda_k Zconst_k >= 0
            %               lambda_k >= 0
            %   Rank relaxation gives the dual form of this SDP and typically
            %   yields the same lower bound:
            %      lb = min trace[Znot Q]
            %         s.t.  Q>=0, Q_{11} = 1
            %               Tr[Zconst_k Q] <= 0   (k=1,...,N)
            %
            %   LB = LOWERBND(CMOPTIMOBJ) produces a lower bound LB on the optimal value
            %   of the optimization stored in NQCQPOBJ. If the problem is infeasible,
            %   LB will be +Inf.
            %
            %   LB = LOWERBND(CMOPTIMOBJ,OPT) uses options taken from the DCOptions
            %   object OPT. This method used OPT.display and OPT.sedumiParEps.
            %
            %   [LB, LAMBDA] = LOWERBND(...) additionally returns a structure array
            %   LAMBDA with fields .lower, .upper, .quadineq, .relaxedeq. Each field
            %   contains the S-procedure multipliers lambda for the corresponding
            %   constraints.
            %
            %   [LB, LAMBDA, Q] = LOWERBND(...) additionally returns the optimal
            %   (n+1)x(n+1) matrix Q from the rank relaxation.
            
            
            %TODO use linprog if everything is linear
            
            error(nargchk(1,2,nargin));
            error(nargoutchk(0,3,nargout));
            
            % Initialize inputs
            error(nargchk(1,2,nargin));
            if nargin==1
                opt = DClab.DCOptions;
            end
            
            % Note this problem is always feasible.
            
            % Problem Data
            n = cmObj.nVars;
            tau = cmObj.tau;
            Zquad = cmObj.ZquadFinal;
            ZrelaxedEq = cmObj.ZrelaxedEq;
            
            %convert variableBnds to quadratic constraints
            xbnds = cmObj.variableBnds;
            finiteXBnds = find(~isinf(xbnds(:,1)));
            
            nn = length(finiteXBnds);
            % Create quadratic representations of the box constraints
            
            if opt.fullBlockHRelaxation
                
                % TODO we assume all bounds are +/- 1. As of 2/12/08, our software only
                % produced such objects.
                
                % full box version.
                
                Zbox = cell(nn^2,1);
                [Zbox{:}] = deal(spalloc(n+1,n+1,7));
                counter = 1;
                for i1 = reshape(finiteXBnds,1,nn)
                    for i2 = finiteXBnds'
                        % affine
                        Zbox{counter}(1,1) = -1;
                        
                        % linear
                        Zbox{counter}(i1+1,1) = 0.5;
                        Zbox{counter}(1,i1+1) = 0.5;
                        Zbox{counter}(i2+1,1) = Zbox{counter}(i2+1,1) - 0.5;
                        Zbox{counter}(1,i2+1) = Zbox{counter}(1,i2+1) - 0.5;
                        
                        % quadratic
                        Zbox{counter}(i1+1,i2+1) = 0.5;
                        Zbox{counter}(i2+1,i1+1) = Zbox{counter}(i2+1,i1+1) + 0.5;
                        
                        counter = counter+1;
                    end
                end
                
            else
                
                % alpha trick version
                Zbox = cell(2*nn,1);
                [Zbox{:}] = deal(spalloc(n+1,n+1,4));
                epsilon = zeros(n,1);
                affineTerm = zeros(n,1);
                linTerm = zeros(n,1);
                
                epsilon(finiteXBnds) = 0.05*(xbnds(finiteXBnds,2)-xbnds(finiteXBnds,1));
                affineTerm(finiteXBnds) = xbnds(finiteXBnds,1).* (xbnds(finiteXBnds,2)+epsilon(finiteXBnds));
                linTerm(finiteXBnds) = -(xbnds(finiteXBnds,1) + xbnds(finiteXBnds,2) + epsilon(finiteXBnds));
                for i2 = reshape(finiteXBnds,1,nn)
                    Zbox{i2}(1,1) = affineTerm(i2);
                    Zbox{i2}(i2+1,1) = linTerm(i2)/2;
                    Zbox{i2}(1,i2+1) = linTerm(i2)/2;
                    Zbox{i2}(i2+1,i2+1) = 1;
                end
                affineTerm(finiteXBnds) = xbnds(finiteXBnds,2).* (xbnds(finiteXBnds,1)-epsilon(finiteXBnds));
                linTerm(finiteXBnds) = -(xbnds(finiteXBnds,1) + xbnds(finiteXBnds,2) - epsilon(finiteXBnds));
                for i2 = reshape(finiteXBnds,1,nn)
                    Zbox{i2+nn}(1,1) = affineTerm(i2);
                    Zbox{i2+nn}(i2+1,1) = linTerm(i2)/2;
                    Zbox{i2+nn}(1,i2+1) = linTerm(i2)/2;
                    Zbox{i2+nn}(i2+1,i2+1) = 1;
                end
            end
            
            %   % Old version
            %  Zbox = cell(nn,1);
            %  [Zbox{:}] = deal(spalloc(n+1,n+1,4));
            %  affineTerm = xbnds(:,1).* xbnds(:,2);
            %  linTerm = -xbnds(:,1) - xbnds(:,2);
            %  for i2 = 1:nn
            %    Zbox{i2}(1,1) = affineTerm(i2);
            %    Zbox{i2}(i2+1,1) = linTerm(i2)/2;
            %    Zbox{i2}(1,i2+1) = linTerm(i2)/2;
            %    Zbox{i2}(i2+1,i2+1) = 1;
            %  end
            
            Zconst = [Zbox; Zquad; ZrelaxedEq];
            
            Nbox = size(Zbox,1);
            Nquad = size(Zquad,1);
            NrelaxedEq = size(ZrelaxedEq,1);
            N = length(Zconst);
            
            %TODO, go through sedumi help and document the sedumi call better. Maybe
            %even in pdf.
            
            % Sedumi (dual-form) notation:
            %     max b'y  s.t. c-A'*y \in K^*
            
            % Define y:=[lambda_1 ... lambda_N gamma]
            b = spalloc(N+1,1,1); b(end) = 1;
            
            % One primal free variable.
            K.f = 1;
            
            % dual constraint: C - sum \lamb_i*c_i = 0
            cc = 1;
            
            Ac = [zeros(1,Nbox) tau' zeros(1,NrelaxedEq) 0];
            
            % N linear constraints:
            
            % dual constraint: C - lamb_i >= 0
            
            K.l = N;
            clin = spalloc(N,1,0);
            Alin = [-speye(N) spalloc(N,1,0)];
            
            % (n+1)x(n+1) LMI constraint:
            %    Z0 - [gamma 0; 0 0] + sum_k lambda_k Z_k >= 0
            
            % dual constraint C - A'\lamb \in PSD
            K.s = n+1;
            clmi = spalloc((n+1)^2,1,0);
            Almi = [];
            if N > 0
                Almi = zeros(numel(Zconst{1}),N);
                for i2 = 1:N
                    Almi(:,i2) = -Zconst{i2}(:);
                end
            end
            
            temp = spalloc(n+1,n+1,1); temp(1,1) = 1;
            Almi = sparse([Almi temp(:)]);
            
            % Concatenate the constraints
            c = [cc; clin; clmi];
            At = [Ac; Alin; Almi];
            
            % Call Sedumi
            if strcmp(opt.display,'ALL')
                pars.fid = 1;
            else
                pars.fid=0;    % fid=0: no output; fid=1: display output
            end
            pars.eps=opt.sedumiParEps; % Desired accuracy, sedumi default = 1e-8;
            
            K.xcomplex = [];
            K.scomplex = [];
            K.ycomplex = [];
            [x,y,info]=sedumi(At,b,c,K,pars);
            
            %TODO, we may want to revert to empty outputs if there are numerical or
            %infeasibility problems.
            
            % Output information
            lb = y(end);
            
            lambda.lower = zeros(n,1);
            lambda.upper = zeros(n,1);
            
            if opt.fullBlockHRelaxation
                disp('No X sensitivities computed for fullBlockHRelaxation method')
            else
                lambda.lower(finiteXBnds) = y(1:Nbox/2);
                lambda.upper(finiteXBnds) = y(Nbox/2+1:Nbox);
            end
            lambda.quadineq = y(Nbox+1:Nbox+Nquad);
            lambda.relaxedeq = y(Nbox+Nquad+1:end-1);
            
            Q = reshape(x(N+2:end),K.s,K.s);
            
            % Parse info
            str = '';
            if info.numerr == 2 && ismember(opt.display,{'notify','all','ALL'})
                str = '  Warning! sedumi numerical error, returning -Inf outer bound';
                lb(i1,1) = -inf;
            elseif info.numerr == 1 && ismember(opt.display,{'notify','all','ALL'})
                str = '  Warning: sedumi numerical problems warning, results are accurate to the level of PARS.bigeps';
            elseif info.dinf==1 && info.pinf==1 && ismember(opt.display,{'notify','all','ALL'})
                str = '  In sedumi call: both primal and dual infeasible, expect infinite answers';
            elseif info.pinf == 1 && ismember(opt.display,{'notify','all','ALL'})
                str = '  In sedumi call: primal infeasible, dual feasible';
            elseif info.dinf==1 && ismember(opt.display,{'notify','all','ALL'})
                str = '  In sedumi call: primal feasible, dual infeasible';
            else
                %do nothing
            end
            
            if ~isempty(str)
                DClab.dcdispstr(str,opt.guiHandle,false)
            end
        end
        
        
        function [ub,xopt,mult] = upperBnd(varargin)
            %UPPERBND method of cmoptim produces an upper bound on the optimal value
            %
            %   An cmoptim object describes an optimization with box constraints,
            %   nonlinear constraints that are quadratic in x, but may have some nasty
            %   dependence on \gamma. If .type = 'simple', the dependence on \gamma is
            %   linear. Additionally there may be some special purpose linXlogX
            %   equality constraints. This function determines a upper bound on the
            %   optimal value by calling FMINCON. Future functionality includes
            %   providing the capability to use additional solvers.
            %
            %   UB = UPPERBND(CMOPTIMOBJ) produces an upper bound UB on the optimal value
            %   of the optimization stored in CMOPTIMOBJ. The solver will be
            %   initialized at a random point drawn from a uniform distribution over
            %   the rectanbe described by the box constraints. If a feasible point
            %   cannot be found, UB will be +Inf.
            %
            %   UB = UPPERBND(NQCQPOBJ,XINIT) initializes the solver at the n-by-1
            %   vector XINIT. Note, XINIT should contain only the x-variable, not
            %   gamma.
            %
            %   UB = UPPERBND(NQCQPOBJ,XINIT,OPT) uses options taken from the DCOptions
            %   object OPT. This method used OPT.display and XXX. Set XINIT=[] to
            %   randomly generate the initial point.
            %
            %   UB = UPPERBND(NQCQPOBJ,XINIT,OPT,RELAXEQ) use the quadratic relaxations
            %   of any special purpose linXlogX equality constraints instead of the
            %   equality constraints when RELAXEQ==TRUE. (FALSE is the default). This
            %   might be used in order to determine the "gap" of the SDP lower bound,
            %   since it uses these quadratic relaxations. Set XINIT=[] to randomly
            %   generate the initial point; set OPT=[] to use the default options.
            %
            %   [UB XOPT] = UPPERBND(...) additionally returns the locally optimal
            %   decision vector XOPT (not including the gamma variable) that produces
            %   UB.
            %
            %   [UB XOPT MULT] = UPPERBND(...) additionally returns a structure of
            %   Lagrange multipliers LAMBDA with fields .lower, .upper, .quadineq,
            %   .relaxedeq. All elements of .relaxedeq will be 0 if RELAXEQ==FALSE
            
            
            if strcmp(varargin{1}.type,'empty')
                warning('Received empty object, returning empty outputs')
                ub = [];
                xopt = [];
                mult = [];
            elseif strcmp(varargin{1}.type,'simple')
                [ub xopt mult] = DClab.cmoptim.ubSimple(varargin{:});
            else
                %TODO
                error('Code not complete')
            end
        end
        
    end %public methods
    
    methods (Access=private,Static)
        
        %[iconval,econval,icongrd,econgrd] = fmconfun(xRem,Zquad,tau,x2XRem,x2XElim,xRem2XEqFcnXElim,t);
        %[iconval,econval,icongrd,econgrd] = fmconfunRelax(x,Zi,tau);
        %[objval,objgrd] = fmobjfun(xRem);
        %[ub,xopt,mult] = ubSimple(obj,xinit,opt,relaxEq);
        
        function [iconval,econval,icongrd,econgrd] = fmconfun(xRem,Zquad,tau,x2XRem,x2XElim,xRem2XEqFcnXElim,t)
            %   This function evaluates each of the quadratic inequality constraints and
            %   computes their gradiants.
            %
            %   Inputs:
            %   xRem: the decision vector with x^g from any equality constraints
            %      eliminated, and gamma pasted onto the end
            %   Zquad: a cell array of (n+1)-by-(n+1) sparse symmetric matricies.
            %   x2XRem: location in x of the elements of the smaller decision vector
            %      xRem.
            %   x2XElim: location in x of the elements that were eliminated through the
            %      equality constraints.
            %   xRem2XEqFncXElim: location in xRem of the variables that are related
            %      through the equality constraint to the variables eliminated from x.
            %   t: Neq-by-4 matrix that contains the coefficients for the
            %      transformations relating x(x2XElim) and xRem(xRem2XEqFcnXElim).
            
            % Compute inequality constraint functions and and their gradients
            
            gamma = xRem(end);
            xRem = xRem(1:end-1);
            
            nxr = length(xRem);
            
            if ~isempty(Zquad)
                m = length(Zquad);
                n = size(Zquad{1},1)-1;
                
                iconval = zeros(m,1);
                icongrd = zeros(nxr+1,m);
                
                % Recover the full decision vector
                x = zeros(n,1);
                x(x2XRem) = xRem;
                
                % Using the simplications described in qpub_FMINCON, the equality
                % constraint is
                %
                % xg = f_c^{-1}(xl) = t1 * log10(t2*(xl-t3) + 1) + t4;
                
                x(x2XElim) = t(:,1).*log10( t(:,2).*(xRem(xRem2XEqFcnXElim)-t(:,3)) + 1 ) + t(:,4);
                
                xx = [1;x];
                % Determine if we need to compute the gradients.
                if nargout <= 2
                    for i1=1:m
                        iconval(i1) = xx'*Zquad{i1}*xx - tau(i1)*gamma;
                    end
                else
                    for i1 = 1:m
                        tmp1 = zeros(nxr,n+1);
                        tmp1(:,x2XRem+1) = eye(nxr);
                        tmp2 = zeros(nxr,n+1);
                        dirFcInv = (1/log(10))*t(:,1).*t(:,2)./(t(:,2).*(xRem(xRem2XEqFcnXElim)-t(:,3))+1);
                        tmp2(xRem2XEqFcnXElim,x2XElim+1) = diag(dirFcInv);
                        
                        tmp3 = Zquad{i1}*xx;
                        
                        icongrd(:,i1) = [2*(tmp1+tmp2)*tmp3; - tau(i1)];
                        iconval(i1) = xx'*tmp3 - tau(i1)*gamma;
                    end
                    
                end
                
            else
                
                %TODO check this case. Will the function even be called? What dimensions
                %does matlab expect for the output? How does it determine this? Basically
                %we are trying to explicitly represent an empty constraint.
                iconval = -1; %anything feasible will do.
                icongrd = zeros(nxr+1,1);
            end
            
            % No Equality constraints
            econval = 0;
            econgrd = zeros(nxr+1,1);
        end
        
        
        function [iconval,econval,icongrd,econgrd] = fmconfunRelax(x,Zi,tau)
            
            % Compute inequality constraint functions and and their gradients
            
            % The constraints are [1 x(1:end-1)']*Zi*[1;x(1:end-1)] <= tau*x(end)
            
            if ~isempty(Zi)
                
                N = length(Zi);
                n = length(x);
                
                iconval = zeros(N,1);
                icongrd = zeros(n,N);
                
                % gamma is the last variable.
                xx = [1;x(1:end-1)];
                for i1=1:N
                    tmp = Zi{i1}*xx;
                    icongrd(:,i1) = [2*tmp(2:end); -tau(i1)]; % = 2*grad(xx)*Zi*xx = 2*[zeros(n,1) eye(n)]*Zi*xx
                    iconval(i1) = xx'*tmp - tau(i1)*x(end);
                end
                
            else
                n = length(x);
                iconval = -1;
                icongrd = zeros(n,1);
            end
            
            % No Equality constraints
            econval = 0;
            econgrd = zeros(n,1);
        end
        
        
        function [objval,objgrd] = fmobjfun(xRem)
            %   This function evaluates the linear objective function and computes
            %   its gradient.
            %
            %   Inputs:
            %   xRem: the decision vector with x^g from any directly enforced special
            %   purpose equality constraints eliminated.
            
            % Compute the objective function and its gradient.
            nxr = length(xRem);
            
            objval = xRem(end);
            objgrd = zeros(nxr,1);
            objgrd(end) = 1;
        end
        
        
        function [ub,xopt,mult] = ubSimple(obj,xinit,opt,relaxEq)
            %UBSIMPLE produces an upper bound on the optimal value for 'simple' case.
            %
            %   See the help on DClab.cmoptim/upperBnd for a description of the input
            %   arguments.
            %
            %   UB = UBSIMPLE(OBJ) is the most simple syntax.
            %
            %   [UB XOPT MULT] = UBSIMPLE(OBJ,XINIT,OPT,RELAXEQ) is the most complicated
            %   syntax.
            
            
            % This problem is always feasible.
            
            % Initialize inputs
            ni = nargin;
            if ni==1
                xinit = [];
                opt = [];
                relaxEq = [];
            elseif ni==2
                opt = [];
                relaxEq = [];
            elseif ni==3
                relaxEq = [];
            else
                error(nargchk(1,4,ni));
            end
            
            % Problem Data
            n = obj.nVars;      % Number of variables
            nQuad = size(obj.ZquadFinal,1);
            nRelaxEq = size(obj.ZrelaxedEq,1);
            Zquad = obj.ZquadFinal;
            tau = obj.tau;
            LB = obj.variableBnds(:,1);
            UB = obj.variableBnds(:,2);
            linXlogX = obj.linXlogX;
            
            % Initialize empty inputs
            if isempty(xinit)
                % Create a random initial point. Use 0 for unbounded variables. By
                % the requirements of the object they will be unbounded above and
                % below
                bnded = ~isinf(LB);
                xinit = zeros(n,1);
                xinit(bnded) = (UB(bnded)-LB(bnded).*rand(sum(bnded),1)) + LB(bnded);
                
                % Make sure xinit at least satisfies the equality constraints.
                % Let [i j c1 c2 c3] denote a row of LINXLOGX.
                % Now tweak it to ensure at least the equality constraints are satisfied.
                %           UB(i)-LB(i)
                %    X(i) = ---------- c2 ( c1^[(X(j)-LB(j))/(UB(j)-LB(j))] - 1 ) + LB(i)
                %              c3
                for i1 = 1:size(linXlogX,1)
                    ii = linXlogX(i1,1);
                    jj = linXlogX(i1,2);
                    c1 = linXlogX(i1,3);
                    c2 = linXlogX(i1,4);
                    c3 = linXlogX(i1,5);
                    xinit(ii) = [(UB(ii)-LB(ii))/c3]*c2*(c1^[(xinit(jj)-LB(jj))/(UB(jj)-LB(jj))]-1 ) + LB(ii); %#ok
                end
                
            end
            if isempty(opt)
                opt = DClab.DCOptions;
            end
            if isempty(relaxEq)
                relaxEq = false;
            end
            
            % Find feasible value for gamma.
            tmp = zeros(nQuad,1);
            for i1 = 1:nQuad
                tmp(i1) = [1 xinit']*Zquad{i1}*[1; xinit];
            end
            gammainit = max(tmp./tau);
            
            % would we ever want the notify or final messages to splash to the screen?
            if strcmp(opt.display,'ALL')
                dispMode = 'iter';
            else
                dispMode = 'off';
            end
            
            %Determine which optimization case we're dealing with
            % case1: there are no linXlogX equality constraints and all quadratic
            %   constraints are actually liner. use linprog.
            % case2: there are no linXlogX equality constraints, but a constraint is
            %   quadratic. use fmincon.
            % case3: we are supposed to use the quadratic relaxation of the linXlogX
            %   equality constraints. use fmincon
            % case4: we are supposed to enforce the linXlogX equality constraints. use
            %   fmincon.
            
            % Define solveCase, A,B,Aeq,Beq.
            A = [];
            B = [];
            Aeq = [];
            Beq = [];
            if isempty(obj.ZrelaxedEq)
                linIdx = false(nQuad,1);
                for i2 = 1:length(Zquad)
                    if isequal(Zquad{i2}(2:n+1,2:n+1),zeros(n))
                        tmp = [2*Zquad{i2}(1,2:end) -tau(i2)];
                        A = [A; tmp];
                        B = [B; -Zquad{i2}(1,1)];
                        linIdx(i2) = true;
                    end
                end
                if sum(linIdx) == length(Zquad)
                    solveCase = 1;
                else
                    solveCase = 2;
                end
                %remove from the quadratic contraint list those that were actually linear
                Zquad(linIdx) = [];
                tau(linIdx) = [];
            else
                if relaxEq
                    solveCase = 3;
                    linIdx = false(nQuad,1);
                    for i2 = 1:length(Zquad)
                        if isequal(Zquad{i2}(2:n+1,2:n+1),zeros(n))
                            tmp = [2*Zquad{i2}(1,2:end) -tau(i2)];
                            A = [A; tmp];
                            B = [B; -Zquad{i2}(1,1)];
                            linIdx(i2) = true;
                        end
                    end
                    %remove from the quadratic contraint list those that were actually linear
                    Zquad(linIdx) = [];
                    tau(linIdx) = [];
                else
                    solveCase = 4;
                end
            end
            
            switch solveCase
                case 1
                    lpoptions = optimset('linprog');
                    lpoptions = optimset(lpoptions,'Display',dispMode,'TolFun',opt.tolFun,'TolX',opt.tolCon,'TolCon',opt.tolCon);
                    
                    % We need to solve
                    % min \gamma s.t. LB <= x <= UB
                    % Ax <= B + c\gamma.
                    
                    % decision variable = [x;\gamma]
                    
                    % First attempt to solve
                    [x,fval,exitflg,output,lambda] = ...
                        linprog([zeros(n,1); 1],A,B,[],[],[LB;-inf],[UB;inf],[],lpoptions);
                    
                    % Second attempt to solve. If this also fails, we give up.
                    if exitflg <= 0
                        if ismember(opt.display,{'notify','all','ALL'})
                            str = '  LINPROG failed to find a feasible point, retrying';
                            DClab.dcdispstr(str,opt.guiHandle,false)
                        end
                        
                        %Don't use largescale method this time.
                        lpoptions = optimset(lpoptions,'largescale','off');
                        [x,fval,exitflg,output,lambda] = ...
                            linprog([zeros(n,1); 1],A,B,[],[],[LB;-inf],[UB;inf],[xinit;gammainit],lpoptions);
                    end
                    % pull gamma off the end of x;
                    x = x(1:end-1);
                    
                case {2,3}
                    % No special purpose equality constraints to directly enforce.
                    
                    confun = @(x) DClab.cmoptim.fmconfunRelax(x,[Zquad;obj.ZrelaxedEq],[tau; zeros(length(obj.ZrelaxedEq),1)]);
                    objfun = @(x) DClab.cmoptim.fmobjfun(x);
                    
                    % FMINCON options (large scale doesn't work with nonlin constraints)
                    fmoptions = optimset('fmincon');
                    fmoptions = optimset(fmoptions,'GradConstr','on','GradObj','on','TolX',opt.tolCon,...
                        'LargeScale','off','MaxIter',500,'MaxFunEval',12500,...
                        'TolFun',opt.tolFun,'MaxSQPIter',10000,'Display',dispMode,'TolCon',opt.tolFun);
                    
                    % First attempt to solve
                    warning('off','optim:fmincon:NLPAlgLargeScaleConflict');
                    [x,fval,exitflg,output,lambda]=...
                        fmincon(objfun,[xinit;gammainit],A,B,Aeq,Beq,[LB;-inf],[UB;inf],confun,fmoptions);
                    warning('on','optim:fmincon:NLPAlgLargeScaleConflict');
                    
                    % Second attempt to solve. If this also fails, we give up.
                    if exitflg <= 0
                        if ismember(opt.display,{'notify','all','ALL'})
                            str = '  FMINCON failed to find a feasible point, retrying';
                            DClab.dcdispstr(str,opt.guiHandle,false)
                        end
                        
                        % Create a random initial point. Use uniform on [-500, 500] for
                        % unbounded variables.
                        idx = ~isinf(obj.variableBnds(:,1));
                        xinit = 1000*rand(size(xinit)) - 500;
                        xinit(idx) = (UB(idx)-LB(idx)).*rand(sum(idx),1) + LB(idx);
                        
                        % Make sure xinit at least satisfies any special purpose equalities
                        for i1 = 1:size(linXlogX,1)
                            ii = linXlogX(i1,1);
                            jj = linXlogX(i1,2);
                            c1 = linXlogX(i1,3);
                            c2 = linXlogX(i1,4);
                            c3 = linXlogX(i1,5);
                            xinit(ii) = [(UB(ii)-LB(ii))/c3]*c2*(c1^[(xinit(jj)-LB(jj))/(UB(jj)-LB(jj))]-1 ) + LB(ii); %#ok
                        end
                        % Find feasible value for gamma.
                        tmp = zeros(nQuad,1);
                        for i1 = 1:nQuad
                            tmp(i1) = [1 xinit']*Zquad{i1}*[1; xinit];
                        end
                        gammainit = max(tmp./tau);
                        
                        [x,fval,exitflg,output,lambda] = ...
                            fmincon(objfun,[xinit;gammainit],A,B,Aeq,Beq,[LB;-inf],[UB;inf],confun,fmoptions);
                    end
                    x = x(1:end-1);
                case 4
                    %Here the optimization solver is only aware of the smaller decision
                    %vector that results from using the special purpose equality
                    %constraints to reduce the number of optimization variables. However,
                    %in the objective and constraint functions, we recover the entire
                    %decision vector in order to evaluate these.
                    
                    % Define some useful indexing vectors
                    x2XElim = obj.linXlogX(:,2);         %xElim = x(x2XElim)
                    x2XEqFcnXElim = obj.linXlogX(:,1);   %xEqFcnXElim = x(x2XEqFcnXElim)
                    x2XRem = setdiff((1:n)',x2XElim);       %xRem = x(x2XRem)
                    
                    % Determine where in xRem it is related by xElim by the equality
                    % constraint.
                    
                    [trash a2sort b2sort] = intersect(x2XRem,x2XEqFcnXElim);
                    xRem2XEqFcnXElim(b2sort,1) = x2XRem(a2sort);
                    
                    nElim = length(x2XElim);
                    c = zeros(nElim,7);
                    c(:,[1 2 3]) = obj.linXlogX(:,[3 4 5]);
                    c(:,4) = LB(x2XEqFcnXElim);
                    c(:,5) = UB(x2XEqFcnXElim) - LB(x2XEqFcnXElim);
                    c(:,6) = LB(x2XElim);
                    c(:,7) = UB(x2XElim) - LB(x2XElim);
                    
                    % The equality constraint is
                    %
                    % f_c^{-1}(xl) = c7/log10(c1) * log10( c3/(c2*c5)*(xl-c4) + 1 ) + c6
                    %
                    % To reduce the computational time, simplify.
                    % Let t1 = c7/log10(c1),
                    %     t2 = c3/(c2*c5),
                    %     t3 = c4;
                    %     t4 = c6;
                    %
                    % Thus the equality constraint is
                    %
                    % f_c^{-1}(xl) = t1 * log10(t2*(xl-t3) + 1) + t4;
                    
                    t = zeros(nElim,4);
                    t(:,1) = c(:,7)./log10(c(:,1));
                    t(:,2) = c(:,3)./(c(:,2).*c(:,5));
                    t(:,3) = c(:,4);
                    t(:,4) = c(:,6);
                    
                    confun = @(xRem) DClab.cmoptim.fmconfun(xRem,Zquad,tau,x2XRem,x2XElim,xRem2XEqFcnXElim,t);
                    objfun = @(xRem) DClab.cmoptim.fmobjfun(xRem);
                    
                    A=[]; B=[]; %No generalized LP constraints
                    Aeq=[]; Beq=[];  %No equality constraints
                    
                    % FMINCON options (large scale doesn't work with nonlin constraints)
                    fmoptions = optimset('fmincon');
                    fmoptions = optimset(fmoptions,'GradConstr','on','GradObj','on','TolX',opt.tolCon,...
                        'LargeScale','off','MaxIter',500,'MaxFunEval',12500, ...
                        'TolFun',opt.tolFun,'MaxSQPIter',10000,'Display',dispMode,'TolCon',opt.tolFun);
                    
                    %First attempt to solve.
                    state = warning('off','optim:fmincon:NLPAlgLargeScaleConflict');
                    [xRem,fval,exitflg,output,lambda]=...
                        fmincon(objfun,[xinit(x2XRem);gammainit],A,B,Aeq,Beq,[LB(x2XRem);-inf],[UB(x2XRem);inf],confun,fmoptions);
                    warning(state);
                    
                    %Second attempt to solve
                    if exitflg <= 0
                        if ismember(opt.display,{'notify','all','ALL'})
                            str = '  fmincon failed to find a feasible point, retrying';
                            DClab.dcdispstr(str,opt.guiHandle,false)
                        end
                        
                        % Create a random initial point. Use uniform on [-500, 500] for
                        % unbounded variables.
                        idx = ~isinf(obj.variableBnds(:,1));
                        xinit = 1000*rand(size(xinit)) - 500;
                        xinit(idx) = (UB(idx)-LB(idx)).*rand(sum(idx),1) + LB(idx);
                        
                        % Make sure xinit at least satisfies any special purpose equalities
                        for i1 = 1:size(linXlogX,1)
                            ii = linXlogX(i1,1);
                            jj = linXlogX(i1,2);
                            c1 = linXlogX(i1,3);
                            c2 = linXlogX(i1,4);
                            c3 = linXlogX(i1,5);
                            xinit(ii) = [(UB(ii)-LB(ii))/c3]*c2*(c1^[(xinit(jj)-LB(jj))/(UB(jj)-LB(jj))]-1 ) + LB(ii); %#ok
                        end
                        % Find feasible value for gamma.
                        tmp = zeros(nQuad,1);
                        for i1 = 1:nQuad
                            tmp(i1) = [1 xinit']*Zquad{i1}*[1; xinit];
                        end
                        gammainit = max(tmp./tau);
                        
                        [xRem,fval,exitflg,output,lambda]=...
                            fmincon(objfun,[xinit(x2XRem);gammainit],A,B,Aeq,Beq,[LB(x2XRem);-inf],[UB(x2XRem);inf],confun,fmoptions);
                    end
                    
                    %Recover full decision vector.
                    x = zeros(n,1);
                    x(x2XRem) = xRem(1:end-1); %cut off gamma variable.
                    x(x2XElim) = t(:,1).*log10( t(:,2).*(xRem(xRem2XEqFcnXElim)-t(:,3)) + 1 ) + t(:,4);
                    
            end %switch solveCase
            
            %if still infeas, return Inf. Set lambda=[] so we can catch this later
            if exitflg <= 0
                if ismember(opt.display,{'notify','all','ALL'})
                    str = '  FMINCON failed to find a feasible point, returning infinite bound';
                    DClab.dcdispstr(str,opt.guiHandle,false)
                end
                ub = inf;
                xopt = x; %still return x because fmincon tries to minimize
                %constraint violation on infeas problems. this x "may" be in the true feasible set
                lambda = [];
            else
                % Output information
                ub = fval;
                xopt = x;
            end
            
            %parse lambda to create the desired multiplier structure.
            
            if isempty(lambda)
                %solution bombed return all zeros
                mult.lower = zeros(n,1);
                mult.upper = zeros(n,1);
                mult.quadineq = zeros(nQuad,1);
                mult.relaxedeq = zeros(nRelaxEq,1);
                
            else
                switch solveCase
                    case 1
                        mult.lower = lambda.lower(1:end-1);
                        mult.upper = lambda.upper(1:end-1);
                        mult.quadineq = lambda.ineqlin;
                        mult.relaxedeq = zeros(0,1);
                    case 2
                        mult.lower = lambda.lower(1:end-1);
                        mult.upper = lambda.upper(1:end-1);
                        mult.quadineq = zeros(nQuad,1);
                        mult.quadineq(linIdx) = lambda.ineqlin;
                        mult.quadineq(~linIdx) = lambda.ineqnonlin;
                        mult.relaxedeq = zeros(0,1);
                    case 3
                        mult.lower = lambda.lower(1:end-1);
                        mult.upper = lambda.upper(1:end-1);
                        mult.quadineq = zeros(nQuad,1);
                        mult.quadineq(linIdx) = lambda.ineqlin;
                        nQuadAsQuad = sum(~linIdx);
                        mult.quadineq(~linIdx) = lambda.ineqnonlin(1:nQuadAsQuad);
                        mult.relaxedeq = lambda.ineqnonlin(nQuadAsQuad+1:end);
                    case 4
                        mult.lower = zeros(n,1);
                        mult.lower(x2XRem) = lambda.lower(1:end-1);
                        mult.upper = zeros(n,1);
                        mult.upper(x2XRem) = lambda.upper(1:end-1);
                        mult.quadineq = lambda.ineqnonlin;
                        mult.relaxedeq = zeros(nRelaxEq,1);
                end
            end
            
            %TODO, for times when exitflg==0, max iterations have occured, but the
            %returned point might be feasible. If so, we should return non inf values.
            
            
            %       %===Solver could not find a solution. Get user assistance===
            %
            %       if exitflg <= 0 && strcmp(solver,'fmincon')
            %         disp(['fmincon failed to find a soln, adjust options and run the following to ' ...
            %               'retry']);
            %         disp(['[x,fval,exitflg,output,multipliers]= fmincon' ...
            %           '(@fmobjfun,xinit{i1},A,B,Aeq,Beq,LB,UB,@fmconfun,options,Z0,Zi)']);
            %         disp('Alternatively, type return to terminate with no solution found')
            %         x_ = x; fval_ = fval; lambda_ = lambda;
            %         fval = -inf; x = []; lambda = [];
            %         keyboard
            %       end
            %
            %       if exitflg <= 0 && strcmp(solver,'linprog')
            %         disp(['linprog failed to find a soln, adjust options and run the following to ' ...
            %               'retry']);
            %         disp(['[x,fval,exitflg,output,lambda]=' ...
            %           'linprog(2*Z0(1,2:end),A,B,[],[],LB,UB,[],options)']);
            %         disp('Alternatively, type return to terminate with no solution found')
            %         x_ = x; fval_ = fval; lambda_ = lambda;
            %         fval = -inf; x = []; lambda = [];
            %         keyboard
            %       end
            %       %===End user assistance===
        end
        
        
    end %private,static methods
end %classdef
