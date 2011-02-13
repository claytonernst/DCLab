classdef nqcqp < DClab.DCObject
    %NQCQP constructs an object that represents an optimization problem
    %
    %   NQCQP sets up an optimization of the form
    %
    %   min [1;X]^T*ZNOT*[1;X] subj. to: [1;X]^T*ZQUAD{i}*[1;X] <= 0 i=1,...,m
    %    x                               LB <= X <= UB
    %                                    equality constraints of a special form
    %
    %   OBJ = NQCQP(ZNOT,ZQUAD) sets up a problem in n decision variables.
    %   ZNOT should be a sparce symmetric (n+1)-by-(n+1) matrix. ZQUAD should
    %   be an m-by-1 cell array of sparce symmetric (n+1)-by-(n+1) matrices.
    %   Set ZQUAD=[] if no quadratic constraints exist.
    %
    %   OBJ = NQCQP(ZNOT,ZQUAD,LB,UB) defines a set of lower and upper
    %   bounds on the design variables, X, so that a solution is found in
    %   the range LB <= X <= UB. Limitations: if one of LB and UB is empty, the
    %   other must be as well. If LB(i) = -Inf, UB(i) must be Inf, and
    %   conversely.
    %
    %   OBJ = NQCQP(ZNOT,ZQUAD,LB,UB,LINXLOGX) provides for additional
    %   equality constraints of a special form. These can only be used to
    %   constrain variables whose range (as specified by LB and UB) is a
    %   bounded interval. LINXLOGX should be p-by-5. Let [i j c1 c2 c3] denote
    %   a row of LINXLOGX. This row imposes the equality constraint
    %           UB(i)-LB(i)
    %    X(i) = ---------- c2 ( c1^[(X(j)-LB(j))/(UB(j)-LB(j))] - 1 ) + LB(i)
    %              c3
    %   This equality constraint in relaxed using quadratic inequalities by the
    %   lowerBnd method. The upperBnd method handles it directly.
    %
    %   OBJ = NQCQP(ZNOT,ZQUAD,LB,UB,LINXLOGX,KNOWNFEAS) KNOWNFEAS==1
    %   indicates the constraint set is known to be nonempty. KNOWNFEAS==0
    %   indicates this is unknown for the present problem. Setting KNOWNFEAS=1
    %   when the constraint set is known to be nonempty reduces computation
    %   time, as certain routines are bypassed. (Set LB=[] and UB=[] if no
    %   bounds exist; set LINXLOGX=[] if no equality constraints exist.)
    
    properties
        Znot;
        Zquad;
        ZrelaxedEq;
        nVars;
        variableBnds;
        linXlogX;
        knownFeas;
        sProcInfeas;
        xfeas;
    end
    
    methods
        function obj = nqcqp(Znot,Zquad,LB,UB,linXlogX,knownFeas)
            
            error(nargoutchk(0,1,nargout))
            ni = nargin;
            
            switch ni
                case 0
                    return
                case 2
                    LB = [];
                    UB = [];
                    linXlogX = [];
                    knownFeas = 0;
                case 4
                    linXlogX = [];
                    knownFeas = 0;
                case 5
                    knownFeas = 0;
                case 6
                    %do nothing
                otherwise
                    error('Incorrect number of input arguments.')
            end
            
            % ===input checking===
            
            % objective function:
            assert(~isempty(Znot),'Inputs: ZNOT cannot be empty.')
            if ~DClab.isrealsymmatrix(Znot)
                Znot = 0.5*(Znot+Znot');
            end
            n = size(Znot,1)-1;
            
            % constraint functions:
            if isempty(Zquad)
                Zquad = cell(0,1);
            end
            assert(iscell(Zquad) && size(Zquad,2)==1,'Inputs: ZQUAD must be a column cell array')
            assert(all(cellfun('prodofsize',Zquad) == (n+1)^2),'Inputs: ZQUAD must contain matrices of the same dimension as ZNOT')
            for i1 = 1:length(Zquad)
                if ~DClab.isrealsymmatrix(Zquad{i1})
                    Zquad{i1} = 0.5*(Zquad{i1}+Zquad{i1}');
                end
            end
            
            % variable bounds
            assert(~xor(isempty(LB),isempty(UB)),'Inputs: when one of LB,UB is empty, the other must be as well.')
            if isempty(UB) && isempty(LB)
                LB = repmat(-inf,n,1);
                UB = repmat(inf,n,1);
            else
                assert(isequal(size(LB),[n 1]) && isequal(size(UB),[n 1]),'Inputs: when nonempty, LB and UB must be n-by-1 vectors.')
                
                %Make sure the range of each variable is (-Inf,Inf) or a nonempty
                %bounded interval.
                assert(~any(LB==Inf) && ~any(UB==-Inf),'Inputs: LB(i) == Inf or UB(i) == -Inf is not permitted.')
                
                infLB = LB==-Inf;
                infUB = UB==Inf;
                
                assert(isequal(infLB,infUB),'Inputs: if LB(i)==-Inf UB(i) must equal Inf, and conversely.')
                assert(all(UB(~infUB) - LB(~infLB) >= 0),'Inputs: box constraint implied by LB and UB is empty.')

            end
            
            % linXlogX bounds
            if isempty(linXlogX)
                linXlogX = [];
            else
                assert(size(linXlogX,2) == 5,'Inputs: size(LINXLOGX,2) must equal 5.')
                assert(isequal(round(linXlogX(:,[1 2])),linXlogX(:,[1 2])),'Inputs: 1st and 2nd columns of LINXLOGX must be integers.')
                assert(min(min(linXlogX(:,[1 2]))) >= 1 && max(max(linXlogX(:,[1 2]))) <= n,'Inputs: 1st and 2nd columns of LINXLOGX must be members of {1,...,n}.')
                assert(~any(isinf(LB(linXlogX(:,1)))) && ~any(isinf(LB(linXlogX(:,2)))),'Inputs: LINXLOGX can only impose constraints on variables whose ranges are bounded by LB and UB')
            end
            
            % knownFeas
            assert(knownFeas == 0 || knownFeas ==1,'Inputs: KNOWNFEAS must be 0 or 1.')

            
            % Create quadratic relaxations of any equality constraints. We employ four
            % quadratic constraints per equality constraint.
            if isempty(linXlogX)
                Z4relaxedEqualities = [];
            else
                p = size(linXlogX,1);
                Z4relaxedEqualities = cell(4*p,1);
                for i1 = 1:size(linXlogX,1)
                    Z4relaxedEqualities(4*i1-3:4*i1) = DClab.nqcqp.buildRelaxedEqualities(linXlogX(i1,:),LB,UB);
                end
            end
            
            obj.Znot = Znot;
            obj.Zquad = Zquad;
            obj.ZrelaxedEq = Z4relaxedEqualities;
            obj.nVars = n;
            obj.variableBnds = [LB UB];
            obj.linXlogX = linXlogX;
            obj.knownFeas = knownFeas;
            if knownFeas == 1 || isempty(Zquad)
                obj.sProcInfeas = 0;
            else
                obj.sProcInfeas = [];
            end
            obj.xfeas = [];
            
        end
        
        function [ub,xopt,mult] = upperBnd(varargin)
            %UPPERBND method of nqcqp produces an upper bound on the optimal value
            %
            %   An nqcqp object describes an optimization of a quadratic objective
            %   function Znot subject to box constraints, quadratic inequality
            %   constraints, and special purpose linXlogX equality constraints. This
            %   function determines an upper bound on the optimal value by calling
            %   FMINCON. Future functionality includes providing the capability to use
            %   additional solvers.
            %
            %   UB = UPPERBND(NQCQPOBJ) produces an upper bound UB on the optimal value
            %   of the optimization stored in NQCQPOBJ. The solver will be initialized
            %   at a random point drawn from a uniform distribution over the rectanbe
            %   described by the box constraints. If a feasible point cannot be found,
            %   UB will be +Inf.
            %
            %   UB = UPPERBND(NQCQPOBJ,XINIT) initializes the solver at the n-by-1
            %   vector XINIT.
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
            %   decision vector XOPT that produces UB.
            %
            %   [UB XOPT MULT] = UPPERBND(...) additionally returns a structure of
            %   Lagrange multipliers LAMBDA with fields .lower, .upper, .quadineq,
            %   .relaxedeq. All elements of .relaxedeq will be 0 if RELAXEQ==FALSE
            
            % This function is just a wrapper. Choose the desired solver
            NONE = 0;
            NPSOL =1;
            FMINCON=2;
            
            SOLVER=FMINCON;
            switch SOLVER
                case NONE
                    ub=Inf; xopt=[]; mult=[];
                case NPSOL
                    [ub,xopt,mult] = qpub_NPSOL(varargin{:});
                case FMINCON
                    [ub,xopt,mult] = qpub_FMINCON(varargin{:});
            end
        end
        
        function [lb,lambda,Q] = lowerBnd(nqcqpObj,opt)
            %LOWERBND method of nqcqp produces a lower bound on the optimal value
            %
            %   An nqcqp object describes an optimization of a quadratic objective
            %   function Znot subject to box constraints, quadratic inequality
            %   constraints, and special purpose linXlogX equality constraints. The
            %   object also includes quadratic inequality constraints ZrelaxedEq that
            %   are a "good" relaxation of the equality constraints. This function
            %   determines a lower bound on the optimal value. First each dimension of
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
            %   LB = LOWERBND(NQCQPOBJ) produces a lower bound LB on the optimal value
            %   of the optimization stored in NQCQPOBJ. If the problem is infeasible,
            %   LB will be +Inf.
            %
            %   LB = LOWERBND(NQCQPOBJ,OPT) uses options taken from the DCOptions
            %   object OPT. This method used OPT.display and OPT.sedumiParEps.
            %
            %   [LB, LAMBDA] = LOWERBND(...) additionally returns a structure array
            %   LAMBDA with fields .lower, .upper, .quadineq, .relaxedeq. Each field
            %   contains the S-procedure multipliers lambda for the corresponding
            %   constraints.
            %
            %   [LB, LAMBDA, Q] = LOWERBND(...) additionally returns the optimal
            %   (n+1)x(n+1) matrix  Q from the rank relaxation.
            
            
            %TODO use linprog if everything is linear
            
            error(nargchk(1,2,nargin));
            error(nargoutchk(0,3,nargout));
            
            % Initialize inputs
            error(nargchk(1,2,nargin));
            if nargin==1
                opt = DClab.DCOptions;
            end
            
            % Check if S-procedure can prove the constraints are infeasible and
            % populate the object with this info.
            nqcqpObj = checkSProcInfeas(nqcqpObj,opt);
            
            if nqcqpObj.sProcInfeas
                lb = inf;
                lambda = [];
                Q = [];
            else
                % Problem Data
                n = nqcqpObj.nVars;
                Znot = nqcqpObj.Znot;
                
                %convert variableBnds to quadratic constraints
                xbnds = nqcqpObj.variableBnds;
                
                % by requirements of the object, they are either both finite or both
                % infinite.
                finiteXBnds = find(~isinf(xbnds(:,1)));
                epsilon = zeros(n,1);
                affineTerm = zeros(n,1);
                linTerm = zeros(n,1);
                
                % Create quadratic representations of the box constraints
                
                % Constraint LB_i <= x_i
                nn = length(finiteXBnds);
                Zbox = cell(2*nn,1);
                [Zbox{:}] = deal(spalloc(n+1,n+1,4));
                epsilon(finiteXBnds) = 0.05*(xbnds(finiteXBnds,2)-xbnds(finiteXBnds,1));
                affineTerm(finiteXBnds) = xbnds(finiteXBnds,1).* (xbnds(finiteXBnds,2)+epsilon(finiteXBnds));
                linTerm(finiteXBnds) = -(xbnds(finiteXBnds,1) + xbnds(finiteXBnds,2) + epsilon(finiteXBnds));
                for i2 = 1:length(finiteXBnds)
                    i3 = finiteXBnds(i2);
                    Zbox{i2}(1,1) = affineTerm(i3);
                    Zbox{i2}(i3+1,1) = linTerm(i3)/2;
                    Zbox{i2}(1,i3+1) = linTerm(i3)/2;
                    Zbox{i2}(i3+1,i3+1) = 1;
                end
                
                % Constraint x_i <= UB_i
                affineTerm(finiteXBnds) = xbnds(finiteXBnds,2).* (xbnds(finiteXBnds,1)-epsilon(finiteXBnds));
                linTerm(finiteXBnds) = -(xbnds(finiteXBnds,1) + xbnds(finiteXBnds,2) - epsilon(finiteXBnds));
                %for i2 = 1:finiteXBnds'
                for i2 = 1:length(finiteXBnds)
                    i3 = finiteXBnds(i2);
                    Zbox{i2+nn}(1,1) = affineTerm(i3);
                    Zbox{i2+nn}(i3+1,1) = linTerm(i3)/2;
                    Zbox{i2+nn}(1,i3+1) = linTerm(i3)/2;
                    Zbox{i2+nn}(i3+1,i3+1) = 1;
                end
                
                % Old version
                %  [Zbox{:}] = deal(spalloc(n+1,n+1,4));
                %  affineTerm = xbnds(:,1).* xbnds(:,2);
                %  linTerm = -xbnds(:,1) - xbnds(:,2);
                %  for i2 = 1:nn
                %    Zbox{i2}(1,1) = affineTerm(i2);
                %    Zbox{i2}(i2+1,1) = linTerm(i2)/2;
                %    Zbox{i2}(1,i2+1) = linTerm(i2)/2;
                %    Zbox{i2}(i2+1,i2+1) = 1;
                %  end
                
                Zconst = [Zbox; nqcqpObj.Zquad; nqcqpObj.ZrelaxedEq];
                
                N = length(Zconst);
                
                %TODO, go through sedumi help and document the sedumi call better. Maybe
                %even in pdf.
                
                % Sedumi (dual-form) notation:
                %     max b'y  s.t. c-A'*y >= 0
                % Define y:=[lambda_1 ... lambda_N gamma]
                b = spalloc(N+1,1,1); b(end) = 1;
                
                % N linear constraints:
                %  Alin*y <= clin ensures lambda >=0
                K.l = N;
                clin = spalloc(N,1,0);
                Alin = [speye(N) spalloc(N,1,0)];
                
                % (n+1)x(n+1) LMI constraint:
                %    Z0 - [gamma 0; 0 0] + sum_k lambda_k Z_k >= 0
                K.s = n+1;
                clmi = sparse(Znot(:));
                Almi = [];
                if N > 0
                    Almi = zeros(numel(Zconst{1}),N);
                    for i2 = 1:N
                        Almi(:,i2) = Zconst{i2}(:);
                    end
                end
                
                temp = spalloc(n+1,n+1,1); temp(1,1) = 1;
                Almi = sparse([Almi -temp(:)]);
                
                % Concatenate the constraints
                c = [clin; clmi];
                At = -[Alin; Almi];
                
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
                Nbox = size(Zbox,1);
                Nquad = size(nqcqpObj.Zquad,1);
                
                lambda.lower = zeros(n,1);
                lambda.lower(finiteXBnds) = y(1:Nbox/2);
                lambda.upper = zeros(n,1);
                lambda.upper(finiteXBnds) = y(Nbox/2+1:Nbox);
                lambda.quadineq = y(Nbox+1:Nbox+Nquad);
                lambda.relaxedeq = y(Nbox+Nquad+1:end-1);
                
                Q = reshape(x(N+1:end),K.s,K.s);
                
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
                
            end %if/else nqcqpObj.sProcInfeas
        end
        
        function [list,sz] = displayProps(obj)
            list = [];
            sz = '';
        end
        
    end %public methods
    
    methods (Access=private,Static)
        function Zcell = buildRelaxedEqualities(linXlogX,LB,UB)
            %BUILDRELAXEDEQUALITIES creates a quad relaxation of a special equality constraint.
            %
            %   Inputs: linXlogX should be 1-by-5
            %           LB and UB should be n-by-1
            %
            %   Outputs: Zcell will be a 4-by-1 cell array of sparse symmetric
            %   (n+1)-by-(n+1) matrices. Each imposes the constraint
            %   [1 x]^T*Zcell{i}*[1;x] <= 0
            %
            
            ii = linXlogX(1);
            jj = linXlogX(2);
            
            % xl = x(ii)
            % xg = x(jj)
            
            c1 = linXlogX(3); % = b/a
            c2 = linXlogX(4); % = a
            c3 = linXlogX(5); % = b - a
            c4 = LB(ii);
            c5 = UB(ii) - LB(ii);
            c6 = LB(jj);
            c7 = UB(jj) - LB(jj);
            
            % The equality constraint we seek to enforce is
            %
            % xl = f_c(xg) = c5*c2/c3 * [c1^((xg-c6)/c7) - 1] + c4
            %
            % or equivalently
            %
            % xg = f_c^{-1}(xl) = c7/log10(c1) * log10( c3/(c2*c5)*(xl-c4) + 1 ) + c6
            
            N = 100; % number of grid points.
            n = size(LB,1); % number of optimization variables.
            
            % === determine q ===
            
            % solve min \gamma subject to: for all xg \in [LB(jj),UB(jj)]
            % f_c(xg) - 0.5\gamma <= q(xg) <= f_c(xg) + 0.5\gamma
            
            % let q(xg) = c0 + c1*xg + c2*xg^2
            
            xg = linspace(LB(jj),UB(jj),N)';
            fxg = c5*c2/c3 * ( c1.^((xg-c6)./c7) - 1 ) + c4;
            
            % X=LINPROG(f,A,b) solves  min f'*X    subject to:   A*X <= b
            %                           X
            
            % let X = [\gamma c0 c1 c2]'
            f = [1 0 0 0]';
            
            % f_c(xg) - 0.5\gamma <= q(xg) constraint
            A1 = [-0.5*ones(N,1) -ones(N,1) -xg -xg.^2];
            b1 = -fxg;
            
            % q(xg) <= f_c(xg) + 0.5\gamma constraint
            A2 = [-0.5*ones(N,1) ones(N,1) xg xg.^2];
            b2 = fxg;
            
            %Although this problem is alway feasible, i've had linprog fail on it when
            %using the largescale method.
            opt = optimset('linprog');
            opt = optimset(opt,'largescale','off','display','off');
            
            [X trash exitflg] = linprog(f,[A1;A2],[b1;b2],[],[],[],[],[],opt);
            
            if exitflg < 1
                disp('failed call to linprog')
                keyboard
            end
            
            gamma = X(1);
            
            % Create q(xl,xg) = q(xg) - xl
            %                 = c0 + c10*xl + c20*xg + c11*xl^2 + c12*xl*xg + c22*xg^2
            coeff = [X(2) -1 X(3) 0 0 X(4)];
            q = DClab.coeff2quadform(coeff,2);
            
            % === determine p ===
            
            % solve min \delta subject to: for all xl \in [LB(ii),UB(ii)]
            % f^{-1}_c(xl) - 0.5\delta <= p(xl) <= f^{-1}_c(xl) + 0.5\delta
            
            % let p(xl) = c0 + c1*xl + c2*xl^2
            
            xl = linspace(LB(ii),UB(ii),N)';
            
            finvxl = c7/log10(c1) * log10( c3/(c2*c5)*(xl-c4) + 1 ) + c6;
            
            % let X = [\delta c0 c1 c2]'
            f = [1 0 0 0]';
            
            % f^{-1}_c(xl) - 0.5\delta <= p(xl) constraint
            A1 = [-0.5*ones(N,1) -ones(N,1) -xl -xl.^2];
            b1 = -finvxl;
            
            % p(xl) <= f^{-1}_c(xl) + 0.5\delta constraint
            A2 = [-0.5*ones(N,1) ones(N,1) xl xl.^2];
            b2 = finvxl;
            
            [X trash exitflg] = linprog(f,[A1;A2],[b1;b2],[],[],[],[],[],opt);
            if exitflg < 1
                disp('failed call to linprog')
                keyboard
            end
            
            delta = X(1);
            
            % Create p(xl,xg) = p(xl) - xg
            %                 = c0 + c10*xl + c20*xg + c11*xl^2 + c12*xl*xg + c22*xg^2
            coeff = [X(2) X(3) -1 X(4) 0 0];
            p = DClab.coeff2quadform(coeff,2);
            
            % === create the four quadratic constraints ===
            Zcell = cell(4,1);
            
            xlIDX = find(1:n == ii);
            xgIDX = find(1:n == jj);
            
            % -0.575 gamma <= q(xl,xg) <= 0.575 gamma (use .575 instead of 0.5 for 15% fudge
            
            tempL = spalloc(n+1,n+1,3^2); %constraint from -gamma/2 <= q(x)
            tempR = tempL; %constraint from  q(x) <= gamma/2
            
            tempL([1;xlIDX+1;xgIDX+1],[1;xlIDX+1;xgIDX+1]) = -q;
            tempL(1,1) = tempL(1,1) - 0.575*gamma;
            
            tempR([1;xlIDX+1;xgIDX+1],[1;xlIDX+1;xgIDX+1]) = q;
            tempR(1,1) = tempR(1,1) - 0.575*gamma;
            
            Zcell{1} = tempL;
            Zcell{2} = tempR;
            
            % -0.575 delta <= p(xl,xg) <= 0.575 delta (use .575 instead of 0.5 for 15% fudge
            
            tempL = spalloc(n+1,n+1,3^2); %constraint from -gamma/2 <= q(x)
            tempR = tempL; %constraint from  q(x) <= gamma/2
            
            tempL([1;xlIDX+1;xgIDX+1],[1;xlIDX+1;xgIDX+1]) = -p;
            tempL(1,1) = tempL(1,1) - 0.575*delta;
            
            tempR([1;xlIDX+1;xgIDX+1],[1;xlIDX+1;xgIDX+1]) = p;
            tempR(1,1) = tempR(1,1) - 0.575*delta;
            
            Zcell{3} = tempL;
            Zcell{4} = tempR;
            
        end
        
        function [iconval,econval,icongrd,econgrd] = fmconfun(xRem,Zquad,x2XRem,x2XElim,xRem2XEqFcnXElim,t)
            %   This function evaluates each of the quadratic inequality constraints and
            %   computes their gradiants.
            %
            %   Inputs:
            %   xRem: the decision vector with x^g from any equality constraints
            %      eliminated.
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
            nxr = length(xRem);
            
            if ~isempty(Zquad)
                m = length(Zquad);
                n = size(Zquad{1},1)-1;
                
                iconval = zeros(m,1);
                icongrd = zeros(nxr,m);
                
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
                        iconval(i1) = xx'*Zquad{i1}*xx;
                    end
                else
                    for i1 = 1:m
                        tmp1 = zeros(nxr,n+1);
                        tmp1(:,x2XRem+1) = eye(nxr);
                        tmp2 = zeros(nxr,n+1);
                        dirFcInv = (1/log(10))*t(:,1).*t(:,2)./(t(:,2).*(xRem(xRem2XEqFcnXElim)-t(:,3))+1);
                        tmp2(xRem2XEqFcnXElim,x2XElim+1) = diag(dirFcInv);
                        
                        tmp3 = Zquad{i1}*xx;
                        
                        icongrd(:,i1) = 2*(tmp1+tmp2)*tmp3;
                        iconval(i1) = xx'*tmp3;
                    end
                    
                end
                
            else
                
                %TODO check this case. Will the function even be called? What dimensions
                %does matlab expect for the output? How does it determine this?
                iconval = -1; %anything feasible will do.
                icongrd = zeros(nxr,1);
            end
            
            % No Equality constraints
            econval = 0;
            econgrd = zeros(nxr,1);
            
            
            
        end
        
        function [iconval,econval,icongrd,econgrd] = fmconfunRelax(x,Zi)
            
            % Compute inequality constraint functions and and their gradients
            
            if ~isempty(Zi)
                N = length(Zi);
                n = length(x);
                
                iconval = zeros(N,1);
                icongrd = zeros(n,N);
                
                xx = [1;x];
                for i1=1:N
                    tmp = Zi{i1}*xx;
                    icongrd(:,i1) = 2*tmp(2:end); % = 2*grad(xx)*Zi*xx = 2*[zeros(n,1) eye(n)]*Zi*xx
                    iconval(i1) = xx'*tmp;
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
        
        function [objval,objgrd] = fmobjfun(xRem,Znot,x2XRem,x2XElim,xRem2XEqFcnXElim,t)
            %   This function evaluates the quadratic objective function and computes
            %   its gradient.
            %
            %   Inputs:
            %   xRem: the decision vector with x^g from any equality constraints
            %      eliminated.
            %   Znot: a sparse symmetric matrix indicating the objective function.
            %   x2XRem: location in x of the elements of the smaller decision vector
            %      xRem.
            %   x2XElim: location in x of the elements that were eliminated through the
            %      equality constraints.
            %   xRem2XEqFncXElim: location in xRem of the variables that are related
            %      through the equality constraint to the variables eliminated from x.
            %   t: Neq-by-4 matrix that contains the coefficients for the
            %      transformations relating x(x2XElim) and xRem(xRem2XEqFcnXElim).
            
            % Compute the objective function and its gradient.
            nxr = length(xRem);
            
            n = size(Znot,1)-1;
            
            % Recover the full decision vector
            x = zeros(n,1);
            x(x2XRem) = xRem;
            
            % Using the simplications described in qpub_FMINCON, the equality
            % constraint is
            %
            % xg = f_c^{-1}(xl) = t1 * log10(t2*(xl-t3) + 1) + t4;
            
            x(x2XElim) = t(:,1)*log10( t(:,2).*(xRem(xRem2XEqFcnXElim)-t(:,3)) + 1 ) + t(:,4);
            
            xx = [1;x];
            % Determine if we need to compute the gradient.
            if nargout <= 1
                objval = xx'*Znot*xx;
                objgrd = zeros(nxr,1);
            else
                tmp1 = zeros(nxr,n+1);
                tmp1(:,x2XRem+1) = eye(nxr);
                tmp2 = zeros(nxr,n+1);
                dirFcInv = (1/log(10))*t(:,1).*t(:,2)./(t(:,2).*(xRem(xRem2XEqFcnXElim)-t(:,3))+1);
                tmp2(xRem2XEqFcnXElim,x2XElim+1) = diag(dirFcInv);
                
                tmp3 = Znot*xx;
                
                objgrd = 2*(tmp1+tmp2)*tmp3;
                objval = xx'*tmp3;
            end
        end
        
        function [objval,objgrd] = fmobjfunRelax(x,Znot)
            
            
            % Compute objective function value and
            % its gradient
            
            xx = [1; x];
            
            tmp = Znot*xx;
            objgrd = 2*tmp(2:end);   % = 2*grad(xx)*Znot*xx = 2*[zeros(n,1) eye(n)]*Znot*xx
            objval = xx'*tmp;
        end
    end
end %classdef




