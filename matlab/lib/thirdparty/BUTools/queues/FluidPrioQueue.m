%  Ret = FluidPrioQueue(Q, Rin, d, ...)
%  
%  Returns various performane measures of a continuous time 
%  fluid priority queue, see [1]_.
%  
%  Parameters
%  ----------
%  Q : matrix of shape (N,N)
%      The generator of the Markov chain modulating the arrival process
%  Rin : matrix of shape (K,N)
%      The matrix defining the input fluid rates in various states of the 
%      background process of all fluid types
%  d : real number
%      The state independent fluid service rate
%  further parameters : 
%      The rest of the function parameters specify the options
%      and the performance measures to be computed.
%  
%      The supported performance measures and options in this 
%      function are:
%  
%      +----------------+--------------------+----------------------------------------+
%      | Parameter name | Input parameters   | Output                                 |
%      +================+====================+========================================+
%      | "flMoms"       | Number of moments  | The moments of the fluid level         |
%      +----------------+--------------------+----------------------------------------+
%      | "flDistr"      | A vector of points | The fluid level distribution at the    |
%      |                |                    | requested points                       |
%      +----------------+--------------------+----------------------------------------+
%      | "stMoms"       | Number of moments  | The sojourn time moments of fluid      | 
%      |                |                    | drops                                  |
%      +----------------+--------------------+----------------------------------------+
%      | "stDistr"      | A vector of points | The sojourn time distribution at the   |
%      |                |                    | requested points (cummulative, cdf)    |
%      +----------------+--------------------+----------------------------------------+
%      | "prec"         | The precision      | Numerical precision used as a stopping |
%      |                |                    | condition when solving the Riccati and |
%      |                |                    | the matrix-quadratic equations         |
%      +----------------+--------------------+----------------------------------------+
%      | "erlMaxOrder"  | Integer number     | The maximal Erlang order used in the   |
%      |                |                    | erlangization procedure. The default   |
%      |                |                    | value is 200.                          |
%      +----------------+--------------------+----------------------------------------+
%      | "classes"      | Vector of integers | Only the performance measures          |
%      |                |                    | belonging to these classes are         |
%      |                |                    | returned. If not given, all classes    |
%      |                |                    | are analyzed.                          |
%      +----------------+--------------------+----------------------------------------+
%      
%  Returns
%  -------
%  Ret : list of the performance measures
%      Each entry of the list corresponds to a performance 
%      measure requested. Each entry is a matrix, where the
%      columns belong to the various job types.
%      If there is just a single item, 
%      then it is not put into a list.
%  
%  References
%  ----------
%  .. [1] G. Horvath, "Efficient analysis of the MMAP[K]/PH[K]/1
%         priority queue", European Journal of Operational 
%         Research, 246(1), 128-139, 2015.

function varargout = FluidPrioQueue(Q, R, d, varargin)
    
    K = size(R,1);

    % parse options
    erlMaxOrder = 200;
    precision = 1e-14;
    classes = 1:K;
    eaten = [];
    for i=1:length(varargin)
        if strcmp(varargin{i},'erlMaxOrder')
            erlMaxOrder = varargin{i+1};
            eaten = [eaten, i, i+1];
        elseif strcmp(varargin{i},'prec')
            precision = varargin{i+1};
            eaten = [eaten, i, i+1];
        elseif strcmp(varargin{i},'classes')
            classes = varargin{i+1};
            eaten = [eaten, i, i+1];
        end
    end

    global BuToolsCheckInput;

    if isempty(BuToolsCheckInput)
        BuToolsCheckInput = true;
    end   

    if BuToolsCheckInput && ~CheckGenerator(Q)
        error('FluidPrioQueue: Matrix Q is not a valid continuous time generator matrix!');
    end
    
    if BuToolsCheckInput && d<precision
        error('FluidPrioQueue: The fluid service rate must be positive!');
    end
    
    if BuToolsCheckInput
        for k=1:K
            if ~all(all(R>=-precision))
                error('FluidPrioQueue: The fluid arrival rate can not be negative!');
            end
        end
    end

    % Auxiliary functions
    % ===================
    
    % calculates the nth derivative of inv(v*R-Q), even if R contains zero elements
    function dr = DReward (Q, R, n)
        NQ = size(Q,1);
        ix = (1:NQ);
        ixz = ix(abs(diag(R))<=precision);
        ixp = [ix(diag(R)>precision), ix(diag(R)<-precision)];
        Nz = length(ixz);
        Np = length(ixp);
        Per = zeros(NQ);
        for iw=1:Nz
            Per(iw,ixz(iw))=1;
        end
        for iw=1:Np
            Per(Nz+iw,ixp(iw))=1;
        end
        iPer = inv(Per);
        Rp = R(ixp,ixp);
        Qpp = Q(ixp,ixp);
        Qpz = Q(ixp,ixz);
        Qzp = Q(ixz,ixp);
        Qzz = Q(ixz,ixz);
        dXvn = (-1)^n * factorial(n) * inv(inv(Rp)*(-Qpp-Qpz*inv(-Qzz)*Qzp))^(n+1) * inv(Rp);
        drpar = [inv(-Qzz)*Qzp*dXvn*Qpz*inv(-Qzz), inv(-Qzz)*Qzp*dXvn; dXvn*Qpz*inv(-Qzz), dXvn];
        dr = iPer*drpar*Per;
    end
    
    function BPM = BusyPeriodRewardMoms (F, C, D, numOfMoms)   
        % block partitioning of the input
        NF = size(F,1);
        ix = (1:NF);
        ixz = ix(abs(diag(C))<=precision);
        ixp = ix(diag(C)>precision);
        ixn = ix(diag(C)<-precision);
        Nz = length(ixz);
        Np = length(ixp);
        Nn = length(ixn);   
        % permutation matrix that converts between the original and the partitioned state ordering
        Per = zeros(NF);
        for i=1:Nz
            Per(i,ixz(i))=1;
        end
        for i=1:Np
            Per(Nz+i,ixp(i))=1;
        end
        for i=1:Nn
            Per(Nz+Np+i,ixn(i))=1;
        end
        iPer = inv(Per);
        Fc = mat2cell(Per*F*iPer, [Nz, Np, Nn], [Nz, Np, Nn]);
        [Fzz,Fpz,Fmz,Fzp,Fpp,Fmp,Fzm,Fpm,Fmm] = Fc{:};
        Cm = C(ixn,ixn);
        Cp = C(ixp,ixp);
        Dm = D(ixn,ixn);
        Dp = D(ixp,ixp);
        Dz = D(ixz,ixz);
        
        % detivatives of F(v)
        Fppd = cell(1,numOfMoms+1);
        Fpmd = cell(1,numOfMoms+1);
        Fmpd = cell(1,numOfMoms+1);
        Fmmd = cell(1,numOfMoms+1);
        Fppd{1} = inv(Cp)*(Fpp+Fpz*inv(-Fzz)*Fzp);
        Fpmd{1} = inv(Cp)*(Fpm+Fpz*inv(-Fzz)*Fzm);
        Fmpd{1} = inv(-Cm)*(Fmp+Fmz*inv(-Fzz)*Fzp);
        Fmmd{1} = inv(-Cm)*(Fmm+Fmz*inv(-Fzz)*Fzm);
        for i=1:numOfMoms
            dr = DReward(Fzz, Dz, i);
            Fppd{i+1} = inv(Cp) * Fpz * dr * Fzp;
            Fpmd{i+1} = inv(Cp) * Fpz * dr * Fzm;
            Fmpd{i+1} = inv(-Cm) * Fmz * dr * Fzp;
            Fmmd{i+1} = inv(-Cm) * Fmz * dr * Fzm;
            if i==1
                Fppd{i+1} = Fppd{i+1} -  inv(Cp)*Dp;
                Fmmd{i+1} = Fmmd{i+1} -  inv(-Cm)*Dm;
            end
        end       
        Psi = FluidFundamentalMatrices(Fppd{1}, Fpmd{1}, Fmpd{1}, Fmmd{1}, 'P', precision);
        BPM = cell(1,numOfMoms+1);       
        BPM{1} = Psi;
        for i=1:numOfMoms
            X = -Psi*Fmpd{i+1}*Psi + Fpmd{i+1};
            for m=0:i-1
                X = X + nchoosek(i,m) * ((Fppd{i-m+1} + Psi*Fmpd{i-m+1})*BPM{m+1} + BPM{m+1}*(Fmmd{i-m+1}+Fmpd{i-m+1}*Psi));
            end
            for l=1:i-1
                for m=1:i-l
                    X = X + nchoosek(i,l)*nchoosek(i-l,m)*BPM{l+1}*Fmpd{i-l-m+1}*BPM{m+1};
                end
            end
            BPM{i+1} = lyap(Fppd{1}+Psi*Fmpd{1}, Fmmd{1}+Fmpd{1}*Psi, X);
        end
        % re-order states
        for i=1:length(BPM)
            BPM{i} = iPer*[zeros(Nz,NF);zeros(Np,Nz+Np),BPM{i};zeros(Nn,NF)]*Per;
        end
    end

    function [pr, Pn] = BusyPeriodRewardDistr (F, C, D, t)
        % block partitioning of the input
        NF = size(F,1);
        ix = (1:NF);
        ixz = ix(abs(diag(C))<=precision);
        ixp = ix(diag(C)>precision);
        ixn = ix(diag(C)<-precision);
        Nz = length(ixz);
        Np = length(ixp);
        Nn = length(ixn);   
        % permutation matrix that converts between the original and the partitioned state ordering
        Per = zeros(NF);
        for i=1:Nz
            Per(i,ixz(i))=1;
        end
        for i=1:Np
            Per(Nz+i,ixp(i))=1;
        end
        for i=1:Nn
            Per(Nz+Np+i,ixn(i))=1;
        end
        iPer = inv(Per);
        Fc = mat2cell(Per*F*iPer, [Nz, Np, Nn], [Nz, Np, Nn]);
        [Fzz,Fpz,Fmz,Fzp,Fpp,Fmp,Fzm,Fpm,Fmm] = Fc{:};
        Cm = C(ixn,ixn);
        Cp = C(ixp,ixp);
        Dm = D(ixn,ixn);
        Dp = D(ixp,ixp);
        Dz = D(ixz,ixz);        
        % start erlangization
            L = erlMaxOrder;
            nu = L/t;
            Z = inv(nu*Dz-Fzz);
            Psie = FluidFundamentalMatrices (inv(Cp)*(Fpp-nu*Dp+Fpz*Z*Fzp), inv(Cp)*(Fpm+Fpz*Z*Fzm), inv(-Cm)*(Fmp+Fmz*Z*Fzp), inv(-Cm)*(Fmm-nu*Dm+Fmz*Z*Fzm), 'P', precision);
            Pn = {Psie};
            pr = Psie;
            AM = inv(Cp)*(Fpp-nu*Dp+Fpz*Z*Fzp) + Psie*inv(-Cm)*(Fmp+Fmz*Z*Fzp);
            BM = inv(-Cm)*(Fmm-nu*Dm+Fmz*Z*Fzm) + inv(-Cm)*(Fmp+Fmz*Z*Fzp)*Psie;
            for n=1:L-1
                CM = inv(Cp)*nu*Dp*Pn{n} + Pn{n}*inv(-Cm)*nu*Dm;
                for i=1:n-1
                    CM = CM + Pn{i+1}*inv(-Cm)*Fmp*Pn{n-i+1};
                end
                CM = CM + inv(Cp)*Fpz*Z*(nu*Dz*Z)^n*Fzm - Psie*inv(-Cm)*Fmz*Z*(nu*Dz*Z)^n*Fzp*Psie;
                if ~isempty(ixz)
                    for i=0:n-1
                        CM = CM + Pn{i+1}*inv(-Cm)*Fmz*Z*(nu*Dz*Z)^(n-i)*(Fzm+Fzp*Psie);
                        CM = CM + (inv(Cp)*Fpz+Psie*inv(-Cm)*Fmz)*Z*(nu*Dz*Z)^(n-i)*Fzp*Pn{i+1};
                    end
                    for i=1:n-1
                        for j=1:n-i
                            CM = CM + Pn{i+1}*inv(-Cm)*Fmz*Z*(nu*Dz*Z)^(n-i-j)*Fzp*Pn{j+1};
                        end
                    end
                end
                PM = lyap(AM, BM, CM);
                Pn{n+1} = PM;               
                % accumulation
                pr = pr + PM;
            end
            % re-order states
            pr = iPer*[zeros(Nz,NF);zeros(Np,Nz+Np),pr;zeros(Nn,NF)]*Per;
            for i=1:length(Pn)
                Pn{i} = iPer*[zeros(Nz,NF);zeros(Np,Nz+Np),Pn{i};zeros(Nn,NF)]*Per;
            end
        % end of erlangization
    end        

    % some preparation
    pi = CTMCSolve(Q);
    lambda = pi*R';
    N = size(Q,1);
    
    % step 2. calculate performance measures
    % ======================================
    Ret = {};
    for k=classes
        % step 1. solution of the workload process for fluid types having
        % the same or higher priority
        % ============================================================
        [mass0, ini, Km, clo] = GeneralFluidSolve (Q, diag(sum(R(k:end,:),1))/d-eye(N), [], precision);
        KN = size(Km,1);
        clok = clo*diag(R(k,:)) / lambda(k);
        
        % similarity transformation
        Delta = diag(inv(-Km)*sum(clok,2));
        K0 = inv(Delta)*Km*Delta;
        K1 = inv(Delta)*clok;
        kappa = ini*Delta;

        Km = K0;
        clok = K1;
        ini = kappa;

        if k<K          
            % step 4.3. calculate the performance measures
            % ==========================================   
            argIx = 1;
            while argIx<=length(varargin)
                if any(ismember(eaten, argIx))
                    argIx = argIx + 1;
                    continue;
                elseif strcmp(varargin{argIx},'stMoms') 
                    % MOMENTS OF THE SOJOURN TIME
                    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    numOfSTMoms = varargin{argIx+1};
                    F = [Km, clok; zeros(N,KN), Q];
                    C = [eye(KN), zeros(KN,N); zeros(N,KN), diag(sum(R(k+1:end,:),1))/d-eye(N)];
                    D = [zeros(KN,KN+N); zeros(N,KN), eye(N)];
                    Tmp = BusyPeriodRewardMoms (F, C, D, numOfSTMoms);
                    inis = [ini, zeros(1,N)];
                    stMoms = zeros(1,numOfSTMoms);
                    for i=1:length(Tmp)-1
                        stMoms(i) = (-1)^i*sum(inis * Tmp{i+1});
                    end
                    Ret{end+1} = stMoms;
                    argIx = argIx + 1;
                elseif strcmp(varargin{argIx},'flMoms') 
                    % MOMENTS OF THE FLUID LEVEL 
                    % ~~~~~~~~~~~~~~~~~~~~~~~~~~                    
                    % first compute moments right after fluid drop
                    % departures
                    numOfFLMoms = varargin{argIx+1};
                    F = [Km, clok; zeros(N,KN), Q];
                    C = [eye(KN), zeros(KN,N); zeros(N,KN), diag(sum(R(k+1:end,:),1))/d-eye(N)];
                    D = [zeros(KN,KN+N); zeros(N,KN), diag(R(k,:))];
                    FLDPn = BusyPeriodRewardMoms (F, C, D, numOfFLMoms);
                    for i=1:length(FLDPn)
                        X = FLDPn{i};
                        FLDPn{i} = X(:,KN+1:end);
                    end
                    inis = [ini, zeros(1,N)];
                    fldMoms = zeros(1,numOfFLMoms);                    
                    for i=1:length(FLDPn)
                        FLDPn{i} = inis*FLDPn{i};
                        if i==1
                            FLDPn{i} = FLDPn{i} + mass0*diag(R(k,:))/lambda(k);
                        end
                        fldMoms(i) = (-1)^(i-1)*sum(FLDPn{i});
                    end
                    % calculate moments in random point of time
                    FLPn = {pi};
                    flMoms = zeros(1,numOfFLMoms);
                    iTerm = inv(ones(N,1)*pi - Q);
                    for n=1:numOfFLMoms
                        sumP = sum(FLDPn{n+1}) + n*(-FLDPn{n} + FLPn{n}*diag(R(k,:))/lambda(k))*iTerm*R(k,:)';
                        P = sumP*pi + n*(-FLPn{n}*diag(R(k,:)) + FLDPn{n}*lambda(k))*iTerm;
                        FLPn{n+1} = P;
                        flMoms(n) = (-1)^n*sum(P);
                    end
                    Ret{end+1} = flMoms;
                    argIx = argIx + 1;                    
                elseif strcmp(varargin{argIx},'stDistr') 
                    % DISTRIBUTION OF THE SOJOURN TIME
                    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    stCdfPoints = varargin{argIx+1};
                    F = [Km, clok; zeros(N,KN), Q];
                    C = [eye(KN), zeros(KN,N); zeros(N,KN), diag(sum(R(k+1:end,:),1))/d-eye(N)];
                    D = [zeros(KN,KN+N); zeros(N,KN), eye(N)];
                    inis = [ini, zeros(1,N)];
                    res = zeros(1,length(stCdfPoints));
                    for x=1:length(stCdfPoints)
                        Tmp = BusyPeriodRewardDistr (F, C, D, stCdfPoints(x));
                        res(x) = sum(mass0*diag(R(k,:))/lambda(k)) + sum(inis*Tmp);
                    end
                    Ret{end+1} = res;
                    argIx = argIx + 1;
                elseif strcmp(varargin{argIx},'flDistr') 
                    % DISTRIBUTION OF THE FLUID LEVEL
                    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    flCdfPoints = varargin{argIx+1};
                    F = [Km, clok; zeros(N,KN), Q];
                    C = [eye(KN), zeros(KN,N); zeros(N,KN), diag(sum(R(k+1:end,:),1))/d-eye(N)];
                    D = [zeros(KN,KN+N); zeros(N,KN), diag(R(k,:))];
                    inis = [ini, zeros(1,N)];
                    res = zeros(1,length(flCdfPoints));
                    % resDep = zeros(1,length(flCdfPoints));
                    for x=1:length(flCdfPoints)
                        nu = erlMaxOrder/flCdfPoints(x);
                        [Tmp, Psix] = BusyPeriodRewardDistr (F, C, D, flCdfPoints(x));
                        Psiy = lambda(k)*nu*(mass0*diag(R(k,:))/lambda(k)+inis*Psix{1}(:,KN+1:end)) * inv(nu*diag(R(k,:))-Q);
                        for i=2:length(Psix)
                            Psiy = nu*(lambda(k)*inis*Psix{i}(:,KN+1:end) + Psiy*diag(R(k,:)))*inv(nu*diag(R(k,:))-Q);
                        end
                        res(x) = sum(Psiy);
                        % resDep(x) = sum(mass0*diag(R(k,:))/lambda(k)) + sum(inis*Tmp);
                    end                  
                    Ret{end+1} = res;
                    argIx = argIx + 1;
                else
                    error (['FluidPrioQueue: Unknown parameter ' varargin{argIx}])
                end
                argIx = argIx + 1;
            end
        elseif k==K
            % step 3. calculate the performance measures
            % ==========================================   
            argIx = 1;
            while argIx<=length(varargin)
                if any(ismember(eaten, argIx))
                    argIx = argIx + 1;
                    continue;
                elseif strcmp(varargin{argIx},'stMoms') 
                    % MOMENTS OF THE SOJOURN TIME
                    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    numOfSTMoms = varargin{argIx+1};
                    stMoms = zeros(1,numOfSTMoms);
                    for i=1:numOfSTMoms
                        stMoms(i) = factorial(i) * sum(ini*inv(-Km)^(i+1)*clok);
                    end
                    Ret{end+1} = stMoms;
                    argIx = argIx + 1;
                elseif strcmp(varargin{argIx},'stDistr') 
                    % DISTRIBUTION OF THE SOJOURN TIME
                    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    stCdfPoints = varargin{argIx+1};
                    res = zeros(1,length(stCdfPoints));
                    for x=1:length(stCdfPoints)
                        res(x) = sum(mass0*diag(R(k,:))/lambda(k)) + sum(ini*inv(-Km)*(eye(size(Km))-expm(Km*stCdfPoints(x)))*clok);
                    end
                    Ret{end+1} = res;
                    argIx = argIx + 1;
                    
                elseif strcmp(varargin{argIx},'flMoms') 
                    % MOMENTS OF THE FLUID LEVEL 
                    % ~~~~~~~~~~~~~~~~~~~~~~~~~~                    
                    % first compute moments right after fluid drop
                    % departures
                    numOfFLMoms = varargin{argIx+1};
                    Ret{end+1} = FluFluQueue(Q,diag(R(k,:)),0,d,false,'flMoms', numOfFLMoms);
                    argIx = argIx + 1;
                elseif strcmp(varargin{argIx},'flDistr') 
                    % DISTRIBUTION OF THE FLUID LEVEL
                    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    % at departures:
                    flCdfPoints = varargin{argIx+1};
                    Ret{end+1} = FluFluQueue(Q,diag(R(k,:)),0,d,false,'flDistr', flCdfPoints);
                    argIx = argIx + 1;
                else
                    error (['FluidPrioQueue: Unknown parameter ' varargin{argIx}])
                end
                argIx = argIx + 1;
            end
        end
    end   
    varargout = Ret;
end