% [masses, iniF, KF, cloF, iniB, KB, cloB] = SecondOrderLevelDependentFluidSolve (Q, R, S, T, boundaryL, boundaryU, Qt, prec)
%
% * Q: cell of generators in the different layers/regimes
% * R: cell of diagonal matrices of fluid rates in the different layers/regimes
% * S: cell of diagonal matrices of fluid variances in the different layers/regimes
% * T: vector of thresholds
% * boundaryL/U: vector defining the lower and the upper boundary behavior
%       one item for each state of the background process
%           =0: reflective boundary (default)
%           =1: absorbing boundary
%
% The performance measures (pdf, cdf, etc.) can be computed by the 
% LevelDependentFluidStationaryDistr function.
% 
function [masses, iniF, KF, cloF, iniB, KB, cloB] = SecondOrderLevelDependentFluidSolve (Q, R, S, T, boundaryL, boundaryU, Qt, prec)

    K = length(T);
    N = size(Q{1},1);

    if ~exist('prec','var')
        prec = 1e-14;
    end
    if ~exist('boundaryL','var')
        % zeros(1,N), not zeros(K,N): boundaryL/U is ONE flag per background
        % state (see the header above, "vector defining the lower and the upper
        % boundary behavior"), and there is a single lower and a single upper
        % boundary rather than one per regime. All four use sites below index
        % ix = 1:N with it, so a K x N default makes the logical mask larger
        % than the array it indexes and the function crashed on its own defaults
        % for every K >= 2 -- that is, for every genuinely multi-regime model.
        % K == 1 survived by accident, zeros(1,N) being the correct shape there.
        boundaryL = zeros(1,N);
    end
    if ~exist('boundaryU','var')
        boundaryU = boundaryL;
    end    
    if ~exist('Qt','var') || isempty(Qt)
        for k=1:K
            Qt{k} = Q{k};
        end
        Qt{K+1} = Q{K};
    elseif length(Qt)==1
        for k=2:K+1
            Qt{k} = Qt{1};
        end
    end
    
    T = [0,T];

    % preparation
    KF = cell(1,K);
    KB = cell(1,K);
    cloF = cell(1,K);
    cloB = cell(1,K);
    Np = zeros(1,K);
    Nn = zeros(1,K);
    Ns = zeros(1,K);
    NbF = zeros(1,K);
    NbB = zeros(1,K);
    vixp = cell(1,K);
    vixn = cell(1,K);
    vix0 = cell(1,K);
    vixs = cell(1,K);
    vixsn = cell(1,K);
    vixsp = cell(1,K);
    vixs0 = cell(1,K);
    for k=1:K
        
        S{k} = S{k}/2;
        % collect zero states
        ix = (1:N);
        ix0 = ix(abs(diag(R{k}))<=prec & diag(S{k})<=prec);
        ixn0 = setdiff(ix,ix0);

        Q00 = Q{k}(ix0,ix0);
        Q0n = Q{k}(ix0,ixn0);
        Qn0 = Q{k}(ixn0,ix0);
        Qnn = Q{k}(ixn0,ixn0);

        Qv = Qnn + Qn0*inv(-Q00)*Q0n;
        Rv = R{k}(ixn0,ixn0);
        Sv = S{k}(ixn0,ixn0);

        Nnz = size(Qv,1);
        
        % partition the state space according to zero, positive and negative fluid rates
        ix = (1:Nnz);
        ixp = ix(diag(Rv)>prec & diag(Sv)<=prec);
        ixn = ix(diag(Rv)<-prec & diag(Sv)<=prec);
        ixs = ix(diag(Sv)>prec);
        Np(k) = length(ixp);
        Nn(k) = length(ixn);
        Ns(k) = length(ixs);

        % FORWARD parameters

        % obtain "c" constant to transform matrix equations to QBD like matrix
        % quadratic equations
        c1 = max(-diag(Qv(ixp,ixp))./diag(Rv(ixp,ixp)));
        discr = (diag(Rv(ixs,ixs)).^2  - 2*diag(2*Sv(ixs,ixs)).*diag(Qv(ixs,ixs)));
        c2 = max((discr>0).*((-diag(Rv(ixs,ixs))+sqrt(discr))./diag(2*Sv(ixs,ixs))));
        c = max([c1,c2,1]);

        ixbF = [ixs,ixp];
        NbF(k) = length(ixbF);

        Bm = blkdiag(c*Sv(ixbF,ixbF),-Rv(ixn,ixn));
        Lm = [-Rv(ixbF,ixbF)-2*c*Sv(ixbF,ixbF),zeros(NbF(k),Nn(k));Qv(ixn,ixbF)/c,Qv(ixn,ixn)/c+Rv(ixn,ixn)];
        Fm = [Qv(ixbF,ixbF)/c+c*Sv(ixbF,ixbF)+Rv(ixbF,ixbF),Qv(ixbF,ixn)/c;zeros(Nn(k),NbF(k)+Nn(k))];

        % Solve QBD for matrix R
        [~, QBDR] = QBD_CR(Bm, Lm, Fm);

        % extract K and Psi from the solution
        KF{k} = (QBDR(1:NbF(k),1:NbF(k)) - eye(NbF(k))) * c;
        PsiF = QBDR(1:NbF(k),NbF(k)+1:end);

        % closing matrix of the stationary density of the fluid level

        clovF = zeros(NbF(k),NbF(k)+Nn(k));
        clovF(:,ixbF) = eye(NbF(k));
        clovF(:,ixn) = PsiF;

        cloF{k} = zeros(NbF(k),N);
        cloF{k}(:,ixn0) = clovF;
        cloF{k}(:,ix0) = clovF*Qn0*inv(-Q00);        
        
        % BACKWARD parameters

        % obtain "c" constant to transform matrix equations to QBD like matrix
        % quadratic equations
        c1 = max(-diag(Qv(ixn,ixn))./diag(-Rv(ixn,ixn)));
        discr = (diag(Rv(ixs,ixs)).^2  - 2*diag(2*Sv(ixs,ixs)).*diag(Qv(ixs,ixs)));
        c2 = max((discr>0).*((diag(Rv(ixs,ixs))+sqrt(discr))./diag(2*Sv(ixs,ixs))));
        c = max([c1,c2,1]);

        ixbB = [ixs,ixn];
        NbB(k) = length(ixbB);

        Bm = blkdiag(c*Sv(ixbB,ixbB),Rv(ixp,ixp));
        Lm = [Rv(ixbB,ixbB)-2*c*Sv(ixbB,ixbB),zeros(NbB(k),Np(k));Qv(ixp,ixbB)/c,Qv(ixp,ixp)/c-Rv(ixp,ixp)];
        Fm = [Qv(ixbB,ixbB)/c+c*Sv(ixbB,ixbB)-Rv(ixbB,ixbB),Qv(ixbB,ixp)/c;zeros(Np(k),NbB(k)+Np(k))];

        % Solve QBD for matrix R
        [~, QBDR] = QBD_CR(Bm, Lm, Fm);

        % extract K and Psi from the solution
        KB{k} = (QBDR(1:NbB(k),1:NbB(k)) - eye(NbB(k))) * c;
        PsiB = QBDR(1:NbB(k),NbB(k)+1:end);

        % closing matrix of the stationary density of the fluid level
        clovB = zeros(NbB(k),NbB(k)+Np(k));
        clovB(:,ixbB) = eye(NbB(k));
        clovB(:,ixp) = PsiB;
        
        cloB{k} = zeros(NbB(k),N);
        cloB{k}(:,ixn0) = clovB;
        cloB{k}(:,ix0) = clovB*Qn0*inv(-Q00);           
        
        vixp{k} = ixn0(ixp);
        vixn{k} = ixn0(ixn);
        vix0{k} = ix0;
        vixs{k} = ixn0(ixs);        
        vixsn{k} = ixn0(diag(Sv)>prec & diag(Rv)<0);
        vixsp{k} = ixn0(diag(Sv)>prec & diag(Rv)>=0);
        vixs0{k} = ixn0(diag(Sv)>prec & diag(Rv)==0);
    end
    
    % Boundary vectors
    % ================
    
    % Construct and solve linear set of equations for the unknows (probability 
    % masses and density parameters iniF and iniB for all regimes)
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Neqns = (K+1)*N + sum(Np) + sum(Nn) + 2*sum(Ns);

    if Neqns>300
        M = sparse(Neqns,Neqns);
    else
        M = zeros(Neqns);
    end
    pos = 1;
    for k=1:K
        pos = [pos, N, NbF(k), NbB(k)];
    end
    pp = cumsum(pos);
    p = 1;
    i = 1;
    % equalities for flux conservation
    for k=0:K
        M(pp(i):pp(i)+N-1,k*N+1:(k+1)*N) = -Qt{k+1};
        if k>0
            M(pp(i-2):pp(i-2)+NbF(k)-1,k*N+1:(k+1)*N) = expm(KF{k}*(T(k+1)-T(k)))*(-cloF{k}*R{k} + KF{k}*cloF{k}*S{k});
            M(pp(i-1):pp(i-1)+NbB(k)-1,k*N+1:(k+1)*N) = -cloB{k}*R{k} - KB{k}*cloB{k}*S{k};
        end
        if k<K
            M(pp(i+1):pp(i+1)+NbF(k+1)-1,k*N+1:(k+1)*N) = cloF{k+1}*R{k+1}-KF{k+1}*cloF{k+1}*S{k+1};
            M(pp(i+2):pp(i+2)+NbB(k+1)-1,k*N+1:(k+1)*N) = expm(KB{k+1}*(T(k+2)-T(k+1)))*(cloB{k+1}*R{k+1} + KB{k+1}*cloB{k+1}*S{k+1});
        end
        i = i + 3;
    end
    
    % further equations
    col = (K+1)*N+1;
    ix = (1:N);
    i = 1;
    for k=0:K
        if k==0
            % mass is zero in positive and in reflective states
            ixr0 = intersect(ix(boundaryL==0), vixs{1});
            Nr0 = length(ixr0);
            ms0 = zeros(N,Np(1)+Nr0);
            ms0([vixp{1},ixr0],:) = eye(Np(1)+Nr0);
            M(pp(i):pp(i)+N-1,col:col+Np(1)+Nr0-1) = ms0;
            col = col + Np(1)+Nr0;
            % density is zero in absorbing states
            ixa0 = intersect(ix(boundaryL==1), vixs{1});
            Na0 = length(ixa0);   
            pdfF = cloF{1};
            pdfB = expm(KB{1}*T(2))*cloB{1};
            M(pp(i+1):pp(i+1)+NbF(k+1)-1,col:col+Na0-1) = pdfF(:,ixa0);
            M(pp(i+2):pp(i+2)+NbB(k+1)-1,col:col+Na0-1) = pdfB(:,ixa0);
            col = col + Na0;
        elseif k==K
            % mass is zero in negative and in reflective states
            ixrB = intersect(ix(boundaryU==0), vixs{K});
            NrB = length(ixrB);
            ms0 = zeros(N,Nn(K)+NrB);
            ms0([vixn{K},ixrB],:) = eye(Nn(K)+NrB);
            M(pp(i):pp(i)+N-1,col:col+Nn(K)+NrB-1) = ms0;
            col = col + Nn(K)+NrB;
            % density is zero in absorbing states
            ixaB = intersect(ix(boundaryU==1), vixs{K});   
            NaB = length(ixaB);
            pdfF = expm(KF{K}*(T(K+1)-T(K)))*cloF{K};
            pdfB = cloB{K};
            M(pp(i-2):pp(i-2)+NbF(k)-1,col:col+NaB-1) = pdfF(:,ixaB);
            M(pp(i-1):pp(i-1)+NbB(k)-1,col:col+NaB-1) = pdfB(:,ixaB);
            col = col + NaB;
        else
            % there is no mass except S+(k) U S-(k+1)
            st0 = setdiff(ix, union(intersect(vixp{k},vixn{k+1}),union(vix0{k},vix0{k+1})));
            N0 = length(st0);
            ms0 = zeros(N,N0);
            ms0(st0,:) = eye(N0);
            M(pp(i):pp(i)+N-1,col:col+N0-1) = ms0;
            col = col + N0;            
            % equations for the second-order states
            sts = setdiff(union(vixs{k},vixs{k+1}), union(vixn{k+1},vixp{k}));
            Nss = length(sts);
            BelowF = expm(KF{k}*(T(k+1)-T(k)))*-cloF{k}*sqrt(S{k});
            BelowB = -cloB{k}*sqrt(S{k});
            AboveF = cloF{k+1}*sqrt(S{k+1});
            AboveB = expm(KB{k+1}*(T(k+2)-T(k+1)))*cloB{k+1}*sqrt(S{k+1});
            M(pp(i-2):pp(i-2)+NbF(k)-1,col:col+Nss-1) = BelowF(:,sts);
            M(pp(i-1):pp(i-1)+NbB(k)-1,col:col+Nss-1) = BelowB(:,sts);
            M(pp(i+1):pp(i+1)+NbF(k+1)-1,col:col+Nss-1) = AboveF(:,sts);
            M(pp(i+2):pp(i+2)+NbB(k+1)-1,col:col+Nss-1) = AboveB(:,sts);
            col = col + Nss;
        end
        i = i + 3;
    end
    
    % this function computes the integral of a matrix exponential from 0 to
    % L if it has a 0 eigenvalue hence can not be inverted
    function KAi = integExp(KA,L)
        l=CRPSolve(KA);
        r=CRPSolve(KA')';  
        l=l/(l*r);
        KAi = inv(-(KA-r*l)) *(eye(size(KA,1))-expm((KA-r*l)*L)) + r*l*(L+exp(-L)-1);       
    end

    % this function receives two matrices and calculates the integrals of
    % the matrix exponentials if exactly one of them has a 0 eigenvalue
    function [KAi,KBi] = integExp2(KA,KB,L)
        if min(abs(eig(KA))) > min(abs(eig(KB)))
            KAi = inv(-KA)*(eye(size(KA,1))-expm(KA*L));
            KBi = integExp(KB, L);
        else
            KAi = integExp(KA, L);
            KBi = inv(-KB)*(eye(size(KB,1))-expm(KB*L));
        end
    end
    
    % normalizing condition
    h = ones(N,1);
    Mi = cell(1,K);
    for k=1:K
        [sumKF, sumKB] = integExp2(KF{k}, KB{k}, T(k+1)-T(k));
        h = [h; sum(sumKF*cloF{k},2); sum(sumKB*cloB{k},2); ones(N,1)];
        Mi{k} = [sum(sumKF*cloF{k},2); sum(sumKB*cloB{k},2)];
    end
    
    % solve linear system
    M(:,1) = h;
    b = [1,zeros(1,length(h)-1)]/M;
        
    % obtain solution
    masses = cell(1,K);
    iniF = cell(1,K);
    iniB = cell(1,K);
    masses{1} = b(1:N);
    i = 2;
    for k=1:K
        iniF{k} = b(pp(i):pp(i)+NbF(k)-1);
        iniB{k} = b(pp(i+1):pp(i+1)+NbB(k)-1);
        masses{k+1} = b(pp(i+2):pp(i+2)+N-1);
        i = i + 3;
    end
end

