%{
%{
 % @file pfqn_clust.m
 % @brief de Souza e Silva-Lavenberg-Muntz Clustering Approximation (CA).
%}
%}

function [XN,QN,UN,RN,it] = pfqn_clust(L,N,Z,subnets,localclasses,inner,tol,maxiter)
%{
%{
 % @brief de Souza e Silva-Lavenberg-Muntz Clustering Approximation (CA).
 %
 % E. de Souza e Silva, S. S. Lavenberg, R. R. Muntz, "A clustering
 % approximation technique for queueing network models with a large number of
 % chains", IEEE Trans. Computers C-35(5), 1986. The network is covered by
 % subnetworks whose union is the whole network but which need not be
 % disjoint. Every class visiting a subnetwork S is either LOCAL to S, and is
 % then solved inside it, or FOREIGN, and is then seen only through the
 % utilization it leaves behind. Each subnetwork is solved by an ordinary
 % approximate MVA algorithm with two replacements: the complement of S is
 % collapsed into a per-class delay P_c (eq. 2.41) and the foreign classes
 % into a per-centre utilization U_k (eq. 2.42),
 %
 %   X_c(N) = N_c / (sum_{k in S} R_ck(N) + Z_c + P_c),                (2.43)
 %   Q_k(N) = [sum_{c in LC(S)} R_ck(N) X_c(N) + U_k] / (1 - U_k).     (2.44)
 %
 % Choosing the PE algorithm for every subnetwork reproduces global PE
 % exactly, so the useful setting is Linearizer inside, PE outside: the cost
 % then sits between pfqn_bs and pfqn_linearizer, which is the point of the
 % method. The answer depends on the decomposition, which is an input.
 %
 % When no decomposition is supplied the criterion of the paper is applied
 % automatically: the cheap PAMB estimate (pfqn_pam) of the centre
 % utilizations is taken, every class is attached to the centre where it
 % loads the most, classes sharing that centre form one cluster, and the
 % subnetwork of a cluster is the set of centres its classes visit.
 %
 % The name avoids pfqn_ca, which is the exact convolution algorithm.
 %
 % @fn pfqn_clust(L, N, Z, subnets, localclasses, inner, tol, maxiter)
 % @param L Service demand matrix (stations x classes).
 % @param N Population vector.
 % @param Z Think time vector.
 % @param subnets Cell array of station index vectors; their union must cover
 %        1..M. Empty for the automatic decomposition above.
 % @param localclasses Cell array of class index vectors, one per subnetwork,
 %        listing the local classes of that subnetwork. Empty for automatic.
 % @param inner 'lin' (default) or 'bs', the algorithm run inside a subnetwork.
 % @param tol Convergence tolerance (default: 1e-6).
 % @param maxiter Maximum number of outer iterations (default: 1000).
 % @return XN System throughput.
 % @return QN Mean queue lengths.
 % @return UN Utilization.
 % @return RN Residence times.
 % @return it Number of outer iterations performed.
%}
%}

[M,R]=size(L);
if nargin<3 || isempty(Z)
    Z = zeros(1,R);
end
Z = sum(Z,1);
if nargin<4
    subnets = [];
end
if nargin<5
    localclasses = [];
end
if nargin<6 || isempty(inner)
    inner = 'lin';
end
if nargin<7 || isempty(tol)
    tol = 1e-6;
end
if nargin<8 || isempty(maxiter)
    maxiter = 1000;
end

%% seed and, if needed, the automatic decomposition
[X0,Q0] = pfqn_pam(L,N,Z,'pamb');
if isempty(subnets) || isempty(localclasses)
    U0 = L .* repmat(X0,M,1);
    bottleneck = ones(1,R);
    for r=1:R
        if any(L(:,r)>0)
            [~,bottleneck(r)] = max(U0(:,r));
        end
    end
    centres = unique(bottleneck);
    subnets = cell(1,numel(centres));
    localclasses = cell(1,numel(centres));
    for g=1:numel(centres)
        cls = find(bottleneck==centres(g));
        localclasses{g} = cls;
        subnets{g} = find(any(L(:,cls)>0,2))';
        if isempty(subnets{g})
            subnets{g} = centres(g);
        end
    end
    % the union of the subnetworks must be the whole network
    missing = setdiff(1:M, unique([subnets{:}]));
    if ~isempty(missing)
        subnets{end+1} = missing;
        localclasses{end+1} = [];
    end
end
G = numel(subnets);
owner = zeros(1,R);
for g=1:G
    owner(localclasses{g}) = g;
end

XN = X0;
QN = Q0;
for it=1:maxiter
    QN_1 = QN;
    Qk = sum(QN,2);
    for g=1:G
        S = subnets{g};
        LC = localclasses{g};
        if isempty(LC) || isempty(S)
            continue
        end
        outside = setdiff(1:M, S);
        FC = setdiff(find(any(L(S,:)>0,1)), LC);

        % (2.41) delay of a local class in the complement of S
        P = zeros(1,numel(LC));
        for a=1:numel(LC)
            c = LC(a);
            if N(c) <= 0
                continue
            end
            for ist=outside
                P(a) = P(a) + L(ist,c)*(1+Qk(ist))/(1 + L(ist,c)*XN(c)/N(c));
            end
        end
        % (2.42) utilization left in S by the foreign classes
        Uk = zeros(numel(S),1);
        for b=1:numel(S)
            ist = S(b);
            for i=FC
                if N(i) > 0
                    Uk(b) = Uk(b) + L(ist,i)*XN(i)/(1 + L(ist,i)*XN(i)/N(i));
                end
            end
        end
        Uk = min(Uk, 1-1e-8);

        [Xs,Qs] = SubnetSolve(L(S,LC),N(LC),Z(LC)+P,Uk,inner,tol,maxiter);
        XN(LC) = Xs;
        QN(S,LC) = Qs;
        % the local classes still hold jobs outside S; charge them the
        % per-station term of (2.41)
        for a=1:numel(LC)
            c = LC(a);
            for ist=outside
                QN(ist,c) = XN(c)*L(ist,c)*(1+Qk(ist))/(1 + L(ist,c)*XN(c)/max(N(c),eps));
            end
        end
    end
    % a class owned by no subnetwork keeps the seed throughput
    for r=find(owner==0)
        QN(:,r) = XN(r)*L(:,r);
    end
    if max(max(abs(QN-QN_1))) < tol
        break
    end
end
UN = repmat(XN,M,1) .* L;
RN = QN ./ repmat(XN,M,1);
RN(:,N==0) = 0;
RN(~isfinite(RN)) = 0;
end

function [X,Q] = SubnetSolve(L,N,Z,Uk,inner,tol,maxiter)
% Approximate MVA restricted to the local classes of one subnetwork, with
% (2.43) for the throughput and (2.44) for the total queue length.
[M,R] = size(L);
X = zeros(1,R);
Q = repmat(N,M,1)/max(M,1);
if M==0 || R==0
    return
end
switch inner
    case 'bs'
        for it=1:maxiter %#ok<NASGU>
            Qold = Q;
            for r=1:R
                if N(r) <= 0
                    continue
                end
                % PE arrival-instant local queue, inflated by (2.44)
                A = (sum(Q,2) - Q(:,r)/N(r) + Uk)./(1-Uk);
                W = L(:,r) .* (1 + A);
                X(r) = N(r)/(Z(r)+sum(W));
                Q(:,r) = X(r)*W;
            end
            if max(max(abs(Q-Qold))) < tol
                break
            end
        end
    otherwise % 'lin'
        Qs = zeros(M,R,1+R);
        for s=0:R
            Qs(:,:,1+s) = Q;
        end
        Delta = zeros(M,R,R);
        for pass=1:3 %#ok<NASGU>
            for s=0:R
                Ns = oner(N,s);
                Qs(:,:,1+s) = SubnetCore(L,Ns,Z,Uk,Qs(:,:,1+s),Delta,tol,maxiter);
            end
            for r=1:R
                for s=1:R
                    Ns = oner(N,s);
                    if N(r) > 0 && Ns(r) > 0
                        Delta(:,r,s) = Qs(:,r,1+s)/Ns(r) - Qs(:,r,1+0)/N(r);
                    elseif N(r) > 0
                        Delta(:,r,s) = -Qs(:,r,1+0)/N(r);
                    else
                        Delta(:,r,s) = 0;
                    end
                end
            end
        end
        Qs(:,:,1+0) = SubnetCore(L,N,Z,Uk,Qs(:,:,1+0),Delta,tol,maxiter);
        [Q,X] = SubnetForward(L,N,Z,Uk,Qs(:,:,1+0),Delta);
end
end

function Q = SubnetCore(L,N,Z,Uk,Q,Delta,tol,maxiter)
for it=1:maxiter %#ok<NASGU>
    Qold = Q;
    Q = SubnetForward(L,N,Z,Uk,Q,Delta);
    if max(max(abs(Q-Qold))) < tol
        break
    end
end
end

function [Qout,X] = SubnetForward(L,N,Z,Uk,Q,Delta)
[M,R] = size(L);
Qout = zeros(M,R);
X = zeros(1,R);
for r=1:R
    if N(r) <= 0
        continue
    end
    % Linearizer estimate of the local queue at N - 1_r, inflated by (2.44)
    Nr = oner(N,r);
    Qm = zeros(M,1);
    for s=1:R
        if N(s) > 0 && Nr(s) > 0
            Qm = Qm + Nr(s)*(Q(:,s)/N(s) + Delta(:,s,r));
        end
    end
    A = (max(Qm,0) + Uk)./(1-Uk);
    W = L(:,r) .* (1 + A);
    X(r) = N(r)/(Z(r)+sum(W));
    Qout(:,r) = X(r)*W;
end
end
