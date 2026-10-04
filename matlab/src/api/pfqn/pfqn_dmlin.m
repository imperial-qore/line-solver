%{
%{
 % @file pfqn_dmlin.m
 % @brief de Souza e Silva-Muntz Improved Linearizer (IL).
%}
%}

function [Q,U,W,C,X,totiter] = pfqn_dmlin(L,N,Z,type,tol,maxiter,QN0,npasses)
%{
%{
 % @brief de Souza e Silva-Muntz Improved Linearizer (IL).
 %
 % E. de Souza e Silva, R. R. Muntz, "A note on the computational cost of the
 % Linearizer algorithm for queueing networks", IEEE Trans. Computers 39(6),
 % 1990. Linearizer evaluates the arrival-instant queue length as
 %
 %   A_k^(c)(n) = sum_i (n_i - delta_c^(i)) [Q_ik(n)/n_i + Delta^(i)_ck],
 %
 % re-summing the C Delta-terms at every Core iteration, at every one of the
 % C+1 populations: O(K C^3) per refresh pass. IL splits that sum into the
 % part that moves with the Core iterate and the part that does not,
 %
 %   A_k^(c)(n) = sum_i (n_i - delta_c^(i)) Q_ik(n)/n_i + xi_ck(n),
 %   xi_ck(N)       = sum_i (N_i - delta_c^(i)) Delta^(i)_ck,
 %   xi_ck(N - 1_j) = xi_ck(N) - Delta^(j)_ck,
 %
 % so the C K aggregates xi are computed ONCE per refresh pass and each Core
 % iteration then costs O(K C) instead of O(K C^2). Time drops to O(K C^2)
 % with the space unchanged at O(K C^2), and, because the split is an
 % identity and not an approximation, the fixed point is the one Linearizer
 % reaches: pfqn_dmlin and pfqn_linearizer agree to round-off. That is what
 % makes IL dominate AQL (pfqn_aql), which buys the same cost by aggregating
 % the queue lengths themselves and does change the answer.
 %
 % @fn pfqn_dmlin(L, N, Z, type, tol, maxiter, QN0, npasses)
 % @param L Service demand matrix (stations x classes).
 % @param N Population vector.
 % @param Z Think time vector.
 % @param type Scheduling strategy type per station (accepted, unused: the
 %        Linearizer family in LINE treats every station as single-server PS).
 % @param tol Convergence tolerance (default: 1e-8).
 % @param maxiter Maximum number of iterations (default: 1000).
 % @param QN0 (M x R) queue lengths that warm-start the Bard-Schweitzer
 %        initialization; empty for the default cold start.
 % @param npasses Number of xi refresh passes (default 3, the Chandy-Neuse rule).
 % @return Q Mean queue lengths.
 % @return U Utilization.
 % @return W Residence times.
 % @return C Cycle times.
 % @return X System throughput.
 % @return totiter Total iterations performed.
%}
%}

if nargin<8 || isempty(npasses)
    npasses = 3;
end
if nargin<7
    QN0 = [];
end
if nargin<6 || isempty(maxiter)
    maxiter = 1000;
end
if nargin<5 || isempty(tol)
    tol = 1e-8;
end
if nargin<4
    type = []; %#ok<NASGU>
end

[M,R]=size(L);
if isempty(Z)
    Z = zeros(1,R);
end
Z = sum(Z,1);
if isempty(L) || all(max(L)==0)
    X = N./Z;
    Q = zeros(M,R);
    U = zeros(M,R);
    W = zeros(M,R);
    C = zeros(1,R);
    for r=1:R
        for i=1:M
            U(i,r) = X(r)*L(i,r);
        end
    end
    totiter = 0;
    return
end

% Initialize, as Linearizer does, from Bard-Schweitzer at every population
Qs = zeros(M,R,1+R);
for s=0:R
    N_1 = oner(N,s);
    if isempty(QN0)
        [~,q] = pfqn_bs(L,N_1,Z);
    else
        [~,q] = pfqn_bs(L,N_1,Z,tol,maxiter,QN0);
    end
    Qs(:,:,1+s) = q;
end
Delta = zeros(M,R,R);   % Delta(i,r,c): the Delta^(r)_c term of station i
xi = zeros(M,R);        % xi(i,c) = xi_ck(N)

totiter = 0;
for I=1:npasses %#ok<NASGU>
    for s=0:R
        N_1 = oner(N,s);
        % xi at population N - 1_s, exactly; s = 0 leaves xi at N
        if s == 0
            xis = xi;
        else
            xis = xi - reshape(Delta(:,s,:),M,R);
        end
        [Qs(:,:,1+s),~,~,iter] = Core(L,M,R,N_1,Z,Qs(:,:,1+s),xis,tol,maxiter-totiter);
        totiter = totiter + iter;
    end
    % Refresh the Delta-terms, then aggregate them into xi once per pass
    for i=1:M
        for r=1:R
            if N(r)==1
                Qs(i,r,1+r) = 0;
            end
            for s=1:R
                Ns = oner(N,s);
                if N(r) > 0 && Ns(r) > 0
                    Delta(i,r,s) = Qs(i,r,1+s)/Ns(r) - Qs(i,r,1+0)/N(r);
                elseif N(r) > 0
                    Delta(i,r,s) = -Qs(i,r,1+0)/N(r);
                else
                    Delta(i,r,s) = 0;
                end
            end
        end
    end
    for c=1:R
        acc = zeros(M,1);
        for r=1:R
            acc = acc + max(N(r)-double(r==c),0)*Delta(:,r,c);
        end
        xi(:,c) = acc;
    end
end

[Q,W,X,iter] = Core(L,M,R,N,Z,Qs(:,:,1+0),xi,tol,maxiter-totiter);
totiter = totiter + iter;
U = zeros(M,R);
for i=1:M
    for r=1:R
        U(i,r)=X(r)*L(i,r);
    end
end
C = N./X-Z;
end

function [Q,W,T,iter] = Core(L,M,R,N_1,Z,Q,xi,tol,maxiter)
% Fixed point of the aggregated arrival-instant estimate with (2.2)-(2.6)
hasConverged = false;
W = L;
T = zeros(1,R);
iter = 0;
while ~hasConverged
    Qlast = Q;
    A = zeros(M,R);
    for c=1:R
        acc = zeros(M,1);
        for r=1:R
            if N_1(r) > 0
                nr = N_1(r)-double(r==c);
                if nr > 0
                    acc = acc + nr*Q(:,r)/N_1(r);
                end
            end
        end
        A(:,c) = acc + xi(:,c);
    end
    W = L .* (1 + A);
    for r=1:R
        if N_1(r) > 0
            T(r) = N_1(r) / (Z(r)+sum(W(:,r)));
        else
            T(r) = 0;
        end
        Q(:,r) = T(r) * W(:,r);
    end
    if enorm(Q-Qlast)<tol || iter > maxiter
        hasConverged = true;
    end
    iter = iter + 1;
end
end
