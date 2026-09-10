%{
%{
 % @file pfqn_conwayms.m
 % @brief Multiserver Linearizer approximation (Conway 1989).
%}
%}

%{
%{
 % @brief Multiserver Linearizer approximation (Conway 1989).
 % @fn pfqn_conwayms(L, N, Z, nservers, type, tol, maxiter)
 % @param L Service demand matrix.
 % @param N Population vector.
 % @param Z Think time vector.
 % @param nservers Number of servers per station.
 % @param type Scheduling strategy type per station (default: FCFS).
 % @param tol Convergence tolerance (default: 1e-8).
 % @param maxiter Maximum number of iterations (default: 1000).
 % @param QN0 (M x R) queue lengths that warm-start the Bard-Schweitzer initialization; empty for the default cold start.
 % @return Q Mean queue lengths.
 % @return U Utilization.
 % @return R Residence times.
 % @return C Cycle times.
 % @return X System throughput.
 % @return totiter Total number of iterations.
%}
%}
function [Q,U,R,C,X,totiter] = pfqn_conwayms(L,N,Z,nservers,type,tol,maxiter,QN0)
% Multiserver version of Linearizer as described in Conway 1989,  Fast
% Approximate Solution of Queueing Networks with Multi-Server Chain-
% Dependent FCFS Queues

[M,R]=size(L);
if nargin<5
    type = SchedStrategy.FCFS * ones(M,1);
end
if nargin<6
    tol = 1e-8;
end
if nargin<7
    maxiter = 1000;
end
if nargin<8
    QN0 = [];
end

if isempty(Z)
    Z = zeros(1,R);
end

Z = sum(Z,1);
% Initialize
Q = zeros(M,R,1+R);
PB = zeros(M,1+R);
P = zeros(M,max(nservers(:)),1+R);
Delta = zeros(M,R,R);
for i=1:M
    for r=1:R
        for s=1:R
            Delta(i,r,s) = 0;
        end
        for s=0:R
            N_1 = oner(N,s);
            if isempty(QN0)
                Q(i,r,1+s) = N_1(r)/M;
            else
                Q(i,r,1+s) = QN0(i,r);   % warm start from supplied queue lengths
            end
        end
    end
end
for i=1:M
    for r=1:R
        for s=0:R
            N_1 = oner(N,s);
            pop = sum(N_1);
            if nservers(i)>1
                if pop == 0
                    % empty network: the station is idle with probability one
                    P(i,1+(1:(nservers(i)-1)),1+s) = 0;
                    PB(i,1+s) = 0;
                    P(i,1+0,1+s) = 1;
                    continue
                end
                for j=1:(nservers(i)-1)
                    P(i,1+j,1+s) = 2*sum(Q(i,:,1+s))/(pop*(pop+1));
                end
                if pop > nservers(i)-1
                    PB(i,1+s) = 2*sum(Q(i,:,1+s))/(pop+1-nservers(i))/(pop*(pop+1));
                else % fewer jobs than servers: they cannot all be busy
                    PB(i,1+s) = 0;
                end
                P(i,1+0,1+s) = 1 - PB(i,1+s) - sum(P(i,1+(1:(nservers(i)-1)),1+s));
            end
        end
    end
end

totiter = 0;
% Main loop
for I=1:2
    for s=0:R
        N_1 = oner(N,s); % for k=0 it just returns N
        % Core(N_1)
        [Q(:,:,1+s),~,~,P(:,:,1+s),PB(:,1+s),iter] = Core(L,M,R,N_1,Z,nservers,Q(:,:,1+s),P(:,:,1+s),PB(:,1+s),Delta,type,tol,maxiter-totiter);
        totiter = totiter + iter;
    end
    % Update_Delta
    for i=1:M
        for r=1:R
            for s=1:R
                Ns = oner(N,s);
                if N(s)>2 && N(r)>0 % an empty class has no F_ir to correct
                    Delta(i,r,s) = Q(i,r,1+s)/Ns(r) - Q(i,r,1+0)/N(r);
                end
            end
        end
    end
end

% Core(N)
[Q,W,X,~,~,iter] = Core(L,M,R,N,Z,nservers,Q(:,:,1+0),P(:,:,1+0),PB(:,1+0),Delta,type,tol,maxiter);
totiter = totiter + iter;
% Compute performance metrics
U = zeros(M,R);
for i=1:M
    for r=1:R
        if nservers(i)==1
            U(i,r)=X(r)*L(i,r);
        else
            U(i,r)=X(r)*L(i,r) / nservers(i);
        end
    end
end

Q = Q(1:M,1:R,1+0);
C = N./X-Z;
R = W;
end

function [Q,W,T,P,PB,iter] = Core(L,M,R,N_1,Z,nservers,Q,P,PB,Delta,type,tol,maxiter)
hasConverged = false;
W = L;
Wlast = [];
T = zeros(1,R);
iter = 1;
while ~hasConverged
    Qlast = Q;
    % Estimate population at
    [Q_1,P_1,PB_1,T_1] = Estimate(M,R,N_1,nservers,Q,P,PB,Delta,W);
    % Forward MVA
    [Q,W,T,P,PB] = ForwardMVA(L,M,R,N_1,Z,nservers,type,Q_1,P_1,PB_1,T_1);
    % W must enter the test: Q alone is satisfied on the FIRST sweep whenever Q
    % cannot move (M=1 seeds Q at its own fixed point), and the residence times
    % returned then are still the seed W=L, so T_1=Q/W is unbounded and the
    % throughput exceeds the station's own service capacity.
    if isempty(Wlast)
        moved = Inf;
    else
        moved = max(norm(Q-Qlast), norm(W-Wlast));
    end
    Wlast = W;
    if moved < tol || iter > maxiter
        hasConverged = true;
    end
    iter = iter + 1;
end % it
end

function [Q_1,P_1,PB_1,T_1] = Estimate(M,R,N_1,nservers,Q,P,PB,Delta,W)
P_1 = zeros(M,max(nservers(:)),1+R);
PB_1 = zeros(M,1+R);
Q_1 = zeros(M,R,1+R);
T_1 = zeros(R,1+R);
for i=1:M
    if nservers(i)>1
        for j=0:(nservers(i)-1)
            for s=0:R
                P_1(i,1+j,1+s) = P(i,1+j);
            end
        end
        for s=0:R
            PB_1(i,1+s) = PB(i,1);
        end
    end
    for r=1:R
        for s=1:R
            Ns = oner(N_1,s);
            if N_1(r) > 0
                Q_1(i,r,1+s) = Ns(r)*(Q(i,r,1+0)/N_1(r) + Delta(i,r,s));
            else % a class with no jobs left has an empty queue everywhere
                Q_1(i,r,1+s) = 0;
            end
        end
    end
end
% T_1 is Little's law over the queueing part of the cycle, sum_i Q_1 / sum_i W,
% and not the ratio at the FIRST station with a positive residence time: the
% per-station estimates disagree, so picking one made the answer depend on the
% station order. The demand matrix carries no order, so a model symmetric under
% permuting classes and stations together must return equal class throughputs,
% and with the single-station pick it did not.
for r=1:R
    for s=1:R
        Nr = oner(N_1,r);
        num = 0; den = 0;
        for i=1:M
            if W(i,s,1+0)>0 && N_1(s)>0 % a class with no jobs left has no throughput
                % Delta is indexed (station, queued class, removed class), as the
                % Q_1 loop above uses it: here class s queues and class r is removed
                num = num + Nr(s)*(Q(i,s,1+0)/N_1(s) + Delta(i,s,r));
                den = den + W(i,s,1+0);
            end
        end
        if den > 0
            T_1(s,1+r) = max(0, num/den);
        end
    end
end
end

function [Q,W,T,P,PB] = ForwardMVA(L,M,R,N_1,Z,nservers,type,Q_1,P_1,PB_1,T_1)
W = zeros(M,R);
T = zeros(1,R);
Q = zeros(M,R);
P = zeros(M,max(nservers(:)));
PB = zeros(M,1);
XR = zeros(M,R);
C = zeros(M,R+1);
XE = zeros(M,R,R);

F = cell(1,R);
for r=1:R
    F{r} = zeros(M,R);
end
for r=1:R
    for ist=1:M
        den = (L(ist,:)*T_1(:,1+r));
        for c=1:R
            F{r}(ist,c) = T_1(c,1+r)*L(ist,c)/den;
        end
    end
end
% Compute XR
mu = 1./L;
for ist=1:M
    for r=1:R
        if nservers(ist) > 1
            XR(ist,r) = 0;
            C(ist,1+r) = 0;
            [s,n,S,D]=sprod(R,nservers(ist));
            while s>=0
                if all(n(:)'<=oner(N_1,r)) % Br set
                    n = n(:)';
                    Ai = exp(multinomialln(n) + n*log(F{r}(ist,:)'));
                    C(ist,1+r) = C(ist,1+r) + Ai;
                    XR(ist,r) = XR(ist,r) + Ai*(mu(ist,:)*n(:))^(-1);
                end
                [s,n]=sprod(s,S,D);
            end
            if C(ist,1+r) > 0
                XR(ist,r) = XR(ist,r) / C(ist,1+r);
            else % Br empty: fewer jobs than servers, so all of them can never be busy
                XR(ist,r) = 0;
            end
        end
    end
end

% Compute XE
Cx = zeros(M,1+R);
for ist=1:M
    for r=1:R
        if nservers(ist) > 1
            for c=1:R
                XE(ist,r,c) = 0;
                Cx(ist,1+r) = 0;
                [s,n,S,D]=sprod(R,nservers(ist));
                while s>=0
                    if all(n(:)'<= oner(N_1,r) & n(c)>=1) % Axr set
                        n = n(:)';
                        Aix = exp(multinomialln(n) + n*log(F{r}(ist,:)'));
                        Cx(ist,1+r) = Cx(ist,1+r) + Aix;
                        XE(ist,r,c) = XE(ist,r,c) + Aix*(mu(ist,:)*n(:))^(-1);
                    end
                    [s,n]=sprod(s,S,D);
                end
                if Cx(ist,1+r) > 0
                    XE(ist,r,c) = XE(ist,r,c) / Cx(ist,1+r);
                else % Axr empty: class c has no job left in N_1-e_r, so its term is 0
                    XE(ist,r,c) = 0;
                end
            end
        end
    end
end

% Compute residence time
for ist=1:M
    for r=1:R
        if nservers(ist) == 1
            if type == SchedStrategy.FCFS
                W(ist,r) = L(ist,r);
                for c=1:R
                    W(ist,r) = W(ist,r) + L(ist,c)*Q_1(ist,c,1+r);
                end
            else
                W(ist,r) = L(ist,r);
                for c=1:R
                    W(ist,r) = W(ist,r) + L(ist,r)*Q_1(ist,c,1+r);
                end
            end
        else
            W(ist,r) = L(ist,r) + PB_1(ist,1+r)*XR(ist,r);
            for c=1:R
                W(ist,r) = W(ist,r) + XE(ist,r,c)*(Q_1(ist,c,1+r)-L(ist,c)*T_1(c,1+r));
            end
        end
    end
end
% Compute throughputs and qlens
for r=1:R
    T(r) = N_1(r) / (Z(r)+sum(W(:,r)));
    for ist=1:M
        Q(ist,r) = T(r) * W(ist,r);
    end
end
% Queue-length marginals. The relations
%   p_j = A*p_{j-1}/j,  pB = A*(pB + p_{ms-1})/ms,  p_0 = 1 - pB - sum_j p_j
% with A = sum_s X_s*L_is the mean number of busy servers are solved in closed
% form rather than iterated. As a Jacobi iteration they amplify by A per sweep,
% and since the convergence test watches Q and W but not P the routine returned
% marginals whose mass had run to 334 behind the p_0 = max(0,1-...) floor.
% Estimate hands the same marginals to every reduced population, so the
% population corrections that pfqn_linearizerms carries here are all zero.
for ist=1:M
    ms = nservers(ist);
    if ms > 1 && isfinite(ms)
        A = 0;
        for s=1:R
            A = A + L(ist,s)*T(s);
        end
        if A >= ms
            % Saturated: the closed form is singular and its limit is the
            % degenerate marginal, every server busy with probability one.
            % N = m with Z = 0 reaches it exactly, so this is a legal input.
            for j=0:(ms-1)
                P(ist,1+j) = 0;
            end
            PB(ist) = 1;
        else
            % p_j = alpha(1+j)*p_0
            alpha = zeros(1,ms);
            alpha(1) = 1;
            for j=1:(ms-1)
                alpha(1+j) = A*alpha(j)/j;
            end
            alphaB = A*alpha(ms)/(ms-A);
            P(ist,1+0) = 1/(1 + sum(alpha(2:ms)) + alphaB);
            for j=1:(ms-1)
                P(ist,1+j) = alpha(1+j)*P(ist,1+0);
            end
            PB(ist) = alphaB*P(ist,1+0);
        end
    end
end
end