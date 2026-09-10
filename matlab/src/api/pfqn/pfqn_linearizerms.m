%{
%{
 % @file pfqn_linearizerms.m
 % @brief Multiserver Linearizer (Krzesinski/Conway/De Souza-Muntz).
%}
%}

%{
%{
 % @brief Multiserver Linearizer (Krzesinski/Conway/De Souza-Muntz).
 % @fn pfqn_linearizerms(L, N, Z, nservers, type, tol, maxiter)
 % @param L Service demand matrix.
 % @param N Population vector.
 % @param Z Think time vector.
 % @param nservers Number of servers per station.
 % @param type Scheduling strategy per station (default: PS).
 % @param tol Convergence tolerance (default: 1e-8).
 % @param maxiter Maximum number of iterations (default: 1000).
 % @param QN0 (M x R) queue lengths that warm-start the Bard-Schweitzer initialization; empty for the default cold start.
 % @return Q Mean queue lengths.
 % @return U Utilization.
 % @return R Residence times.
 % @return C Cycle times.
 % @return X System throughput.
 % @return totiter Total iterations performed.
%}
%}
function [Q,U,R,C,X,totiter] = pfqn_linearizerms(L,N,Z,nservers,type,tol,maxiter,QN0)
% Multiserver version of Krzesinski's Linearizer as described in Conway
% 1989,  Fast Approximate Solution of Queueing Networks with Multi-Server
% Chain- Dependent FCFS Queues.
% Some minor adjustments based on De Souza-Muntz's description of the
% algorithm.

[M,R]=size(L);
if nargin<5
    type = SchedStrategy.PS * ones(M,1);
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

% Initialize
Q = zeros(M,R,1+R);
PB = zeros(M,1+R);
P = zeros(M,max(nservers(:)),1+R);
Delta = zeros(M,R,R);
% Linearizer corrections for the queue-length marginals. Without them the
% marginals stay at population N while the queue lengths are reduced to N-e_s,
% which breaks Q + sum_j (m-1-j) p_j >= m-1 and lets W fall below the mean
% service time.
DeltaP = zeros(M,max(nservers(:)),R);
DeltaPB = zeros(M,R);
for i=1:M
    for r=1:R
        for s=1:R
            Delta(i,r,s) = 0;
        end
        for s=0:R

        N_1 = oner(N,s);
        if isempty(QN0)
            [~,q] = pfqn_bs(L,N_1,Z);
        else
            [~,q] = pfqn_bs(L,N_1,Z,tol,maxiter,QN0);   % warm start
        end
        Q(:,r,1+s) = q(:,r);
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
        [Q(:,:,1+s),~,~,P(:,:,1+s),PB(:,1+s),iter] = Core(L,M,R,N_1,Z,nservers,Q(:,:,1+s),P(:,:,1+s),PB(:,1+s),Delta,DeltaP,DeltaPB,type,tol,maxiter-totiter);
        totiter = totiter + iter;
    end
    % Update_Delta
    for i=1:M
        for r=1:R
            for s=1:R
                Ns = oner(N,s);
                if Ns(r) > 0
                    Delta(i,r,s) = Q(i,r,1+s)/Ns(r) - Q(i,r,1+0)/N(r);
                else % Chandy-Neuse 0/0 convention: F_ir(N-e_s) = 0
                    Delta(i,r,s) = -Q(i,r,1+0)/N(r);
                end
            end
        end
    end
    % Update_DeltaP: probabilities do not scale with the population, so the
    % analogue of Delta is a plain difference
    for i=1:M
        if nservers(i) > 1 && isfinite(nservers(i))
            for s=1:R
                for j=0:(nservers(i)-1)
                    DeltaP(i,1+j,s) = P(i,1+j,1+s) - P(i,1+j,1+0);
                end
                DeltaPB(i,s) = PB(i,1+s) - PB(i,1+0);
            end
        end
    end
end

% Core(N)
[Q,W,X,~,~,iter] = Core(L,M,R,N,Z,nservers,Q(:,:,1+0),P(:,:,1+0),PB(:,1+0),Delta,DeltaP,DeltaPB,type,tol,maxiter-totiter);
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

function [Q,W,T,P,PB,iter] = Core(L,M,R,N_1,Z,nservers,Q,P,PB,Delta,DeltaP,DeltaPB,type,tol,maxiter)
iter = 0;
W = zeros(M,R);
T = zeros(1,R);
hasConverged = false;
while ~hasConverged
    iter = iter + 1;
    Qlast = Q;
    % Estimate population at
    [Q_1,P_1,PB_1] = Estimate(M,R,N_1,nservers,Q,P,PB,Delta,DeltaP,DeltaPB);
    % Forward MVA
    [Q,W,T,P,PB] = ForwardMVA(L,M,R,N_1,Z,nservers,type,Q_1,P_1,PB_1);
    if norm(Q-Qlast)<tol || iter > maxiter
        hasConverged = true;
    end
end % it
end

function [Q_1,P_1,PB_1] = Estimate(M,R,N_1,nservers,Q,P,PB,Delta,DeltaP,DeltaPB)
P_1 = zeros(M,max(nservers(:)),1+R);
PB_1 = zeros(M,1+R);
% (M,R,1+R), not (M,R): the writes below are 3-D. The interpreter grew the array
% on first assignment, codegen rejects that, and the MEX build was broken on it.
Q_1 = zeros(M,R,1+R);
for i=1:M
    if nservers(i)>1
        for j=0:(nservers(i)-1)
            P_1(i,1+j,1+0) = P(i,1+j);
            for s=1:R
                P_1(i,1+j,1+s) = P(i,1+j) + DeltaP(i,1+j,s);
            end
        end
        PB_1(i,1+0) = PB(i,1);
        for s=1:R
            PB_1(i,1+s) = PB(i,1) + DeltaPB(i,s);
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
end

function [Q,W,T,P,PB] = ForwardMVA(L,M,R,N_1,Z,nservers,type,Q_1,P_1,PB_1)
W = zeros(M,R);
T = zeros(1,R);
Q = zeros(M,R);
P = zeros(M,max(nservers(:)));
PB = zeros(M,1);
for ist=1:M
    for r=1:R
        W(ist,r) = L(ist,r)/nservers(ist);
        if L(ist,r) == 0
            % 0 service demand at this station => this class does not visit the current node
            continue;
        end
        if type == SchedStrategy.FCFS
            for s=1:R
                W(ist,r) = W(ist,r) + (L(ist,s)/nservers(ist))*Q_1(ist,s,1+r);
            end
        else
            for s=1:R
                W(ist,r) = W(ist,r) + (L(ist,r)/nservers(ist))*Q_1(ist,s,1+r);
            end
        end
        % Partially-idle-server correction. It compensates the 1/m scaling of the
        % arriving job's OWN service, so it carries L(ist,r)/m and no sum over the
        % other classes: at N=e_r the terms must collapse to W = L(ist,r).
        if nservers(ist) > 1
            for j=0:(nservers(ist)-2)
                W(ist,r) = W(ist,r) + (L(ist,r)/nservers(ist))*(nservers(ist)-1-j)*P_1(ist,1+j,1+r);
            end
        end
    end
end
for r=1:R
    T(r) = N_1(r) / (Z(r)+sum(W(:,r)));
    for ist=1:M
        Q(ist,r) = T(r) * W(ist,r);
    end
end
% Queue-length marginals. The relations
%   p_j = (A*p_{j-1} + d_{j-1})/j,  pB = (A*(pB + p_{ms-1}) + dB)/ms,
%   p_0 = 1 - pB - sum_j p_j
% with A = sum_s X_s*L_is the mean number of busy servers and d the population
% corrections, are solved in closed form rather than iterated: as a Jacobi
% iteration they amplify by A per sweep and diverge once A approaches ms.
for ist=1:M
    ms = nservers(ist);
    if ms > 1 && isfinite(ms)
        A = 0; d = zeros(1,ms); dB = 0;
        for s=1:R
            a_s = L(ist,s)*T(s);
            A = A + a_s;
            for j=0:(ms-1)
                d(1+j) = d(1+j) + a_s*(P_1(ist,1+j,1+s) - P_1(ist,1+j,1+0));
            end
            dB = dB + a_s*(PB_1(ist,1+s) - PB_1(ist,1+0));
        end
        if A >= ms
            % int32 casts: codegen rejects %d on a double
            line_error(mfilename,sprintf(['Station %d offers %g busy servers out of %d: the model is ' ...
                'saturated and its queue-length marginals do not exist.'],int32(ist),A,int32(ms)));
        end
        % p_j = alpha(1+j)*p_0 + beta(1+j)
        alpha = zeros(1,ms); beta = zeros(1,ms);
        alpha(1) = 1;
        for j=1:(ms-1)
            alpha(1+j) = A*alpha(j)/j;
            beta(1+j) = (A*beta(j) + d(j))/j;
        end
        alphaB = A*alpha(ms)/(ms-A);
        betaB = (A*beta(ms) + dB + d(ms))/(ms-A);
        P(ist,1+0) = (1 - sum(beta(2:ms)) - betaB)/(1 + sum(alpha(2:ms)) + alphaB);
        for j=1:(ms-1)
            P(ist,1+j) = alpha(1+j)*P(ist,1+0) + beta(1+j);
        end
        PB(ist) = alphaB*P(ist,1+0) + betaB;
    end
end
end
