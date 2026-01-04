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
 % @return Q Mean queue lengths.
 % @return U Utilization.
 % @return R Residence times.
 % @return C Cycle times.
 % @return X System throughput.
 % @return totiter Total number of iterations.
%}
%}
function [Q,U,R,C,X,totiter] = pfqn_conwayms(L,N,Z,nservers,type,tol,maxiter)
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

if isempty(Z)
    Z = zeros(1,R);
end

Z = sum(Z,1);
% Initialize
Q = zeros(M,R,1+R);
PB = zeros(M,1+R);
P = zeros(M,max(nservers),1+R);
Delta = zeros(M,R,R);
for i=1:M
    for r=1:R
        for s=1:R
            Delta(i,r,s) = 0;
        end
        for s=0:R
            N_1 = oner(N,s);
            Q(i,r,1+s) = N_1(r)/M;
        end
    end
end
for i=1:M
    for r=1:R
        for s=0:R
            N_1 = oner(N,s);
            pop = sum(N_1);
            if nservers(i)>1
                for j=1:(nservers(i)-1)
                    P(i,1+j,1+s) = 2*sum(Q(i,:,1+s))/(pop*(pop+1));
                end
                PB(i,1+s) = 2*sum(Q(i,:,1+s))/(pop+1-nservers(i))/(pop*(pop+1));
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
                if N(s)>2
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
iter = 1;
while ~hasConverged
    Qlast = Q;
    % Estimate population at
    [Q_1,P_1,PB_1,T_1] = Estimate(M,R,N_1,nservers,Q,P,PB,Delta,W);
    % Forward MVA
    [Q,W,T,P,PB] = ForwardMVA(L,M,R,N_1,Z,nservers,type,Q_1,P_1,PB_1,T_1);
    if norm(Q-Qlast)<tol || iter > maxiter
        hasConverged = true;
    end
    iter = iter + 1;
end % it
end

function [Q_1,P_1,PB_1,T_1] = Estimate(M,R,N_1,nservers,Q,P,PB,Delta,W)
P_1 = zeros(M,max(nservers),1+R);
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
            Q_1(i,r,1+s) = Ns(r)*(Q(i,r,1+0)/N_1(r) + Delta(i,r,s));
        end
    end
end
for r=1:R
    for s=1:R
        Nr = oner(N_1,r);
        for i=1:M
            if W(i,s,1+0)>0
                T_1(s,1+r) = Nr(s)*(Q(i,s,1+0)/N_1(s) + Delta(i,r,s))/W(i,s,1+0);
                break;
            end
        end
    end
end
end

function [Q,W,T,P,PB] = ForwardMVA(L,M,R,N_1,Z,nservers,type,Q_1,P_1,PB_1,T_1)
W = zeros(M,R);
T = zeros(1,R);
Q = zeros(M,R);
P = zeros(M,max(nservers));
PB = zeros(M,1);
XR = zeros(M,R);
C = zeros(M,R+1);
XE = zeros(M,R,R);

F = cell(1,R);
for r=1:R
    F{r} = zeros(M,R);    
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
            XR(ist,r) = XR(ist,r) / C(ist,1+r);
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
                XE(ist,r,c) = XE(ist,r,c) / Cx(ist,1+r);
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
% Compute marginal probabilities
for ist=1:M
    if nservers(ist) > 1
        P(ist,:) = 0;
        for j=1:(nservers(ist)-1)
            for c=1:R
                P(ist,1+j) = P(ist,1+j) + L(ist,c)*T(c)*P_1(ist,1+(j-1),1+c)/j;
            end
        end
    end
end
for ist=1:M
    if nservers(ist) > 1
        PB(ist) = 0;
        for c=1:R
            PB(ist) = PB(ist) + L(ist,c)*T(c)*(PB_1(ist,1+c)+P_1(ist,1+nservers(ist)-1,1+c))/nservers(ist);
        end
    end
end
for ist=1:M
    if nservers(ist) > 1
        P(ist,1+0) = max(0,1 - PB(ist));
        for j=1:(nservers(ist)-1)
            P(ist,1+0) = max(0,P(ist,1+0) - P(ist,1+j));
        end
    end
end
end