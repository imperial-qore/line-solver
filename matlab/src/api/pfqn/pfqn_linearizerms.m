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
 % @return Q Mean queue lengths.
 % @return U Utilization.
 % @return R Residence times.
 % @return C Cycle times.
 % @return X System throughput.
 % @return totiter Total iterations performed.
%}
%}
function [Q,U,R,C,X,totiter] = pfqn_linearizerms(L,N,Z,nservers,type,tol,maxiter)
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

if isempty(Z)
    Z = zeros(1,R);
end

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
        [~,q] = pfqn_bs(L,N_1,Z);
        Q(:,r,1+s) = q(:);
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
                Delta(i,r,s) = Q(i,r,1+s)/Ns(r) - Q(i,r,1+0)/N(r);
            end
        end
    end
end

% Core(N)
[Q,W,X,~,~,iter] = Core(L,M,R,N,Z,nservers,Q(:,:,1+0),P(:,:,1+0),PB(:,1+0),Delta,type,tol,maxiter-totiter);
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
iter = 0;
hasConverged = false;
while ~hasConverged
    iter = iter + 1;
    Qlast = Q;
    % Estimate population at
    [Q_1,P_1,PB_1] = Estimate(M,R,N_1,nservers,Q,P,PB,Delta);
    % Forward MVA
    [Q,W,T,P,PB] = ForwardMVA(L,M,R,N_1,Z,nservers,type,Q_1,P_1,PB_1);
    if norm(Q-Qlast)<tol || iter > maxiter
        hasConverged = true;
    end
end % it
end

function [Q_1,P_1,PB_1] = Estimate(M,R,N_1,nservers,Q,P,PB,Delta)
P_1 = zeros(M,max(nservers),1+R);
PB_1 = zeros(M,1+R);
Q_1 = zeros(M,R);
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
end

function [Q,W,T,P,PB] = ForwardMVA(L,M,R,N_1,Z,nservers,type,Q_1,P_1,PB_1)
W = zeros(M,R);
T = zeros(1,R);
Q = zeros(M,R);
P = zeros(M,max(nservers));
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
        if nservers(ist) > 1
            for j=0:(nservers(ist)-2)
                if type == SchedStrategy.FCFS
                    for s=1:R
                        W(ist,r) = W(ist,r) + L(ist,s)*(nservers(ist)-1-j)*P_1(ist,1+j,1+r);
                    end
                else
                    for s=1:R
                        W(ist,r) = W(ist,r) + L(ist,r)*(nservers(ist)-1-j)*P_1(ist,1+j,1+r);
                    end
                end

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
for ist=1:M
    if nservers(ist) > 1
        P(ist,:) = 0;
        for j=1:(nservers(ist)-1)
            for s=1:R
                P(ist,1+j) = P(ist,1+j) + L(ist,s)*T(s)*P_1(ist,1+(j-1),1+s)/j;
            end
        end
    end
end
for ist=1:M
    if nservers(ist) > 1
        PB(ist) = 0;
        for s=1:R
            PB(ist) = PB(ist) + L(ist,s)*T(s)*(PB_1(ist,1+s)+P_1(ist,1+nservers(ist)-1,1+s))/nservers(ist);
        end
    end
end
for ist=1:M
    if nservers(ist) > 1
        P(ist,1+0) = 1 - PB(ist);
        for j=1:(nservers(ist)-1)
            P(ist,1+0) = P(ist,1+0) - P(ist,1+j);
        end
    end
end
end
