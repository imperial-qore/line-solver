%{
%{
 % @file pfqn_egflinearizer.m
 % @brief Extended generalized fixed-point Linearizer approximation.
%}
%}

%{
%{
 % @brief Extended generalized fixed-point Linearizer approximation.
 % @fn pfqn_egflinearizer(L, N, Z, type, tol, maxiter, alpha)
 % @param L Service demand matrix.
 % @param N Population vector.
 % @param Z Think time vector.
 % @param type Scheduling strategy type per station.
 % @param tol Convergence tolerance (default: 1e-8).
 % @param maxiter Maximum number of iterations (default: 1000).
 % @param alpha Per-class scaling exponent vector.
 % @return Q Mean queue lengths.
 % @return U Utilization.
 % @return W Waiting times.
 % @return C Cycle times.
 % @return X System throughput.
 % @return totiter Total iterations performed.
%}
%}
function [Q,U,W,C,X,totiter] = pfqn_egflinearizer(L,N,Z,type,tol,maxiter,alpha)
% Single-server version of linearizer

if nargin<5
    maxiter = 1000;
end
if nargin<4
    tol = 1e-8;
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

% Initialize
Q = zeros(M,R,1+R);
Delta = zeros(M,R,R);
for s=0:R
    N_1 = oner(N,s);
    [~,q] = pfqn_bs(L,N_1,Z);
    for r=1:R
        Q(:,r,1+s) = q(:,r);
    end
end

totiter = 0;
% Main loop
for I=1:3
    for s=0:R
        N_1 = oner(N,s); % for k=0 it just returns N
        % Core(N_1)
        [Q(:,:,1+s),~,~,iter] = Core(L,M,R,N_1,Z,Q(:,:,1+s),Delta,type,tol,maxiter-totiter,alpha);
        totiter = totiter + iter;
    end
    % Update_Delta
    for i=1:M
        for r=1:R
            if N(r)==1
                Q(i,:,1+r) = 0;
            end
            for s=1:R
                if N(s)>1
                    Ns = oner(N,s);
                    Delta(i,r,s) = Q(i,r,1+s)/Ns(r)^alpha(r) - Q(i,r,1+0)/N(r)^alpha(r);
                end
            end
        end
    end
end


% Core(N)
[Q,W,X,iter] = Core(L,M,R,N,Z,Q(:,:,1+0),Delta,type,tol,maxiter-totiter,alpha);
totiter = totiter + iter;
% Compute performance metrics
U = zeros(M,R);
for i=1:M
    for r=1:R
        U(i,r)=X(r)*L(i,r);
    end
end
Q = Q(1:M,1:R,1+0);
C = N./X-Z;
end

function [Q,W,T,iter] = Core(L,M,R,N_1,Z,Q,Delta,type,tol,maxiter,alpha)
hasConverged = false;
W = L;
iter = 0;
while ~hasConverged
    Qlast = Q;
    % Estimate population at
    Q_1 = Estimate(L,M,R,N_1,Z,Q,Delta,W,alpha);
    % Forward MVA
    [Q,W,T] = ForwardMVA(L,M,R,type,N_1,Z,Q_1);
    if enorm(Q-Qlast)<tol || iter > maxiter
        hasConverged = true;
    end
    iter = iter + 1;
end % it
end

function [Q_1,T_1] = Estimate(~,M,R,N_1,~,Q,Delta,~,alpha)
Q_1 = zeros(M,R);
T_1 = zeros(R,1+R);
for i=1:M
    for r=1:R
        for s=1:R
            Ns = oner(N_1,s);
            Q_1(i,r,1+s) = Ns(r)^alpha(r)*(Q(i,r,1+0)/N_1(r)^alpha(r) + Delta(i,r,s));
        end
    end
end

% This part is not used in Core so commented out
% for r=1:R
%     Nr = oner(N_1,r);
%     for s=1:R
%         % initial guess based on balanced job bound
%         % helpful in case no stations with positive demand exists
%         T_1(s,1+r) = Nr(s) / (Z(s) + max(L(:,s))*(sum(Nr)-1));
%         for i=1:M
%             if W(i,s,1+0)>0
%                 T_1(s,1+r) = Nr(s)*(Q(i,s)/N_1(s) + Delta(i,r,s))/W(i,s,1+0);
%                 break;
%             end
%         end
%     end
% end
end

function [Q,W,T] = ForwardMVA(L,M,R,type,N_1,Z,Q_1)
W = zeros(M,R);
T = zeros(1,R);
Q = zeros(M,R);

% Compute residence time
for ist=1:M
    for r=1:R
        if type(ist) == SchedStrategy.FCFS
            W(ist,r) = L(ist,r);
            if L(ist,r) ~= 0
                for s=1:R
                    W(ist,r) = W(ist,r) + L(ist,s)*Q_1(ist,s,1+r);
                end
            end
        else
            W(ist,r) = L(ist,r)*(1+sum(Q_1(ist,:,1+r)));
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
end


