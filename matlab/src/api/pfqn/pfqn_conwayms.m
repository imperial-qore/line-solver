function [Q,U,R,C,X,totiter] = pfqn_conwayms(L,N,Z,nservers,type,tol,maxiter)
% Multiserver version of Linearizer as described in Conway 1989,  Fast
% Approximate Solution of Queueing Networks with Multi-Server Chain-
% Dependent FCFS Queues

[M,R]=size(L);
if nargin<5
    type = SchedStrategy.ID_FCFS * ones(M,1);
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
Q_1 = zeros(M,R);
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
        T_1(s,1+r) = Nr(s)*(Q(i,s)/N_1(s) + Delta(i,r,s))/W(i,s,1+0);
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
XE = zeros(M,R,R);

% Compute XR
mu = 1./L;
for i=1:M
    for r=1:R
        if nservers(i) > 1
            XR(i,r) = 0;
            C(i,1+r) = 0;
            F = T_1(:,1+r).*L(i,:)/sum(T_1(:,1+r).*L(i,:));
            [k,n,S,D]=sprod(R,nservers(i));
            while k>=0
                if all(n<= oner(N_1,r)) % Br set
                    n = n(:)';
                    Ai = exp(multinomialln(n) + n*log(F(:)));
                    C(i,1+r) = C(i,1+r) + Ai;
                    XR(i,r) = XR(i,r) + Ai*(mu(i,:)*n(:))^(-1);
                end
                [k,n]=sprod(k,S,D);
            end
            XR(i,r) = XR(i,r) / C(i,1+r);
        end
    end
end

% Compute XE
Cx = zeros(M,1+R);
for i=1:M
    for r=1:R
        if nservers(i) > 1
            for s=1:R
                XE(i,r,s) = 0;
                Cx(i,1+r) = 0;
                F = T_1(:,1+r).*L(i,:)/sum(T_1(:,1+r).*L(i,:));
                [k,n,S,D]=sprod(R,nservers(i));
                while k>=0
                    if all(n<= oner(N_1,r) & n(s)>=1) % Axr set
                        n = n(:)';
                        Aix = exp(multinomialln(n) + n*log(F(:)));
                        Cx(i,1+r) = Cx(i,1+r) + Aix;
                        XE(i,r,s) = XE(i,r,s) + Aix*(mu(i,:)*n(:))^(-1);
                    end
                    [k,n]=sprod(k,S,D);
                end
                XE(i,r,s) = XE(i,r,s) / Cx(i,1+r);
            end
        end
    end
end

% Compute residence time
for i=1:M
    for r=1:R
        if nservers(i) == 1
            if type == SchedStrategy.ID_FCFS
                W(i,r) = L(i,r);
                for s=1:R
                    W(i,r) = W(i,r) + L(i,s)*Q_1(i,s,1+r);
                end
            else
                W(i,r) = L(i,r);
                for s=1:R
                    W(i,r) = W(i,r) + L(i,r)*Q_1(i,s,1+r);
                end
            end
        else
            W(i,r) = L(i,r) + PB_1(i,1+r)*XR(i,r);
            for s=1:R
                W(i,r) = W(i,r) + XE(i,r,s)*(Q_1(i,s,1+r)-L(i,s)*T_1(s,1+r));
            end
        end
    end
end
% Compute throughputs and qlens
for r=1:R
    T(r) = N_1(r) / (Z(r)+sum(W(:,r)));
    for i=1:M
        Q(i,r) = T(r) * W(i,r);
    end
end
% Compute marginal probabilities
for i=1:M
    if nservers(i) > 1
        P(i,:) = 0;
        for j=1:(nservers(i)-1)
            for s=1:R
                P(i,1+j) = P(i,1+j) + L(i,s)*T(s)*P_1(i,1+(j-1),1+s)/j;
            end
        end
    end
end
for i=1:M
    if nservers(i) > 1
        PB(i) = 0;
        for s=1:R
            PB(i) = PB(i) + L(i,s)*T(s)*(PB_1(i,1+s)+P_1(i,1+nservers(i)-1,1+s))/nservers(i);
        end
    end
end
for i=1:M
    if nservers(i) > 1
        P(i,1+0) = max(0,1 - PB(i));
        for j=1:(nservers(i)-1)
            P(i,1+0) = max(0,P(i,1+0) - P(i,1+j));
        end
    end
end
end