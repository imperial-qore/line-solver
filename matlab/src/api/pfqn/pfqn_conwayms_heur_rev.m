function [Q,U,W,C,X,totiter] = pfqn_conwayms_heur_rev(L,N,Z,nservers,type,tol,maxiter)
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
if isempty(L) || all(max(L)==0)
    X = N./Z;
    Q = zeros(1,R);
    U = zeros(M,R);
    W = zeros(M,R);
    for r=1:R
        for i=1:M
            if nservers(i)==1
                U(i,r)=X(r)*L(i,r);
            else
                U(i,r)=X(r)*L(i,r) / nservers(i);
            end
            W(i,r) = Q(i,r) / X(r);
        end
    end
    totiter = 0;
    return
end

% Initialize
Q = zeros(M,R,1+R);
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

totiter = 0;
% Main loop
for I=1:3
    for s=0:R
        N_1 = oner(N,s); % for k=0 it just returns N
        % Core(N_1)
        if I==1
            [Q(:,:,1+s),~,~,iter] = Core(L,M,R,N_1,Z,ones(1,M),Q(:,:,1+s),Delta,type,tol,maxiter-totiter);
        else
            [Q(:,:,1+s),~,~,iter] = Core(L,M,R,N_1,Z,nservers,Q(:,:,1+s),Delta,type,tol,maxiter-totiter);
        end
        totiter = totiter + iter;
    end
    % Update_Delta
    for i=1:M
        for r=1:R
            if N(r)==1
                Q(i,:,1+r) = 0;
            end
            for s=1:R
                if N(s)>2
                    Ns = oner(N,s);
                    Delta(i,r,s) = Q(i,r,1+s)/Ns(r) - Q(i,r,1+0)/N(r);
                end
            end
        end
    end
end


% Core(N)
[Q,W,X,iter] = Core(L,M,R,N,Z,nservers,Q(:,:,1+0),Delta,type,tol,maxiter);
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
end

function [Q,W,T,iter] = Core(L,M,R,N_1,Z,nservers,Q,Delta,type,tol,maxiter)
hasConverged = false;
W = L;
iter = 0;
while ~hasConverged
    Qlast = Q;
    % Estimate population at
    [Q_1,T_1] = Estimate(L,M,R,N_1,Z,nservers,Q,Delta,W);
    % Forward MVA
    [Q,W,T] = ForwardMVA(L,M,R,N_1,Z,nservers,type,Q_1,T_1);
    if norm(Q-Qlast)<tol || iter > maxiter
        hasConverged = true;
    end
    iter = iter + 1;
end % it
end

function [Q_1,T_1] = Estimate(L,M,R,N_1,Z,~,Q,Delta,W)
Q_1 = zeros(M,R);
T_1 = zeros(R,1+R);
for i=1:M
    for r=1:R
        for s=1:R
            Ns = oner(N_1,s);
            Q_1(i,r,1+s) = Ns(r)*(Q(i,r,1+0)/N_1(r) + Delta(i,r,s));
        end
    end
end

for r=1:R
    Nr = oner(N_1,r);
    for s=1:R
        % initial guess based on balanced job bound
        % helpful in case no stations with positive demand exists
        T_1(s,1+r) = Nr(s) / (Z(s) + max(L(:,s))*(sum(Nr)-1));
        for i=1:M
            if W(i,s,1+0)>0
                T_1(s,1+r) = Nr(s)*(Q(i,s)/N_1(s) + Delta(i,r,s))/W(i,s,1+0);
                break;
            end
        end
    end
end
end

function [Q,W,T] = ForwardMVA(L,M,R,N_1,Z,nservers,type,Q_1,T_1)
W = zeros(M,R);
T = zeros(1,R);
Q = zeros(M,R);
XR = zeros(M,R);
XE = zeros(M,R,R);

% Compute XR
mu = 1./L;
for i=1:M
    for r=1:R
        if nservers(i) > 1
            F = T_1(:,1+r).*L(i,:)/sum(T_1(:,1+r).*L(i,:));
            N_r = oner(N_1,r)';
            QLen = F' .* N_r' ./ (ones(1,length(F))+F');
            XR(i,r)=1 / (mu(i,:)*QLen(:));
        end
    end
end

% Compute XE
for i=1:M
    for r=1:R
        if nservers(i) > 1
            F = T_1(:,1+r).*L(i,:)/sum(T_1(:,1+r).*L(i,:));
            for s=1:R
                N_sr = oner(oner(N_1,s),r)';
                QLen_s = F' .* N_sr' ./ (ones(1,length(F))+F');
                XE(i,r,s)=1 / (mu(i,s)+mu(i,:)*QLen_s(:));
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
            Ui = min(1,L(i,:)*T_1(:,1+r)/nservers(i));
            PBir = Ui.^nservers(i); % blocking probability as in Rolia-Sevcik
            W(i,r) = L(i,r) + 0*PBir*XR(i,r);
            for s=1:R
                W(i,r) = W(i,r) + XE(i,r,s)*max(0,(Q_1(i,s,1+r)-0*L(i,s)*T_1(s,1+r))/nservers(i));
            end
        end
    end
end

% Compute throughputs and qlens
for r=1:R
    T(r) = N_1(r) / ( Z(r) + sum(W(:,r)) );
    for i=1:M
        Q(i,r) = T(r) * W(i,r);
    end
end
end

% % Compute XE
% for i=1:M
%     for r=1:R
%         if nservers(i) > 1
%             F = T_1(:,1+r).*L(i,:)/sum(T_1(:,1+r).*L(i,:));
%             for s=1:R
%                     N_sr = oner(oner(N_1,s),r)';
%                     nnzClasses = find(N_sr>0);
%                     QLen_s = lossn_erlangfp(F(nnzClasses)',eye(length(nnzClasses)),N_sr(nnzClasses));
%                 if s~=r
%                     XE(i,r,s)=1/(mu(i,s)+mu(i,r)+mu(i,nnzClasses)*QLen_s(:));
%                 else
%                     XE(i,r,r)=0; % 1/(mu(i,s)+mu(i,r)+mu(i,nnzClasses)*QLen_s(:));
%                 end
%             end
%         end
%     end
% end
% end