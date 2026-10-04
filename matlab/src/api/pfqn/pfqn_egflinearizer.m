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
 % @param tol Convergence tolerance (default: 1e-8); 'cn' or NaN selects the Chandy-Neuse (1982) population-scaled termination test, see pfqn_cntol.
 % @param maxiter Maximum number of iterations (default: 1000).
 % @param alpha Per-class scaling exponent vector.
 % @param QN0 (M x R) queue lengths that warm-start the Bard-Schweitzer initialization; empty for the default cold start.
 % @param npasses Number of Delta refresh passes (default 3, the Chandy-Neuse rule; pfqn_scat sets 1).
 % @return Q Mean queue lengths.
 % @return U Utilization.
 % @return W Waiting times.
 % @return C Cycle times.
 % @return X System throughput.
 % @return totiter Total iterations performed.
%}
%}
function [Q,U,W,C,X,totiter] = pfqn_egflinearizer(L,N,Z,type,tol,maxiter,alpha,QN0,npasses)
% Single-server version of linearizer

if nargin<9 || isempty(npasses)
    npasses = 3; % Chandy-Neuse (1982) fixed three-iteration rule; pfqn_scat passes 1
end
if nargin<8
    QN0 = [];
end
if nargin<5
    maxiter = 1000;
end
if nargin<4
    tol = 1e-8;
end

% tol = 'cn' (or NaN) selects the published Linearizer termination test of
% Chandy and Neuse, Commun. ACM 25(2), 1982, p.129: Core stops when
% max_{i,r}|dQ(i,r)|/N_r < pfqn_cntol(N), where N is the population Core is
% being run at, so each of the R+1 Core calls gets its own cutoff. The default
% instead stops on enorm(dQ) < tol. Carried as NaN so that the pfqn_bs
% initialization below inherits the same test.
cntest = false;
if ischar(tol) || isstring(tol)
    if strcmpi(tol,'cn')
        cntest = true;
    else
        line_error(mfilename,sprintf('pfqn_egflinearizer: unknown tolerance specifier ''%s''.',char(tol)));
    end
elseif isnan(tol)
    cntest = true;
end
if cntest
    tol = NaN;
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
    if isempty(QN0)
        [~,q] = pfqn_bs(L,N_1,Z);
    else
        % warm-start the Bard-Schweitzer initialization from the supplied Q
        [~,q] = pfqn_bs(L,N_1,Z,tol,maxiter,QN0);
    end
    for r=1:R
        Q(:,r,1+s) = q(:,r);
    end
end

totiter = 0;
% Main loop
for I=1:npasses
    for s=0:R
        N_1 = oner(N,s); % for k=0 it just returns N
        % Core(N_1)
        [Q(:,:,1+s),~,~,iter] = Core(L,M,R,N_1,Z,Q(:,:,1+s),Delta,type,tol,maxiter-totiter,alpha,cntest);
        totiter = totiter + iter;
    end
    % Update_Delta
    for i=1:M
        for r=1:R
            if N(r)==1
                % At population N-e_r only class r itself vanishes; the other
                % classes keep the queue lengths Core just computed.
                Q(i,r,1+r) = 0;
            end
            for s=1:R
                Ns = oner(N,s);
                if Ns(r) > 0
                    Delta(i,r,s) = Q(i,r,1+s)/Ns(r)^alpha(r) - Q(i,r,1+0)/N(r)^alpha(r);
                else
                    % (N-e_s)_r = 0, i.e. r==s and N(r)==1: class r is absent at
                    % N-e_s, so F_ir(N-e_s) = 0 by the 0/0 convention of Chandy
                    % and Neuse (1982) eq (10). Their worked trace confirms it:
                    % D222 = 0 - L22(8,1)/1 = -1.
                    Delta(i,r,s) = -Q(i,r,1+0)/N(r)^alpha(r);
                end
            end
        end
    end
end


% Core(N)
[Q,W,X,iter] = Core(L,M,R,N,Z,Q(:,:,1+0),Delta,type,tol,maxiter-totiter,alpha,cntest);
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

function [Q,W,T,iter] = Core(L,M,R,N_1,Z,Q,Delta,type,tol,maxiter,alpha,cntest)
hasConverged = false;
W = L;
T = zeros(1,R);
iter = 0;
if cntest
    % Chandy and Neuse (1982), p.129 and appendix: the cutoff is a function of
    % the population Core is running at, so it is recomputed here rather than
    % once for the whole Linearizer.
    tol = pfqn_cntol(N_1);
    nz = N_1 > 0;
end
while ~hasConverged
    Qlast = Q;
    % Estimate population at
    Q_1 = Estimate(L,M,R,N_1,Z,Q,Delta,W,alpha);
    % Forward MVA
    [Q,W,T] = ForwardMVA(L,M,R,type,N_1,Z,Q_1);
    if cntest
        % max_{i,r} |dQ(i,r)| / N_r over the non-empty classes; an empty class
        % would divide by zero and it carries no jobs to converge.
        if isempty(find(nz,1))
            dev = 0;
        else
            dev = max(max(abs(Q(:,nz)-Qlast(:,nz))./repmat(N_1(nz),M,1)));
        end
    else
        dev = enorm(Q-Qlast);
    end
    if dev<tol || iter > maxiter
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
            % A class with no jobs at N_1 (or at N_1-e_s) has queue length 0
            % there. Guarding this is required, not cosmetic: without it
            % N_1(r)=0 divides by zero and oner returns a negative population
            % (e.g. oner([0 2],1) = [-1 2]), so Q_1 becomes NaN for the empty
            % class and poisons any Delta computed from this slice.
            if N_1(r) <= 0 || Ns(r) <= 0
                Q_1(i,r,1+s) = 0;
            else
                Q_1(i,r,1+s) = Ns(r)^alpha(r)*(Q(i,r,1+0)/N_1(r)^alpha(r) + Delta(i,r,s));
            end
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
% see _kb/03-api-layer.md (pfqn/ family: scaling, log-domain switches, dispatch gates)
for ist=1:M
    for r=1:R
        W(ist,r) = L(ist,r)*(1+sum(Q_1(ist,:,1+r)));
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


