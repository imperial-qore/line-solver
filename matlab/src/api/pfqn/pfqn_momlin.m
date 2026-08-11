%{
%{
 % @file pfqn_momlin.m
 % @brief Moment linearizer: approximate mean queue lengths and their second
 %        moments (variance / covariance) for large closed product-form
 %        queueing networks.
%}
%}

function [Q,X,U,R,QVar,QCov,dQ] = pfqn_momlin(L,N,Z,tol,maxiter)
%{
%{
 % @brief Approximate first and second queue-length moments of a closed
 %        product-form network, scalable to large populations and many classes
 %        where exact MVA (exponential in the number of classes) and CoMoM
 %        (single-station) are infeasible.
 %
 %        Means are obtained from the Schweitzer-Bard AMVA fixed point. Second
 %        moments use the exact product-form identity
 %             Cov[n_{i,r},n_{j,s}] = D_{j,s} dQ_{i,r}/dD_{j,s},
 %        with the demand derivatives obtained by analytically linearizing the
 %        AMVA fixed point (a "moment linearizer" in the sense of Strelen and
 %        Akyildiz, here developed per class). Both moments carry the AMVA
 %        approximation error and become exact only in the limits where
 %        Schweitzer-Bard is exact; for exact moments on tractable models use
 %        pfqn_sens (differentiated MVA / CoMoM).
 %
 % @fn pfqn_momlin(L, N, Z, tol, maxiter)
 % @param L  Service demand matrix (M x R).
 % @param N  Closed population vector (1 x R).
 % @param Z  Think time vector (1 x R). Default: zeros.
 % @param tol Convergence tolerance on the queue-length fixed point. Default 1e-8.
 % @param maxiter Maximum iterations. Default 1000.
 % @return Q    Mean queue length (M x R).
 % @return X    Throughput per class (1 x R).
 % @return U    Utilization (M x R).
 % @return R    Residence time (M x R).
 % @return QVar Queue-length variance (M x R).
 % @return QCov Queue-length covariance tensor (M x R x M x R).
 % @return dQ   Demand-derivative tensor, dQ(i,r,j,s) = dQ_{i,r}/dD_{j,s} (M x R x M x R).
%}
%}
[M,R] = size(L);
N = ceil(N(:)');
if nargin < 3 || isempty(Z), Z = zeros(1,R); end
Z = Z(:)';
if nargin < 4 || isempty(tol), tol = 1e-8; end
if nargin < 5 || isempty(maxiter), maxiter = 1000; end
if any(isinf(N))
    line_error(mfilename,'pfqn_momlin supports closed classes only.');
end

% Schweitzer population-scaling coefficients c_s^{(r)} = (N_s-delta_{rs})/N_s
% approximating Q_{i,s}(N-1_r) ~ c_s^{(r)} Q_{i,s}(N).
c = ones(R,R);   % c(r,s)
for r = 1:R
    for s = 1:R
        if N(s) > 0
            c(r,s) = (N(s) - (r==s)) / N(s);
        else
            c(r,s) = 0;
        end
    end
end

% ---- Schweitzer-Bard AMVA fixed point for the means -------------------------
Q = zeros(M,R);
for r = 1:R
    if N(r) > 0, Q(:,r) = N(r)/M; end   % uniform initial guess
end
X = zeros(1,R); Rmat = zeros(M,R);
for it = 1:maxiter
    Qold = Q;
    for r = 1:R
        if N(r) == 0, X(r) = 0; Rmat(:,r) = 0; continue; end
        for i = 1:M
            Rmat(i,r) = L(i,r) * (1 + c(r,:) * Q(i,:)');
        end
        X(r) = N(r) / (Z(r) + sum(Rmat(:,r)));
        Q(:,r) = X(r) * Rmat(:,r);
    end
    if max(abs(Q(:)-Qold(:))) < tol, break; end
end
U = zeros(M,R);
for r = 1:R
    U(:,r) = X(r) * L(:,r);
end

% ---- analytic linearization of the fixed point ------------------------------
% see _kb/03-api-layer.md (pfqn/ family: scaling, log-domain switches, dispatch gates)
dQ = zeros(M,R,M,R);
denom = zeros(1,R);
for r = 1:R
    denom(r) = Z(r) + sum(Rmat(:,r));
end
for j = 1:M
    for s0 = 1:R
        if N(s0) == 0, continue; end   % empty class: derivative 0
        dq = zeros(M,R);
        for it = 1:maxiter
            dqold = dq;
            for r = 1:R
                if N(r) == 0, continue; end
                dR = zeros(M,1);
                for i = 1:M
                    dDir = (i==j) && (r==s0);
                    dR(i) = dDir * (1 + c(r,:)*Q(i,:)') + L(i,r) * (c(r,:) * dq(i,:)');
                end
                dXr = -(X(r)^2/N(r)) * sum(dR);   % dZ/dtheta = 0
                dq(:,r) = dXr * Rmat(:,r) + X(r) * dR;
            end
            if max(abs(dq(:)-dqold(:))) < tol, break; end
        end
        dQ(:,:,j,s0) = dq;
    end
end

% ---- second moments via the product-form covariance identity ----------------
QCov = zeros(M,R,M,R);
for i = 1:M
    for r = 1:R
        for j = 1:M
            for s = 1:R
                QCov(i,r,j,s) = L(j,s) * dQ(i,r,j,s);
            end
        end
    end
end
QVar = zeros(M,R);
for i = 1:M
    for r = 1:R
        QVar(i,r) = QCov(i,r,i,r);
    end
end

R = Rmat;   % assign residence-time output last (R held the class count above)
end
