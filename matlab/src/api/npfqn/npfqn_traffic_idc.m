function ctx = npfqn_traffic_idc(lambda0, P, c2a0, a0IdcFun, mu, cs2, sIdcFun, corrections)
% ctx = NPFQN_TRAFFIC_IDC(lambda0,P,c2a0,a0IdcFun,mu,cs2,sIdcFun)
%
% Traffic variability equations for the Robust Queueing Network Analyzer
% (RQNA) of W. Whitt and W. You (2018), "A Robust Queueing Network Analyzer
% Based on Indices of Dispersion". Assembles and solves:
%   - the limiting variability equations (eq. 42/44) for the asymptotic
%     total-arrival variability parameters c2_{a,i} = I_{a,i}(Inf);
%   - a solver (handle ctx.IaFun) for the time-dependent index of dispersion
%     for counts (IDC) equations (eq. 40/43), returning I_{a,i}(t) for all
%     internal arrival flows, using the default correction terms alpha_{i,j}
%     (eq. 34) and beta_i (eqs. 38-39) and tuning function h(rho)=rho^2.
%
% This models a single-class open queueing network of K single-server FCFS
% queues with Markovian routing P (P(i,j)=p_{i,j}).
%
% Inputs (all K-vectors are column vectors, K = number of queues):
%   lambda0   : external arrival rate into each queue
%   P         : KxK routing matrix among queues
%   c2a0      : asymptotic IDC (SCV) of each external arrival process
%   a0IdcFun  : handle, a0IdcFun(t) -> Kx1 external arrival IDC I_{a,0,i}(t)
%   mu        : service rate at each queue
%   cs2       : service SCV c2_{s,i}
%   sIdcFun   : handle, sIdcFun(t) -> Kx1 service IDC I_{s,i}(t)
%
% Output ctx: struct with fields
%   lambda,rho,Xi,c2a,c2d,c2aij,c2x and function handle IaFun(t) returning the
%   Kx1 vector of total arrival IDCs I_{a,i}(t).

K = numel(mu);
lambda0 = lambda0(:); mu = mu(:); cs2 = cs2(:); c2a0 = c2a0(:);
if nargin < 8 || isempty(corrections)
    corrections = struct();
end
if ~isfield(corrections,'alpha'), corrections.alpha = true; end
if ~isfield(corrections,'beta'),  corrections.beta  = true; end
useAlpha = corrections.alpha; useBeta = corrections.beta;

% ----- traffic rate equations (eq. 20-21) -----
Xi = inv(eye(K) - P');           % fundamental matrix (I-P')^{-1}
lambda = Xi * lambda0;           % total arrival rate at each queue
rho = lambda ./ mu;
lam_ji = (lambda(:) * ones(1,K)) .* P;   % lam_ji(j,i) = lambda_j p_{j,i}

% ----- correction terms (asymptotic, w*(Inf)=1) -----
% alpha: c2alpha_{i,j} = 2 Xi_{i,j} p_{i,j} (1-p_{i,j})
c2alpha = 2 * Xi .* P .* (1 - P);        % KxK
if ~useAlpha, c2alpha = zeros(K,K); end

% beta: zeta_{j,i;k,i} from eq. (39), then c2beta_i = (2/lambda_i) sum_{j<k} zeta
Sigma = cell(K,1);               % splitting covariance at each station l
for l = 1:K
    pl = P(l,:);                 % 1xK routing out of l
    Sl = -(pl' * pl) * lambda(l);
    Sl(1:K+1:end) = pl .* (1 - pl) * lambda(l);
    Sigma{l} = Sl;
end
Amat = diag(c2a0 .* lambda0);    % external arrival Brownian variance-rate
for l = 1:K
    Amat = Amat + Sigma{l};
end
% zetaAll{i} is KxK with entry (j,k) = zeta_{j,i;k,i}
zetaAll = cell(K,1);
c2beta = zeros(K,1);
for i = 1:K
    nu = (P(:,i) * ones(1,K)) .* Xi;   % nu(l,:) = p_{l,i} * Xi(l,:)   (row l)
    Z = nu * Amat * nu';               % first term for all (j,k)
    % add cross terms nu_k*Sigma_j*e_i + nu_j*Sigma_k*e_i
    for j = 1:K
        for k = 1:K
            Z(j,k) = Z(j,k) + nu(k,:)*Sigma{j}(:,i) + nu(j,:)*Sigma{k}(:,i);
        end
    end
    zetaAll{i} = Z;
    s = 0;
    for j = 1:K
        for k = j+1:K
            s = s + Z(j,k);
        end
    end
    if lambda(i) > 0
        c2beta(i) = (2/lambda(i)) * s;
    end
end
if ~useBeta, c2beta = zeros(K,1); zetaAll = repmat({zeros(K,K)},K,1); end

% ----- limiting variability equations (eq. 44): (E - Minf) c = binf -----
% variable ordering: [c2a(1..K), c2aij(i,j) row-major, c2d(1..K)]
Na = K; Naij = K*K; Nd = K;
N = Na + Naij + Nd;
ia  = @(i) i;
iaij= @(i,j) Na + (i-1)*K + j;
id  = @(i) Na + Naij + i;

Minf = zeros(N,N);
binf = zeros(N,1);
for i = 1:K
    % c2a_i = sum_j (lam_ji/lambda_i) c2aij_{j,i} + (lambda0_i/lambda_i) c2a0_i + c2beta_i
    if lambda(i) > 0
        for j = 1:K
            Minf(ia(i), iaij(j,i)) = lam_ji(j,i) / lambda(i);
        end
        binf(ia(i)) = (lambda0(i)/lambda(i)) * c2a0(i) + c2beta(i);
    end
    for j = 1:K
        % c2aij_{i,j} = p_{i,j} c2d_i + (1-p_{i,j}) + c2alpha_{i,j}
        Minf(iaij(i,j), id(i)) = P(i,j);
        binf(iaij(i,j)) = (1 - P(i,j)) + c2alpha(i,j);
    end
    % c2d_i = c2a_i
    Minf(id(i), ia(i)) = 1;
end
csol = (eye(N) - Minf) \ binf;
c2a  = csol(1:K);
c2aij= reshape(csol(Na+1:Na+Naij), K, K)';   % c2aij(i,j)
c2d  = csol(Na+Naij+1:end);
c2x  = c2a + cs2;                              % c2_{x,i} = I_{a,i}(Inf)+c2_{s,i}

% ----- assemble context and time-dependent solver -----
ctx = struct();
ctx.K = K; ctx.lambda = lambda; ctx.lambda0 = lambda0; ctx.mu = mu;
ctx.rho = rho; ctx.Xi = Xi; ctx.P = P; ctx.cs2 = cs2;
ctx.c2a = c2a; ctx.c2aij = c2aij; ctx.c2d = c2d; ctx.c2x = c2x;
ctx.c2a0 = c2a0; ctx.lam_ji = lam_ji;
ctx.zetaAll = zetaAll; ctx.c2alpha = c2alpha;
ctx.a0IdcFun = a0IdcFun; ctx.sIdcFun = sIdcFun;
ctx.IaFun = @(t) local_idc_at(ctx, t);
end

function Ia = local_idc_at(ctx, t)
% Solve the time-dependent IDC equations (43) at a single time t; return the
% Kx1 vector of total arrival IDCs I_{a,i}(t).
K = ctx.K; P = ctx.P; lambda = ctx.lambda; rho = ctx.rho;
c2x = ctx.c2x; lam_ji = ctx.lam_ji; Xi = ctx.Xi;
h = rho.^2;                              % tuning function h(rho)=rho^2

% departure weights w_i(t) = w*((1-rho_i)^2 lambda_i t /(h_i c2x_i))
warg = zeros(K,1);
for i = 1:K
    if h(i) > 0 && c2x(i) > 0
        warg(i) = (1-rho(i))^2 * lambda(i) * t / (h(i) * c2x(i));
    else
        warg(i) = Inf;
    end
end
w = npfqn_rqna_weight(warg);

Ia0 = ctx.a0IdcFun(t); Ia0 = Ia0(:);
Is  = ctx.sIdcFun(rho .* t); Is = Is(:);   % service IDC at scaled time rho*t

% time-dependent correction terms
% alpha_{i,j}(t) = c2alpha_{i,j} * w_i(t)
alpha_t = ctx.c2alpha .* (w * ones(1,K));  % row i scaled by w(i)
% beta_i(t) = (1/lambda_i) sum_{j~=k} zeta_{j,i;k,i} w*(arg_j), where the weight
% for source station j feeding i uses arg_j = (1-rho_j)^2 p_{j,i} lambda_j t/(h_j c2x_j)
beta_t = zeros(K,1);
for i = 1:K
    Z = ctx.zetaAll{i};
    % weight for source station j feeding i:
    wj = zeros(K,1);
    for j = 1:K
        if h(j) > 0 && c2x(j) > 0 && P(j,i) > 0
            aj = (1-rho(j))^2 * P(j,i) * lambda(j) * t / (h(j) * c2x(j));
            wj(j) = npfqn_rqna_weight(aj);
        end
    end
    s = 0;
    for j = 1:K
        for k = 1:K
            if j ~= k
                s = s + Z(j,k) * wj(j);
            end
        end
    end
    if lambda(i) > 0
        beta_t(i) = s / lambda(i);
    end
end

% assemble (E - M(t)) I = b(t)
Na = K; Naij = K*K; Nd = K; N = Na+Naij+Nd;
ia  = @(i) i; iaij= @(i,j) Na+(i-1)*K+j; idd = @(i) Na+Naij+i;
M = zeros(N,N); b = zeros(N,1);
for i = 1:K
    if lambda(i) > 0
        for j = 1:K
            M(ia(i), iaij(j,i)) = lam_ji(j,i)/lambda(i);
        end
        b(ia(i)) = (ctx.lambda0(i)/lambda(i))*Ia0(i) + beta_t(i);
    end
    for j = 1:K
        M(iaij(i,j), idd(i)) = P(i,j);
        b(iaij(i,j)) = (1 - P(i,j)) + alpha_t(i,j);
    end
    M(idd(i), ia(i)) = w(i);
    b(idd(i)) = (1 - w(i)) * Is(i);
end
sol = (eye(N) - M) \ b;
Ia = sol(1:K);
end
