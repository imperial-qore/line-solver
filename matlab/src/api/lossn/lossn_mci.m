%{ @file lossn_mci.m
 %  @brief Monte Carlo importance-sampling summation for loss networks
 %
 %  @author LINE Development Team
%}

%{
 % @brief Estimates the product-form normalization constant and class
 %        blocking probabilities of a loss network by Monte Carlo summation.
 %
 % @details
 % Implements the importance-sampling Monte Carlo summation method of
 % Ross and Wang, "Monte Carlo Summation Applied to Product-Form Loss
 % Networks", Probability in the Engineering and Informational Sciences,
 % 6 (1992), 323-348.
 %
 % A loss network has links j=1..J with capacity C(j) and classes (routes)
 % r=1..R with offered load nu(r) and per-link circuit requirement A(j,r).
 % The state n=(n_1,..,n_R) is feasible iff A*n <= C (set Omega). The
 % equilibrium distribution is product form with normalization constant
 %      g(C) = sum_{n in Omega} prod_r nu(r)^n_r / n_r! .
 % The class-r acceptance probability is 1-beta_r = g(C-A(:,r))/g(C).
 %
 % States are sampled from the importance distribution (Eq. 6)
 %      p(n) = (1/c) prod_r gamma_r^n_r / n_r!,  n in {0..N_1}x..x{0..N_R},
 % with N_r = min_j floor(C_j/A_jr) over A_jr>0. Ratio estimators (Eq. 8)
 % give unbiased/consistent g and consistent blocking with delta-method
 % confidence intervals.
 %
 % @par Syntax:
 % @code
 % [QLen, Loss, lG, ci, nsamples] = lossn_mci(nu, A, C)
 % [QLen, Loss, lG, ci, nsamples] = lossn_mci(nu, A, C, options)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>nu<td>Offered load of route (class) r (1xR vector)
 % <tr><td>A<td>Circuit requirement of link j for route r (JxR matrix)
 % <tr><td>C<td>Available capacity of link j (Jx1 vector)
 % <tr><td>options<td>struct with optional fields: samples (default 1e5),
 %                    gamma (1xR importance params, default Sec. 3.4
 %                    heuristic), seed (rng seed), alpha (CI level,
 %                    default 0.05)
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>QLen<td>Mean carried load (E[n_r]) for route r (1xR)
 % <tr><td>Loss<td>Blocking probability beta_r for route r (1xR)
 % <tr><td>lG<td>Log of the estimated normalization constant g(C)
 % <tr><td>ci<td>struct with confidence intervals and point estimates:
 %              .accept (Rx2), .loss (Rx2), .acceptPoint (1xR),
 %              .lossPoint (1xR), .level (1-alpha)
 % <tr><td>nsamples<td>Number of Monte Carlo samples used
 % </table>
%}
function [QLen, Loss, lG, ci, nsamples] = lossn_mci(nu, A, C, options)
if nargin < 4 || isempty(options)
    options = struct();
end
nu = nu(:)';         % 1xR
C  = C(:);           % Jx1
R = numel(nu);
J = numel(C);
if size(A,1) ~= J || size(A,2) ~= R
    line_error(mfilename, sprintf('A must be %dx%d (J x R).', J, R));
end

nsamples = getfielddef(options, 'samples', 1e5);
alpha    = getfielddef(options, 'alpha', 0.05);
seed     = getfielddef(options, 'seed', []);
gamma    = getfielddef(options, 'gamma', []);
if ~isempty(seed)
    rng(seed);
end
S = nsamples;

% Per-class maximum feasible occupancy N_r = min_j floor(C_j/A_jr)
N = zeros(1, R);
for k = 1:R
    pos = A(:,k) > 0;
    if any(pos)
        N(k) = floor(min(C(pos) ./ A(pos, k)));
    else
        N(k) = 0;  % class uses no link: cannot admit connections
    end
end

% Importance-sampling parameters gamma (Section 3.4 heuristic)
if isempty(gamma)
    load_j = zeros(J, 1);
    for j = 1:J
        load_j(j) = sum(A(j,:) .* nu) / C(j);
    end
    delta = max(load_j);
    b = max(A, [], 1);                 % 1xR, b_k = max_j A_jk
    base = 1 - 0.15 * (1 - delta);
    base = max(base, 1e-6);
    gamma = nu .* (base .^ b);
end
gamma = gamma(:)';
gamma = max(gamma, 1e-300);

% Normalization constant c of the importance distribution (log space)
log_c = 0;
for k = 1:R
    l = 0:N(k);
    logterms = l * log(gamma(k)) - gammaln(l + 1);
    log_c = log_c + logsumexp(logterms);
end

% Draw S i.i.d. samples V (S x R), each column truncated Poisson(gamma_k)
V = zeros(S, R);
for k = 1:R
    l = 0:N(k);
    logpmf = l * log(gamma(k)) - gammaln(l + 1);
    pmf = exp(logpmf - logsumexp(logpmf));
    cdf = cumsum(pmf);
    cdf(end) = 1;                       % guard rounding
    u = rand(S, 1);
    idx = sum(u > cdf(:)', 2);          % # of thresholds exceeded = value
    V(:, k) = idx;
end

% Feasibility indicators
AV = V * A';                            % S x J, row s = A*n_s
inOmega = all(AV <= C(:)', 2);          % S x 1
inOmegaK = false(S, R);
for k = 1:R
    Ck = (C - A(:,k))';                 % 1 x J capacity for Omega(C-A_.k)
    inOmegaK(:, k) = all(AV <= Ck, 2);
end

% Likelihood ratio alpha^i = prod_k (nu_k/gamma_k)^V_k  (log space)
logratio = log(nu) - log(gamma);        % 1xR
log_alpha = V * logratio';              % S x 1

% Normalization constant estimate g(C) = (c/S) sum alpha 1(Omega)
laO = log_alpha(inOmega);
if isempty(laO)
    lG = -Inf;
else
    m = max(laO);
    lG = log_c + m + log(sum(exp(laO - m))) - log(S);
end

% Ratio estimators for acceptance (Eq. 8) with shifted weights for stability
if any(inOmega)
    M = max(log_alpha(inOmega));
else
    M = 0;
end
w = exp(log_alpha - M);                 % S x 1 scaled likelihood ratios
Z = w .* inOmega;                       % denominator summand
meanZ = mean(Z);

crit = sqrt(2) * erfinv(1 - alpha);     % standard-normal 1-alpha/2 quantile
accept = zeros(1, R);
acceptCI = zeros(R, 2);
for k = 1:R
    Y = w .* inOmegaK(:, k);
    meanY = mean(Y);
    if meanZ <= 0
        phi = NaN; half = NaN;
    else
        phi = meanY / meanZ;
        varY = var(Y);
        varZ = var(Z);
        covYZ = sum((Y - meanY) .* (Z - meanZ)) / (S - 1);
        sig2 = (varY - 2*phi*covYZ + phi^2*varZ) / (S * meanZ^2);
        sig2 = max(sig2, 0);
        half = crit * sqrt(sig2);
    end
    accept(k) = phi;
    acceptCI(k, :) = [phi - half, phi + half];
end

Loss = 1 - accept;
QLen = nu .* accept;

ci = struct();
ci.accept = acceptCI;
ci.loss = [1 - acceptCI(:,2), 1 - acceptCI(:,1)];
ci.acceptPoint = accept;
ci.lossPoint = Loss;
ci.level = 1 - alpha;
end

function s = logsumexp(x)
x = x(:)';
m = max(x);
if isinf(m)
    s = m;
else
    s = m + log(sum(exp(x - m)));
end
end

function v = getfielddef(s, f, d)
if isfield(s, f) && ~isempty(s.(f))
    v = s.(f);
else
    v = d;
end
end
