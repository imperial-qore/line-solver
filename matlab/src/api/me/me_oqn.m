function [L, W, Ca, Cd, lambda, rho, iter] = me_oqn(M, R, lambda0, Ca0, mu, Cs, P, c, insens, options)
%ME_OQN Maximum Entropy algorithm for Open Queueing Networks
%
% Implements the ME algorithm from Kouvatsos (1994) "Entropy Maximisation
% and Queueing Network Models", Section 3.2, with the GE/GE/c building
% block of Section 3.4 (eq. 3.9) and the GE/GE/inf building block.
%
% INPUTS:
%   M       - Number of queues (stations)
%   R       - Number of job classes
%   lambda0 - External arrival rates [M x R matrix], lambda0(i,r) = lambda_oi,r
%   Ca0     - External arrival scv [M x R matrix], Ca0(i,r) = Caoi,r
%   mu      - Service rates [M x R matrix], mu(i,r) = mu_i,r
%   Cs      - Service scv [M x R matrix], Cs(i,r) = Cs_i,r
%   P       - Routing probability matrix [M x M x R], P(j,i,r) = p_ji,r
%             (probability class r goes from queue j to queue i)
%   c       - (optional) Servers per queue [M x 1 vector]; Inf marks an
%             infinite-server (IS) queue (default: ones(M,1))
%   insens  - (optional) Logical vector [M x 1]; true marks a station with
%             an insensitive scheduling discipline (PS, LCFS-PR), solved
%             with the product-form mean queue length L_r = rho_r/(1-rho)
%             instead of the FCFS GE formula (default: false(M,1))
%   options - (optional) struct with fields:
%             .tol     - convergence tolerance (default: 1e-6)
%             .maxiter - maximum iterations (default: 1000)
%             .verbose - print iteration info (default: false)
%
% OUTPUTS:
%   L      - Mean queue lengths [M x R matrix]
%   W      - Mean response times [M x R matrix], W = L ./ lambda
%   Ca     - Arrival scv at each queue [M x R matrix]
%   Cd     - Departure scv at each queue [M x R matrix]
%   lambda - Total arrival rates [M x R matrix], inclusive of self-loop
%            revisits (visit-based throughput)
%   rho    - Utilizations [M x R matrix]: per-server utilization for finite
%            c, mean number of busy servers for IS queues
%   iter   - Number of iterations until convergence
%
% Reference:
%   D.D. Kouvatsos, "Entropy Maximisation and Queueing Network Models",
%   Annals of Operations Research, 48:63-126, 1994. Equations (3.3), (3.6),
%   (3.7), (3.9) and the multiclass GE/GE/1/FCFS mql of Section 3.1.1.

% Handle optional arguments
if nargin < 8 || isempty(c)
    c = ones(M, 1);
end
if nargin < 9 || isempty(insens)
    insens = false(M, 1);
end
if nargin < 10
    options = struct();
end
insens = logical(insens(:));
if ~isfield(options, 'tol')
    options.tol = 1e-6;
end
if ~isfield(options, 'maxiter')
    options.maxiter = 1000;
end
if ~isfield(options, 'verbose')
    options.verbose = false;
end
c = c(:);

% Step 1: Feedback correction. A job that returns immediately to queue i
% with probability pii receives a geometric number of service passes, so
% the composite service time has rate mu*(1-pii) and scv pii+(1-pii)*Cs.
% The self-loop is removed and the residual routing renormalized.
P_eff = P;
mu_eff = mu;
Cs_eff = Cs;
for i = 1:M
    for r = 1:R
        pii = P(i, i, r);
        if pii > 0
            mu_eff(i, r) = mu(i, r) * (1 - pii);
            Cs_eff(i, r) = pii + (1 - pii) * Cs(i, r);
            P_eff(i, :, r) = P(i, :, r) / (1 - pii);
            P_eff(i, i, r) = 0;
        end
    end
end

% Step 3: Solve the job flow balance equations lambda = lambda0 + P'*lambda
% on the original routing; lambda counts self-loop revisits (visit-based
% throughput), while lambda_eff excludes them and is the arrival rate seen
% by the feedback-corrected queue.
lambda = zeros(M, R);
lambda_eff = zeros(M, R);
for r = 1:R
    Pr = P(:, :, r);
    A = eye(M) - Pr';
    lambda(:, r) = A \ lambda0(:, r);
    lambda_eff(:, r) = lambda(:, r) .* (1 - diag(Pr));
end

% Utilizations: per-server utilization for finite-server queues (invariant
% to the feedback correction since lambda_eff/mu_eff = lambda/mu); mean
% number of busy servers for IS queues.
rho = zeros(M, R);
for i = 1:M
    for r = 1:R
        if mu(i, r) > 0
            if isinf(c(i))
                rho(i, r) = lambda_eff(i, r) / mu_eff(i, r);
            else
                rho(i, r) = lambda(i, r) / (c(i) * mu(i, r));
            end
        end
    end
end

% Stability check (finite-server queues only)
unstable = false(M, 1);
for i = 1:M
    if ~isinf(c(i)) && sum(rho(i, :)) >= 1
        unstable(i) = true;
    end
end
if any(unstable)
    warning('me_oqn:unstable', 'Network is unstable (utilization >= 1 at some queues)');
end

% Step 2: Initialize arrival scvs
Ca = ones(M, R);
Cd = ones(M, R);
L = zeros(M, R);

% Steps 4-5: fixed-point iteration on the arrival scvs
delta = Inf;
iter = 0;
for iter = 1:options.maxiter
    Ca_old = Ca;

    % Step 4: GE-type mean queue length formulae
    for i = 1:M
        rho_i = sum(rho(i, :));
        if isinf(c(i))
            % GE/GE/inf queue: L = lambda/mu
            for r = 1:R
                if lambda_eff(i, r) > 0 && mu_eff(i, r) > 0
                    L(i, r) = lambda_eff(i, r) / mu_eff(i, r);
                end
            end
        elseif unstable(i)
            for r = 1:R
                if lambda_eff(i, r) > 0
                    L(i, r) = Inf;
                end
            end
        elseif c(i) == 1
            if insens(i)
                % Insensitive disciplines (PS, LCFS-PR): product-form mql,
                % exact irrespective of the service distribution
                for r = 1:R
                    if lambda_eff(i, r) > 0 && mu_eff(i, r) > 0
                        L(i, r) = rho(i, r) / (1 - rho_i);
                    end
                end
            else
                % Multiclass GE/GE/1/FCFS mql (Section 3.1.1):
                % L_r = rho_r*(Ca_r+1)/2
                %       + lambda_r*sum_u lambda_u*(Cs_u+Ca_u)/mu_u^2/(2*(1-rho))
                resid = 0;
                for u = 1:R
                    if lambda_eff(i, u) > 0 && mu_eff(i, u) > 0
                        resid = resid + lambda_eff(i, u) * (Cs_eff(i, u) + Ca(i, u)) / mu_eff(i, u)^2;
                    end
                end
                for r = 1:R
                    if lambda_eff(i, r) > 0 && mu_eff(i, r) > 0
                        L(i, r) = rho(i, r) * (Ca(i, r) + 1) / 2 + lambda_eff(i, r) * resid / (2 * (1 - rho_i));
                    end
                end
            end
        else
            % GE/GE/c/FCFS (eq. 3.9) on the class-aggregated stream; the
            % per-class disaggregation assigns each class its mean number
            % in service plus a share of the common FCFS waiting line
            % proportional to its arrival rate.
            lam_a = 0;
            for u = 1:R
                if lambda_eff(i, u) > 0 && mu_eff(i, u) > 0
                    lam_a = lam_a + lambda_eff(i, u);
                end
            end
            if lam_a > 0
                inv_a = 0;
                ES = 0;
                ES2 = 0;
                for u = 1:R
                    if lambda_eff(i, u) > 0 && mu_eff(i, u) > 0
                        wu = lambda_eff(i, u) / lam_a;
                        inv_a = inv_a + wu / (Ca(i, u) + 1);
                        ES = ES + wu / mu_eff(i, u);
                        ES2 = ES2 + wu * (Cs_eff(i, u) + 1) / mu_eff(i, u)^2;
                    end
                end
                Ca_a = -1 + 1 / inv_a;
                Cs_a = ES2 / ES^2 - 1;
                L_a = ge_gec_mql(lam_a, Ca_a, 1 / ES, Cs_a, c(i));
                Lq_a = L_a - lam_a * ES; % mean waiting-line length
                for r = 1:R
                    if lambda_eff(i, r) > 0 && mu_eff(i, r) > 0
                        L(i, r) = c(i) * rho(i, r) + (lambda_eff(i, r) / lam_a) * Lq_a;
                    end
                end
            end
        end
    end

    % Step 5a: departure scvs
    for j = 1:M
        rho_j = sum(rho(j, :));
        for r = 1:R
            if lambda_eff(j, r) > 0
                if isinf(c(j))
                    % GE/GE/inf queue: interdeparture scv = interarrival scv
                    Cd(j, r) = Ca(j, r);
                elseif unstable(j)
                    % Saturated server: departures follow the service process
                    Cd(j, r) = Cs_eff(j, r);
                elseif c(j) == 1
                    % Eq. (3.6) on the class-r virtual queue, with the
                    % marginal utilization rhohat_r of eq. (3.3)
                    rhohat = rho(j, r) * L(j, r) / (L(j, r) + rho_j - rho(j, r));
                    Cd(j, r) = 2 * L(j, r) * (1 - rhohat) + Ca(j, r) * (1 - 2 * rhohat);
                else
                    % GE/GE/c interdeparture scv (Section 4.2)
                    Cd(j, r) = rho_j * (1 - rho_j) + (1 - rho_j) * Ca(j, r) + rho_j^2 * Cs_eff(j, r);
                end
            end
        end
    end

    % Step 5b: arrival scvs by GE-type merging, eq. (3.7), with the
    % splitting (thinning) formula Cdji = 1 + pji*(Cd - 1) applied to the
    % feedback-corrected flows
    for i = 1:M
        for r = 1:R
            if lambda_eff(i, r) > 0
                sum_inv = 0;
                for j = 1:M
                    pji = P_eff(j, i, r);
                    if pji > 0 && lambda_eff(j, r) > 0
                        Cdji = 1 + pji * (Cd(j, r) - 1);
                        sum_inv = sum_inv + (lambda_eff(j, r) * pji / lambda_eff(i, r)) / (Cdji + 1);
                    end
                end
                if lambda0(i, r) > 0
                    sum_inv = sum_inv + (lambda0(i, r) / lambda_eff(i, r)) / (Ca0(i, r) + 1);
                end
                if sum_inv > 0
                    Ca(i, r) = -1 + 1 / sum_inv;
                end
            end
        end
    end

    % Check convergence
    delta = max(abs(Ca(:) - Ca_old(:)));
    if options.verbose
        fprintf('Iteration %d: max delta = %e\n', int32(iter), delta);
    end
    if delta < options.tol
        if options.verbose
            fprintf('Converged after %d iterations\n', int32(iter));
        end
        break;
    end
end

if iter == options.maxiter && delta >= options.tol
    warning('me_oqn:noconverge', 'Did not converge within %d iterations (delta=%e)', int32(options.maxiter), delta);
end

% Step 6: response times by Little's law on the reported arrival rates
W = zeros(M, R);
for i = 1:M
    for r = 1:R
        if lambda(i, r) > 0
            W(i, r) = L(i, r) / lambda(i, r);
        end
    end
end

end

function L = ge_gec_mql(lambda, Ca, mu, Cs, c)
% Mean queue length of a stable GE/GE/c/FCFS queue via the exact ME
% solution of Kouvatsos (1994), eq. (3.9).
alpha2 = 2 / (Cs + 1);
alpha1 = 1 - alpha2;
beta2 = 2 / (Ca + 1);
beta1 = 1 - beta2;
lambda2 = beta2 * lambda;
mu2 = alpha2 * mu;
g = zeros(c, 1);
for j = 1:(c - 1)
    g(j) = (lambda2 + (j - 1) * mu2 * beta1) * alpha2 / (j * mu2 * (1 - alpha1 * beta1));
end
g(c) = (lambda2 + (c - 1) * mu2 * beta1) * alpha2 / (lambda2 * alpha1 + c * mu2);
x = (lambda2 + c * mu2 * beta1) / (lambda2 * alpha1 + c * mu2);
Gn = cumprod(g);
Z = 1 + sum(Gn(1:c-1)) + Gn(c) / (1 - x);
S1 = 0;
for n = 1:(c - 1)
    S1 = S1 + n * Gn(n);
end
S2 = Gn(c) * (c / (1 - x) + x / (1 - x)^2);
L = (S1 + S2) / Z;
end
