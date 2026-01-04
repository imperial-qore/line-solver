function [L, W, Ca, Cd, lambda, rho] = me_oqn(M, R, lambda0, Ca0, mu, Cs, P, options)
%ME_OQN Maximum Entropy algorithm for Open Queueing Networks
%
% Implements the ME algorithm from Kouvatsos (1994) "Entropy Maximisation
% and Queueing Network Models", Section 3.2.
%
% INPUTS:
%   M       - Number of queues (stations)
%   R       - Number of job classes
%   lambda0 - External arrival rates [M x R matrix], lambda0(i,r) = lambda_oi,r
%   Ca0     - External arrival scv [M x R matrix], Ca0(i,r) = Caoi,r
%   mu      - Service rates [M x R matrix], mu(i,r) = μi,r
%   Cs      - Service scv [M x R matrix], Cs(i,r) = Csi,r
%   P       - Routing probability matrix [M x M x R], P(j,i,r) = pji,r
%             (probability class r goes from queue j to queue i)
%   options - (optional) struct with fields:
%             .tol     - convergence tolerance (default: 1e-6)
%             .maxiter - maximum iterations (default: 1000)
%             .verbose - print iteration info (default: false)
%
% OUTPUTS:
%   L      - Mean queue lengths [M x R matrix]
%   W      - Mean waiting times [M x R matrix]
%   Ca     - Arrival scv at each queue [M x R matrix]
%   Cd     - Departure scv at each queue [M x R matrix]
%   lambda - Total arrival rates [M x R matrix]
%   rho    - Utilizations [M x R matrix]
%
% Reference: Kouvatsos (1994), Equations 3.6 and 3.7

% Handle optional arguments
if nargin < 8
    options = struct();
end
if ~isfield(options, 'tol')
    options.tol = 1e-6;
end
if ~isfield(options, 'maxiter')
    options.maxiter = 1000;
end
if ~isfield(options, 'verbose')
    options.verbose = false;
end

% Step 1: Feedback correction
% If pii,r > 0, apply feedback elimination transformation
P_eff = P;
mu_eff = mu;
Cs_eff = Cs;
for i = 1:M
    for r = 1:R
        pii = P(i, i, r);
        if pii > 0
            % Feedback correction: adjust service rate and scv
            mu_eff(i, r) = mu(i, r) * (1 - pii);
            % Adjusted service scv accounting for feedback
            % Cs_eff = (Cs + pii*(1-pii)) / (1-pii)^2 simplified form
            Cs_eff(i, r) = Cs(i, r) / (1 - pii) + pii / (1 - pii);
            % Remove self-loop from routing
            P_eff(i, i, r) = 0;
            % Renormalize routing probabilities
            row_sum = sum(P_eff(i, :, r));
            if row_sum > 0
                P_eff(i, :, r) = P_eff(i, :, r) / row_sum * (1 - pii);
            end
        end
    end
end

% Step 2: Initialize arrival scv
Ca = ones(M, R);

% Step 3: Solve job flow balance equations using ORIGINAL P (including feedback)
% lambda_i,r = lambda_oi,r + Σj lambda_j,r * pji,r
% This is: lambda = lambda0 + P' * lambda (for each class)
lambda = zeros(M, R);
for r = 1:R
    % Build the system (I - P') * lambda = lambda0
    % where P' is the transpose of the routing matrix for class r
    % Use ORIGINAL P to get correct effective arrival rates with feedback
    Pr = squeeze(P(:, :, r))';  % P(j,i,r) -> need sum over j
    A = eye(M) - Pr;
    lambda(:, r) = A \ lambda0(:, r);
end

% Compute utilizations using ORIGINAL mu (not mu_eff)
% mu_eff is only for variability calculations after feedback elimination
rho = zeros(M, R);
for i = 1:M
    for r = 1:R
        if mu(i, r) > 0
            rho(i, r) = lambda(i, r) / mu(i, r);
        end
    end
end

% Check stability
rho_total = sum(rho, 2);
if any(rho_total >= 1)
    warning('me_oqn:unstable', 'Network is unstable (utilization >= 1 at some queues)');
end

% Initialize outputs
L = zeros(M, R);
Cd = ones(M, R);

% Step 4-5: Iterative computation of Ca, Cd, and L
for iter = 1:options.maxiter
    Ca_old = Ca;

    % Step 4: Apply GE-type formulae for mean queue length
    % For GE/GE/1 queue under ME principle:
    % L = ρhat + ρhat^2(Ca + Cs) / (2(1 - ρhat))
    for i = 1:M
        rho_i = sum(rho(i, :));  % Total utilization at queue i
        if rho_i > 0 && rho_i < 1
            for r = 1:R
                if rho(i, r) > 0
                    % Proportion of class r traffic
                    prop_r = rho(i, r) / rho_i;
                    % GE/GE/1 mean queue length formula (ME approximation)
                    % Aggregate Ca and Cs weighted by traffic
                    Ca_agg = Ca(i, r);
                    Cs_agg = Cs_eff(i, r);
                    L_total = rho_i + (rho_i^2 * (Ca_agg + Cs_agg)) / (2 * (1 - rho_i));
                    L(i, r) = prop_r * L_total;
                end
            end
        end
    end

    % Step 5a: Compute departure scv using equation (3.6)
    % Cdj,r = 2*Lj,r*(1 - ρhatj) + Caj,r*(1 - 2*ρhatj)
    for j = 1:M
        rho_j = sum(rho(j, :));
        for r = 1:R
            if lambda(j, r) > 0
                Cd(j, r) = 2 * L(j, r) * (1 - rho_j) + Ca(j, r) * (1 - 2 * rho_j);
                Cd(j, r) = max(0, Cd(j, r));  % Ensure non-negative
            end
        end
    end

    % Step 5b: Compute arrival scv using equation (3.7)
    % Cai,r = -1 + {Σj (lambda_j,r*pji,r/lambda_i,r)*[Cdji,r + 1]^-1 + (lambda_oi,r/lambda_i,r)*[Caoi,r + 1]^-1}^-1
    % where Cdji,r = 1 + pji,r*(Cdj,r - 1) (thinning)
    for i = 1:M
        for r = 1:R
            if lambda(i, r) > 0
                sum_inv = 0;

                % Contribution from other queues
                for j = 1:M
                    pji = P_eff(j, i, r);
                    if pji > 0 && lambda(j, r) > 0
                        % Thinning formula for departure scv after splitting
                        Cdji = 1 + pji * (Cd(j, r) - 1);
                        weight = (lambda(j, r) * pji) / lambda(i, r);
                        sum_inv = sum_inv + weight / (Cdji + 1);
                    end
                end

                % Contribution from external arrivals
                if lambda0(i, r) > 0
                    weight0 = lambda0(i, r) / lambda(i, r);
                    sum_inv = sum_inv + weight0 / (Ca0(i, r) + 1);
                end

                % Compute new arrival scv
                if sum_inv > 0
                    Ca(i, r) = -1 + 1 / sum_inv;
                    Ca(i, r) = max(0, Ca(i, r));  % Ensure non-negative
                end
            end
        end
    end

    % Check convergence
    delta = max(abs(Ca(:) - Ca_old(:)));
    if options.verbose
        fprintf('Iteration %d: max delta = %e\n', iter, delta);
    end
    if delta < options.tol
        if options.verbose
            fprintf('Converged after %d iterations\n', iter);
        end
        break;
    end
end

if iter == options.maxiter && delta >= options.tol
    warning('me_oqn:noconverge', 'Did not converge within %d iterations (delta=%e)', options.maxiter, delta);
end

% Step 6: Compute final statistics
% Mean waiting time via Little's law: W = L / lambda_
W = zeros(M, R);
for i = 1:M
    for r = 1:R
        if lambda(i, r) > 0
            W(i, r) = L(i, r) / lambda(i, r);
        end
    end
end

end
