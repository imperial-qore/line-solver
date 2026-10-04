function [L, W, Ca, Cd, lambda, rho, X, iter] = me_mqn(M, R, openClasses, lambda0, Ca0, N, mu, Cs, P, c, refstat, insens, options)
%ME_MQN Maximum Entropy algorithm for Mixed Queueing Networks
%
% Extension of the Kouvatsos (1994) Maximum Entropy Method to mixed
% open/closed multiclass networks. The 1994 survey notes that the closed
% two-stage treatment carries over to mixed networks (Section 2.3) but
% gives no algorithm; this implementation composes the open (Section 3.2)
% and closed (Section 3.3) ME algorithms by product-form-style
% conditioning:
%   1. the open classes are solved by the open GE-type fixed point on
%      the station set, ignoring the closed classes;
%   2. the closed classes are solved by the two-stage pseudo-open plus
%      convolution algorithm on servers whose capacity is reduced by the
%      open-class utilization, mu_c(i,r) = mu(i,r)*(1-rho_o(i));
%   3. the open mean queue lengths are inflated by the closed occupancy,
%      L_o(i,r) <- L_o(i,r)*(1 + Lc(i)) at single-server stations.
% Steps 2-3 are exact in the BCMP product-form limit (exponential
% services, where they reduce to the classical mixed MVA treatment) and
% GE-type approximations otherwise. Only single-server and infinite-
% server stations are supported.
%
% INPUTS:
%   M           - Number of queues (stations)
%   R           - Number of job classes
%   openClasses - Logical vector [1 x R], true for open classes
%   lambda0     - External arrival rates [M x R], zero columns for closed
%   Ca0         - External arrival scv [M x R]
%   N           - Class populations [1 x R], Inf for open classes
%   mu          - Service rates [M x R matrix]
%   Cs          - Service scv [M x R matrix]
%   P           - Routing probability matrix [M x M x R], P(j,i,r) = p_ji,r
%   c           - (optional) Servers per queue [M x 1]; Inf marks an IS
%                 queue; finite values must be 1 (default: ones(M,1))
%   refstat     - (optional) Reference station per closed class [1 x R]
%   options     - (optional) struct with tol / maxiter / verbose fields
%
% OUTPUTS:
%   L      - Mean queue lengths [M x R matrix]
%   W      - Mean response times [M x R matrix], W = L ./ lambda
%   Ca     - Arrival scv at each queue [M x R matrix]
%   Cd     - Departure scv at each queue [M x R matrix]
%   lambda - Throughputs at each queue [M x R matrix]
%   rho    - Utilizations [M x R matrix]
%   X      - Class throughputs [1 x R] (arrival rates for open classes,
%            reference-station throughputs for closed classes)
%   iter   - Total number of fixed-point iterations
%
% Reference:
%   D.D. Kouvatsos, "Entropy Maximisation and Queueing Network Models",
%   Annals of Operations Research, 48:63-126, 1994. Sections 3.2-3.3.

if nargin < 10 || isempty(c)
    c = ones(M, 1);
end
if nargin < 11 || isempty(refstat)
    refstat = zeros(1, R);
end
if nargin < 12 || isempty(insens)
    insens = false(M, 1);
end
if nargin < 13
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
openClasses = logical(openClasses(:)');

oc = find(openClasses);
cc = find(~openClasses);
Ro = length(oc);
Rc = length(cc);

L = zeros(M, R);
W = zeros(M, R);
Ca = ones(M, R);
Cd = ones(M, R);
lambda = zeros(M, R);
rho = zeros(M, R);
X = zeros(1, R);
iter = 0;

% Step 1: open classes by the Section 3.2 GE-type fixed point
rho_o = zeros(M, 1); % per-station aggregate open utilization
if Ro > 0
    [Lo, ~, Cao, Cdo, lamo, rhoo, itero] = me_oqn(M, Ro, lambda0(:, oc), Ca0(:, oc), ...
        mu(:, oc), Cs(:, oc), P(:, :, oc), c, insens, options);
    iter = iter + itero;
    L(:, oc) = Lo;
    Ca(:, oc) = Cao;
    Cd(:, oc) = Cdo;
    lambda(:, oc) = lamo;
    rho(:, oc) = rhoo;
    for i = 1:M
        if ~isinf(c(i))
            rho_o(i) = sum(rhoo(i, :));
        end
    end
    for k = 1:Ro
        X(oc(k)) = sum(lambda0(:, oc(k)));
    end
end

% Step 2: closed classes by the Section 3.3 algorithm on servers with
% capacity reduced by the open-class utilization
if Rc > 0
    mu_c = mu(:, cc);
    for i = 1:M
        if ~isinf(c(i))
            mu_c(i, :) = mu_c(i, :) * max(1 - rho_o(i), 0);
        end
    end
    [Lc, Wc, Cac, Cdc, lamc, rhoc, Xc, iterc] = me_cqn(M, Rc, N(cc), mu_c, Cs(:, cc), ...
        P(:, :, cc), c, refstat(cc), insens, options);
    iter = iter + iterc;
    L(:, cc) = Lc;
    W(:, cc) = Wc;
    Ca(:, cc) = Cac;
    Cd(:, cc) = Cdc;
    lambda(:, cc) = lamc;
    X(cc) = Xc;
    % Closed utilizations are relative to the reduced capacity; rescale to
    % the busy fraction of the physical server
    for i = 1:M
        for k = 1:Rc
            if isinf(c(i))
                rho(i, cc(k)) = rhoc(i, k);
            else
                rho(i, cc(k)) = rhoc(i, k) * max(1 - rho_o(i), 0);
            end
        end
    end
end

% Step 3: inflate the open queue lengths by the closed occupancy at
% single-server stations (exact in the product-form limit) and derive
% response times by Little's law
if Ro > 0 && Rc > 0
    for i = 1:M
        if ~isinf(c(i))
            Lc_i = sum(L(i, cc));
            L(i, oc) = L(i, oc) * (1 + Lc_i);
        end
    end
end
for i = 1:M
    for k = 1:Ro
        if lambda(i, oc(k)) > 0
            W(i, oc(k)) = L(i, oc(k)) / lambda(i, oc(k));
        end
    end
end

end
