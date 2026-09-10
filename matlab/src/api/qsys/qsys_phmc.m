function result = qsys_phmc(alpha, T, mu, c, varargin)
% QSYS_PHMC Exact PH/M/c (Neuts' matrix-geometric).
%
% Solves the GI/M/c QBD via R^2*A2 + R*A1 + A0 = 0, where
%   A0 = (-T*1)*alpha,  A1 = T - c*mu*I,  A2 = c*mu*I.
% pi_n = pi_c * R^(n-c) for n >= c. Boundary states pi_0..pi_c are obtained
% from balance equations + normalization sum_{n<c} pi_n*1 + pi_c*(I-R)\1 = 1.
%
% Inputs:
%   ALPHA - PH entry probability row vector (1 x k)
%   T     - PH sub-generator (k x k)
%   MU    - Exponential service rate per server (positive scalar)
%   C     - Number of servers (positive integer)
%
% Optional name-value:
%   'maxIter' - max iterations for R fixed-point (default 50000)
%   'tol'     - convergence tolerance (default 1e-14)
%
% RESULT struct fields:
%   meanQueueLength, meanWaitingQueue, meanWaitingTime,
%   meanSojournTime, utilization, analyzer

p = inputParser;
addParameter(p, 'maxIter', 50000);
addParameter(p, 'tol', 1e-14);
parse(p, varargin{:});
maxIter = p.Results.maxIter;
tol = p.Results.tol;

if mu <= 0
    line_error(mfilename, 'Service rate MU must be positive');
end
if c < 1 || floor(c) ~= c
    line_error(mfilename, 'C must be a positive integer');
end
alpha = alpha(:)';
k = size(T, 1);
if size(T, 2) ~= k || numel(alpha) ~= k
    line_error(mfilename, 'ALPHA / T dimensions inconsistent');
end

ones_k = ones(k, 1);
t_vec = -T * ones_k;
D1 = t_vec * alpha;

mean_ia = alpha * ((-T) \ ones_k);
if mean_ia <= 0
    line_error(mfilename, sprintf('Non-positive mean inter-arrival: %g', mean_ia));
end
lambda = 1.0 / mean_ia;
rho = lambda / (c * mu);
if rho >= 1 - 1e-12
    line_error(mfilename, sprintf('Load rho=%g must be strictly less than 1', rho));
end

A0 = D1;
A1 = T - c*mu*eye(k);
A2 = c*mu*eye(k);
R = zeros(k);
for it = 1:maxIter
    Mtmp = A1 + R*A2;
    if abs(det(Mtmp)) < 1e-30
        break;
    end
    R_new = -A0 / Mtmp;
    if max(max(abs(R_new - R))) < tol
        R = R_new;
        break;
    end
    R = R_new;
end

% Boundary linear system
n_var = (c+1)*k;
M = zeros(n_var, n_var);
block = @(n) (n*k);

% Level 0: pi_0 T + pi_1 (mu I) = 0
for j = 1:k
    for i = 1:k
        M(j, block(0)+i) = M(j, block(0)+i) + T(i, j);
    end
    M(j, block(1)+j) = M(j, block(1)+j) + mu;
end

% Levels 1..c-1: pi_{n-1} D1 + pi_n (T - n mu I) + pi_{n+1} ((n+1) mu I) = 0
for n = 1:c-1
    eqrow_base = n*k;
    A1n = T - n*mu*eye(k);
    for j = 1:k
        row = eqrow_base + j;
        for i = 1:k
            M(row, block(n-1)+i) = M(row, block(n-1)+i) + D1(i, j);
            M(row, block(n)+i) = M(row, block(n)+i) + A1n(i, j);
        end
        M(row, block(n+1)+j) = M(row, block(n+1)+j) + (n+1)*mu;
    end
end

% Level c: pi_{c-1} D1 + pi_c (T - c mu I + R*(c mu I)) = 0
A1c_plus_RA2 = T - c*mu*eye(k) + R*(c*mu*eye(k));
eqrow_base = c*k;
for j = 1:k
    row = eqrow_base + j;
    for i = 1:k
        M(row, block(c-1)+i) = M(row, block(c-1)+i) + D1(i, j);
        M(row, block(c)+i) = M(row, block(c)+i) + A1c_plus_RA2(i, j);
    end
end

% Replace last row with normalization
IR = eye(k) - R;
sum_geom = IR \ ones_k;
norm_row = zeros(1, n_var);
for n = 0:c-1
    for i = 1:k
        norm_row(block(n)+i) = 1.0;
    end
end
for i = 1:k
    norm_row(block(c)+i) = sum_geom(i);
end
M(end, :) = norm_row;
b = zeros(n_var, 1); b(end) = 1.0;

x = M \ b;
pis = cell(c+1, 1);
for n = 0:c
    pis{n+1} = x(block(n)+1 : block(n)+k);
end
pi_c = pis{c+1};

IR_inv = inv(IR);
IR_inv2 = IR_inv * IR_inv;
Lq = pi_c' * R * IR_inv2 * ones_k;
L_bulk = pi_c' * (c*IR_inv + R*IR_inv2) * ones_k;
L = 0;
for n = 0:c-1
    L = L + n * sum(pis{n+1});
end
L = L + L_bulk;
Wq = Lq / lambda;
W = Wq + 1.0/mu;

result = struct();
result.meanQueueLength = L;
result.meanWaitingQueue = Lq;
result.meanWaitingTime = Wq;
result.meanSojournTime = W;
result.utilization = rho;
result.analyzer = sprintf('PHMC:matrix-geom:c=%d', c);

end
