function [G,lG] = pfqn_conv(L, N, Z, cdscaling, options)
% [G,LG] = PFQN_CONV(L, N, Z, CDSCALING, OPTIONS)
%
% Multichain convolution algorithm for closed queueing networks with
% class-dependent service rates.
%
% Implements the convolution algorithm of Sauer (1983), Section 5.2,
% "Computational Algorithms for State-Dependent Queueing Networks",
% ACM TOCS, Vol. 1, No. 1, pp. 67-92.
%
% The algorithm computes G(N) = (X_1 * X_2 * ... * X_M)(N) where X_m(n)
% is the station factor at population vector n, and * denotes the
% multivariate discrete convolution:
%   A(n) = sum_{i: 0<=i<=n} B(i) * C(n-i)
%
% For class-dependent stations, X_m(n) is computed recursively via Sauer
% eq. (40):
%   X_m(n) = (u_km / mu_km(n)) * X_m(n - e_k)
% where mu_km(n) = (n_k/|n|) * beta_{m,k}(n) and beta is the DIMENSIONLESS
% class-dependent scaling of the service demand supplied by CDSCALING{m}: a
% handle of the per-class population vector n at station m, returning either a
% scalar (shared by every class) or a length-R vector. Equivalently
%   X_m(n) = (|n|/n_k) * (L(m,k)/beta_{m,k}(n)) * X_m(n - e_k),
% which at beta = 1 is exactly the load-independent multinomial form, so a unit
% scaling means "no correction". This is the same convention AMVA and CTMC use
% (effective service time ST/beta). Any saturation/cutoff is applied inside the
% handle.
%
% For standard (load-independent) stations, X_m(n) reduces to the
% multinomial form and the convolution uses the efficient recurrence:
%   G_m(n) = G_{m-1}(n) + sum_r L(m,r) * G_m(n - e_r)
%
% Parameters:
%   L  - Service demand matrix (M x R)
%   N  - Population vector (1 x R), must be finite (closed network)
%   Z  - Think time vector (1 x R), or empty
%   cdscaling   - Cell array {M,1} of class-dependence handles beta_m(n);
%                 empty entries denote load-independent stations
%   options - Solver options (optional)
%
% Returns:
%   G  - Normalizing constant G(N)
%   lG - log(G(N))

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

[M, R] = size(L);

if nargin < 3 || isempty(Z)
    Z = zeros(1, R);
end
if nargin < 4 || isempty(cdscaling)
    cdscaling = cell(M, 1);
end

if any(~isfinite(N))
    line_error(mfilename, 'Convolution algorithm requires finite (closed) populations.');
end

% Total state space size
stateSpaceSize = prod(N + 1);

% Identify which stations carry a class-dependence function
isCd = false(M, 1);
for ist = 1:M
    isCd(ist) = ~isempty(cdscaling{ist});
end

% --- Precompute X_m(n) tables for class-dependent stations ---
% see _kb/03-api-layer.md (pfqn/ family: scaling, log-domain switches, dispatch gates)

Xm = cell(M, 1);
for ist = 1:M
    if isCd(ist)
        Xm{ist} = zeros(stateSpaceSize, 1);
        Xm{ist}(1) = 1; % X_m(0) = 1

        % Enumerate all population vectors and build X_m(n) recursively
        n = pprod(N);
        while n(1) >= 0
            idx = hashpop(n, N);
            if sum(n) == 0
                Xm{ist}(idx) = 1;
            else
                % Use eq. (40): X_m(n) = (u_km / mu_km(n)) * X_m(n-e_k)
                % Pick first class k with n(k) > 0
                for r = 1:R
                    if n(r) > 0
                        % Get service rate mu_km(n) from the class-dependence handle
                        % see _kb/03-api-layer.md (pfqn/ family: scaling, log-domain switches, dispatch gates)
                        bval = cdscaling{ist}(n);
                        if numel(bval) > 1
                            beta = bval(r);
                        else
                            beta = bval;
                        end

                        % see _kb/03-api-layer.md (pfqn/ family: scaling, log-domain switches, dispatch gates)
                        tot = sum(n);
                        nr = n(r);
                        n(r) = n(r) - 1;
                        idx_prev = hashpop(n, N);
                        n(r) = n(r) + 1;
                        if beta > 0
                            Xm{ist}(idx) = (tot / nr) * (L(ist, r) / beta) * Xm{ist}(idx_prev);
                        else
                            Xm{ist}(idx) = 0;
                        end
                        break
                    end
                end
            end
            n = pprod(n, N);
        end
    end
end

% --- Convolution ---
% G_0(n) = F_Z(n) (delay contribution)
% G_m(n) = (G_{m-1} * X_m)(n) for class-dependent stations
% G_m(n) = G_{m-1}(n) + sum_r L(m,r) * G_m(n-e_r) for LI stations

G_curr = zeros(stateSpaceSize, 1);

% Initialize G_0(n) = F_Z(n): delay server contribution
n = pprod(N);
while n(1) >= 0
    idx = hashpop(n, N);
    G_curr(idx) = Fz(Z, n);
    n = pprod(n, N);
end

% Convolve one station at a time
for ist = 1:M
    if isCd(ist)
        % class-dependent station: direct convolution sum
        % G_new(n) = sum_{i: 0<=i<=n} X_m(i) * G_old(n-i)
        G_old = G_curr;
        G_curr = zeros(stateSpaceSize, 1);

        n = pprod(N);
        while n(1) >= 0
            idx_n = hashpop(n, N);
            conv_sum = 0;

            % Inner loop: enumerate all i from 0 to n
            i = pprod(n);
            while i(1) >= 0
                idx_i = hashpop(i, N);
                nmi = n - i; % n - i (component-wise)
                idx_nmi = hashpop(nmi, N);
                conv_sum = conv_sum + Xm{ist}(idx_i) * G_old(idx_nmi);
                i = pprod(i, n);
            end

            G_curr(idx_n) = conv_sum;
            n = pprod(n, N);
        end
    else
        % Load-independent station: efficient recurrence
        % G_m(n) = G_{m-1}(n) + sum_r L(m,r) * G_m(n - e_r)
        n = pprod(N);
        while n(1) >= 0
            idx_n = hashpop(n, N);
            % G_curr(idx_n) already has G_{m-1}(n) from previous iteration
            for r = 1:R
                if n(r) >= 1
                    n(r) = n(r) - 1;
                    idx_n1r = hashpop(n, N);
                    n(r) = n(r) + 1;
                    G_curr(idx_n) = G_curr(idx_n) + L(ist, r) * G_curr(idx_n1r);
                end
            end
            n = pprod(n, N);
        end
    end
end

G = G_curr(end); % G(N) is at the last index (hashpop(N,N) = prod(N+1))
lG = log(G);
end

%% --- Local functions ---

function idx = hashpop(n, N)
% HASHPOP Map population vector to linear index (1-based)
% idx = 1 + n(1) + n(2)*(N(1)+1) + n(3)*(N(1)+1)*(N(2)+1) + ...
idx = 1;
R = length(N);
for r = 1:R
    idx = idx + prod(N(1:r-1) + 1) * n(r);
end
end

function [n] = pprod(n, N)
% PPROD Sequentially generate all vectors n: 0 <= n <= N
% n = pprod(N)    - initialize to zeros
% n = pprod(n, N) - advance to next vector, returns n(1)=-1 when done
if nargin == 1
    N = n;
    n = zeros(size(N));
    return
end

R = length(N);
if sum(n == N) == R
    n = -1 * ones(1, R);
    return
end

s = R;
while s > 0 && n(s) == N(s)
    n(s) = 0;
    s = s - 1;
end
if s > 0
    n(s) = n(s) + 1;
end
end

function f = Fz(Z, n)
% FZ Delay server unnormalized probability factor
% F = (Z(1)^n(1) / n(1)!) * ... * (Z(R)^n(R) / n(R)!)
R = length(n);
if sum(n) == 0
    f = 1;
    return
end
f = 0;
for r = 1:R
    if Z(r) > 0
        f = f + log(Z(r)) * n(r);
        f = f - gammaln(1 + n(r));
    elseif n(r) > 0
        f = 0;
        return
    end
end
f = exp(f);
end
