function [muf, Psi, mu] = pfqn_oi_fnc(Phi, N, f, options)
% [MUF, PSI, MU] = PFQN_OI_FNC(PHI, N, F, OPTIONS)
%
% Order-independent (OI) generalization of the load-dependent functional
% server (FNC) of Casale, "On Single-Class Load-Dependent Normalizing
% Constant Equations", QEST 2006 (Theorem 3, Corollary 1). It is the OI
% counterpart of PFQN_FNC.
%
% Given the balance function Phi of an existing OI station in the model and a
% queue-dependent target f(n) (default f(n) = sum(n), the total occupancy),
% the routine builds a *functional server*: an auxiliary OI station whose
% balance function Psi satisfies the convolution identity
%
%   (Psi * Phi)(n) = (1 + f(n)) Phi(n),        (*)
%
% where * is the OI (population-lattice) convolution
% (Psi * Phi)(n) = sum_{0<=k<=n} Psi(k) Phi(n-k). By construction Psi(0)=1
% (since f(0)=0). Inserting this server (with the same per-class demand as
% the existing station) into the network and letting G, G^{+} be the OI
% normalizing constants without and with it, (*) yields the FNC identity
%
%   E[f(n)] = G^{+}/G - 1,
%
% i.e. the mean of the queue-dependent function is read off a ratio of
% normalizing constants, with no probabilities and no Little's law. For
% f(n)=sum(n) this returns the exact total mean queue length of the station.
%
% Construction (two steps):
%   1. Deconvolution of (*) for the FNC balance function, solved triangularly
%      in column-major order (k<n precedes n):
%        Psi(n) = (1+f(n)) Phi(n) - sum_{0<=k<n} Psi(k) Phi(n-k).
%   2. Balanced-fairness inversion of Psi to the FNC rate function:
%        mu_f(n) = ( sum_{r: n_r>0} Psi(n-e_r) ) / Psi(n).
%   Step 2 alone is the scalar PFQN_FNC when R=1 and Phi is the trivial LI
%   balance (Phi==1): there Psi==1 and mu_f==1, i.e. the FNC degenerates to
%   an identical LI copy of the station, as in the QEST 2006 example.
%
% As noted in that reference (Sec. 6), the FNC balance/rate may be signed or
% non-physical; this is immaterial because only the final normalizing-constant
% ratio is used, and it stays positive with its usual interpretation. Hence
% the functional server is best convolved through its balance function Psi
% rather than through the positivity-guarded rate peeling of PFQN_NCOI.
%
% Parameters:
%   Phi - balance function of the existing OI station over the lattice, an
%         R-dimensional array of size (N_1+1) x ... x (N_R+1) (Phi(n) at
%         subscript n+1), or a flat column-major vector. Obtain it by the
%         forward balanced-fairness recursion from the station rate function.
%   N   - (1 x R) closed population vector, finite. Optional when Phi is a
%         full R-dimensional array (then N = size(Phi) - 1).
%   f   - target queue-dependent function handle f(n), n a (1 x R) count
%         vector, with f(0)=0 (default f = @(n) sum(n)).
%   options - solver options (optional, accepted for signature parity).
%
% Returns:
%   muf - function handle muf(n) giving the FNC rate mu_f(n) (out-of-lattice
%         states return Inf). mu_f(0)=0; states with Psi(n)=0 return Inf.
%   Psi - R-dimensional array (lattice shape) with the FNC balance function.
%   mu  - R-dimensional array with the tabulated FNC rate mu_f.
%
% See also PFQN_FNC, PFQN_NCOI.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 4
    options = struct(); %#ok<NASGU>
end
if nargin < 3 || isempty(f)
    f = @(n) sum(n);
end
if nargin < 2 || isempty(N)
    N = size(Phi) - 1;
end
N = round(N(:)');
R = numel(N);
shp = N + 1;

Phiv = Phi(:);
total = prod(shp);
if numel(Phiv) ~= total
    line_error(mfilename, 'numel(Phi) must equal prod(N+1).');
end

% Column-major strides and decoded subscripts.
stride = ones(1, R);
for d = 2:R, stride(d) = stride(d-1) * shp(d-1); end
subs = zeros(total, R);
for i = 1:total
    li = i - 1;
    for d = 1:R, subs(i, d) = mod(li, shp(d)); li = floor(li / shp(d)); end
end

% ---- Step 1: deconvolve (Psi * Phi)(n) = (1 + f(n)) Phi(n) for Psi --------
Psiv = zeros(total, 1);
for i = 1:total
    n = subs(i, :);
    acc = (1 + f(n)) * Phiv(i);
    % subtract sum over proper sub-lattice k < n of Psi(k) Phi(n-k)
    for j = 1:(i-1)
        k = subs(j, :);
        if all(k <= n)
            idx = 1 + sum((n - k) .* stride);   % linear index of n-k
            acc = acc - Psiv(j) * Phiv(idx);
        end
    end
    Psiv(i) = acc;
end

% ---- Step 2: balanced-fairness inversion of Psi to the FNC rate -----------
muv = inf(total, 1);
for i = 1:total
    n = subs(i, :);
    if all(n == 0)
        muv(i) = 0;                 % empty state
        continue
    end
    denom = Psiv(i);
    if denom == 0
        muv(i) = Inf;               % non-physical / undefined rate
        continue
    end
    num = 0;
    for r = 1:R
        if n(r) > 0, num = num + Psiv(i - stride(r)); end
    end
    muv(i) = num / denom;
end

if R == 1
    Psi = Psiv; mu = muv;
else
    Psi = reshape(Psiv, shp); mu = reshape(muv, shp);
end
muf = @(nq) oi_fnc_eval(nq, muv, shp, stride);
end

function val = oi_fnc_eval(n, muv, shp, stride)
n = round(n(:)');
if numel(n) ~= numel(shp) || any(n < 0) || any(n > shp - 1)
    val = Inf;                      % outside the tabulated lattice
    return
end
val = muv(1 + sum(n .* stride));
end
