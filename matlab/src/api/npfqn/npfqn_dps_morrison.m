function [Q, R, X, aux] = npfqn_dps_morrison(N, Z, S, w)
% [Q,R,X,AUX] = NPFQN_DPS_MORRISON(N, Z, S, W)
%
% Two-term heavy-usage asymptotic approximation for a closed queueing network
% with one infinite-server (think) station and one discriminatory
% processor-sharing (DPS) station, after
%
%   J.A. Morrison, "Asymptotic analysis of a large closed queueing network
%   with discriminatory processor sharing", Queueing Systems 9 (1991) 191-214.
%
% The network is NOT product-form, so there is no normalizing constant here:
% the method expands the GENERATING FUNCTION of the balance equations. The
% substitution P(n) = <w,n> f(n) clears the DPS denominator and turns the
% balance recursion into a linear PDE with affine coefficients (Morrison eq.
% 2.5); rescaling z = 1 - xi/sqrt(N) and expanding in powers of N^(-1/2) leaves
% a degenerate leading operator whose kernel is the functions of the similarity
% variable eta, and the solvability condition along its characteristic gives an
% ODE for the amplitude (eq. 2.20). RESULT 1 (eq. 4.11) and RESULT 2 (eq. 4.17)
% are the two-term approximations returned here.
%
% Input:
%   N - per-class populations (1,K), positive
%   Z - per-class mean think time at the delay station (1,K), positive
%   S - per-class mean service time at the DPS station (1,K), positive
%   w - per-class DPS weights (1,K), positive
%
% Output:
%   Q   - mean number of class-k jobs at the DPS station (1,K)
%   R   - mean class-k sojourn time per visit to the DPS station (1,K)
%   X   - per-class throughput (1,K), by Little's law on the think station
%   aux - struct of intermediate quantities: the usage rho, the heavy-usage
%         parameter a, Morrison's constants, the vector sigma of eq. (4.12),
%         the W_m of eq. (3.23), and Qlead/Rlead, the leading-order (one-term)
%         values, whose gap to Q/R measures the size of the correction.
%
% Scaling. Morrison writes K_j = N b_j and lambda_j = N r_j g_j with N large
% and usage rho = sum_j b_j/g_j = 1 - a/sqrt(N). N is bookkeeping only: b, g
% and a all move with it and the approximation is invariant (verified
% numerically over four decades of N), so this function fixes N = 1, i.e.
% b = N_pop, g = Z./S, r = 1./Z and a = 1 - rho. The accuracy is governed by
% the PHYSICAL regime -- large populations with rho near 1 -- and not by any
% choice made here. rho > 1 (a < 0) is admissible: it is the saturated regime
% of Morrison's appendix A, where the leading term agrees with the Mitra-Weiss
% fluid approximation.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

N = N(:).'; Z = Z(:).'; S = S(:).'; w = w(:).';
p = numel(N);
if numel(Z) ~= p || numel(S) ~= p || numel(w) ~= p
    line_error(mfilename, 'N, Z, S and w must have the same number of classes.');
end
if any(~isfinite(N)) || any(N <= 0)
    line_error(mfilename, 'The Morrison approximation requires finite positive class populations (closed classes only).');
end
if any(~isfinite(Z)) || any(Z <= 0) || any(~isfinite(S)) || any(S <= 0)
    line_error(mfilename, 'Think times Z and DPS service times S must be finite and positive.');
end
if any(~isfinite(w)) || any(w <= 0)
    line_error(mfilename, 'DPS weights must be finite and positive.');
end

% ---- Morrison's parameters at the bookkeeping scale N = 1 ------------------
b = N;                       % K_j = N b_j
r = 1 ./ Z;                  % think rate
g = Z ./ S;                  % lambda_j = N r_j g_j
rho = sum(b ./ g);           % eq. (1.2): usage
a = 1 - rho;                 % eq. (1.2): rho = 1 - a/sqrt(N)

% ---- constants, eqs. (2.18), (2.19), (3.11), (3.13)-(3.15) -----------------
cB = sum(b ./ (r .* g.^2 .* w));
cC = sum(b ./ (g.^2 .* w));
cD = sum(b ./ (r .* g.^2));
cH = sum(b ./ (r.^2 .* g.^3 .* w));
cI = sum(b ./ (r.^2 .* g.^3 .* w.^2));
cJ = sum(b ./ (r .* g.^3 .* w.^2));
cK = sum(b ./ (g.^3 .* w.^2));
cL = sum(b ./ (r .* g.^3 .* w));
cM = sum(b ./ (r.^2 .* g.^3));
cQ = sum(b ./ g.^2);

% ---- sigma: eq. (4.12) with the normalization (4.13) -----------------------
% The p equations have rank p-1 (sum_i sigma_i/(r_i g_i) annihilates both
% sides, Morrison p.197), so (4.13) is appended and the system solved in the
% least-squares sense, which is exact because the system is consistent.
Asys = zeros(p, p);
rhs = zeros(p, 1);
for i = 1:p
    Asys(i, i) = Asys(i, i) + rho;
    for j = 1:p
        den = r(i) * g(i) * w(i) + r(j) * g(j) * w(j);
        Asys(i, i) = Asys(i, i) - w(j) * b(j) * r(j) / den;
        Asys(i, j) = Asys(i, j) - w(j) * b(i) * r(i) / den;
    end
    rhs(i) = rho * (b(i) / g(i)) * (cD / (cB * w(i)) - 1);
end
sigma = ([Asys; 1 ./ (r .* g)] \ [rhs; 0]).';

% alpha from eq. (4.9), then delta of eq. (3.15)
alpha = sigma - (b ./ g) .* (cD ./ (cB * w) - 1);
delta = sum(alpha ./ g);

% ---- eqs. (3.19)-(3.21) ---------------------------------------------------
cR = 3 * (cB * cL - cD * cJ) / (cB * cD);
cS = (2 * cB * (cD * cH - cB * cM) - cD * (cD * cI - cB * cH)) / (2 * cB^2 * cD^2);
cU = (cQ - cC * cD / cB - delta) / rho - cD * cR / cB + (a^2 - cC * cD / cB) * cS;
cA = cS * cC^2 + cR * cC - cK;
cV = cR + 2 * cS * cC;

% ---- W_m of eq. (3.23) ----------------------------------------------------
W = local_W(cB, cC, cD, a);
W0 = W(1); W1 = W(2); W2 = W(3); W3 = W(4); W4 = W(5);

% ---- RESULT 1 (4.11) and RESULT 2 (4.17), at sqrt(N) = 1 ------------------
% NOTE the numerator bracket carries U*W2: eq. (4.10) of the paper misprints it
% as U*W1, but (4.7), (4.11), (A6) and (B2) all agree on U*W2, and it is what
% the derivation from (4.4)-(4.9) gives.
eps = cB / cD;
num = W1 - eps * (cA / 3 * W4 + a / 2 * cV * W3 + cU * W2);
den = W0 - eps * (cA / 3 * W3 + a / 2 * cV * W2 + cU * W1 + cS);
if den == 0 || ~isfinite(den)
    line_error(mfilename, 'The Morrison expansion is degenerate for this model (vanishing denominator); the usage is too far from the moderately-heavy regime.');
end

Qlead = b .* W1 ./ (g .* w * W0);
Q = b .* num ./ (g .* w * den) ...
    - b * W2 ./ (g.^2 .* w.^2 * W0) ...
    - sigma / rho;

Rlead = W1 ./ (r .* g .* w * W0);
R = num ./ (r .* g .* w * den) ...
    + ((W1 / W0)^2 - W2 / W0) ./ (r .* g.^2 .* w.^2) ...
    - sigma ./ (rho * r .* b);

X = r .* (b - Q);            % Little's law at the think station

% Morrison's constants keep their paper names under a 'c' prefix: aux.cQ, aux.cR
% and aux.cS are the constants of (3.15)/(3.19), NOT the queue lengths, response
% times or service times this function otherwise deals in.
aux = struct('rho', rho, 'a', a, 'cB', cB, 'cC', cC, 'cD', cD, 'cH', cH, ...
    'cI', cI, 'cJ', cJ, 'cK', cK, 'cL', cL, 'cM', cM, 'cQ', cQ, 'delta', delta, ...
    'cR', cR, 'cS', cS, 'cU', cU, 'cA', cA, 'cV', cV, 'sigma', sigma, ...
    'alpha', alpha, 'W', W, 'Qlead', Qlead, 'Rlead', Rlead);
end

% ==========================================================================
function W = local_W(cB, cC, cD, y)
% W(m+1) = W_m(y) of eq. (3.23),
%   W_m(y) = (B/D)^2 int_0^inf exp(-(BC/2D) z^2) exp(-(B/D) y z) z^m dz.
% Substituting z = sigma s with sigma = sqrt(D/(BC)) normalizes the Gaussian,
%   W_m(y) = (B/D)^2 sigma^(m+1) I_m(yh),  yh = y sqrt(B/(CD)),
%   I_m(yh) = int_0^inf s^m exp(-s^2/2 - yh s) ds,
% and the I_m follow from I_0 = sqrt(pi/2) erfcx(yh/sqrt(2)) -- the scaled
% complementary error function, which is what keeps large yh from overflowing --
% with I_1 = 1 - yh I_0 and I_m = (m-1) I_{m-2} - yh I_{m-1}. That recursion
% subtracts nearly equal terms once yh is large, so a loss of positivity (the
% I_m are integrals of positive integrands) triggers a quadrature fallback.
mmax = 4;
sig = sqrt(cD / (cB * cC));
yh = y * sqrt(cB / (cC * cD));

Iv = zeros(1, mmax + 1);
Iv(1) = sqrt(pi / 2) * erfcx(yh / sqrt(2));
if ~isfinite(Iv(1))
    line_error('npfqn_dps_morrison', ['The usage is so far above saturation (rho = %g) that the ' ...
        'Morrison expansion overflows. This model is outside the moderately-heavy regime the ' ...
        'approximation is derived for; use SolverFLD, SolverMVA or SolverCTMC.'], 1 - y);
end
Iv(2) = 1 - yh * Iv(1);
for m = 2:mmax
    Iv(m + 1) = (m - 1) * Iv(m - 1) - yh * Iv(m);
end
if any(Iv <= 0)
    for m = 0:mmax
        Iv(m + 1) = integral(@(s) s.^m .* exp(-s.^2 / 2 - yh * s), 0, Inf, ...
            'AbsTol', 1e-14, 'RelTol', 1e-12);
    end
end
W = (cB / cD)^2 * sig.^(1:(mmax + 1)) .* Iv;
end
