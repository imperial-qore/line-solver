function ilt = matlab_ilt(fun, T, maxFnEvals, method)
% ILT = MATLAB_ILT(FUN, T, MAXFNEVALS, METHOD)
% Numerical inverse Laplace transform in the Abate-Whitt framework.
%
%   FUN         handle evaluating the Laplace-domain transform at a complex s.
%               It is called once per (time point, node) pair and need not be
%               vectorized.
%   T           vector of positive real time points.
%   MAXFNEVALS  budget of Laplace-domain evaluations. For 'cme' it selects the
%               tabulated parameter set; for 'euler'/'gaver' it sets the number
%               of summation terms.
%   METHOD      'cme' (default), 'euler' or 'gaver'.
%
%   Every method reduces to the same weighted sum
%
%       f(t) ~ (1/t) * sum_k Re[ eta_k * F(beta_k / t) ]
%
%   and differs only in how (eta, beta) are built. ILT has the shape of T.
%
%   THIS FILE WAS MISSING FROM THE MATLAB TREE UNTIL 2026-09-11, while
%   `solver_mam_transient_qbd.m` and `laplace_invert_cme.m` both called it and
%   the C++ and python ports both carried it (`cpp/include/line/api/mam/
%   matlab_ilt.h`, `python/line_solver/lib/thirdparty/iltcme/matlab_ilt.py`).
%   Only `iltcme.json` had been vendored here. Nothing caught it because
%   SolverMAM's 'ldqbd' gate refused every model that would have reached
%   `solver_mam_transient_qbd`, so the one route into this function was
%   unreachable; see the note in @SolverMAM/SolverMAM.m. The port is taken from
%   the python twin, whose comments quote the original formulas.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 4 || isempty(method)
    method = 'cme';
end
if nargin < 3 || isempty(maxFnEvals)
    line_error(mfilename, 'matlab_ilt requires an evaluation budget maxFnEvals.');
end

switch lower(method)
    case 'cme'
        [eta, beta] = iltcme_cme_weights(maxFnEvals);
    case 'euler'
        [eta, beta] = iltcme_euler_weights(maxFnEvals);
    case 'gaver'
        [eta, beta] = iltcme_gaver_weights(maxFnEvals);
    otherwise
        line_error(mfilename, sprintf(['Unknown inverse Laplace transform ' ...
            'method ''%s''. Supported: cme, euler, gaver.'], method));
end

Tin = T;
t = double(T(:)).';
ilt = zeros(1, numel(t));
for i = 1:numel(t)
    % Summed in complex arithmetic and realified once, which is the python
    % twin's sum(real(...)) up to floating-point associativity.
    acc = 0;
    for k = 1:numel(eta)
        acc = acc + eta(k) * fun(beta(k) / t(i));
    end
    ilt(i) = real(acc) / t(i);
end

% A caller writing into a logical slice (EN(posMask) = matlab_ilt(...)) needs
% the orientation it passed in, not this function's working row.
if size(Tin, 2) == 1 && size(Tin, 1) > 1
    ilt = ilt(:);
end
end


function [eta, beta] = iltcme_cme_weights(maxFnEvals)
% The tabulated Talbot-style contour: among the entries that fit the budget,
% the most concentrated one (smallest cv2). The first entry seeds the search
% whether or not it fits, which is what the python twin and the original do --
% it is the fallback when no entry fits at all.
global cmeParams;
if isempty(cmeParams)
    % Located on the path, as CME.table does; both share this one decode.
    cmeParams = jsondecode(fileread('iltcme.json'));
end

best = iltcme_entry(cmeParams, 1);
for i = 2:numel(cmeParams)
    cand = iltcme_entry(cmeParams, i);
    if cand.cv2 < best.cv2 && cand.n + 1 <= maxFnEvals
        best = cand;
    end
end

a = double(best.a(:)).';
b = double(best.b(:)).';
c = double(best.c);
mu1 = double(best.mu1);
omega = double(best.omega);
n = double(best.n);

eta = [c * mu1, (a + 1i * b) * mu1];
k = 1:n;
beta = [mu1, (1 + 1i * k * omega) * mu1];
end


function entry = iltcme_entry(params, i)
% jsondecode returns a struct array when every entry has the same field shapes
% and a cell array otherwise; iltcme.json's a/b run from length 1 to 1000, so
% it is the cell case here. CME.tableEntry reads the same table the same way.
if iscell(params)
    entry = params{i};
else
    entry = params(i);
end
end


function [eta, beta] = iltcme_euler_weights(maxFnEvals)
% Euler summation of the Bromwich integral: binomial weights accumulated in
% log space so the factorials do not overflow, on geometric abscissae.
nEuler = floor((maxFnEvals - 1) / 2);

eta = zeros(1, 2 * nEuler + 1);
eta(1) = 0.5;
eta(2:nEuler+1) = 1.0;
eta(2 * nEuler + 1) = 2^(-nEuler);
for k = 1:nEuler-1
    eta(2 * nEuler - k + 1) = eta(2 * nEuler - k + 2) + ...
        exp(sum(log(1:nEuler)) - nEuler * log(2) ...
            - sum(log(1:k)) - sum(log(1:(nEuler - k))));
end

k = 0:2 * nEuler;
beta = nEuler * log(10) / 3 + 1i * pi * k;
eta = 10^(nEuler / 3) * (1 - 2 * mod(k, 2)) .* eta;
end


function [eta, beta] = iltcme_gaver_weights(maxFnEvals)
% Gaver functional with Stehfest's weights: real logarithmic abscissae, and
% the summand formed in log space for the same overflow reason.
if mod(maxFnEvals, 2) == 1
    maxFnEvals = maxFnEvals - 1;
end
ndiv2 = maxFnEvals / 2;

eta = zeros(1, maxFnEvals);
beta = zeros(1, maxFnEvals);
ln2 = log(2);

for k = 1:maxFnEvals
    insideSum = 0;
    for j = floor((k + 1) / 2):min(k, ndiv2)
        insideSum = insideSum + ...
            exp((ndiv2 + 1) * log(j) ...
                - sum(log(1:(ndiv2 - j))) ...
                + sum(log(1:2 * j)) ...
                - 2 * sum(log(1:j)) ...
                - sum(log(1:(k - j))) ...
                - sum(log(1:(2 * j - k))));
    end
    eta(k) = ln2 * (-1)^(k + ndiv2) * insideSum;
    beta(k) = k * ln2;
end
end
