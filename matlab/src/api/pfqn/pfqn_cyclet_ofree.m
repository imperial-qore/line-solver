function [f, F, mom, out] = pfqn_cyclet_ofree(v, mu, N, path, tset, options)
% [f, F, MOM, OUT] = PFQN_CYCLET_OFREE(V, MU, N, PATH, TSET, OPTIONS)
%
% Exact passage-time density, CDF and moments along an OVERTAKE-FREE PATH of a
% closed single-chain tree-like product-form network with population N.
%
% V (1,M) visit ratios, MU (1,M) service rates, N the population, PATH the node
% list z = (z_1,...,z_m) of the overtake-free path with z_1 the root, TSET the
% time grid. PATH may instead be a cell array of paths, in which case
% OPTIONS.pathprob weights them and the outputs are the mixture; that is how a
% cycle time is assembled when the root branches.
%
% Reference: P. G. Harrison and W. J. Knottenbelt, "Passage Time Distributions
% in Large Markov Chains", 2002, Sec. 7.1, Theorems 1 and 2, after
% P. G. Harrison, J. Appl. Prob. 27, 1990 and H. Duduna, Adv. Appl. Prob. 14,
% 1982. The underlying sojourn-time result for overtake-free paths is
% F. Kelly and P. Pollett, Adv. Appl. Prob. 15, 1983.
%
% THE ONE FACT THAT MAKES ALL THREE ROUTES WORK. Conditional on the path,
%
%     T | z  =  sum_{j in z} Erlang(u_{z_j} + 1, mu_{z_j})
%
% with u distributed as the network's equilibrium population vector AT N-1 (the
% arrival theorem). Hence the transform of Theorem 1 collapses to
%
%     L(s|z) = prod_{j in z} mu_j/(s+mu_j) * G(y(s), N-1) / G(x, N-1)
%
% where x_i = v_i/mu_i and y_i(s) = x_i mu_i/(s+mu_i) on the path, x_i off it.
% One Buzen convolution per value of s.
%
% OPTIONS.method selects the density route:
%   'auto'  (default) 'exact' when the path rates are separated, else 'lt'
%   'exact' Theorem 2 in closed form; REQUIRES DISTINCT RATES on the path,
%           since its partial fractions divide by prod_{i~=j}(mu_i - mu_j)
%   'lt'    the transform above inverted through api/lti
%           (OPTIONS.lti_method, default 'euler')
%
% MOMENTS ARE NEVER TAKEN FROM THE DENSITY. They come from running the same
% Buzen convolution in the ring of truncated power series in s, so they are
% exact to machine precision, are unaffected by the time grid, and stay valid
% when the rates coincide and Theorem 2 does not apply.
%
% NOTE ON THE PAPER. The inner sum of Theorem 2 reads (v_j t)^(c-i)/(c-i)! and
% that is CORRECT as printed, however odd the visit ratio looks against a time:
% substituting the service rate instead returns negative densities. Verified
% against a direct mixture-of-Erlangs oracle to 1e-15, and at the paper's own
% N = 18 example against the transform route to 1e-11.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 6 || isempty(options)
    options = struct();
end
if ~isfield(options,'method') || isempty(options.method)
    options.method = 'auto';
end
if ~isfield(options,'nmom') || isempty(options.nmom)
    options.nmom = 3;
end
if ~isfield(options,'lti_method') || isempty(options.lti_method)
    options.lti_method = 'euler';
end
if ~isfield(options,'tol') || isempty(options.tol)
    options.tol = 1e-8;
end

v = reshape(v, 1, []);
mu = reshape(mu, 1, []);
M = numel(v);
if numel(mu) ~= M
    line_error(mfilename, 'V and MU must name the same number of nodes.');
end
if any(mu <= 0)
    line_error(mfilename, 'Every service rate must be positive.');
end
if N < 1 || N ~= round(N)
    line_error(mfilename, 'The population N must be a positive integer.');
end

if ~iscell(path)
    paths = {reshape(path,1,[])};
    pathprob = 1;
else
    paths = cellfun(@(z) reshape(z,1,[]), path, 'UniformOutput', false);
    if isfield(options,'pathprob') && ~isempty(options.pathprob)
        pathprob = reshape(options.pathprob,1,[]);
    else
        pathprob = ones(1,numel(paths))/numel(paths);
    end
    if numel(pathprob) ~= numel(paths)
        line_error(mfilename, 'OPTIONS.pathprob must carry one probability per path.');
    end
end

tset = reshape(tset, 1, []);
f = zeros(size(tset));
F = zeros(size(tset));
mom = zeros(1, options.nmom);
out = struct('method', cell(1,numel(paths)), 'lG', [], 'path', []);

x = v ./ mu;
lGfull = local_conv(x, N-1);
Gn1 = lGfull(end);
if Gn1 <= 0
    line_error(mfilename, 'The network normalizing constant at population N-1 vanished; check V and MU.');
end

for ip = 1:numel(paths)
    z = paths{ip};
    if isempty(z)
        line_error(mfilename, 'An overtake-free path must contain at least the root node.');
    end
    if any(z < 1) || any(z > M) || numel(unique(z)) ~= numel(z)
        line_error(mfilename, 'A path must be a set of distinct node indices within the network.');
    end

    method = options.method;
    if strcmpi(method,'auto')
        mp = mu(z);
        if numel(mp) == 1
            method = 'exact';
        else
            d = abs(mp(:) - mp(:).');
            d(logical(eye(numel(mp)))) = Inf;
            if min(d(:)) > options.tol * max(mp)
                method = 'exact';
            else
                method = 'lt';
            end
        end
    end

    switch lower(method)
        case 'exact'
            [fi, Fi] = local_thm2(v, mu, N, z, tset, Gn1, x);
        case 'lt'
            Lfun = @(s) local_lst(v, mu, N, z, s, Gn1);
            fi = laplace_invert_pdf(Lfun, tset, options.lti_method);
            Fi = laplace_invert_cdf(Lfun, tset, options.lti_method);
        otherwise
            line_error(mfilename, sprintf('Unknown method: %s. Supported: auto, exact, lt.', method));
    end

    momi = local_moments(v, mu, N, z, Gn1, x, options.nmom);

    f = f + pathprob(ip) * fi;
    F = F + pathprob(ip) * Fi;
    mom = mom + pathprob(ip) * momi;
    out(ip).method = method;
    out(ip).lG = log(Gn1);
    out(ip).path = z;
end
end

% -------------------------------------------------------------------------

function g = local_conv(y, n)
% Buzen's convolution: g(k+1) = G at population k for the node set y, k = 0..n.
% This is the k(y,a,b) recursion of Sec. 7.1 with the node index rolled up,
% k(y,a,b) = k(y,a-1,b) + y_a k(y,a,b-1), k(y,a,0) = 1, k(y,0,b>0) = 0.
g = zeros(1, n+1);
g(1) = 1;
for i = 1:numel(y)
    for k = 1:n
        g(k+1) = g(k+1) + y(i) * g(k);
    end
end
end

function [f, F] = local_thm2(v, mu, N, z, tset, Gn1, x)
% Theorem 2 in closed form. The density is a finite sum of terms
% t^k exp(-mu_j t), so its integral is an incomplete gamma and the CDF comes
% out in closed form too rather than by quadrature.
M = numel(v);
m = numel(z);
off = setdiff(1:M, z);
mup = mu(z);
vp = v(z);

Gm = local_conv(x(off), N-1);   % constants of the network minus the path

coef = zeros(m, N);             % coef(j,k+1) multiplies t^k exp(-mu_j t)
for j = 1:m
    den = 1;
    for i = 1:m
        if i ~= j
            den = den * (mup(i) - mup(j));
        end
    end
    if den == 0
        line_error(mfilename, 'Theorem 2 needs distinct service rates on the path; two coincide. Use OPTIONS.method = ''lt''.');
    end
    idx = setdiff(1:m, j);
    w = (vp(idx) - vp(j)) ./ (mup(idx) - mup(j));
    K = local_conv(w, N-1);     % K^m(j,l), l = 0..N-1
    for c = 0:(N-1)
        Gmc = Gm(N-c);          % G_m(N-c-1)
        if Gmc == 0
            continue
        end
        for i = 0:c
            k = c - i;
            coef(j,k+1) = coef(j,k+1) + Gmc * K(i+1) / den;
        end
    end
end

pref = prod(mup) / Gn1;
f = zeros(size(tset));
F = zeros(size(tset));
for j = 1:m
    kk = find(coef(j,:) ~= 0) - 1;
    for k = kk
        cjk = coef(j,k+1);
        f = f + pref * cjk * (vp(j).^k) .* (tset.^k) ./ factorial(k) .* exp(-mup(j)*tset);
        % int_0^t s^k exp(-mu s) ds = k!/mu^(k+1) * gammainc(mu t, k+1)
        F = F + pref * cjk * (vp(j).^k) / (mup(j)^(k+1)) .* gammainc(mup(j)*tset, k+1);
    end
end
f = max(f, 0);
F = min(max(F, 0), 1);
end

function L = local_lst(v, mu, N, z, s, Gn1)
% L(s|z) = prod_j mu_j/(s+mu_j) * G(y(s), N-1) / G(x, N-1).
M = numel(v);
y = complex(v ./ mu);
for j = z
    y(j) = y(j) * (mu(j) / (s + mu(j)));
end
gy = local_conv(y, N-1);
L = gy(end) / Gn1;
for j = z
    L = L * (mu(j) / (s + mu(j)));
end
end

function mom = local_moments(v, mu, N, z, Gn1, x, nmom)
% The same Buzen convolution run in the ring of truncated power series in s.
% Every operation in the recursion is an addition or a multiplication, so the
% series ring carries it unchanged, and E[T^q] = (-1)^q q! [s^q] L(s).
K = nmom;
M = numel(v);
Y = zeros(M, K+1);
for i = 1:M
    if any(z == i)
        Y(i,:) = x(i) * ((-1./mu(i)) .^ (0:K));
    else
        Y(i,1) = x(i);
    end
end

G = zeros(N, K+1);          % G(n+1,:) is the series of the constant at pop n
G(1,1) = 1;
for i = 1:M
    for n = 1:(N-1)
        G(n+1,:) = G(n+1,:) + local_series_mul(Y(i,:), G(n,:), K);
    end
end

L = G(N,:) / Gn1;
for j = z
    e = (-1./mu(j)) .^ (0:K);
    L = local_series_mul(L, e, K);
end

mom = zeros(1, nmom);
for q = 1:nmom
    mom(q) = ((-1)^q) * factorial(q) * L(q+1);
end
end

function c = local_series_mul(a, b, K)
c = zeros(1, K+1);
for i = 0:K
    if a(i+1) == 0
        continue
    end
    for j = 0:(K-i)
        c(i+j+1) = c(i+j+1) + a(i+1) * b(j+1);
    end
end
end
