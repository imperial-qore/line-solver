function [p, logp, theta, info] = ctmc_saddlepoint(D0, D1, t, k, method, pi0)
% [P, LOGP, THETA, INFO] = CTMC_SADDLEPOINT(D0, D1, T, K, METHOD, PI0)
%
% Saddlepoint approximation of Pr{N(t)=k}, the probability that the counting
% process of the Markovian arrival process (D0,D1) records exactly k events
% in (0,t], starting from the phase distribution PI0. Throughout, K is the
% number of phases and k the event count.
%
% The counting generating function of a MAP is the matrix exponential
%
%   sum_{k>=0} P(k,t) z^k = exp(t*(D0 + z*D1)),
%
% where P(k,t) is the (i,j) probability of k events and phase j at time t. Its
% cumulant generating function is therefore
%
%   eta(theta) = spectral abscissa of A(theta) = D0 + exp(theta)*D1,
%
% the Perron root of an irreducible Metzler matrix, which is real, simple and
% strictly convex in theta with eta(0)=0 and eta'(0)=lambda, the arrival rate.
% Inverting the generating function by the method of steepest descent gives
% Daniels (1954),
%
%   Pr{N(t)=k} ~ g(theta*) * exp(t*eta(theta*) - k*theta*) / sqrt(2*pi*t*eta''(theta*)),
%
% where the saddle theta* solves eta'(theta*) = k/t and g is the amplitude of
% the Perron projection, g(theta) = (pi0*v(theta))*(u(theta)*1) with u,v the
% left and right Perron vectors normalised by u*v = 1.
%
% THE EXPANSION PARAMETER IS K2 = t*eta''(theta*), THE VARIANCE OF THE COUNT,
% not its mean and not t. Measured across the process family, the two error
% laws are
%
%   err('daniels')  = 0.083 / K2,      err('daniels2') = 0.017 / K2^2,
%
% with the constants flat to two digits over Erlang orders 1..8 and horizons
% 10..160. For a renewal Erlang(r) the count variance rate is lambda/r, so
% K2 = lambda*t/r and an Erlang-4 at t=50 is as accurate as a Poisson at
% t=12.5: low variability SHRINKS the parameter, it does not break the method.
% Below K2 = 5 the expansion is outside its regime and the call warns; INFO.K2
% carries the value per point.
%
% The approximation therefore IMPROVES with the horizon, which is where the
% block-by-block expansions of exp(t*X) lose their digits to cancellation, and
% it keeps relative accuracy on probabilities far below the double-precision
% floor when read through LOGP. Its cost is one eigenproblem of order K per
% Newton step and is independent of both t and k.
%
% This is an asymptotic method, not a quadrature: use it for rare-event and
% large-deviation coefficients, where k/t is away from lambda or where the
% probability underflows. For the bulk of the transient distribution, i.e.
% every block k=0..N-1 at once at moderate t, uniformization
% (CTMC_UNIFORMIZATION, CTMC_FOXGLYNN) is both exact and faster.
%
% @param D0 Generator of the phase process with the counted transitions
%        removed (K x K). May instead be a MAP cell {D0,D1}, in which case the
%        remaining arguments shift left by one.
% @param D1 Rates of the counted transitions (K x K, nonnegative). D0+D1 must
%        be an irreducible generator.
% @param t Time horizon; scalar, or an array broadcast against the count
% @param k Event count, a nonnegative integer; scalar, or an array broadcast
%        against the horizon
% @param method 'daniels2' (default) second-order saddlepoint, the Daniels
%        bracket plus the amplitude's own curvature along the contour, error
%        O(1/K2^2); 'daniels' is the first-order form with the Perron
%        amplitude, error O(1/K2); 'plain' is the bare first-order form with
%        the amplitude set to 1
% @param pi0 Initial phase distribution (1 x K); the stationary distribution
%        of D0+D1 if empty or omitted
% @return p Approximation of Pr{N(t)=k}, of the broadcast size of t and k
% @return logp Its natural logarithm, evaluated without forming P, so that it
%         stays accurate below the smallest positive double
% @return theta The saddle theta*, -Inf where k=0, the one point at which it
%         runs off to minus infinity
% @return info Struct with the per-point fields eta, deta, d2eta, d3eta,
%         d4eta, ampl (the Perron amplitude g), corr (the bracket multiplying
%         the leading term), k2 (t*eta''(theta*), the expansion parameter --
%         read this to know how far the answer can be trusted), iter (Newton
%         steps), exact (true where the value
%         was computed exactly rather than approximated), and the scalar
%         field lambda (the stationary event rate eta'(0))
%
% Examples:
%   MAP = {[-1 0.2 0.1; 0.05 -2 0.3; 0.1 0.1 -0.9], ...
%          [0.4 0.2 0.1; 0.5 0.15 1.0; 0.3 0.2 0.2]};
%   p = ctmc_saddlepoint(MAP, 200, 66)                   % one coefficient
%   [p, logp] = ctmc_saddlepoint(MAP, 1e6, 330298);      % logp = -1.95e5
%   p = ctmc_saddlepoint(MAP, 100, 120:140)              % a stretch of counts
%   [p, ~, ~, info] = ctmc_saddlepoint(MAP, 40, 34);     % info.k2 = 36.06
%
% Measured relative error on that MAP, against Pr{N(t)=k} read off the block
% chain by dense expm (independently cross-checked against uniformization
% summed from n=0):
%
%      t     k      K2    Pr{N(t)=k}    daniels    daniels2   plain
%     20    26    27.2    1.109e-02     3.1e-03    6.0e-06    5.0e-03
%    100   132   137.8    9.749e-07     6.1e-04    2.4e-07    2.7e-03
%    200    66    71.4    5.354e-19     1.1e-03    6.2e-07    1.1e-02
%    500   165   178.6    1.316e-44     4.6e-04    1.0e-07    1.1e-02
%    500   660   688.8    2.984e-25     1.2e-04    9.6e-09    2.2e-03
%
% The three columns are the three regimes of the expansion: 'plain' carries an
% amplitude error that does not vanish, 'daniels' is O(1/K2), 'daniels2' is
% O(1/K2^2). Over the point processes and horizons of LINE's own transient
% examples (MAP, MMPP2, Erlang, APH and Exp sources at t = 5..1000; 25 cases,
% 222 sampled counts at k/lambda*t in [0.15, 2.5]), 'daniels2' was never worse
% than 'daniels' at any point, and better by one to three orders almost
% everywhere. Where both fail is small K2, and there they fail together: the
% second-order correction goes nonpositive, the call falls back to first order
% and warns.
%
%
% ATTRIBUTION. The first-order form is Daniels (1954). The AMPLITUDE g and the
% whole 'daniels2' bracket are NOT a rederivation: they are Jensen, "Saddlepoint
% Expansions for Sums of Markov Dependent Variables on a Continuous State
% Space", Probab. Th. Rel. Fields 89, 1991, Eq. (4.4) with the coefficients on
% p.191. His gamma_0(s) = (sum_i c_i)(sum_i r_i P(Y_0=i)) is exactly g, under
% his own normalisation sum_i r_i c_i = 1; expanding his
% alpha_0 + (1/n){-alpha_3/2 + alpha_4/8 - 5*alpha_5/24} reproduces
% g*(1 + lam4/8 - 5*lam3^2/24) - g''/(2*K2) + g'*K3/(2*K2^2) term for term, and
% his Theorem 4.1 gives the O(n^-2) error this file measures as 0.017/K2^2.
% Jensen works with discrete-n sums over a Markov chain; the continuous-time MAP
% counting process here is that result transcribed, n -> t and the kernel
% eigenvalue -> the Perron root of D0+exp(theta)*D1. He also states the reason
% the amplitude cannot be dropped (p.183): the large-deviation limit discards
% any factor of order one, but a local approximation must carry the projection
% onto the eigenspace, not just the maximal eigenvalue.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

nargs = nargin;
if iscell(D0)
    % MAP cell form {D0,D1}: every later argument sits one position to the
    % left, so the argument count is restated in terms of the split form
    if nargin < 3
        line_error(mfilename, 'The MAP cell form still requires the horizon t and the count k.');
    end
    if nargin > 5
        line_error(mfilename, 'Too many arguments: the MAP cell form takes at most {D0,D1}, t, k, method, pi0.');
    end
    if nargin >= 5, pi0 = method; else, pi0 = []; end
    if nargin >= 4, method = k; else, method = ''; end
    k = t;
    t = D1;
    D1 = D0{2};
    D0 = D0{1};
    nargs = nargin + 1;
end
if nargs < 4
    line_error(mfilename, 'Both the horizon t and the count k are required.');
end
if nargs < 5 || isempty(method)
    method = 'daniels2';
end
if nargs < 6
    pi0 = [];
end

D0 = full(D0);
D1 = full(D1);
nph = size(D0, 1);
if size(D0, 2) ~= nph || any(size(D1) ~= [nph, nph])
    line_error(mfilename, 'D0 and D1 must be square matrices of the same order.');
end
if any(D1(:) < 0)
    line_error(mfilename, 'D1 must be nonnegative.');
end
Q = D0 + D1;
if max(abs(sum(Q, 2))) > 1e-8 * max(1, max(abs(Q(:))))
    line_error(mfilename, 'D0+D1 must be an infinitesimal generator (zero row sums).');
end
if max(abs(D1(:))) <= 0
    line_error(mfilename, 'D1 has no counted transitions, the counting process is identically zero.');
end

method = lower(strtrim(char(method)));
switch method
    case {'daniels', 'sp1'}
        order = 1; useampl = true;
    case {'daniels2', 'sp2'}
        order = 2; useampl = true;
    case {'plain', 'bare'}
        order = 1; useampl = false;
    otherwise
        line_error(mfilename, sprintf('Unknown method ''%s'', expected daniels, daniels2 or plain.', method));
end

if isempty(pi0)
    pi0 = ctmc_solve(Q);
end
pi0 = reshape(full(pi0), 1, []);
if numel(pi0) ~= nph
    line_error(mfilename, 'pi0 must have one entry per phase.');
end
if abs(sum(pi0) - 1) > 1e-8
    line_error(mfilename, 'pi0 must sum to one.');
end

% Broadcast the horizons against the counts
t = double(t);
k = double(k);
if isscalar(t) && ~isscalar(k)
    t = repmat(t, size(k));
elseif isscalar(k) && ~isscalar(t)
    k = repmat(k, size(t));
elseif ~isequal(size(t), size(k))
    line_error(mfilename, 't and k must be scalars or arrays of the same size.');
end
if any(t(:) < 0)
    line_error(mfilename, 'The horizon t must be nonnegative.');
end
if any(k(:) < 0) || any(k(:) ~= round(k(:)))
    line_error(mfilename, 'The count k must be a nonnegative integer.');
end

n = numel(t);
p = zeros(size(t));
logp = -inf(size(t));
theta = -inf(size(t));
info = struct();
info.eta = nan(size(t));
info.deta = nan(size(t));
info.d2eta = nan(size(t));
info.d3eta = nan(size(t));
info.d4eta = nan(size(t));
info.ampl = nan(size(t));
info.corr = nan(size(t));
info.k2 = nan(size(t));
info.iter = zeros(size(t));
info.exact = false(size(t));

% Below this value of K2 = t*eta''(theta*) the expansion is out of its regime.
% K2 is the VARIANCE of the count, not its mean: the measured error laws are
% err(daniels) = 0.083/K2 and err(daniels2) = 0.017/K2^2, uniformly over the
% process family, so K2 = 5 is where the second order stops holding 1e-3. Do
% NOT threshold on lambda*t: for Erlang(r) the count variance rate is lambda/r,
% so K2 = lambda*t/r, and lambda*t over-warns on Poisson-like processes while
% under-warning on low-variability ones.
K2_MIN = 5;
worstK2 = inf;
worstAt = [0, 0];

onesvec = ones(nph, 1);
maxrate = max(abs(D1(:)));
% exp(theta) multiplies D1, so the saddle is confined to the range over which
% A(theta) is representable; the bounds are never active for a feasible k/t
THMAX = log(realmax / 1e6) - log(maxrate);
THMIN = log(realmin * 1e6) - log(maxrate);

s0 = perronstate(0);
info.lambda = s0.deta;

% Sorting by the rate k/t lets each Newton solve warm-start from the previous
% saddle, the saddle being a monotone function of that rate alone
rate = zeros(n, 1);
for i = 1:n
    if t(i) > 0
        rate(i) = k(i) / t(i);
    end
end
[~, ord] = sort(rate);
thprev = 0;

for idx = 1:n
    i = ord(idx);
    ti = t(i);
    ki = k(i);
    if ti == 0
        % No time has elapsed, so the count is zero with probability one
        info.exact(i) = true;
        if ki == 0
            p(i) = 1;
            logp(i) = 0;
        end
        continue
    end
    if ki == 0
        % The saddle runs off to -Inf; the exact value is one matrix
        % exponential of the taboo generator and costs no more than a step of
        % the approximation itself
        info.exact(i) = true;
        p(i) = pi0 * expm(ti * D0) * onesvec;
        logp(i) = log(p(i));
        continue
    end

    [th, iters] = solvesaddle(ki / ti, thprev);
    thprev = th;
    theta(i) = th;
    info.iter(i) = iters;

    s = perronstate(th);
    info.eta(i) = s.eta;
    info.deta(i) = s.deta;
    info.d2eta(i) = s.d2eta;

    K2 = ti * s.d2eta;
    info.k2(i) = K2;
    if K2 < worstK2
        worstK2 = K2;
        worstAt = [ti, ki];
    end
    if ~(K2 > 0)
        line_error(mfilename, sprintf(['The cumulant generating function is not strictly convex at the saddle ' ...
            '(t=%g, k=%d): eta''''=%g. D0+D1 is probably reducible.'], ti, ki, s.d2eta));
    end
    base = ti * s.eta - ki * th - 0.5 * log(2 * pi * K2);

    if useampl
        ampl = s.ampl;
    else
        ampl = 1;
    end
    info.ampl(i) = s.ampl;

    if order == 1
        corr = ampl;
    else
        % The higher cumulants and the derivatives of the amplitude come from
        % central differences of the analytic eta'' and g, both of which carry
        % full precision at each evaluation point
        h = 1e-3 * max(1, abs(th));
        sp = perronstate(th + h);
        sm = perronstate(th - h);
        d3 = (sp.d2eta - sm.d2eta) / (2 * h);
        d4 = (sp.d2eta - 2 * s.d2eta + sm.d2eta) / h^2;
        info.d3eta(i) = d3;
        info.d4eta(i) = d4;
        K3 = ti * d3;
        K4 = ti * d4;
        lam3sq = K3^2 / K2^3;
        lam4 = K4 / K2^2;
        if useampl
            gp = (sp.ampl - sm.ampl) / (2 * h);
            gpp = (sp.ampl - 2 * s.ampl + sm.ampl) / h^2;
        else
            gp = 0;
            gpp = 0;
        end
        % Steepest descent to O(1/K2), Jensen (1991) Eq. (4.4): the Daniels
        % bracket on the amplitude, plus the two terms the amplitude
        % contributes through its own curvature along the contour
        corr = ampl * (1 + lam4 / 8 - 5 * lam3sq / 24) - gpp / (2 * K2) + gp * K3 / (2 * K2^2);
        if corr <= 0
            line_warning(mfilename, ...
                'The second-order correction is nonpositive at t=%g, k=%d; the expansion has broken down, returning the first-order value.\n', ...
                ti, ki);
            corr = ampl;
        end
    end
    info.corr(i) = corr;

    logp(i) = base + log(corr);
    p(i) = exp(logp(i));
end

if worstK2 < K2_MIN
    % Once per call, not once per point: a vectorised call spans hundreds of
    % counts and the caller needs the worst one, not a page of repetitions
    line_warning(mfilename, ...
        ['K2 = t*eta''''(theta*) = %.2f at t=%g, k=%d is below %g, so the saddlepoint expansion is ' ...
        'outside its asymptotic regime there and the result is unreliable (expect a relative error ' ...
        'near %.0e). K2 is the variance of the count, not its mean: a low-variability process needs ' ...
        'a longer horizon than its rate suggests. Take the exact value from the block chain instead.\n'], ...
        worstK2, worstAt(1), worstAt(2), K2_MIN, 0.017 / worstK2^2);
end

    function st = perronstate(th)
        % ST = PERRONSTATE(TH)
        %
        % Perron root of A(th)=D0+exp(th)*D1 with its first two derivatives in
        % th and the amplitude of the Perron projection between pi0 and 1.
        %
        % Only EIGENVALUES are taken from the eigensolver; the Perron vectors
        % come from bordered solves, the idiom CTMC_SOLVE already uses. That
        % keeps all four codebases on ONE algorithm: neither the JAR
        % (commons-math hands back Schur blocks, not eigenvectors, as soon as a
        % complex pair appears) nor C++ (eig.h exposes values only) can supply a
        % left eigenvector. Agreement with the inv(V) form is 5e-15.
        W = exp(th) * D1;               % A'(th) = A''(th) = exp(th)*D1
        A = D0 + W;
        st.eta = max(real(eig(A)));
        Ashift = A - st.eta * eye(nph);
        rhs = zeros(nph, 1);
        rhs(nph) = 1;
        % (A-eta*I)v = 0 with the last row replaced by sum(v)=1. A row may be
        % dropped because A-eta*I is a singular irreducible M-matrix, every
        % proper principal submatrix of which is nonsingular
        M = Ashift;
        M(nph, :) = 1;
        v = M \ rhs;
        % u(A-eta*I) = 0 by the same construction on the transpose
        Mt = Ashift.';
        Mt(nph, :) = 1;
        u = (Mt \ rhs).';
        u = u / (u * v);                % u*v = 1 fixes the residual scale freedom
        st.deta = u * W * v;
        % First-order eigenvector perturbation (A-eta*I)v' = (eta'*I-W)v, taken
        % with u*v'=0; the bordered system is nonsingular because the Perron
        % root of an irreducible Metzler matrix is simple
        M2 = [Ashift, v; u, 0];
        sol = M2 \ [(st.deta * eye(nph) - W) * v; 0];
        vp = sol(1:nph);
        st.d2eta = st.deta + 2 * (u * W * vp);
        st.ampl = (pi0 * v) * (u * onesvec);
    end

    function [th, iters] = solvesaddle(r, th0)
        % [TH, ITERS] = SOLVESADDLE(R, TH0)
        %
        % Saddle of the counting cumulant generating function at rate R, i.e.
        % the root of eta'(th)=R. eta' is continuous and strictly increasing
        % from 0 to +Inf, so the root exists and is unique for every R>0; it is
        % bracketed by geometric expansion from TH0 and then refined by Newton
        % on log(eta'), safeguarded by bisection.
        TOL = 1e-13;
        MAXIT = 200;
        th = min(max(th0, THMIN), THMAX);
        d1 = getderiv(th);
        lo = th; dlo = d1;
        hi = th; dhi = d1;
        step = 1;
        while dlo > r
            hi = lo; dhi = dlo;
            lo = lo - step;
            if lo <= THMIN
                lo = THMIN;
                dlo = getderiv(lo);
                if dlo > r
                    line_error(mfilename, sprintf('The rate k/t=%g is below the representable range of eta''.', r));
                end
                break
            end
            dlo = getderiv(lo);
            step = 2 * step;
        end
        step = 1;
        while dhi < r
            lo = hi; dlo = dhi;
            hi = hi + step;
            if hi >= THMAX
                hi = THMAX;
                dhi = getderiv(hi);
                if dhi < r
                    line_error(mfilename, sprintf('The rate k/t=%g is above the representable range of eta''.', r));
                end
                break
            end
            dhi = getderiv(hi);
            step = 2 * step;
        end
        th = min(max(th, lo), hi);
        logr = log(r);
        iters = 0;
        for it = 1:MAXIT
            iters = it;
            si = perronstate(th);
            f = log(si.deta) - logr;
            if abs(f) <= TOL
                break
            end
            if f > 0
                hi = th;
            else
                lo = th;
            end
            thn = th - f * si.deta / si.d2eta;
            if ~isfinite(thn) || thn <= lo || thn >= hi
                thn = 0.5 * (lo + hi);
            end
            if abs(thn - th) <= TOL * max(1, abs(th))
                th = thn;
                break
            end
            th = thn;
        end
    end

    function d1 = getderiv(th)
        % D1 = GETDERIV(TH)
        %
        % eta'(TH) alone, used while the saddle is being bracketed. The local
        % name differs from the one SOLVESADDLE uses because two sibling
        % nested functions sharing a name share the variable itself.
        sd = perronstate(th);
        d1 = sd.deta;
    end

end
