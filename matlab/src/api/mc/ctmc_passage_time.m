function [F, f, out] = ctmc_passage_time(Q, pi0, target, tset, options)
% [F, f, OUT] = CTMC_PASSAGE_TIME(Q, PI0, TARGET, TSET, OPTIONS)
%
% Cumulative distribution F and density f of the first passage time from the
% initial law PI0 into the target state set, on the grid TSET.
%
%     F(t) = 1 - alpha exp(St) 1        f(t) = alpha exp(St) s0
%
% OPTIONS.method selects the route:
%   'expm' (default)  exact, one matrix exponential reused along a uniform
%                     grid, in the manner of the CTMC response-time getters
%   'lt'              the transform of Eqs. 1-2 inverted through api/lti;
%                     OPTIONS.lti_method picks the inverter ('euler' default)
%
% 'lt' exists for chains whose non-target block is too large for a dense
% exp(St), not because it needs fewer time points: see CTMC_PASSAGE_LST. On a
% small chain 'expm' is both faster and more accurate, which is why it is the
% default.
%
% OUT carries .atom (the mass of PI0 already inside the target, i.e. F(0)),
% .alpha, .S, .s0 and .keep.
%
% Reference: P. G. Harrison and W. J. Knottenbelt, "Passage Time Distributions
% in Large Markov Chains", 2002.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 5 || isempty(options)
    options = struct();
end
if ~isfield(options,'method') || isempty(options.method)
    options.method = 'expm';
end
if ~isfield(options,'lti_method') || isempty(options.lti_method)
    options.lti_method = 'euler';
end

[alpha, S, s0, keep, atom] = ctmc_passage_ph(Q, pi0, target);
out = struct('atom', atom, 'alpha', alpha, 'S', S, 's0', s0, 'keep', keep, ...
    'method', options.method);

tset = reshape(tset, 1, []);
F = zeros(size(tset));
f = zeros(size(tset));

switch lower(options.method)
    case 'expm'
        nA = size(S,1);
        e = ones(nA,1);
        Sf = full(S);
        dt = diff(tset);
        uniform = numel(tset) > 2 && all(abs(dt - dt(1)) < 1e-12 * max(1,abs(dt(1)))) && dt(1) > 0;
        if uniform
            % One exponential, then propagate: the reference recomputes
            % expm(S*t) at every grid point, which is the same answer at a
            % cost linear in the grid rather than constant.
            E = expm(Sf * dt(1));
            v = alpha * expm(Sf * tset(1));
            for i = 1:numel(tset)
                if i > 1
                    v = v * E;
                end
                F(i) = 1 - v * e;
                f(i) = v * s0;
            end
        else
            for i = 1:numel(tset)
                if tset(i) < 0
                    F(i) = 0;
                    continue
                end
                v = alpha * expm(Sf * tset(i));
                F(i) = 1 - v * e;
                f(i) = v * s0;
            end
        end
    case 'lt'
        Lfun = @(s) local_lst(alpha, S, s0, atom, s);
        F = laplace_invert_cdf(Lfun, tset, options.lti_method);
        f = laplace_invert_pdf(@(s) Lfun(s) - atom, tset, options.lti_method);
    otherwise
        line_error(mfilename, sprintf('Unknown passage-time method: %s. Supported: expm, lt.', options.method));
end

F = min(max(F, 0), 1);
f = max(f, 0);
end

function L = local_lst(alpha, S, s0, atom, s)
I = speye(size(S,1));
L = alpha * ((s*I - S) \ s0) + atom;
end
