function [F, f, out] = smp_passage_time(P, hlst, pi0, target, tset, options)
% [F, f, OUT] = SMP_PASSAGE_TIME(P, HLST, PI0, TARGET, TSET, OPTIONS)
%
% Cumulative distribution and density of the first passage time into the target
% state set for a semi-Markov chain, by inverting SMP_PASSAGE_LST through
% api/lti.
%
% There is no matrix-exponential route here: a semi-Markov chain has no
% generator to exponentiate, which is exactly the case uniformization does not
% reach and the transform does. This is the point the paper makes for
% preferring transform inversion over uniformization.
%
% OPTIONS.lti_method defaults to 'euler' RATHER THAN 'weeks'. Semi-Markov
% passage densities are the case Sec. 4.2 singles out as slow-converging for a
% Laguerre series: a kernel with a deterministic or discontinuous holding time
% gives a density whose derivatives jump, and LAPLACE_WEEKS_SCALING then
% refuses by name rather than returning noise.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 6 || isempty(options)
    options = struct();
end
if ~isfield(options,'lti_method') || isempty(options.lti_method)
    options.lti_method = 'euler';
end

n = size(P,1);
target = unique(reshape(target,1,[]));
if isempty(pi0)
    pi0 = ones(1,n)/n;
end
pi0 = reshape(pi0,1,[]);
atom = sum(pi0(target));

Lfun = @(sv) smp_passage_lst(P, hlst, pi0, target, sv);
tset = reshape(tset,1,[]);
F = laplace_invert_cdf(Lfun, tset, options.lti_method);
f = laplace_invert_pdf(@(sv) Lfun(sv) - atom, tset, options.lti_method);
out = struct('atom', atom, 'lti_method', options.lti_method);
end
