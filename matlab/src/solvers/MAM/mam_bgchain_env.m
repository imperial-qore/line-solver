function [A, phi, esup] = mam_bgchain_env(bg, i)
% [A, PHI, ESUP] = MAM_BGCHAIN_ENV(BG, I)
%
% Lumps the background modulating chain BG onto the number of closed jobs held
% by station I, giving the Markovian environment that the open classes at that
% station see.
%
% Station i does not observe the whole closed population vector, only how many
% closed jobs compete with the open ones for its server. The lumped generator
% is the stationary-weighted aggregation of BG.Q over the level sets
% {s : totocc(s,i) = e},
%
%   A(e,e') = sum_{s in e} pi(s) * sum_{s' in e'} Q(s,s') / sum_{s in e} pi(s),
%
% which is exact when the partition is lumpable in Kemeny-Snell's sense and is
% the standard exact-aggregation approximation otherwise. The diagonal is set
% from the off-diagonal row sums, so A is a proper generator whatever the
% lumping error is, and PHI (the aggregated stationary vector) is by
% construction the stationary vector of the aggregated chain when the partition
% is lumpable.
%
% Outputs
%   A     (me x me)  lumped environment generator
%   PHI   (1 x me)   stationary probability of each environment state
%   ESUP  (1 x me)   number of closed jobs each environment state stands for
%
% Environment states of zero stationary probability are unreachable and are
% dropped, so ESUP need not be 0:N.
%
% The aggregation being exact here is not luck: the background chain is
% product-form by construction, so this is Norton's flow-equivalent rather than
% an approximation. MAM_BGCHAIN_ENVFULL is the oracle that re-establishes it.
%
% See also MAM_BGCHAIN_CTMC, MAM_BGCHAIN_ENVFULL, SOLVER_MAM_BGCHAIN.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

occ = bg.totocc(:, i);
esup = unique(occ(:))';
me = numel(esup);

w = zeros(1, me);
for e = 1:me
    w(e) = sum(bg.pi(occ == esup(e)));
end
keep = w > GlobalConstants.Zero;
if ~any(keep)
    % degenerate chain: the station never holds a closed job
    A = 0;
    phi = 1;
    esup = 0;
    return;
end
esup = esup(keep);
w = w(keep);
me = numel(esup);

A = zeros(me, me);
if me > 1
    Qfull = bg.Q;
    for e = 1:me
        rows = find(occ == esup(e));
        if isempty(rows)
            continue;
        end
        flow = bg.pi(rows) * Qfull(rows, :);   % (1 x nstates)
        for ep = 1:me
            if ep == e
                continue;
            end
            A(e, ep) = sum(flow(occ == esup(ep)));
        end
        A(e, :) = A(e, :) / w(e);
    end
    A = A - diag(sum(A, 2));
end

phi = w / sum(w);
end
