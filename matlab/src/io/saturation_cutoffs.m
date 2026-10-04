function [tau, bounded] = saturation_cutoffs(fAt, K, cmax)
% [TAU, BOUNDED] = SATURATION_CUTOFFS(FAT, K, CMAX)
% Per-class saturation cutoffs of a count-indexed function f(n), the search
% behind the OI/PAS and class/joint-dependence table extents of linemodel_save.
%
% FAT(n) returns the row of values the writer would tabulate at the 0-based
% count vector n (1xK), ALREADY SANITIZED (non-finite -> 0), so the search
% compares what the table holds. Starting from tau = 1, a class whose extra job
% still changes f anywhere on the current box has its cutoff raised, until no
% class moves. Only pairs (m, m + e_r) with sum(m) + 1 <= CMAX are compared: a
% state above the buffer is never visited, so whatever a reloaded table answers
% there cannot move the search, which is what makes save -> load -> save a fixed
% point. Passing 256, or a cutoff reaching CMAX (the box would be no smaller than
% the buffer), gives BOUNDED = false. CMAX defaults to Inf.
% Twins: _saturation_cutoffs (Python linemodel_io), io::saturation_cutoffs (C++),
% LineModelIO.saturationCutoffs (JAR).
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 3 || isempty(cmax)
    cmax = Inf;
end
maxCutoff = 256;
tol = 1e-9;
tau = ones(1, K);
changed = true;
while changed
    changed = false;
    for r = 1:K
        lim = tau + 1;      % box sizes per class, class r pinned at tau(r)
        lim(r) = 1;
        total = prod(lim);
        sat = true;
        for li = 0:(total - 1)
            m = zeros(1, K);
            rem = li;
            for j = 1:K
                if j == r
                    m(j) = tau(r);
                else
                    m(j) = mod(rem, lim(j));
                    rem = floor(rem / lim(j));
                end
            end
            if sum(m) + 1 > cmax
                continue
            end
            mp = m;
            mp(r) = mp(r) + 1;
            if max(abs(fAt(mp) - fAt(m))) > tol
                sat = false;
                break
            end
        end
        if ~sat
            tau(r) = tau(r) + 1;
            changed = true;
            if tau(r) > maxCutoff
                bounded = false;
                return
            end
        end
    end
end
bounded = all(tau < cmax);
end
