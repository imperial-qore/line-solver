function visited = fluid_visited_pairs(sn, M, K)
% VISITED = FLUID_VISITED_PAIRS(SN, M, K)
%
% Mark the (station, class) pairs the model actually routes a job into, as an
% (M x K) logical read off the per-chain visit ratios SN.VISITS.
%
% A fluid result cannot decide that question from the SIZE of QN or TN. Both
% carry a decaying remnant of the initial state, which is spread over pairs the
% class never reaches, and the remnant is whatever the integrator left behind
% when it stopped: measured at QN = 1.3e-12 and TN = 1.3e-13 on picard05 for
% test_CQN_Cox_CS_7, i.e. ABOVE GlobalConstants.Zero, so a threshold on them
% divides one remnant by the other and reports the station's own service time
% as a response time. The visit ratios come from the routing solve instead,
% where an unrouted pair is zero to the last bits (2.7e-17 on that pair).
%
% SN.VISITS{c} is indexed by STATEFUL node, hence the SN.STATIONTOSTATEFUL
% lookup; see _kb/04-networkstruct.md. A struct carrying no visit information
% at all decides nothing and every pair is reported visited.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

visited = false(M, K);
haveVisits = false;
maxCols = 0;
for c = 1:numel(sn.visits)
    Vc = sn.visits{c};
    if isempty(Vc)
        continue
    end
    haveVisits = true;
    maxCols = max(maxCols, size(Vc,2));
    kk = 1:min(K, size(Vc,2));
    for ist = 1:M
        if ist > numel(sn.stationToStateful)
            visited(ist,:) = true;
            continue
        end
        isf = sn.stationToStateful(ist);
        if isf < 1 || isf > size(Vc,1)
            visited(ist,:) = true;
            continue
        end
        visited(ist,kk) = visited(ist,kk) | (abs(Vc(isf,kk)) > GlobalConstants.Zero);
    end
end
if ~haveVisits
    visited = true(M, K);
elseif maxCols < K
    % A class NO visit matrix reaches is not evidence of a non-visit, only of a
    % struct whose visits were refreshed against fewer classes. Leave it alone.
    visited(:, maxCols+1:K) = true;
end
end
