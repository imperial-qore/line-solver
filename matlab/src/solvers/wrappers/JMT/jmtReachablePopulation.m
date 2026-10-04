function n = jmtReachablePopulation(sn, ist)
% N = JMTREACHABLEPOPULATION(SN, IST)
%
% The most jobs that can be present at station IST, read the way
% REFRESHCAPACITY derives the capacity itself: per CHAIN, because a chain's
% whole population can reach a station that serves any one of its classes
% (class switching moves jobs between them), and a chain none of whose classes
% is served there cannot put a single job on it.
%
% Deliberately NOT read off sn.classcap, which refreshCapacity has already
% clamped by the station's own cap: comparing a capacity against a quantity
% derived from it would make every user-declared buffer look non-binding.
% Inf when an open chain is served here, which is what sum(sn.njobs) gave
% before and which sends the station to JMTSTATIONCAPREFUSAL, where the open
% classes are skipped by name.
%
% Lifted out of @JMTIO/saveBufferCapacity.m so that JMTMETHODREFUSAL can apply
% the same binding test the writer applies: sn.cap is DERIVED for a station the
% user never capped, so "cap is finite" is not the question -- "cap is below
% what can reach the station" is.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

n = 0;
for c = 1:sn.nchains
    inchain = sn.inchain{c};
    if any(~isnan(sn.rates(ist, inchain)))
        n = n + sum(sn.njobs(inchain));
    end
end
end
