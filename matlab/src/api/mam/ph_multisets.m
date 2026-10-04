function M = ph_multisets(p, k)
% PH_MULTISETS Configurations of k identical servers over p service phases.
%
% M = PH_MULTISETS(P, K) returns the compositions of K into P nonnegative
% parts, one per row: M(r,i) is the number of the K busy servers sitting in
% phase i. There are NCHOOSEK(K+P-1, P-1) rows, the multiset count of Asmussen
% and Moller (2001): identical servers are exchangeable, so only the phase
% COUNTS carry information and the ordered space of size P^K collapses onto
% this one.
%
% The order is fixed and shared by every caller, so a configuration index means
% the same thing in each of them: the first part descends. K = 1 therefore
% yields the identity rows e_1 ... e_p in phase order, which is what makes the
% c = 1 case of LDQBD_MPHC coincide with the plain phase indexing.
%
% Examples:
%   ph_multisets(2, 0) -> [0 0]
%   ph_multisets(3, 1) -> [1 0 0; 0 1 0; 0 0 1]
%   ph_multisets(2, 2) -> [2 0; 1 1; 0 2]
%
% References:
%   S. Asmussen and J.R. Moller, "Calculation of the steady state waiting time
%   distribution in GI/PH/c and MAP/PH/c queues", Queueing Systems 37(1):9-29,
%   2001.
%
% See also ldqbd_mphc, qsys_mapphc
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if k == 0
    M = zeros(1, p);
    return;
end
if p == 1
    M = k;
    return;
end
M = [];
for first = k:-1:0
    rest = ph_multisets(p - 1, k - first);
    M = [M; [repmat(first, size(rest, 1), 1), rest]]; %#ok<AGROW>
end
end
