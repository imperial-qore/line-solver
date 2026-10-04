function grp = mam_bgchain_groups(D, G)
% GRP = MAM_BGCHAIN_GROUPS(D, G)
%
% Groups the columns of D into G clusters by similarity of their per-station
% SERVICE DEMAND, for the aggregation SOLVER_MAM_BGCHAIN applies to the closed
% chains it is not carrying exactly.
%
% Inputs
%   D  (Mc x n)  one column per chain; D(i,c) is the demand chain c places on
%                station i (visits times mean service time, Lchain)
%   G  scalar    number of groups wanted, clamped to [1, n]
%
% Output
%   GRP (1 x n)  group index in 1..G of each chain
%
% WHY DEMAND IS THE RIGHT CRITERION. The aggregate that replaces a group carries
% the flow-weighted mean of its members' service times and routing, so the group
% aggregates EXACTLY when its members place the same demand at every station and
% distorts both quantities in proportion to how far apart they are. The distance
% is therefore the symmetric relative L1 gap between the demand vectors,
%
%   dist(a,b) = sum_i |D(i,a) - D(i,b)| / ((sum_i D(i,a) + sum_i D(i,b))/2),
%
% which is scale-relative rather than absolute: it separates two chains whose
% demand PROFILE across the stations differs and two chains whose profile agrees
% but whose magnitude does not, and being dimensionless it groups a model the
% same way whatever its time unit.
%
% WHY COMPLETE LINKAGE. The clustering is agglomerative from singletons, merging
% at each step the pair of clusters whose WORST member-to-member distance is
% smallest. The aggregation error inside a group is driven by its worst mismatch
% and not by its average one, so complete linkage is the criterion that bounds
% what the aggregation actually costs; average or single linkage would let one
% distant chain ride along inside an otherwise tight group.
%
% DETERMINISM. Ties are broken by the lexicographically smallest pair of cluster
% indices and the groups are relabelled by their smallest member, so the same D
% gives the same grouping in MATLAB, the JAR, Python and C++.
%
% See also SOLVER_MAM_BGCHAIN, MAM_BGCHAIN_CTMC.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

n = size(D, 2);
grp = zeros(1, n);
if n == 0
    return;
end
G = min(max(round(G), 1), n);

%% Pairwise symmetric relative L1 distance between the demand vectors
tot = sum(D, 1);
dist = zeros(n, n);
for a = 1:n
    for b = a+1:n
        den = (tot(a) + tot(b)) / 2;
        if den > GlobalConstants.Zero
            dist(a, b) = sum(abs(D(:, a) - D(:, b))) / den;
        end
        dist(b, a) = dist(a, b);
    end
end

%% Complete-linkage agglomeration down to G clusters
clusters = cell(1, n);
for a = 1:n
    clusters{a} = a;
end
active = true(1, n);
while sum(active) > G
    ids = find(active);
    best = Inf; bp = 0; bq = 0;
    for x = 1:numel(ids)-1
        for y = x+1:numel(ids)
            d = max(max(dist(clusters{ids(x)}, clusters{ids(y)})));
            if d < best - GlobalConstants.Zero
                best = d; bp = ids(x); bq = ids(y);
            end
        end
    end
    if bp == 0
        break;  % every remaining pair is at distance Inf/NaN: stop merging
    end
    clusters{bp} = sort([clusters{bp}, clusters{bq}]);
    clusters{bq} = [];
    active(bq) = false;
end

%% Relabel by smallest member, so the group numbering is canonical
ids = find(active);
firsts = zeros(1, numel(ids));
for x = 1:numel(ids)
    firsts(x) = clusters{ids(x)}(1);
end
[~, ord] = sort(firsts);
for g = 1:numel(ord)
    grp(clusters{ids(ord(g))}) = g;
end
end
