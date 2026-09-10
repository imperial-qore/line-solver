function tf = sn_has_classdep_routing(sn)
% SN_HAS_CLASSDEP_ROUTING True when classes are not routed alike.
%
% TF = SN_HAS_CLASSDEP_ROUTING(SN) returns true when the routing probabilities
% differ between job classes, either because a class switches class on a hop or
% because two classes leave the same station with different probabilities. It
% is false for a single-class model and for a multiclass model in which every
% class traverses the network identically.
%
% This is the condition under which per-class visit ratios diverge, so a method
% that aggregates classes into a per-chain demand vector stops being exact.
% Used to gate Marie's aggregation-decomposition in SolverMVA: with all classes
% routed alike that method reproduces the exact solution, and it degrades as the
% per-class demand vectors separate.
%
% sn.rt is indexed station-major, (i-1)*K+r, matching sn_refresh_visits.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

tf = false;
K = sn.nclasses;
M = sn.nstations;
if K <= 1
    return;
end

tol = GlobalConstants.FineTol;
for i = 1:M
    for j = 1:M
        shared = [];
        for r = 1:K
            % Class switching: leaving station i as class r and arriving at j
            % as some other class s makes the routing class-dependent outright.
            for s = 1:K
                if r ~= s && sn.rt((i-1)*K+r, (j-1)*K+s) > tol
                    tf = true;
                    return;
                end
            end
            p = sn.rt((i-1)*K+r, (j-1)*K+r);
            if isempty(shared)
                shared = p;
            elseif abs(p - shared) > tol
                tf = true;
                return;
            end
        end
    end
end
end
