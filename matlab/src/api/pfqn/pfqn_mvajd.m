function [X, Qjd, Qli, Qdelay, Sjd] = pfqn_mvajd(Z, N, mu, Dli, visits, options)
% [X, QJD, QLI, QDELAY, SJD] = PFQN_MVAJD(Z, N, MU, DLI, VISITS, OPTIONS)
%
% Joint-dependent name of PFQN_MVAOI: the mean-value analysis of a closed
% network whose station rates read the whole per-class occupancy vector.
%
% The two names denote the SAME routine because the recursion evaluates the
% rate handle at a full occupancy vector, mu_i(s_i + e_r) with s_i the shift
% (the occupancy already committed at the bottom of station i), and never
% inspects the structure of mu_i. This is precisely the "third form" of the
% Conditional MVA of Casale, "A Note on Stable Flow-Equivalent Aggregation in
% Closed Networks" (QUESTA 2009): a rate depending on the full per-class
% occupancy vector. Order independence (mu_i constant on each support) is a
% modelling restriction, not an algorithmic one, so any joint-dependent scaling
% eta_i(n) (sn.jdscaling) is admissible.
%
% Unlike the AMVA joint-dependence route (solver_amvald with sn.jdscaling,
% which evaluates eta at the MEAN arrival-instant vector 1 + E[Q] and therefore
% collapses a support indicator to 1), this routine evaluates the rate at exact
% integer occupancies and is exact for the balanced-fair station, at the cost of
% walking prod_r C(N_r+K+1,K+1) states with K joint-dependent stations.
%
% Arguments and returns are exactly those of PFQN_MVAOI.
%
% See also PFQN_MVAOI, PFQN_NCJD, PFQN_CLWJD, PFQN_MVAOI_MARG.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 6
    options = struct();
end
if nargin < 5
    visits = [];
end
if nargin < 4
    Dli = [];
end
if nargin < 3
    mu = {};
end
[X, Qjd, Qli, Qdelay, Sjd] = pfqn_mvaoi(Z, N, mu, Dli, visits, options);
end
