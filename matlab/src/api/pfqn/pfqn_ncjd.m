function [G, lG, Gtab] = pfqn_ncjd(Z, N, mu, visits, options)
% [G, LG, GTAB] = PFQN_NCJD(Z, N, MU, VISITS, OPTIONS)
%
% Joint-dependent name of PFQN_NCOI: the balance-function convolution of a
% closed network whose station rates read the whole per-class occupancy vector.
%
% The two names denote the SAME routine because the balanced-fairness recursion
%   Phi_i(0) = 1,   mu_i(n) Phi_i(n) = sum_{r: n_r>0} v_{i,r} Phi_i(n - e_r)
% never inspects the structure of mu_i: it evaluates the handle at the full
% count vector n. Order independence (mu_i constant on each support) is a
% modelling restriction that buys insensitivity and a physical reading of Phi,
% not something the convolution uses. Any joint-dependent scaling eta_i(n)
% (sn.jdscaling) is therefore admissible here, with the usual proviso that the
% product form it induces is the balanced-fair one matched to that rate: Phi
% must stay positive for the result to be a distribution, which the routine
% enforces by zeroing the balance value of a non-positive rate.
%
% Use PFQN_NCOI when the model is genuinely order independent and the name
% should say so; use PFQN_NCJD when the rate is a general joint dependence. See
% PFQN_CLWJD for the transform route, which needs the rate to saturate at a
% finite cutoff and is not merely a renaming.
%
% Arguments and returns are exactly those of PFQN_NCOI.
%
% See also PFQN_NCOI, PFQN_MVAJD, PFQN_CLWJD, PFQN_PAS_NC.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 5
    options = struct();
end
if nargin < 4
    visits = [];
end
if nargin < 3
    mu = {};
end
[G, lG, Gtab] = pfqn_ncoi(Z, N, mu, visits, options);
end
