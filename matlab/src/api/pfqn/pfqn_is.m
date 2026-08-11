function [G, lG] = pfqn_is(L, N, Z, options)
% [G, LG] = PFQN_IS(L, N, Z, OPTIONS)
%
% Importance-sampling (IS) estimate of the normalizing constant of a closed
% LOAD-INDEPENDENT product-form queueing network with M single-server queues of
% per-class demand L and an aggregated delay of think time Z.
%
% This is the load-independent case of PFQN_LD_IS (capacities mu_i(k) = 1), and
% the ordinary-network counterpart of the order-independent PFQN_OI_IS and the
% pass-and-swap PFQN_PAS_IS: all four are the same sample-an-ordering estimator,
% differing only in the per-position factor of each station's balance function.
% For a single-server queue that factor is the demand of the class at that
% position, L(i,q_p); for the delay it is Z(q_p)/p; for an OI/P&S station it is
% the reciprocal rank rate 1/mu_i(supp(q_1..q_p)).
%
% Writing ell = sum(N), an ordering c of all ell jobs is drawn by placing a
% uniformly random present class at each step (probability p(c) = product of the
% reciprocal branching factors), and the sum over ALL ways of cutting c into
% contiguous per-station segments is computed exactly by dynamic programming:
%   G(N) = E_{C~p}[ S(C)/p(C) ],
%   S(c)  = sum_{cuts} prod_m prod_{p} L(m, seg_m(p)),
% which is unbiased for the exact constant of PFQN_NC.
%
% Parameters:
%   L       - (M x R) per-class service demands at the M single-server queues.
%   N       - (1 x R) closed population vector, finite.
%   Z       - (1 x R) aggregated think time (delay) demand; [] or zeros if none.
%   options - solver options (optional). Fields used:
%               .samples  number of IS samples (default 1e4);
%               .seed     RNG seed for reproducibility (optional).
%
% Returns:
%   G  - IS estimate of the normalizing constant G(N).
%   lG - log(G).
%
% Example (2 queues + delay):
%   L = [0.5 0.3; 0.2 0.4];  N = [3 2];  Z = [1 1];
%   [G,lG] = pfqn_is(L, N, Z, struct('samples',1e5));
%
% See also PFQN_LD_IS, PFQN_OI_IS, PFQN_PAS_IS, PFQN_NC.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 4, options = struct(); end
if nargin < 3, Z = []; end

[G, lG] = pfqn_ld_is(L, N, Z, [], options);
end
