function [QN,UN,RN,TN,CN,XN,totiter] = solver_mam_basic_mmap(sn, options)
% [QN,UN,RN,TN,CN,XN,TOTITER] = SOLVER_MAM_BASIC_MMAP(SN, OPTIONS)
%
% Top-level dispatcher for the MAM/MMAP fork-join decomposition.
% Open networks call solver_mam_basic_mmap_inner directly with arrival
% rates derived from the source/refstat. Closed networks go through
% solver_mam_basic_mmap_closed, which wraps the inner algorithm in an
% MNA-style bisection on per-class throughput.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if sn_is_open_model(sn)
    K = sn.nclasses;
    C = sn.nchains;
    lambda = zeros(1,K);
    for c=1:C
        inchain = sn.inchain{c};
        lambdas_inchain = sn.rates(sn.refstat(inchain(1)), inchain);
        lambdas_inchain = lambdas_inchain(isfinite(lambdas_inchain));
        lambda(inchain) = sum(lambdas_inchain);
    end
    [QN,UN,RN,TN,CN,XN,totiter] = solver_mam_basic_mmap_inner(sn, options, lambda);
else
    [QN,UN,RN,TN,CN,XN,totiter] = solver_mam_basic_mmap_closed(sn, options);
end
end
