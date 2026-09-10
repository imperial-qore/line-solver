function [Q,U,R,T,C,X,lG,iter] = solver_sqd(sn, options)
% [Q,U,R,T,C,X,LG,ITER] = SOLVER_BAS(SN, OPTIONS)
%
% Blocking-After-Service (BAS) approximate MVA handler.
%
% Wraps npfqn_sqd (chain-aggregated approximation): solves the single closed chain,
% then disaggregates the chain-level throughput, queue length and utilization back to
% per-class results via sn_deaggregate_chain_results. Supports single-chain closed
% networks (one or more classes connected by class switching); multi-chain models are
% rejected, since the BAS approximation models a single circulating population.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

iter = 0;
lG = NaN;

if sn.nchains ~= 1
    line_warning(mfilename, sprintf(['SQD (Smith Queue Decomposition) supports single-chain closed ' ...
        'networks only; this model is multichain (nchains=%d) - returning empty results.'], sn.nchains));
    M0 = sn.nstations; K0 = sn.nclasses;
    Q = nan(M0,K0); U = nan(M0,K0); R = nan(M0,K0); T = nan(M0,K0); C = nan(1,K0); X = nan(1,K0);
    return;
end

[Lchain,STchain,Vchain,alpha,~,~,refstatchain] = sn_get_demands_chain(sn);

M = sn.nstations;
[Xst,Qst,Ust,Rst] = npfqn_sqd(sn, sn.nclosedjobs);

% Assemble chain-level (M x 1) matrices for the single closed chain.
% Vchain is normalized to 1 at the reference station, so the per-station
% throughput there equals the chain reference throughput.
refstat = refstatchain(1);
Xchain = Xst(refstat);
Tchain = Xst(:);
Qchain = Qst(:);
Uchain = Ust(:);
Rchain = Rst(:);

[Q,U,R,T,C,X] = sn_deaggregate_chain_results(sn, Lchain, [], STchain, Vchain, alpha, ...
    Qchain, Uchain, Rchain, Tchain, [], Xchain);

end
