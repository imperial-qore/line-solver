function [QN,UN,RN,TN,CN,XN,lG,hitprob,missprob,delayedprob,hitproblist,itemprob,latency,runtime,method] = solver_mva_cacheqn_retrieval_analyzer(sn, options)
% [...] = SOLVER_MVA_CACHEQN_RETRIEVAL_ANALYZER(SN, OPTIONS)
%
% Analyzer for a CLOSED integrated cache-queueing model whose Cache node has a
% delayed-hit retrieval system. Delegates to da_cacheqn_retrieval, which relabels
% the cache as a class switch and lets the finite-population coalescing emerge from
% the closed AMVA via a load-dependent fetch station. Returns the true cache
% hit/miss probabilities (hit = P(item cached), miss = 1 - hit); the delayed-hit
% fraction is folded into miss (delayedprob = 0).
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

T0 = tic;
[res, hitprob, missprob, delayedprob, ~, ~] = da_cacheqn_retrieval(sn, @(snit) netsolve(snit, options), options);

QN = res.QN; UN = res.UN; RN = res.RN; TN = res.TN; XN = res.XN;
if isfield(res,'CN'), CN = res.CN; else, CN = NaN(1, sn.nclasses); end
if isfield(res,'lG'), lG = res.lG; else, lG = NaN; end

K = sn.nclasses;
hitproblist = NaN(K, numel(sn.nodeparam{find(sn.nodetype==NodeType.Cache,1)}.itemcap));
itemprob = [];
latency = NaN(1, K);
method = 'fpi';
runtime = toc(T0);

    function res = netsolve(snit, options)
        res = struct();
        if ~isempty(snit.lldscaling) || ~isempty(snit.cdscaling) || ~isempty(snit.jdscaling)
            [res.Q,res.U,res.R,res.T,res.C,res.X,res.lG,res.runtime] = solver_mvald_analyzer(snit, options);
        else
            [res.Q,res.U,res.R,res.T,res.C,res.X,res.lG,res.runtime] = solver_mva_analyzer(snit, options);
        end
        res.QN = res.Q; res.UN = res.U; res.RN = res.R; res.TN = res.T; res.XN = res.X; res.CN = res.C;
    end
end
