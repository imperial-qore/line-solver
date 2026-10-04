function varargout = solver_tr_chains_analyzer(self, phase, varargin)
% VARARGOUT = SOLVER_TR_CHAINS_ANALYZER(SELF, PHASE, VARARGIN)
%
% Chain-aggregation strategy for TRANSFORMSOLVE. Collapses every chain onto a
% single class, solves that model with the CALLER'S OWN solver, and maps the
% chain metrics back onto the classes.
%
% WHY THIS EXISTS. The state space of a multiclass model grows with the
% per-class populations, so a model with several classes in one chain can be
% intractable while the same model with one class per chain is not.
% ModelAdapter.aggregateChains builds the collapsed model and
% sn_deaggregate_chain_results maps its metrics back through alpha, the
% per-station share of the chain's visits each class carries.
%
% WHAT IS TRADED. The aggregation is EXACT on a product-form model: the chain
% is the unit MVA and convolution already solve in, and the deaggregation is the
% same alpha-weighted split those solvers apply. It is an APPROXIMATION
% otherwise, because one aggregate service law, fitted to the alpha-weighted
% first two moments, replaces the per-class ones. A caller who needs the exact
% multiclass answer must leave the transform off and pay the state space.
%
% Phases: 'expand' builds the single aggregated submodel, 'lift' deaggregates.
% The transformation is SINGLE PASS: one solve of the aggregate determines the
% answer, so ctx.iterated stays false and no coupling or convergence phase is
% implemented.
%
% See also TRANSFORMSOLVE, TRANSFORM_METHOD, SN_DEAGGREGATE_CHAIN_RESULTS.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

switch phase
    case 'expand'
        [varargout{1}, varargout{2}] = expand_(self, varargin{:});
    case 'lift'
        [varargout{1}, varargout{2}, varargout{3}, ...
         varargout{4}, varargout{5}, varargout{6}] = lift_(self, varargin{:});
    otherwise
        line_error(mfilename, sprintf('Unknown chain-transformation phase: %s', phase));
end
end

function [submodels, ctx] = expand_(self, sn, options)
% With one class per chain the aggregation is the IDENTITY, and aggregateChains
% then returns no deaggregation tables at all. Refusing by name beats failing
% later on a missing field; SolverCTMC's own chain_aggregation entry applies the
% same nchains < nclasses guard before it ever reaches here.
if sn.nchains >= sn.nclasses
    line_error(mfilename, sprintf(['chain aggregation needs more classes than chains: this ' ...
        'model has %d classes in %d chains, so the transform is the identity.'], ...
        sn.nclasses, sn.nchains));
end
[chainModel, alpha, deagg] = ModelAdapter.aggregateChains(self.model);
ctx = struct('sn_orig', sn, 'alpha', alpha, 'deagg', deagg, 'iterated', false);
submodels = {chainModel};
line_debug(options, '%s: chain aggregation, %d classes collapsed onto %d chains.', ...
    self.getName(), sn.nclasses, sn.nchains);
end

function [QN,UN,RN,TN,CN,XN] = lift_(~, ctx, res)
r = res{1};
% sn_deaggregate_chain_results reads the ORIGINAL struct and the chain tables;
% ST is left empty so it recovers the per-class service times from sn.rates,
% which is the documented call in @MNetwork/aggregateChains. XN is the per-chain
% SYSTEM throughput, which is why transformSolve collects it on its own channel
% rather than taking the sixth output of getAvg.
[QN,UN,RN,TN,CN,XN] = sn_deaggregate_chain_results(ctx.sn_orig, ctx.deagg.Lchain, [], ...
    ctx.deagg.STchain, ctx.deagg.Vchain, ctx.alpha, r.QN, r.UN, r.RN, r.TN, [], r.XN);

if isempty(CN)
    CN = sum(RN,1);
end
end
