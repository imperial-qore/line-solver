function [chainModel, alpha, deaggInfo] = aggregateChains(self, suffix)
% [CHAINMODEL, ALPHA, DEAGGINFO] = AGGREGATECHAINS(SELF)
% [CHAINMODEL, ALPHA, DEAGGINFO] = AGGREGATECHAINS(SELF, SUFFIX)
%
% Return a copy of this model in which all classes belonging to the same
% chain are merged into a single aggregate class, so that the aggregated
% model has one class per chain of the original. Class switching is
% eliminated in the process. SUFFIX is appended to the names of the
% aggregate classes (default: '').
%
% ALPHA is the (nstations x nclasses) matrix of aggregation factors and
% DEAGGINFO carries everything needed to map chain-level metrics back to
% class-level metrics with sn_deaggregate_chain_results:
%
%   [Q,U,R,T,C,X] = sn_deaggregate_chain_results(self.getStruct(), ...
%       deaggInfo.Lchain, [], deaggInfo.STchain, deaggInfo.Vchain, ...
%       deaggInfo.alpha, Qchain, Uchain, Rchain, Tchain, [], Xchain);
%
% The aggregation is exact for product-form models and approximate
% otherwise, since a single aggregate service process replaces the
% per-class ones weighted by ALPHA.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 2
    suffix = '';
end
[chainModel, alpha, deaggInfo] = ModelAdapter.aggregateChains(self, suffix);
end
