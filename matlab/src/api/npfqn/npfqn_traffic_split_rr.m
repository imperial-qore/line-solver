function kRR = npfqn_traffic_split_rr(sn)
% KRR = NPFQN_TRAFFIC_SPLIT_RR(SN)
%
% Deterministic (round-robin) split degree of the departure stream of each
% station-class. KRR(i,r)=k>1 means that the class-r departures of station i
% are dispatched one-in-k by a round-robin node, so that a downstream flow
% carrying a fraction p of them is the k-fold convolution thinned with
% probability q=k*p and has SCV 1+p*(d2-k), against the Markovian 1+p*(d2-1).
% KRR(i,r)=1 marks an ordinary probabilistic split.
%
% Only two topologies admit the deterministic rule: the station dispatches
% round-robin itself, or it feeds with probability one a router that does and
% whose pointer no other flow advances. Anything else falls back to k=1.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

M = sn.nstations;
K = sn.nclasses;
kRR = ones(M,K);

if ~isfield(sn,'routing') || isempty(sn.routing)
    return
end
if ~any(sn.routing(:) == RoutingStrategy.RROBIN)
    return
end

for ist=1:M
    ind = sn.stationToNode(ist);
    for r=1:K
        if sn.routing(ind,r) == RoutingStrategy.RROBIN
            kRR(ist,r) = rr_degree(sn, ind, r);
            continue
        end
        % a sure transition into a router that dispatches round-robin
        dest = find(sn.rtnodes((ind-1)*K+r, :) > 0);
        if ~isscalar(dest)
            continue
        end
        jnd = ceil(dest/K);
        s = dest - (jnd-1)*K;
        if sn.isstation(jnd) || sn.routing(jnd,s) ~= RoutingStrategy.RROBIN
            continue
        end
        if sn.rtnodes((ind-1)*K+r, dest) < 1 - GlobalConstants.FineTol
            continue
        end
        % the round-robin pointer must be advanced by this stream alone
        if ~isscalar(find(sn.rtnodes(:, dest) > 0))
            continue
        end
        kRR(ist,r) = rr_degree(sn, jnd, s);
    end
end
end

function k = rr_degree(sn, ind, r)
% number of destinations the round-robin pointer of (ind,r) cycles through
if ~isempty(sn.nodeparam) && numel(sn.nodeparam) >= ind && ...
        ~isempty(sn.nodeparam{ind}) && numel(sn.nodeparam{ind}) >= r && ...
        isfield(sn.nodeparam{ind}{r}, 'outlinks')
    k = numel(sn.nodeparam{ind}{r}.outlinks);
else
    k = nnz(sn.connmatrix(ind,:));
end
k = max(1, k);
end
