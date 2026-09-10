function pools = serverPools(self, ind)
% POOLS = SERVERPOOLS(IND)
% Heterogeneous server pools and job parallelism of node IND, or [] when the
% station declares neither.
%
% JMT's Server section takes classParallelism, serverNames,
% serversPerServerType, serverCompatibilities and schedulingPolicy as one
% positional block of its constructor (jmt.engine.NodeSections.Server), so the
% five are emitted together or not at all, and always after the service
% strategies. A station declaring parallelism alone is therefore given one
% synthetic pool holding all of its servers, since the pool counts, not
% numberOfServers, size the server pool once any pool is declared.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

pools = [];
sn = self.getStruct;
np = sn.nodeparam{ind};
if ~isstruct(np)
    np = struct();
end
hasTypes = isfield(np, 'nservertypes') && np.nservertypes > 0 && ...
    isfield(np, 'servertypenames') && ~isempty(np.servertypenames);
hasParallelism = isfield(np, 'serverparallelism') && any(np.serverparallelism > 1);
if ~hasTypes && ~hasParallelism
    return
end

K = sn.nclasses;
pools.parallelism = ones(1, K);
if isfield(np, 'serverparallelism')
    declared = np.serverparallelism;
    ncopy = min(numel(declared), K);
    pools.parallelism(1:ncopy) = max(1, declared(1:ncopy));
end

if hasTypes
    pools.names = np.servertypenames;
    pools.counts = np.serverspertype;
    pools.compat = np.servercompat > 0;
    if isfield(np, 'heteroschedpolicy') && ~isempty(np.heteroschedpolicy)
        pools.policy = np.heteroschedpolicy;
    else
        pools.policy = HeteroSchedPolicy.ORDER;
    end
else
    ist = sn.nodeToStation(ind);
    pools.names = {sprintf('%s - Server Type 1', sn.nodenames{ind})};
    pools.counts = sn.nservers(ist);
    pools.compat = true(1, K);
    pools.policy = HeteroSchedPolicy.ORDER;
end
end
