function reason = nc_lcfs_refusal(sn)
% REASON = NC_LCFS_REFUSAL(SN)
%
% May the paired LCFS / LCFS-PR route of SOLVER_NC run on this model? '' when
% it may, otherwise the reason it may not, in the words SOLVER_NC refuses
% with. A model with no LCFS station is not that route's business and gets ''.
%
% A non-preemptive LCFS station is served by SOLVER_NC_LCFSQN alone, the
% convolution of Casale (QUESTA 2026) for the closed two-station network of one
% LCFS and one LCFS-PR station; PFQN_LCFSQN_CA reads the two stations' per-class
% rates and the per-class populations and nothing else, so anything the shape
% does not name is not ignored by accident but has no term in the product form:
% a third station, a self-loop, a server count above one, a rate lattice, a
% class that switches, an open chain (the fork-join image carries the
% parallelism as OPEN auxiliary classes, so a Fork counts as one).
%
% ONE PREDICATE, TWO CALLERS. SOLVER_NC asks it ahead of the LCFS arm and turns
% a non-empty answer into an error; NC_METHOD_REFUSAL asks it so that
% MODEL.HELP never offers a pair the run refuses. A new rule goes here.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

reason = '';
lcfsStat = find(sn.sched == SchedStrategy.LCFS);
if isempty(lcfsStat)
    return
end
lcfsprStat = find(sn.sched == SchedStrategy.LCFSPR);
if isempty(lcfsprStat)
    reason = 'LCFS scheduling requires a paired LCFS-PR station.';
    return
end
if numel(lcfsStat) ~= 1 || numel(lcfsprStat) ~= 1
    reason = 'LCFS NC requires exactly one LCFS and one LCFS-PR station.';
    return
end
if sn.nstations ~= 2
    reason = ['LCFS NC requires a two-station network: pfqn_lcfsqn_ca convolves the LCFS and ' ...
        'LCFS-PR stations only, so any other station would be dropped from the answer.'];
    return
end
if any(isinf(sn.njobs)) || any(sn.nodetype == NodeType.Fork)
    reason = 'LCFS NC requires a closed queueing network.';
    return
end
if sn.nchains ~= sn.nclasses
    reason = ['LCFS NC requires one class per chain: the convolution is over per-class ' ...
        'populations with one rate per class at each station, which a switching chain has not.'];
    return
end
nclasses = sn.nclasses;
for ist = [lcfsStat, lcfsprStat]
    isf = sn.stationToStateful(ist);
    for r = 1:nclasses
        if sn.rt((isf-1)*nclasses+r, (isf-1)*nclasses+r) > 0
            reason = 'LCFS NC does not support self-loops at stations.';
            return
        end
    end
end
if any(sn.nservers([lcfsStat, lcfsprStat]) ~= 1)
    reason = ['LCFS NC requires single-server LCFS and LCFS-PR stations: the order-dependent ' ...
        'product form is stated for one server at each, and the convolution reads no server count.'];
    return
end
if ~isempty(sn.lldscaling) || ~isempty(sn.cdscaling) || ~isempty(sn.jdscaling)
    reason = ['LCFS NC does not support load-, class- or joint-dependent rates: those routes ' ...
        '(solver_ncld, solver_nc_conv) never reach the LCFS arm and would serve the station as ' ...
        'an ordinary product-form queue.'];
    return
end
end
