function reason = nc_oi_refusal(sn, route)
% REASON = NC_OI_REFUSAL(SN, ROUTE)
%
% May the order-independent route ROUTE run on this model? '' when it may,
% otherwise the reason it may not, in the words the analyzer refuses with.
%
%   ROUTE = 'oi'   the exact OI convolution, SOLVER_NC_OI_ANALYZER, reached on
%                  the shape NC_IS_OI_MODEL tests for;
%   ROUTE = 'pas'  the importance sampler of the two-station pass-and-swap
%                  tandem, SOLVER_NC_PAS_IS_ANALYZER, reached on the shape
%                  NC_IS_PAS_MODEL tests for (and on an OI model asked for by
%                  'is' or 'sampling').
%
% The shape predicates say which analyzer a model reaches; this one says what
% that analyzer then refuses. The rank rate mu(supp n) of an OI station is a
% function of the RAW classes present, so a chain that switches class has no
% rank rate; the sampler places one job of each class per position, so it
% needs unit per-class visits; and both routes are intercepted by
% @SolverNC/runAnalyzer AHEAD of the load-dependent and fork-join dispatches
% and build every BCMP station table from the server count alone, so a rate
% lattice would be dropped and a Fork would send the model to an image with
% open auxiliary classes, on which neither shape holds.
%
% ONE PREDICATE, THREE CALLERS. Each analyzer asks it first and turns a
% non-empty answer into an error; NC_METHOD_REFUSAL asks it so that MODEL.HELP
% never offers a pair the run refuses. A new rule goes here.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

reason = '';
if nargin < 2 || isempty(route)
    route = 'oi';
end
switch lower(route)
    case 'pas'
        who = 'solver_nc_pas_is';
    otherwise
        who = 'solver_nc_oi';
end

for c = 1:sn.nchains
    if numel(sn.inchain{c}) > 1
        reason = sprintf('%s requires one class per chain (no class switching).', who);
        return
    end
end
if any(isinf(sn.njobs))
    reason = sprintf('%s requires a closed queueing network.', who);
    return
end
if any(sn.nodetype == NodeType.Fork)
    reason = sprintf(['%s cannot serve a fork-join model: the fork-join transformation carries the ' ...
        'parallelism as open auxiliary classes, and the order-independent routes are stated for a ' ...
        'closed network.'], who);
    return
end
if ~isempty(sn.lldscaling) || ~isempty(sn.cdscaling) || ~isempty(sn.jdscaling)
    reason = sprintf(['%s builds every BCMP station table from the server count and reads no load-, ' ...
        'class- or joint-dependent rate; declare the dependence through the OI rank rate instead.'], who);
    return
end
M = sn.nstations;
for ist = 1:M
    if sn.sched(ist) ~= SchedStrategy.PAS && sn.sched(ist) ~= SchedStrategy.OI
        continue
    end
    ind = sn.stationToNode(ist);
    if ind < 1 || ind > numel(sn.nodeparam) || ~isstruct(sn.nodeparam{ind}) ...
            || ~isfield(sn.nodeparam{ind}, 'svcRateFun') || isempty(sn.nodeparam{ind}.svcRateFun)
        reason = sprintf('OI station %d has no service rate function; set it via setService(@(c) ...).', ist);
        return
    end
end
if strcmpi(route, 'pas')
    if M ~= 2
        reason = sprintf('%s models a two-station pass-and-swap tandem (got %d stations).', who, M);
        return
    end
    K = sn.nclasses;
    N = round(sn.njobs(:)');
    V = zeros(M, K);
    for r = 1:K
        c = find(sn.chains(:, r));
        vis = sn.visits{c};
        for ist = 1:M
            V(ist, r) = vis(sn.stationToStateful(ist), r);
        end
        vref = V(sn.refstat(r), r);
        if vref > 0
            V(:, r) = V(:, r) / vref;
        end
    end
    for ist = 1:M
        for r = 1:K
            if N(r) > 0 && abs(V(ist, r) - 1) > 1e-9
                reason = sprintf('%s requires unit per-class visits (station %d, class %d, V=%g).', who, ist, r, V(ist, r));
                return
            end
        end
    end
end
end
