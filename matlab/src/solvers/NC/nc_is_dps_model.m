function isdps = nc_is_dps_model(sn)
% ISDPS = NC_IS_DPS_MODEL(SN)
%
% True when the model is the closed two-station network that Morrison's
% heavy-usage asymptotic expansion is derived for: one infinite-server (think)
% station and one single-server discriminatory-processor-sharing station, with
% exponential service everywhere and every class alternating between the two.
% Such models are solved by SOLVER_NC_DPS_ANALYZER (npfqn_dps_morrison).
%
% The shape is checked exactly, not approximately: outside it the expansion has
% no derivation behind it, so a model that misses any clause here is left to
% the ordinary NC routes (which refuse DPS) rather than answered wrongly.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
isdps = false;

% closed, exactly two stations
if any(isinf(sn.njobs)) || ~any(sn.njobs > 0)
    return
end
if sn.nstations ~= 2
    return
end

% one INF station and one single-server DPS station
iInf = find(sn.sched == SchedStrategy.INF);
iDps = find(sn.sched == SchedStrategy.DPS);
if numel(iInf) ~= 1 || numel(iDps) ~= 1
    return
end
if isfinite(sn.nservers(iDps)) && sn.nservers(iDps) ~= 1
    return    % multi-server DPS: the min(n,c) share is not Morrison's
end

% no class switching: chain == class
if sn.nchains ~= sn.nclasses
    return
end

% no load-, class- or joint-dependent scaling
if ~isempty(sn.lldscaling) || ~isempty(sn.cdscaling) || ~isempty(sn.jdscaling)
    return
end

K = sn.nclasses;
for r = 1:K
    if sn.njobs(r) <= 0
        return    % an empty class leaves the expansion's b_j = 0
    end
    % exponential service at both stations
    for ist = [iInf, iDps]
        if ~isfinite(sn.rates(ist, r)) || sn.rates(ist, r) <= 0
            return
        end
        if isfinite(sn.scv(ist, r)) && abs(sn.scv(ist, r) - 1) > 1e-6
            return
        end
    end
    % positive DPS weight
    if ~isfinite(sn.schedparam(iDps, r)) || sn.schedparam(iDps, r) <= 0
        return
    end
end

% every class visits the two stations equally often (think -> DPS -> think)
V = zeros(sn.nstations, K);
for r = 1:K
    c = find(sn.chains(:, r));
    if numel(c) ~= 1
        return
    end
    vis = sn.visits{c};
    for ist = 1:sn.nstations
        V(ist, r) = vis(sn.stationToStateful(ist), r);
    end
    vref = V(sn.refstat(r), r);
    if vref <= 0
        return
    end
    V(:, r) = V(:, r) / vref;
    if abs(V(iInf, r) - V(iDps, r)) > 1e-9 * max(1, V(iInf, r))
        return    % unequal visits: not the alternating cycle of the paper
    end
end

isdps = true;
end
