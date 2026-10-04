function [b, lG, lH] = getAvgBusyPeriod(self, stations, n)
% [B,LG,LH] = GETAVGBUSYPERIOD(STATIONS, N)
%
% Mean duration of the busy period of order N for the subnetwork made of
% STATIONS, that is the time from the instant a job entering the subnetwork
% finds N-1 jobs in it up to the next instant when fewer than N remain.
%
% Input:
%   stations - stations forming the subnetwork, as objects, names, or station
%              indexes; must be a non-empty proper subset of the stations
%   n        - busy period order(s), 1 <= n <= population for a closed model
%
% Output:
%   b  - mean busy period duration(s), same size as n
%   lG - log normalizing constants of the subnetwork, orders 0,1,...
%   lH - log normalizing constants of the complement (closed models only)
%
% H. Daduna, "Busy Periods for Subnetworks in Stochastic Networks: Mean Value
% Analysis", J. ACM 35(3), 1988. The result is exact on the single-chain
% product-form class and, by the insensitivity of Section 5 of that paper,
% depends on the service processes only through their mean rates.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if GlobalConstants.DummyMode
    b = NaN; lG = []; lH = [];
    return
end

sn = self.model.getStruct(false);
if sn.nchains > 1
    line_error(mfilename, 'The busy period of a subnetwork is defined for single-chain models only. Section 5 of Daduna (1988) sketches the multichain extension, which is not implemented.');
end

% station indexes of the subnetwork
subnet = zeros(1, numel(stations));
for t = 1:numel(stations)
    if isnumeric(stations)
        subnet(t) = stations(t);
    else
        if iscell(stations)
            st = stations{t};
        else
            st = stations(t);
        end
        if ischar(st) || isstring(st)
            subnet(t) = sn.nodeToStation(self.model.getNodeIndex(char(st)));
        else
            subnet(t) = sn.nodeToStation(st.index);
        end
    end
end

[~, STchain, Vchain, ~, Nchain] = sn_get_demands_chain(sn);
[rtst, Vst] = sn_rt_stations(sn);
M = sn.nstations;
K = sn.nclasses;

% station-to-station routing of the chain: the class-level probabilities
% weighted by the class visits, which is exact because it is a flow balance
Pst = zeros(M, M);
for i = 1:M
    for j = 1:M
        blk = rtst((i-1)*K + (1:K), (j-1)*K + (1:K));
        Pst(i,j) = Vst(i,:) * sum(blk, 2);
    end
end
Vtot = sum(Vst, 2)';
nz = Vtot > 0;
Pst(nz,:) = Pst(nz,:) ./ Vtot(nz)';

% mu(j,k): the load-dependent rate of station j holding k jobs. STchain is the
% chain mean service time, lldscaling (or min(k,c), k for an infinite server)
% is the dimensionless speedup.
lld = sn.lldscaling;
nservers = sn.nservers;
    function r = ratefun(j, kvec)
        % same precedence as solver_ncld: the infinite server first, then the
        % declared load-dependent scaling, then the multiserver staircase. A
        % scaling table shorter than kvec keeps its last entry.
        if isinf(nservers(j))
            scale = kvec;
        elseif ~isempty(lld)
            scale = lld(j, min(kvec, size(lld,2)));
        else
            scale = min(kvec, nservers(j));
        end
        r = scale / STchain(j);
    end

isOpen = ~isfinite(Nchain(1));
if isOpen
    % the Source is not a node of the Jackson network of the paper: its outflow
    % is the external stream gamma
    source = find(sn.nodetype(sn.stationToNode) == NodeType.Source, 1);
    if isempty(source)
        line_error(mfilename, 'An open model must own a Source station.');
    end
    if any(subnet == source)
        line_error(mfilename, 'The Source cannot belong to the subnetwork.');
    end
    lambda = sum(sn.rates(source, ~isnan(sn.rates(source,:))));
    keep = setdiff(1:M, source);
    remap = zeros(1, M);
    remap(keep) = 1:numel(keep);
    alpha = lambda * Vchain(keep)' / Vchain(source);
    gamma = lambda * Pst(source, keep);
    P = Pst(keep, keep);
    murows = @(j, kvec) ratefun(keep(j), kvec);
    [b, lG, lH] = pfqn_busyp(alpha, murows, P, Inf, remap(subnet), n, gamma);
else
    alpha = Vchain(:)';
    murows = @(j, kvec) ratefun(j, kvec);
    [b, lG, lH] = pfqn_busyp(alpha, murows, Pst, Nchain(1), subnet, n);
end
end
