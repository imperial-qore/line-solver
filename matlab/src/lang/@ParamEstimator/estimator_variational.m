function estVal = estimator_variational(self, nodes)
% ESTIMATOR_VARIATIONAL Variational inference for Markovian queueing networks
%
% Estimates service rates from noisy queue-length readings taken over time,
% with the method of I. Perez, G. Casale, "Variational Inference for Markovian
% Queueing Networks", Advances in Applied Probability 53(3), 2021. The data are
% QLen timeseries, one per (node, class); the estimator treats each reading as
% exact with probability 1-epsilon and uniform over the remaining feasible
% values otherwise, and returns the mean service time under the conjugate Gamma
% posterior of each estimated rate.
%
% The network is translated into the transition set eta=(i,j,c) of the paper,
% with the transition rate lambda_eta = mu_{i,c} p^c_{i,j}. Routing
% probabilities are taken as known from the model; only the station rates of
% the requested nodes are estimated, every other rate is held at its model
% value.
%
% Options (self.options): epsilon, prior_shape, iter_max, tol, ngrid, ymax,
% nsamples, delta, verbose. The posterior of each estimated rate is returned in
% self.options.posterior as [alpha beta] rows.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
% This code is released under the 3-Clause BSD License.

sn = self.model.getStruct;

if ~iscell(nodes)
    nodes = {nodes};
end

M = sn.nstations;
R = sn.nclasses;
if any(sn.nodetype == NodeType.Cache) || any(any(sn.csmask & ~eye(R)))
    line_error(mfilename, 'The variational estimator does not support class switching or caches.');
end

% station-to-station routing, with the pseudo-closed feedback of an open
% class stripped: a departure to the sink is a transition out of the network
rtst = sn_rt_stations(sn);
sourceIdx = find(sn.sched == SchedStrategy.EXT);

% transitions eta = (i,j,c), one per positive routing probability
arcs = zeros(0,3);
routeprob = zeros(0,1);
for c = 1:R
    for i = 1:M
        for j = 1:M
            p = rtst((i-1)*R + c, (j-1)*R + c);
            if p <= 0 || i == j
                continue
            end
            if any(sourceIdx == j)
                continue % the sink-to-source feedback is not a job transition
            end
            arcs(end+1,:) = [i j c]; %#ok<AGROW>
            routeprob(end+1,1) = p; %#ok<AGROW>
        end
    end
end
if isempty(arcs)
    line_error(mfilename, 'The model has no job transitions to infer from.');
end

% station disciplines: 0 = infinite server, 1 = shared server, 2 = external
sched = ones(M,1);
for i = 1:M
    switch sn.sched(i)
        case SchedStrategy.INF
            sched(i) = 0;
        case SchedStrategy.EXT
            sched(i) = 2;
        case {SchedStrategy.PS, SchedStrategy.FCFS, SchedStrategy.DPS, ...
              SchedStrategy.GPS, SchedStrategy.SIRO, SchedStrategy.LCFS}
            sched(i) = 1;
        otherwise
            line_error(mfilename, sprintf('The variational estimator does not support scheduling %s at station %d.', SchedStrategy.toText(sn.sched(i)), i));
    end
end
nservers = sn.nservers(:);
nservers(~isfinite(nservers)) = 1;

% which station-class rates are being estimated
estimated = zeros(M, R);
P = 0;
for n = 1:numel(nodes)
    i = sn.nodeToStation(nodes{n}.index);
    if i <= 0 || i > M
        line_error(mfilename, 'A node handed to the estimator is not a station.');
    end
    for r = 1:R
        if sn.rates(i,r) > 0 && isfinite(sn.rates(i,r))
            P = P + 1;
            estimated(i,r) = P;
        end
    end
end
if P == 0
    line_error(mfilename, 'No station-class pair with a positive service rate was selected.');
end

narcs = size(arcs,1);
arcparam = zeros(narcs,1);
arcrate = NaN(narcs,1);
for e = 1:narcs
    i = arcs(e,1);
    c = arcs(e,3);
    if estimated(i,c) > 0
        arcparam(e) = estimated(i,c);
    else
        arcrate(e) = sn.rates(i,c);
        if ~(arcrate(e) > 0) || ~isfinite(arcrate(e))
            line_error(mfilename, sprintf('Station %d class %d has no usable rate to hold fixed.', i, c));
        end
    end
end

% observations: QLen timeseries, one column per (station, class) pair
obsTimes = [];
for n = 1:numel(nodes)
    for r = 1:R
        qlData = self.getQLen(nodes{n}, self.model.classes{r});
        if isempty(qlData)
            continue
        end
        if iscell(qlData)
            qlData = qlData{1};
        end
        obsTimes = union(obsTimes, qlData.t(:));
    end
end
allNodes = self.model.getNodes;
for n = 1:numel(allNodes)
    i = sn.nodeToStation(allNodes{n}.index);
    if i <= 0 || i > M
        continue
    end
    for r = 1:R
        qlData = self.getQLen(allNodes{n}, self.model.classes{r});
        if isempty(qlData)
            continue
        end
        if iscell(qlData)
            qlData = qlData{1};
        end
        obsTimes = union(obsTimes, qlData.t(:));
    end
end
obsTimes = sort(obsTimes(:));
if isempty(obsTimes)
    line_error(mfilename, 'The variational estimator needs QLen timeseries data.');
end
K = numel(obsTimes);
obsData = NaN(K, M*R);
for n = 1:numel(allNodes)
    i = sn.nodeToStation(allNodes{n}.index);
    if i <= 0 || i > M
        continue
    end
    for r = 1:R
        qlData = self.getQLen(allNodes{n}, self.model.classes{r});
        if isempty(qlData)
            continue
        end
        if iscell(qlData)
            qlData = qlData{1};
        end
        for k = 1:numel(qlData.t)
            [tf, pos] = ismember(qlData.t(k), obsTimes);
            if tf
                obsData(pos, (r-1)*M + i) = round(qlData.data(k));
            end
        end
    end
end

% population per class bounds both the contamination support and the load
popr = zeros(1,R);
for r = 1:R
    if isfinite(sn.njobs(r))
        popr(r) = sn.njobs(r);
    else
        col = obsData(:, (r-1)*M + (1:M));
        popr(r) = max(1, 2*max(col(~isnan(col)))); % open class: twice the peak observed
    end
end
obsRange = zeros(M,R);
capacity = Inf(M,R);
for r = 1:R
    for i = 1:M
        obsRange(i,r) = popr(r);
        if isfinite(sn.njobs(r))
            capacity(i,r) = popr(r);
        end
    end
end

% initial state: the model's own if set, otherwise the closed population at
% the reference station
x0 = zeros(M,R);
for r = 1:R
    if isfinite(sn.njobs(r)) && sn.njobs(r) > 0
        ref = sn.refstat(r);
        if ref >= 1 && ref <= M
            x0(ref,r) = sn.njobs(r);
        else
            x0(1,r) = sn.njobs(r);
        end
    end
end

% Gamma priors centred on the model's current rates
shape0 = 1;
if isfield(self.options, 'prior_shape') && ~isempty(self.options.prior_shape)
    shape0 = self.options.prior_shape;
end
alpha0 = zeros(P,1);
beta0 = zeros(P,1);
for i = 1:M
    for r = 1:R
        p = estimated(i,r);
        if p > 0
            alpha0(p) = shape0;
            beta0(p) = shape0 / sn.rates(i,r);
        end
    end
end

epsilon = 0.05;
if isfield(self.options, 'epsilon') && ~isempty(self.options.epsilon)
    epsilon = self.options.epsilon;
end

spec = struct('arcs', arcs, 'x0', x0, 'sched', sched, 'nservers', nservers, ...
    'routeprob', routeprob, 'arcparam', arcparam, 'arcrate', arcrate, ...
    'alpha0', alpha0, 'beta0', beta0, 'obsTimes', obsTimes, 'obsData', obsData, ...
    'obsRange', obsRange, 'epsilon', epsilon, 'capacity', capacity);

vopt = struct();
for fn = {'iter_max','tol','ngrid','dt','ymax','nsamples','delta','rate_max','verbose','tmax'}
    if isfield(self.options, fn{1}) && ~isempty(self.options.(fn{1}))
        vopt.(fn{1}) = self.options.(fn{1});
    end
end
if ~isfield(vopt, 'iter_max')
    vopt.iter_max = 20;
end

out = infer_variational(spec, vopt);
self.options.posterior = [out.alpha(:), out.beta(:)];
self.options.bound = out.bound;

estVal = zeros(numel(nodes), R);
for n = 1:numel(nodes)
    i = sn.nodeToStation(nodes{n}.index);
    for r = 1:R
        p = estimated(i,r);
        if p > 0
            estVal(n,r) = out.meanServiceTime(p);
        elseif sn.rates(i,r) > 0
            estVal(n,r) = 1 / sn.rates(i,r);
        end
    end
end
end
