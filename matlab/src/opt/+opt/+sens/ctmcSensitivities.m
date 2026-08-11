function sens = ctmcSensitivities(model, options)
% ctmcSensitivities  d(metric)/d(rate) for an arbitrary Markovian model, via
% the generator-derivative equation of Trivedi and Bobbio (2017), Eq. (9.81).
%
% This is the fallback for the cases the differentiated-MVA primitive cannot
% reach: open, multiserver, and non-unit-visit networks, for which
% opt.sens.openSensitivities and opt.sens.closedSensitivities return []. It
% is exact wherever SolverCTMC is exact, and correspondingly it is limited by
% the state-space size rather than by the product-form assumptions.
%
% Only queue-length sensitivities are produced. Utilization, response time,
% and throughput are reward rates whose definition involves the solved
% metrics themselves, so they need the reward-derivative term of Eq. (9.83)
% and are not covered here.
%
% Returns [] if the model has no CTMC representation of tractable size.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

sens = [];
if nargin < 2 || isempty(options)
    options = SolverCTMC.defaultOptions();
end

sn = model.getStruct();
M = sn.nstations;
K = sn.nclasses;

solver = SolverCTMC(model, options);
try
    Q = full(solver.getGenerator());
    spaceAggr = solver.getStateSpaceAggr();
catch
    return;   % state space unavailable or over the memory gate
end
if isempty(spaceAggr)
    return;
end

pi = ctmc_solve(Q);

% One parameter per finite positive exponential service rate.
%
% The rate setter substitutes an Exp of the perturbed rate, which is a
% perturbation of theta only where the nominal process is itself exponential.
% For an Erlang or Coxian service the substitution would change the
% distribution family, and the resulting difference quotient would not be
% dQ/dtheta at all, so those stations are skipped rather than reported wrong.
params = {};
for ist = 1:M
    for k = 1:K
        rate = sn.rates(ist, k);
        if ~isfinite(rate) || rate <= 0
            continue;
        end
        if isempty(sn.procid) || size(sn.procid, 1) < ist || size(sn.procid, 2) < k
            continue;
        end
        if sn.procid(ist, k) ~= ProcessType.EXP
            continue;
        end
        nodeIdx = sn.stationToNode(ist);
        params{end+1} = struct(...
            'value', rate, ...
            'node', nodeIdx, ...
            'class', k, ...
            'pkey', opt.SensitivityData.paramKey(sn.nodenames{nodeIdx}, sn.classnames{k})); %#ok<AGROW>
    end
end

sens = opt.SensitivityData();
for p = 1:numel(params)
    pr = params{p};
    param = struct();
    param.value = pr.value;
    param.set = makeRateSetter(pr.node, pr.class);

    try
        [~, ~, dpi] = solver.getSensitivity(param);
    catch
        continue;   % perturbation changed the state space; skip this parameter
    end

    for ist = 1:M
        nameI = sn.nodenames{sn.stationToNode(ist)};
        for r = 1:K
            col = (ist - 1) * K + r;
            if col > size(spaceAggr, 2)
                continue;
            end
            rvec = spaceAggr(:, col);
            if any(~isfinite(rvec))
                continue;   % Source stations carry an infinite population
            end
            mkey = opt.SensitivityData.metricKey(nameI, sn.classnames{r});
            sens.add('QLen', mkey, pr.pkey, dpi * rvec);
        end
    end
end
end

function h = makeRateSetter(nodeIdx, k)
% h = makeRateSetter(nodeIdx, k)
% Handle setting the service rate of node NODEIDX class K on a model copy.
% The distribution is replaced by an exponential of the requested rate, so
% this is only meaningful where the nominal service is itself exponential.
% Node and class are resolved inside the copy, since the objects of the
% original model do not belong to it.

h = @(m, value) setRate(m, nodeIdx, k, value);
end

function setRate(m, nodeIdx, k, value)
nodes = m.getNodes();
classes = m.getClasses();
nodes{nodeIdx}.setService(classes{k}, Exp(value));
end
