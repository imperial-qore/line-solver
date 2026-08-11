function sens = openSensitivities(model)
% openSensitivities  Analytic d(metric)/d(rate) for open product-form networks
% (single-server queueing and infinite-server delay stations, which decouple),
% mirroring native-Python _open_sensitivities. Returns [] for closed, mixed, or
% multiserver networks (finite-difference fallback).
%
% Reference: Z. Liu and P. Nain, INRIA RR-1144 (1989), Thm 3.2 (open BCMP).

sens = [];
sn = model.getStruct();
R = sn.nclasses;
njobs = sn.njobs(:);
if isempty(njobs) || ~all(isinf(njobs))
    return;   % open only
end

nstations = sn.nstations;
rates = sn.rates;
nservers = sn.nservers(:);

% visits per (station,class): sum over chains (sn.visits is a cell per chain)
visits = zeros(nstations, R);
gotVisits = false;
if iscell(sn.visits) && ~isempty(sn.visits)
    for c = 1:numel(sn.visits)
        vm = sn.visits{c};
        if ~isempty(vm) && size(vm, 1) == nstations && size(vm, 2) == R
            visits = visits + full(vm);
            gotVisits = true;
        end
    end
end
if ~gotVisits || all(visits(:) == 0), visits(:) = 1.0; end

extId = SchedStrategy.toId(SchedStrategy.EXT);
infId = SchedStrategy.toId(SchedStrategy.INF);

lam = zeros(1, R);
for i = 1:nstations
    if sn.sched(i) == extId
        for r = 1:R
            if isfinite(rates(i, r)), lam(r) = lam(r) + rates(i, r); end
        end
    end
end

sens = opt.SensitivityData();
for i = 1:nstations
    sc = sn.sched(i);
    if sc == extId, continue; end
    isDelay = (sc == infId);
    if ~isDelay && nservers(i) > 1
        sens = []; return;   % multiserver M/M/c not handled here
    end
    stName = sn.nodenames{sn.stationToNode(i)};

    D = zeros(1, R); rho = zeros(1, R);
    for r = 1:R
        mu = rates(i, r);
        if isfinite(mu) && mu > 0 && visits(i, r) > 0
            D(r) = visits(i, r) / mu;
            rho(r) = lam(r) * D(r);
        end
    end
    U = 0.0; if ~isDelay, U = sum(rho); end
    denom = 1.0 - U;
    if ~isDelay && denom <= 0, sens = []; return; end

    for s = 1:R
        muS = rates(i, s);
        if isfinite(muS) && muS > 0 && visits(i, s) > 0
            sens.add('Util', stName, opt.SensitivityData.paramKey(stName, sn.classnames{s}), -rho(s) / muS);
        end
    end

    for r = 1:R
        mu = rates(i, r);
        if visits(i, r) <= 0 || ~(isfinite(mu) && mu > 0), continue; end
        cl = sn.classnames{r};
        mkey = opt.SensitivityData.metricKey(stName, cl);
        for s = 1:R
            muS = rates(i, s);
            if ~(isfinite(muS) && muS > 0 && visits(i, s) > 0), continue; end
            pkey = opt.SensitivityData.paramKey(stName, sn.classnames{s});
            dU = 0.0; if ~isDelay, dU = -rho(s) / muS; end
            if s == r, dDr = -D(r) / muS; dRhoR = -rho(r) / muS; else, dDr = 0.0; dRhoR = 0.0; end
            if isDelay
                dResp = dDr; dQ = dRhoR;
            else
                dResp = (dDr * denom + D(r) * dU) / (denom * denom);
                dQ = (dRhoR * denom + rho(r) * dU) / (denom * denom);
            end
            sens.add('RespT', mkey, pkey, dResp);
            sens.add('QLen', mkey, pkey, dQ);
        end
        sens.add('Tput', mkey, opt.SensitivityData.paramKey(stName, cl), 0.0);
    end
end
end
