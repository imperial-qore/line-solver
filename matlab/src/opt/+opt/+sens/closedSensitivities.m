function sens = closedSensitivities(model)
% closedSensitivities  Analytic d(metric)/d(rate) for closed single-server
% product-form networks with unit visit ratios, via the differentiated-MVA
% primitive pfqn_sens. Mirrors native-Python _closed_sensitivities. Returns []
% for open, multiserver, or non-unit-visit networks.

sens = [];
sn = model.getStruct();
R = sn.nclasses;
N = sn.njobs(:);
if isempty(N) || any(~isfinite(N))
    return;   % closed only
end

[~, D, ~, Zmat, ~, S] = sn_get_product_form_params(sn);   % D: Mq x R
S = S(:);                 % server counts per queue
if any(S(isfinite(S)) > 1), return; end   % single-server only
Mq = size(D, 1);

Z = zeros(1, R);
if ~isempty(Zmat)
    zsum = sum(Zmat, 1);
    if numel(zsum) == R, Z = zsum(:).'; end
end

res = pfqn_sens(D, N(:).', Z);

rates = sn.rates;
nodeToStation = sn.nodeToStation;
queueNodes = find(sn.nodetype(:).' == double(NodeType.Queue));
if numel(queueNodes) ~= Mq, return; end

params = {};   % each {j, s, p, chain, pkey}
for j = 1:numel(queueNodes)
    nodeJ = queueNodes(j);
    sj = nodeToStation(nodeJ);
    for s = 1:R
        ratej = rates(sj, s);
        if ~isfinite(ratej) || ratej <= 0 || D(j, s) <= 0, continue; end
        p = (j - 1) * R + s;   % 1-based param index into pfqn_sens outputs
        params{end+1} = {j, s, p, -D(j, s) / ratej, ...
            opt.SensitivityData.paramKey(sn.nodenames{nodeJ}, sn.classnames{s})}; %#ok<AGROW>
    end
end

sens = opt.SensitivityData();
for i = 1:numel(queueNodes)
    nodeI = queueNodes(i);
    nameI = sn.nodenames{nodeI};
    for r = 1:R
        if D(i, r) <= 0, continue; end
        mkey = opt.SensitivityData.metricKey(nameI, sn.classnames{r});
        for pk = 1:numel(params)
            pr = params{pk};
            p = pr{3}; chain = pr{4}; pkey = pr{5};
            sens.add('RespT', mkey, pkey, res.dR(i, r, p) * chain);
            sens.add('QLen', mkey, pkey, res.dQ(i, r, p) * chain);
            sens.add('Tput', mkey, pkey, res.dX(r, p) * chain);
            sens.add('Util', nameI, pkey, res.dU(i, r, p) * chain);
        end
    end
end
end
