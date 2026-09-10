function demandEst = infer_mlps(model, node, rt, class, ql)
% INFER_MLPS MLPS demand estimation using sn struct-level operations.
%
% Estimates service demands at a PS queue using the Maximum Likelihood
% for Processor Sharing method. Pre-builds augmented CTMC models for
% each unique (tagClass, aQueue) combination once, then uses
% sn_set_service + solver_ctmc directly inside the optimization loop.
%
% Inputs:
%   model  - LINE Network model with delay rates set and queue rates to estimate
%   node   - PS queue node (Station object)
%   rt     - response time samples (column vector, n x 1)
%   class  - class of each sample (column vector, n x 1)
%   ql     - queue lengths at arrival (n x R matrix, per-class)
%
% Returns:
%   demandEst - 1 x R vector of estimated mean service demands
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
% This code is released under the 3-Clause BSD License.

sn = model.getStruct();
R = sn.nclasses;
nCores = node.getNumberOfServers();

% Get delay rates from model
delayIdx = 0;
for i = 1:sn.nstations
    if sn.sched(i) == SchedStrategy.INF
        delayIdx = i;
        break;
    end
end
muZ = zeros(1, R);
for k = 1:R
    muZ(k) = sn.mu{delayIdx}{k}(1);
end

% Initial point estimate
meanQL = mean(sum(ql, 2));
Vtilde = min(meanQL, nCores);
x0 = zeros(1, R);
for j = 1:R
    if ~isempty(rt(class == j))
        x0(j) = Vtilde * mean(rt(class == j)) / meanQL;
    else
        x0(j) = 1e-3;
    end
end

xLB = zeros(1, R);
xUB = max(rt) * ones(1, R);

newR = R + 1;
augMuZ = [muZ, 0]; % placeholder for tagged class rate

% see _kb/03-api-layer.md for rationale
uniqueTC = unique(class);
uniqueQL = unique(ql, 'rows');
ctmcOpts = []; % will be set from first solver instance

% Cached augmented models stored as parallel arrays keyed by state string:
% prebuiltKeys{i} is the key, prebuiltVals{i} the corresponding struct.
prebuiltKeys = {};
prebuiltVals = {};
for u = 1:length(uniqueTC)
    tc = uniqueTC(u);
    augMuZ(newR) = muZ(tc);
    for v = 1:size(uniqueQL, 1)
        aq = uniqueQL(v, :);

        % Population: move 1 job from tagClass to tagged class
        N = aq; N(tc) = N(tc) - 1; N(newR) = 1;

        % Build augmented Network model (once per combo)
        augModel = Network('mlps_aug');
        augNode = cell(1, 2);
        augNode{1} = Delay(augModel, 'Think');
        augNode{2} = Queue(augModel, 'Queue1', SchedStrategy.PS);
        augNode{2}.setNumberOfServers(nCores);
        for r = 1:newR
            augClass = ClosedClass(augModel, sprintf('Class%d', r), N(r), augNode{1}, 0);
            augNode{1}.setService(augClass, Exp(augMuZ(r)));
            augNode{2}.setService(augClass, Exp(1)); % placeholder
        end
        P = augModel.initRoutingMatrix;
        for r = 1:newR
            P{r} = Network.serialRouting(augNode);
        end
        augModel.link(P);

        % Use SolverCTMC to get proper options and run initial solve
        solver = SolverCTMC(augModel);
        if isempty(ctmcOpts)
            ctmcOpts = solver.getOptions();
        end
        [infGen, eventFilt, ev] = solver.getGenerator();
        stateSpaceAggr = solver.getStateSpaceAggr();
        augSn = augModel.getStruct();

        % Find tagged departure event indices (invariant across rate changes)
        queueNodeIdx = augNode{2}.index;
        taggedDepIdx = [];
        for e = 1:length(ev)
            if ev{e}.active{1}.node == queueNodeIdx && ...
               ev{e}.active{1}.class == newR && ...
               ev{e}.active{1}.event == EventType.DEP
                taggedDepIdx(end+1) = e;
            end
        end

        % Find subset where tagged job is at Queue (invariant)
        queueStIdx = augNode{2}.stationIndex;
        taggedColAtQueue = (queueStIdx - 1) * newR + newR;
        subset = find(stateSpaceAggr(:, taggedColAtQueue) == 1);

        % Extract queue state space for subset (invariant)
        queueCols = (queueStIdx - 1) * newR + (1:newR);
        SSqueue = stateSpaceAggr(subset, queueCols);

        cacheKey = mat2str([tc, aq]);
        prebuiltKeys{end+1} = cacheKey; %#ok<AGROW>
        prebuiltVals{end+1} = struct(...
            'sn', augSn, ...
            'queueStIdx', queueStIdx, ...
            'taggedDepIdx', taggedDepIdx, ...
            'subset', subset, ...
            'SSqueue', SSqueue, ...
            'N', N, ...
            'tagClass', tc); %#ok<AGROW>
    end
end

%% Optimization options
options = optimset();
options.Display = 'iter';
options.LargeScale = 'on';
options.MaxIter = 1e10;
options.MaxFunEvals = 1e10;
options.MaxSQPIter = 5000;
options.TolCon = 1e-6;
options.Algorithm = 'interior-point';

[demandEst, ~] = fmincon(@objfun, x0, [], [], [], [], xLB, xUB, [], options);

    function f = objfun(x)
        TOL = 1e-6;
        rates = 1 ./ x; % x = mean demands

        % Build cache: update rates via sn_set_service, call solver_ctmc.
        % Stored as parallel arrays keyed by state string.
        cacheKeys = {};
        cacheVals = {};

        for kk = 1:length(prebuiltKeys)
            key = prebuiltKeys{kk};
            pb = prebuiltVals{kk};

            % Set augmented rates: base classes + tagged class
            augRates = [rates, rates(pb.tagClass)];

            % Update service rates in cached sn struct (cell-of-cells format)
            sn_upd = pb.sn;
            for cc = 1:newR
                sn_upd = sn_set_service_coc(sn_upd, pb.queueStIdx, cc, augRates(cc));
            end

            % Run solver_ctmc with updated rates (reuses sn topology/state structure)
            [infGen, ~, ~, Dfilt] = solver_ctmc(sn_upd, ctmcOpts);

            % Build D1 from cached departure event indices
            D1 = sparse(size(infGen, 1), size(infGen, 2));
            for di = 1:length(pb.taggedDepIdx)
                D1 = D1 + Dfilt{pb.taggedDepIdx(di)};
            end

            % Extract sub-generator using cached subset
            MAPQ1 = infGen - D1;
            A = full(MAPQ1(pb.subset, pb.subset));

            cacheKeys{end+1} = key; %#ok<AGROW>
            cacheVals{end+1} = struct('A', A, 'SSqueue', pb.SSqueue, 'N', pb.N); %#ok<AGROW>
        end

        % Compute likelihoods using cached CTMCs
        ftemp = zeros(size(rt));
        for r = 1:length(rt)
            cacheKey = mat2str([class(r), ql(r, :)]);
            cached = cacheVals{find(strcmp(cacheKeys, cacheKey), 1)};
            ftemp(r) = log(TOL + eval_mlps_likelihood(cached.A, cached.SSqueue, cached.N, rt(r)));
        end
        f = -sum(ftemp);
    end

end


function LIKE = eval_mlps_likelihood(A, SSqueue, N, Rsam)
% Compute MLPS likelihood from pre-built CTMC components.

pie = zeros(1, length(A));
idx = matchrow(SSqueue, N);
if ~isempty(idx) && idx > 0
    pie(idx) = 1;
end

MAP = {A, -A * ones(size(pie))' * pie};
LIKE = map_pdf(MAP, Rsam);

end
