% Debug script for tut10 - compare layer results with Python
clear;
model = LayeredNetwork('ClientDBSystem');

% Create processors
P1 = Processor(model, 'ClientProcessor', 1, SchedStrategy.PS);
P2 = Processor(model, 'DBProcessor', 1, SchedStrategy.PS);

% Create tasks
T1 = Task(model, 'ClientTask', 10, SchedStrategy.REF).on(P1);
T1.setThinkTime(Exp.fitMean(5.0));
T2 = Task(model, 'DBTask', Inf, SchedStrategy.INF).on(P2);

% Create entries
E1 = Entry(model, 'ClientEntry').on(T1);
E2 = Entry(model, 'DBEntry').on(T2);

% Define activities
A1 = Activity(model, 'ClientActivity', Exp.fitMean(1.0)).on(T1);
A1.boundTo(E1).synchCall(E2, 2.5);

A2 = Activity(model, 'DBActivity', Exp.fitMean(0.8)).on(T2);
A2.boundTo(E2).repliesTo(E2);

% Create solver and build layers
solver = LN(model, @(m) MVA(m));
solver.buildLayers();

fprintf('\n=== Number of layers: %d ===\n', solver.nlayers);
fprintf('=== Ensemble size: %d ===\n\n', length(solver.ensemble));

for e = 1:length(solver.ensemble)
    layer = solver.ensemble{e};
    if isempty(layer)
        continue;
    end
    fprintf('\n--- Layer %d: %s ---\n', e-1, layer.name);
    fprintf('  clientIdx: %d\n', layer.attribute.clientIdx);
    fprintf('  serverIdx: %d\n', layer.attribute.serverIdx);
    fprintf('  ishost: %d\n', isfield(layer.attribute, 'ishost') && layer.attribute.ishost);

    nodes = layer.nodes;
    classes = layer.classes;

    fprintf('\n  Nodes (%d):\n', length(nodes));
    for i = 1:length(nodes)
        fprintf('    [%d] %s (%s)\n', i, nodes{i}.name, class(nodes{i}));
    end

    fprintf('\n  Classes (%d):\n', length(classes));
    for j = 1:length(classes)
        cls = classes{j};
        if isa(cls, 'ClosedClass')
            pop = cls.population;
        else
            pop = NaN;
        end
        fprintf('    %s: pop=%g\n', cls.name, pop);
    end
end

% Run one iteration
fprintf('\n\n=== Running solver iteration 1 ===\n');
solver.iterate();

% Print results for each layer
for e = 1:length(solver.ensemble)
    layer = solver.ensemble{e};
    if isempty(layer)
        continue;
    end
    result = solver.results{end, e};
    if isempty(result)
        continue;
    end

    fprintf('\n--- Layer %d (%s) results ---\n', e-1, layer.name);

    RN = result.RN;
    WN = result.WN;

    nodes = layer.nodes;
    classes = layer.classes;

    fprintf('  RN (response times per visit):\n');
    for i = 1:size(RN, 1)
        for j = 1:size(RN, 2)
            if RN(i,j) > 0
                fprintf('    %s, %s: RN=%.6f\n', nodes{i}.name, classes{j}.name, RN(i,j));
            end
        end
    end

    fprintf('  WN (residence times):\n');
    for i = 1:size(WN, 1)
        for j = 1:size(WN, 2)
            if WN(i,j) > 0
                fprintf('    %s, %s: WN=%.6f\n', nodes{i}.name, classes{j}.name, WN(i,j));
            end
        end
    end
end

% Print servt/residt
fprintf('\n\n=== Activity servt/residt ===\n');
lqn = solver.lqn;
for a = 1:lqn.nacts
    aidx = lqn.ashift + a;
    if solver.servt(aidx) > 0 || solver.residt(aidx) > 0
        fprintf('  Activity %d (aidx=%d): servt=%.6f, residt=%.6f\n', a, aidx, solver.servt(aidx), solver.residt(aidx));
    end
end

fprintf('\n=== Call servt/residt ===\n');
for cidx = 1:lqn.ncalls
    if solver.callservt(cidx) > 0 || solver.callresidt(cidx) > 0
        fprintf('  Call %d: callservt=%.6f, callresidt=%.6f\n', cidx, solver.callservt(cidx), solver.callresidt(cidx));
    end
end

fprintf('\n=== Entry servt ===\n');
for e = 1:lqn.nentries
    eidx = lqn.eshift + e;
    if solver.entry_servt(eidx) > 0
        fprintf('  Entry %d (eidx=%d): entry_servt=%.6f\n', e, eidx, solver.entry_servt(eidx));
    end
end
