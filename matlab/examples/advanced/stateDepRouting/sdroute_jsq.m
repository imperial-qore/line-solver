clear node jobclass solver AvgTable

model = Network('myModel');

N = 3;    % number of queues
rho = 0.7; % per-queue load
mu = 1;    % service rate per queue
lambda = N * rho * mu; % total arrival rate

% Block 1: nodes
source = Source(model, 'Source');
router = Router(model, 'Router');
queue = cell(1, N);
for i = 1:N
    queue{i} = Queue(model, sprintf('Queue%d', i), SchedStrategy.FCFS);
end
sink = Sink(model, 'Sink');

% Block 2: classes
oclass = OpenClass(model, 'Class1');
source.setArrival(oclass, Exp(lambda));
for i = 1:N
    queue{i}.setService(oclass, Exp(mu));
end

% Block 3: topology
model.addLink(source, router);
for i = 1:N
    model.addLink(router, queue{i});
    model.addLink(queue{i}, sink);
end

router.setRouting(oclass, RoutingStrategy.JSQ);

solver = {};
solver{end+1} = JMT(model, 'seed', 23000);
solver{end+1} = LDES(model, 'seed', 23000);

AvgTable = {};
for s = 1:length(solver)
    fprintf(1, 'SOLVER: %s\n', strrep(solver{s}.getName(), 'Solver', ''));
    AvgTable{s} = solver{s}.getAvgNodeTable();
    AvgTable{s}
end
