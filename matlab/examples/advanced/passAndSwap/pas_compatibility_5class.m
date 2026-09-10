% Pass-and-swap (PAS) queue: an order-independent queue with a swapping graph.
% Reproduces the five-class, three-server compatibility example of
% Dorsman & Gardner (2024), "New directions in pass-and-swap queues",
% Queueing Systems 107:205-256 (Figs 1-2). A PAS station is parameterized by
% the total service-rate function mu(c) of the ordered state vector c and by a
% swapping graph G; both are properties of the Queue object. The stationary
% distribution is the order-independent product form and is invariant to G.
clear node jobclass

% Server compatibility (Fig 1): three unit-rate servers over five classes.
% mu(c) = number of servers compatible with at least one class present in c.
comp = [1 0 0 1 0;    % server 1: classes {1,4}
        0 1 0 1 0;    % server 2: classes {2,4}
        0 0 1 0 1];   % server 3: classes {3,5}
muFun = @(c) sum( any(comp(:, c(c>0)), 2) );

model = Network('PAScompatibility');
source = Source(model, 'Source');
queue  = Queue(model, 'PASQueue', SchedStrategy.PAS);
sink   = Sink(model, 'Sink');

lambda = [0.5 0.4 0.3 0.2 0.1];
jobclass = cell(1,5);
for r = 1:5
    jobclass{r} = OpenClass(model, sprintf('Class%d', r));
    source.setArrival(jobclass{r}, Exp(lambda(r)));
end

% PAS service: specify the rate function mu(c) as a whole (no per-class rates)
queue.setService(muFun);
% Swapping graph (Fig 2a): E = {(1,3),(1,5),(2,4),(3,4),(4,5)}
G = zeros(5);
for e = [1 3; 1 5; 2 4; 3 4; 4 5]'
    G(e(1), e(2)) = 1; G(e(2), e(1)) = 1;
end
queue.setSwapGraph(G);
queue.setNumberOfServers(3); % three servers (used for utilization reporting)
queue.setCap(3);             % finite buffer (order-independent loss model)

P = model.initRoutingMatrix;
for r = 1:5
    P{jobclass{r}} = Network.serialRouting(source, queue, sink);
end
model.link(P);

AvgTable = CTMC(model, 'exact', 'cutoff', 3).getAvgTable
