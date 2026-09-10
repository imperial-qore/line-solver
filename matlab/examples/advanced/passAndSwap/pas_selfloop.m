% Pass-and-swap (PAS) queue with a swapping graph that includes a self-loop.
% Per Dorsman & Gardner (2024), Sect. 2.3, self-loops in the swapping graph are
% permissible: a completing class-i job may take the place of another class-i
% job further back in the queue. The order-independent product form (Theorem 2)
% still holds and remains invariant to the swapping graph.
clear node jobclass

% Three unit-rate servers, one per class (an M/M/1-per-class compatibility):
% mu(c) = number of distinct classes present in c.
comp = eye(3);
muFun = @(c) sum( any(comp(:, c(c>0)), 2) );

model = Network('PASselfloop');
source = Source(model, 'Source');
queue  = Queue(model, 'PASQueue', SchedStrategy.PAS);
sink   = Sink(model, 'Sink');

lambda = [0.6 0.4 0.3];
jobclass = cell(1,3);
for r = 1:3
    jobclass{r} = OpenClass(model, sprintf('Class%d', r));
    source.setArrival(jobclass{r}, Exp(lambda(r)));
end

queue.setService(muFun);
% Swapping graph with a self-loop on class 1 and an edge (2,3)
G = zeros(3);
G(1,1) = 1;            % self-loop: class 1 swaps with class 1
G(2,3) = 1; G(3,2) = 1;
queue.setSwapGraph(G);
queue.setNumberOfServers(3); % three servers (used for utilization reporting)
queue.setCap(3);

P = model.initRoutingMatrix;
for r = 1:3
    P{jobclass{r}} = Network.serialRouting(source, queue, sink);
end
model.link(P);

AvgTable = CTMC(model, 'exact', 'cutoff', 3).getAvgTable
