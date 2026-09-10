% Pass-and-swap (PAS) station as a multiserver order-independent (M/M/K) queue.
% The M/M/K queue satisfies the order-independence conditions (Dorsman &
% Gardner 2024, Sect. 2.2): with K unit-rate servers the total service rate in
% state c is mu(c) = min(n, K), where n is the number of jobs. With an empty
% swapping graph the PAS station reduces to a plain order-independent queue.
clear node jobclass

K = 2;                          % number of servers
muFun = @(c) min(numel(c(c>0)), K);

model = Network('PASmmk');
source = Source(model, 'Source');
queue  = Queue(model, 'PASQueue', SchedStrategy.PAS);
sink   = Sink(model, 'Sink');

lambda = [0.7 0.5];
jobclass = cell(1,2);
for r = 1:2
    jobclass{r} = OpenClass(model, sprintf('Class%d', r));
    source.setArrival(jobclass{r}, Exp(lambda(r)));
end

queue.setService(muFun);
queue.setSwapGraph(zeros(2));   % empty graph: plain order-independent queue
queue.setNumberOfServers(K);    % K servers (used for utilization reporting)
queue.setCap(4);

P = model.initRoutingMatrix;
for r = 1:2
    P{jobclass{r}} = Network.serialRouting(source, queue, sink);
end
model.link(P);

AvgTable = CTMC(model, 'exact', 'cutoff', 4).getAvgTable
