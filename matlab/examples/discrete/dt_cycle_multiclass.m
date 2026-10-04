% Multichain closed cycle of Bernoulli servers.
%
% Section 3.2 of Daduna (2001) splits the circulating jobs into chains that
% never mix. Service in the cycle is type independent and FCFS forbids
% overtaking, so the cyclic order of the jobs is frozen for all time and the
% joint queue length law is the unichain one at the aggregate population. Each
% chain then holds a share of every station equal to its share of the
% population, which is what SolverNC reports below: the per-station totals
% reproduce dt_cycle at N = N_1 + N_2, and the per-class rows split them in the
% ratio N_1 : N_2.

clear solver AvgTable;

p = [0.5 0.25 0.7];
N1 = 3;    % jobs in chain 1
N2 = 2;    % jobs in chain 2

model = Network('BernoulliCycleMC');

station = cell(1, numel(p));
for j = 1:numel(p)
    station{j} = Queue(model, sprintf('Queue%d', j), SchedStrategy.FCFS);
end

class1 = ClosedClass(model, 'Chain1', N1, station{1});
class2 = ClosedClass(model, 'Chain2', N2, station{1});
for j = 1:numel(p)
    station{j}.setService(class1, Geometric(p(j)));
    station{j}.setService(class2, Geometric(p(j)));
end

% Both chains follow the same cycle and never switch class.
J = numel(p);
P = model.initRoutingMatrix;
P{class1, class1} = circul(J);
P{class1, class2} = zeros(J);
P{class2, class1} = zeros(J);
P{class2, class2} = circul(J);
model.link(P);

%%
options = SolverNC.defaultOptions;
options.config.slotted = true;
solver = SolverNC(model, options);
AvgTable = solver.getAvgTable; AvgTable

[lG, G, G1] = dpfqn_nc(p, N1 + N2);
X = G1(N1 + N2 + 1) / G;
fprintf('aggregate throughput  = %g jobs per slot\n', X);
fprintf('per-chain throughput  = %g and %g\n', X*N1/(N1+N2), X*N2/(N1+N2));
