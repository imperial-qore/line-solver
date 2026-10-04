% Closed cycle of Bernoulli servers on a discrete time scale.
%
% N jobs circulate through J FCFS single-server stations; station j completes
% a service with probability p_j in each slot. The stationary queue length
% vector has the product form of Daduna (2001), corollary 3.4,
%
%   pi(n_1,...,n_J) = prod_j (q_j/p_j)^n_j (1/q_j)^{1{n_j>0}} / G(N,J),
%
% whose extra factor on the busy nodes is what separates it from the
% continuous-time Gordon-Newell form: a homogeneous cycle is uniform on the
% state space in continuous time and is not here. SolverNC evaluates G with
% the three-term recursion of proposition 3.18 and the arrival constants with
% proposition 3.19, and reads the marginals off corollary 3.20.

clear solver AvgTable;

p = [0.5 0.25 0.7];   % per-slot service completion probabilities
N = 5;                % jobs cycling

model = Network('BernoulliCycle');

station = cell(1, numel(p));
for j = 1:numel(p)
    station{j} = Queue(model, sprintf('Queue%d', j), SchedStrategy.FCFS);
end

jobClass = ClosedClass(model, 'Jobs', N, station{1});
for j = 1:numel(p)
    station{j}.setService(jobClass, Geometric(p(j)));
end

model.link(Network.serialRouting(station{:}));

%%
options = SolverNC.defaultOptions;
options.config.slotted = true;
solver = SolverNC(model, options);
AvgTable = solver.getAvgTable; AvgTable

% The same quantities from the normalizing constants directly. Throughput is
% G_1(N,J)/G(N,J), the discrete-time counterpart of G(N-1)/G(N), and it is the
% same at every node of the cycle.
[lG, G, G1] = dpfqn_nc(p, N);
fprintf('log G(N,J)   = %g\n', lG);
fprintf('throughput   = %g jobs per slot\n', G1(N+1)/G);
fprintf('utilizations = %s\n', mat2str(round((G1(N+1)/G) ./ p, 6)));
