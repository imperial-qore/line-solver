% Closed cycle of state dependent Bernoulli servers.
%
% Theorem 3.2 of Daduna (2001) keeps the product form when the service
% probability of station j depends on its own queue length, with node weight
%
%   w_j(n) = prod_{h=1}^{n-1} q_j(h) / prod_{h=1}^{n} p_j(h).
%
% Note where the state dependence sits: the missing q_j in the numerator is
% tied to the actual queue length, not to the node being non-empty, so the
% tidy (1/q_j)^{1{n>0}} of the state independent case cannot be factored out.
%
% Station 2 below runs at p_2(n) = p_2 min(n,2), the discrete-time analogue of
% adding a second server. That is admissible only because the dependence is
% expressed as a state dependent single server: a genuine multiserver node
% inside a cycle of geometrical queues destroys the product form for every
% finite server count (Pestien and Ramakrishnan, cited before example 2.10),
% and SolverNC rejects one rather than approximating it.

clear solver AvgTable;

p = [0.5 0.25 0.7];
N = 5;

model = Network('BernoulliCycleLD');

station = cell(1, numel(p));
for j = 1:numel(p)
    station{j} = Queue(model, sprintf('Queue%d', j), SchedStrategy.FCFS);
end

jobClass = ClosedClass(model, 'Jobs', N, station{1});
for j = 1:numel(p)
    station{j}.setService(jobClass, Geometric(p(j)));
end
station{2}.setLoadDependence(min(1:N, 2));

model.link(Network.serialRouting(station{:}));

%%
options = SolverNC.defaultOptions;
options.config.slotted = true;
solver = SolverNC(model, options);
AvgTable = solver.getAvgTable; AvgTable

% Marginal queue length law of station 2 from the convolution constants.
P = repmat(p(:), 1, N);
P(2, :) = p(2) * min(1:N, 2);
[lG, G, W, Gc] = dpfqn_ncld(P, N);
marg = W(2, :) .* Gc(2, N+1:-1:1) / G(N+1);
fprintf('P(X_2 = 0..%d) = %s\n', N, mat2str(round(marg, 6)));
