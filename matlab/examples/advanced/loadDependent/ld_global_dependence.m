clear node jobclass

% Globally state-dependent (Whittle) model: setGlobalDependence declares a rate
% scaling phi(n) over the FULL (nstations x nclasses) population matrix, not
% just the population local to one station. Here two PS stations share one unit
% of capacity, phi_s(n) = n_s/|n|, which is the single-link allocation every
% alpha-fair rule collapses to.
%
% phi satisfies the Whittle balance property
%   phi_s(n) phi_t(n-e_s) = phi_t(n) phi_s(n-e_t),
% so the chain is reversible, has the product form pi(n) ~ Phi(n) prod rho^n and
% is INSENSITIVE: the means below do not change when the exponential service is
% replaced by an Erlang or a hyperexponential of the same mean.
N = 3; % number of jobs
%%
gdmodel = Network('model');
node{1} = Queue(gdmodel, 'Queue1', SchedStrategy.PS);
node{2} = Queue(gdmodel, 'Queue2', SchedStrategy.PS);
jobclass{1} = ClosedClass(gdmodel, 'Class1', N, node{1}, 0);
node{1}.setService(jobclass{1}, Exp.fitMean(1.0));
node{2}.setService(jobclass{1}, Exp.fitMean(0.5));
gdmodel.link(gdmodel.serialRouting(node));

sn = gdmodel.getStruct();
i1 = sn.nodeToStation(gdmodel.getNodeIndex('Queue1'));
i2 = sn.nodeToStation(gdmodel.getNodeIndex('Queue2'));
M = sn.nstations; K = sn.nclasses;
% peak scaling is 1: no station ever receives more than the whole link
gdmodel.setGlobalDependence(@(n) share(n, i1, i2, M, K), 1.0);

% Only SolverCTMC and SolverSSA plumb the handle; every other solver rejects
% the model rather than silently solving it unscaled (feature
% 'GlobalDependence'). SolverSSA carries the SAME factorization on the sample
% path -- phi(n) is a constant within a state, so it is evaluated once per state
% and multiplies every station service rate there. Its NRM engine cannot (its
% propensity closures see one station's population slice), so a model declaring
% a global dependence is routed to the serial engine.
gdAvgTableCTMC = CTMC(gdmodel,'exact').getAvgTable
gdAvgTableSSA = SSA(gdmodel,'seed',23000,'samples',2e5).getAvgTable

% The balance property is checkable, and is what separates a Whittle network
% from an arbitrary state-dependent rate.
viol = sn_gd_balance(@(n) shareVec(n), [4;4]);
line_printf('\nworst relative balance violation: %.2e (balanced when ~0)\n', viol);

model = gdmodel; % for test compatibility

function v = share(n, i1, i2, M, K)
% phi as SolverCTMC sees it: an (nstations x nclasses) matrix of scalings
v = ones(M,K);
tot = sum(n(i1,:)) + sum(n(i2,:));
if tot > 0
    v(i1,:) = sum(n(i1,:)) / tot;
    v(i2,:) = sum(n(i2,:)) / tot;
end
end

function x = shareVec(n)
% the same allocation in the per-station form sn_gd_balance expects
tot = sum(n);
if tot == 0
    x = zeros(numel(n),1);
else
    x = n(:) / tot;
end
end
