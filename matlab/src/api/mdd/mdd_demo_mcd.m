function mdd_demo_mcd()
% MDD_DEMO_MCD
% Demonstration and validation of the Miner-Ciardo-Donatelli approximate
% aggregation (mdd_mcd) on single-class closed exponential queueing networks.
% For product-form models the method is EXACT (SIGMETRICS 2000, Sec. 5), so it
% must reproduce the exact solution (mdd_closedqn / SolverCTMC) to machine
% precision while solving only K small level-CTMCs instead of the full chain.
%
% See also: mdd_mcd, mdd_descriptor, mdd_closedqn.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

% each case: {name, mu, servers, N, routing}
Pgen4 = [0 0.6 0.4 0; 0 0 0.5 0.5; 0.5 0 0 0.5; 0.7 0.3 0 0];  % general, irreducible
cases = {
    'K=3 cyclic, 1 delay', [0.8 4 3],        [Inf 1 1],   8,  cyclic(3)
    'K=5 cyclic',          [4 3 2 5 3],       [1 1 1 1 1], 10, cyclic(5)
    'K=4 general routing',  [3 4 2 5],         [1 1 1 1],   12, Pgen4
    'K=4 general, 2 delay', [1 4 2 5],         [Inf 1 1 1], 12, Pgen4
    };

line_printf('\n%-24s %8s %8s %6s %6s %10s %10s\n', ...
    'case', '|S|', 'maxMk', 'iters', 'match', 'errQLen', 'errTput');
for c = 1:size(cases, 1)
    nm = cases{c, 1}; mu = cases{c, 2}; servers = cases{c, 3};
    N = cases{c, 4}; P = cases{c, 5};

    oute = mdd_closedqn(mu, P, servers, N);          % exact reference
    mdds = oute.mdd.toStruct();
    desc = mdd_descriptor(mu, P, servers, N);
    outm = mdd_mcd(mdds, desc, struct());            % approximate aggregation

    eQ = max(abs(oute.QLen - outm.QLen));
    eX = max(abs(oute.X - outm.X));
    ok = eQ < 1e-8 && eX < 1e-8;
    line_printf('%-24s %8d %8d %6d %6s %10.2e %10.2e\n', ...
        nm, oute.stats.numStates, max(outm.levelSizes), outm.iters, ...
        ternary(ok, 'OK', 'FAIL'), eQ, eX);
end

% extra cross-check of case 1 against exact SolverMVA
line_printf('\ncross-check (case 1) mean queue lengths:\n');
mu = cases{1, 2}; servers = cases{1, 3}; N = cases{1, 4}; K = numel(mu);
model = Network('mcd_check');
st = cell(1, K); st{1} = Delay(model, 'Think');
for i = 2:K, st{i} = Queue(model, sprintf('Q%d', i), SchedStrategy.PS); st{i}.setNumberOfServers(1); end
job = ClosedClass(model, 'Jobs', N, st{1});
for i = 1:K, st{i}.setService(job, Exp(mu(i))); end
model.link(Network.serialRouting(st{:}));
mvaQ = SolverMVA(model).avgTable().QLen(:)';
oute = mdd_closedqn(mu, cyclic(K), servers, N);
outm = mdd_mcd(oute.mdd.toStruct(), mdd_descriptor(mu, cyclic(K), servers, N), struct());
line_printf('  MDD-MCD : %s\n', mat2str(outm.QLen, 5));
line_printf('  MVA     : %s\n', mat2str(mvaQ, 5));
line_printf('  max |diff vs MVA| = %.2e\n', max(abs(outm.QLen - mvaQ)));
end

function P = cyclic(K)
P = zeros(K);
for i = 1:K, P(i, mod(i, K) + 1) = 1; end
end

function s = ternary(c, a, b)
if c, s = a; else, s = b; end
end
