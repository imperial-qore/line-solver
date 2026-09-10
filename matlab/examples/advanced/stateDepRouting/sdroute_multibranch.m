% Multi-center branches under Krzesinski (1987) state-dependent routing.
%
% Checks that the product form of eq. (16) extends to a branch holding several
% centers, and that the coefficients xi are the branch traffic equations rather
% than the paper Section 3.2 shorthand xi = xi_e, which is exact only when the
% branch departure center is visited once. See _kb/16-state-dependent-routing.md
clear; lineStart;

fprintf('\n==== SDR with multi-center branches ====\n');
run_case('A: branch 2 = 2a -> 2b (series)', 0.0);
run_case('B: branch 2 = 2a -> 2b, 2b -> 2a w.p. 0.5 (feedback onto the branch departure)', 0.5);

function run_case(label, pback)
mu = [1 0.9 0.7 0.5];   % centers 1, 2a, 2b, 3
N  = 3;
model = Network('sdr_multi');
nd{1} = Queue(model, 'CPU',  SchedStrategy.FCFS);
nd{2} = Queue(model, 'B2a',  SchedStrategy.FCFS);
nd{3} = Queue(model, 'B2b',  SchedStrategy.FCFS);
nd{4} = Queue(model, 'B3',   SchedStrategy.FCFS);
cls = ClosedClass(model, 'Class1', N, nd{1}, 0);
for i = 1:4, nd{i}.setService(cls, Exp(mu(i))); end
model.addLink(nd{1}, nd{1});
model.addLink(nd{1}, nd{2});
model.addLink(nd{1}, nd{4});
model.addLink(nd{2}, nd{3});
model.addLink(nd{3}, nd{1});
model.addLink(nd{4}, nd{1});
nd{2}.setProbRouting(cls, nd{3}, 1.0);
if pback > 0
    model.addLink(nd{3}, nd{2});
    nd{3}.setProbRouting(cls, nd{2}, pback);
    nd{3}.setProbRouting(cls, nd{1}, 1-pback);
else
    nd{3}.setProbRouting(cls, nd{1}, 1.0);
end
nd{4}.setProbRouting(cls, nd{1}, 1.0);

d = zeros(2,3); d(1,2) = 2; d(1,3) = 2; d(2,3) = 2;
nd{1}.setStateDepRouting(cls, nd{1}, {[], {nd{2}, nd{3}}, {nd{4}}}, [0 1 2], [-1 -1], d);

sn = model.getStruct();
P = zeros(4,4);
P(1,1) = 1;                       % complement M-V = {1}
P(2,3) = 1;                       % inside branch 2
if pback > 0, P(3,2) = pback; end
xi = pfqn_sdrvisits(sn.sdr, P);
S = (1./mu)';
[Q, X] = pfqn_sdr(S, xi, N, sn.sdr);

ctmc = SolverCTMC(model);
Qc = ctmc.getAvgQLen(); Xc = ctmc.getAvgTput();

fprintf('\n  %s\n', label);
fprintf('    xi          = %s\n', mat2str(xi',5));
fprintf('    pfqn_sdr Q  = %s\n', mat2str(Q',6));
fprintf('    CTMC     Q  = %s\n', mat2str(Qc',6));
fprintf('    pfqn_sdr X  = %s\n', mat2str(X',6));
fprintf('    CTMC     X  = %s\n', mat2str(Xc',6));
fprintf('    max|dQ| = %.3e   max|dX| = %.3e\n', max(abs(Qc-Q)), max(abs(Xc-X)));

% the reading asserted verbatim by the paper for E+D: xi = xi_e everywhere
xiflat = ones(4,1);
[Q2, X2] = pfqn_sdr(S, xiflat, N, sn.sdr);
fprintf('    xi==1 everywhere: max|dQ| = %.3e   max|dX| = %.3e\n', max(abs(Qc-Q2)), max(abs(Xc-X2)));
end
