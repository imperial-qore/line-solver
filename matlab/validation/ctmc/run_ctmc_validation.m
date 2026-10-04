%% LINE against the CTMC case studies of an external model checker.
%
% The reference values are the external checker's, computed from the models it
% distributes, through its explicit engine at a termination tolerance of 1e-12
% and confirmed against its sparse engine with backwards Gauss-Seidel at 1e-14.
% They are goldens here rather than recomputed, so that this script needs only
% LINE; the reference query is recorded beside each block.
%
% What this validates is LINE's stochastic Petri net semantics and its CTMC
% solver against an independent tool, on models neither was written for.
%
% HISTORICAL NOTE, superseded 2026-07-21. This harness originally reported a
% systematic deviation, not noise: LINE latched a transition's servers with an
% ENABLE event emitted at GlobalConstants.Immediate, those latch states were
% never eliminated, and every tangible marginal therefore carried a relative
% error of order (largest rate)/GlobalConstants.Immediate. Raising the constant
% moved the deviation in exact proportion, which is what identified the
% mechanism:
%
%   Immediate = 1e8 (default)   E[s1] = 0.401595740018   rel dev 8.06e-07
%   Immediate = 1e10            E[s1] = 0.401595419517   rel dev 8.06e-09
%   Immediate = 1e12            E[s1] = 0.401595416311   rel dev 8.06e-11
%
% Purging the timed arcs out of vanishing rows removed the bias rather than
% shrinking it. Agreement is now at machine precision and no longer depends on
% GlobalConstants.Immediate:
%
%   tandem, c = 2    worst relative deviation 3.92e-14
%   polling, N = 2   worst relative deviation 1.28e-15
%
% The tolerances were correspondingly tightened on 2026-07-21 from the pre-fix
% 1e-7 and 1e-5, which had been sized to the bias and would no longer have caught
% its reintroduction, to 1e-12 and 1e-13. Each is roughly one to two orders of
% magnitude above the residual measured above, which absorbs the pivoting-order
% variation of the linear solve while still failing on any systematic error.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

addpath(fileparts(mfilename('fullpath')));

fprintf('\n=== Tandem queueing network, c = 2 ===\n');
% Reference: tandem case study at c = 2, with rewards
%   "customers" true : sc+sm;  "E_sc" true : sc;  "E_sm" true : sm;
%   "tput" [route] true : 1;   all queried as R{...}=? [ S ]
ref = struct('E_sc', 1.773522148312772, 'E_sm', 0.5126824318317403, ...
             'customers', 2.286204580144513, 'tput', 1.5739808687725572);
t = SolverCTMC(ctmc_tandem(2)).getAvgTable();
qlen = @(n) t.QLen(strcmp(string(t.Station), n));
tput = @(n) t.Tput(strcmp(string(t.Station), n));
got = struct('E_sc', qlen('sc'), 'E_sm', qlen('sm'), ...
             'customers', qlen('sc') + qlen('sm'), 'tput', tput('sm'));
report(ref, got, 1e-12);   % measured worst residual 3.92e-14

fprintf('\n=== Cyclic server polling, N = 2 ===\n');
% Reference: polling case study at N = 2, with rewards "E_s1" true : s1;
%   "E_s2" true : s2;
%   and properties S=? [ s1=1 ], S=? [ s=1 & a=1 ].
% P(a=1) is twice P(s=1 & a=1) by the symmetry of the two stations.
ref = struct('E_s1', 0.4015954162791172, 'E_s2', 0.40159541627911877, ...
             'P_serving', 2 * 0.29920229186044095, 'P_at_station1', 0.5);
t = SolverCTMC(ctmc_polling()).getAvgTable();
qlen = @(n) t.QLen(strcmp(string(t.Station), n));
got = struct('E_s1', qlen('Q1'), 'E_s2', qlen('Q2'), ...
             'P_serving', qlen('A1'), 'P_at_station1', qlen('S1'));
report(ref, got, 1e-13);   % measured worst residual 1.28e-15

function report(ref, got, tol)
names = fieldnames(ref);
worst = 0;
for i = 1:numel(names)
    a = ref.(names{i});
    b = got.(names{i});
    dev = abs(a - b) / max(1, abs(a));
    worst = max(worst, dev);
    fprintf('  %-14s REF %.12f   LINE %.12f   rel dev %.2e\n', names{i}, a, b, dev);
end
fprintf('  worst relative deviation %.2e\n', worst);
assert(worst < tol, sprintf('worst relative deviation %.2e exceeds %.1e', worst, tol));
end
