function pas_nc_sampling()
% PAS_NC_SAMPLING  Importance-sampling normalizing-constant analysis of a closed
% pass-and-swap (P&S) tandem via SolverNC method='sampling' (Casale, Comte &
% Dorsman, 2026). The model is the Figure 6 closed tandem of Comte & Dorsman
% (2021, arXiv:2009.12299) also used in PAS_CLOSED_TANDEM_FIG6.
%
% A non-empty swap graph makes the ordered-state chain reducible; the recurrent
% communicating class carries the per-class product form pi(c)=Phi_1 Phi_2/G_C.
% SolverNC importance sampling estimates G_C and the mean queue lengths by auto-normalized
% importance sampling (API PFQN_PAS_IS), scaling to populations where the exact
% ordered-state CTMC is intractable. Here we validate it against the exact CTMC.

mu1 = 1.0;  mu2 = 1.3;
edges = [1 3; 1 4; 2 4; 2 5; 3 6; 4 6; 5 6];
G = zeros(6);
for k = 1:size(edges,1)
    G(edges(k,1),edges(k,2)) = 1;
    G(edges(k,2),edges(k,1)) = 1;
end

model = Network('PASsampling');
q1 = Queue(model, 'PASQueue1', SchedStrategy.PAS);
q2 = Queue(model, 'PASQueue2', SchedStrategy.PAS);
jobclass = cell(1,6);
for r = 1:6
    jobclass{r} = ClosedClass(model, sprintf('Class%d', r), 1, q1);
end
q1.setService(@(c) mu1);   q2.setService(@(c) mu2);   % head-only single server
q1.setSwapGraph(G);  q1.setNumberOfServers(1);  q1.setCap(6);
q2.setSwapGraph(G);  q2.setNumberOfServers(1);  q2.setCap(6);
P = model.initRoutingMatrix;
for r = 1:6
    P{jobclass{r}}(q1, q2) = 1.0;
    P{jobclass{r}}(q2, q1) = 1.0;
end
model.link(P);
q1.setState([1 2 3 4 5 6]);   % Fig. 6a initial placement (reducible model input)

fprintf('=== CTMC (exact) ===\n');
Tc = CTMC(model, 'exact', 'cutoff', 6).getAvgTable;
fprintf('=== SolverNC method=''sampling'' (importance sampling, 5e5 samples) ===\n');
Tn = NC(model, 'method', 'sampling', 'samples', 5e5, 'seed', 777, 'verbose', false).getAvgTable;

qc = Tc.QLen;  qn = Tn.QLen;
fprintf('\n  station     class    CTMC      NC-samp\n');
for i = 1:height(Tc)
    fprintf('  %-9s  %-6s %9.5f %9.5f\n', string(Tc.Station(i)), ...
        string(Tc.JobClass(i)), qc(i), qn(i));
end
err = max(abs(qn - qc));
fprintf('\nNC-samp vs CTMC:  max|dQ| = %.3e (importance-sampling noise)\n', err);
assert(err <= 5e-2, 'NC sampling does not match the exact CTMC within IS tolerance');
fprintf('PASS: SolverNC ''sampling'' matches the exact CTMC within importance-sampling noise.\n');
end
