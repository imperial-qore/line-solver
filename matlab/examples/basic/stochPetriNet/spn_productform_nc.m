% SPN_PRODUCTFORM_NC  Solve a stochastic Petri net analytically with SolverNC.
%
% SolverNC's 'rec' method is the first ANALYTICAL route LINE offers for a Petri
% net: SolverCTMC builds the explicit generator, SolverSSA and SolverLDES
% simulate, SolverFLD fluidises. It works in three steps, each with its own
% reference:
%
%   SPN_PF      decides whether the net has a product form and derives the
%               per-place factors g_l, by complex balance
%               (Coleman-Henderson-Taylor, Perform. Eval. 26(3), 1996)
%   MDD_REC     evaluates G = sum over the reachable set of prod_l g_l(m_l) by
%               ONE memoised walk of the decision diagram holding that set
%               (Balsamo-Marin-Stojic, FGCS 111 (2020) 475-490)
%   SPN_METRICS reads the mean tokens, the utilisations and the throughputs off
%               masked walks of the same diagram
%
% Two nets are solved here. The first is a closed cycle, which a queueing
% network could also express. The second FORKS: its transition Tf consumes one
% token and produces two, and Tj consumes two and produces one, so the marking
% is not a conserved job population and there is no queueing-network
% counterpart -- which is exactly the limitation the MDD-rec paper opens with.

%% A closed cycle, where the exact CTMC gives the reference
model = i_cyclic(4);
nc = SolverNC(model);
line_printf('\n--- 3-place cyclic net, N = 4 ---\n');
nc.getAvgTable()

line_printf('The same net through the explicit generator, for comparison:\n');
SolverCTMC(i_cyclic(4)).getAvgTable()

% The certificate the derivation produced. Deficiency zero plus weak
% reversibility is what Feinberg's theorem needs for a positive complex-balanced
% point to exist at ANY choice of rates.
pf = spn_pf(i_cyclic(4), struct('verbose', true));
line_printf('product form: %s, deficiency %d, %d linkage classes, rank %d\n', ...
    pf.kind, pf.deficiency, pf.linkage, pf.srank);

%% A fork-join net, which has no queueing-network form at all
fj = i_forkjoin(3);
line_printf('\n--- fork-join net, P0 -> P1+P2 -> P3 -> P0, 3 tokens at P0 ---\n');
SolverNC(fj).getAvgTable()

pfj = spn_pf(i_forkjoin(3));
met = spn_metrics(pfj.mdds, pfj.g, pfj.info);
% Every token that forks must later join and return, so the three modes share
% one throughput -- a flow-conservation law nothing in the derivation was told.
line_printf('mode throughputs: %s (they must all agree)\n', mat2str(met.modeTput, 8));
% The place invariant of this net is 2*m0 + m1 + m2 + 2*m3 = 6, not the token
% count, which is what "not a conserved population" means concretely.
line_printf('place invariant 2*m0 + m1 + m2 + 2*m3 = %.6f\n', ...
    [2 1 1 2] * met.tokens(:));

function model = i_cyclic(njobs)
rates = [1, 1.5, 2];
model = Network('spn');
P = cell(1,3); T = cell(1,3);
for i = 1:3, P{i} = Place(model, sprintf('P%d', i-1)); end
for i = 1:3, T{i} = Transition(model, sprintf('T%d', i-1)); end
cls = ClosedClass(model, 'Class1', njobs, P{1});
for i = 1:3
    m = T{i}.addMode('fire');
    T{i}.setDistribution(m, Exp(rates(i)));
    T{i}.setNumberOfServers(m, 1);
    T{i}.setEnablingConditions(m, cls, P{i}, 1);
    T{i}.setFiringOutcome(m, cls, P{mod(i,3)+1}, 1);
end
R = model.initRoutingMatrix();
for i = 1:3
    R{1,1}(P{i}, T{i}) = 1;
    R{1,1}(T{i}, P{mod(i,3)+1}) = 1;
end
model.link(R);
P{1}.setState(njobs); P{2}.setState(0); P{3}.setState(0);
end

function model = i_forkjoin(njobs)
model = Network('fj');
P = cell(1,4);
for i = 1:4, P{i} = Place(model, sprintf('P%d', i-1)); end
Tf = Transition(model, 'Tf'); Tj = Transition(model, 'Tj'); Tb = Transition(model, 'Tb');
cls = ClosedClass(model, 'C', njobs, P{1});
m = Tf.addMode('f'); Tf.setDistribution(m, Exp(1.3)); Tf.setNumberOfServers(m, 1);
Tf.setEnablingConditions(m, cls, P{1}, 1);
Tf.setFiringOutcome(m, cls, P{2}, 1); Tf.setFiringOutcome(m, cls, P{3}, 1);
m = Tj.addMode('j'); Tj.setDistribution(m, Exp(0.7)); Tj.setNumberOfServers(m, 1);
Tj.setEnablingConditions(m, cls, P{2}, 1); Tj.setEnablingConditions(m, cls, P{3}, 1);
Tj.setFiringOutcome(m, cls, P{4}, 1);
m = Tb.addMode('b'); Tb.setDistribution(m, Exp(1.9)); Tb.setNumberOfServers(m, 1);
Tb.setEnablingConditions(m, cls, P{4}, 1);
Tb.setFiringOutcome(m, cls, P{1}, 1);
R = model.initRoutingMatrix();
R{1,1}(P{1}, Tf) = 1; R{1,1}(Tf, P{2}) = 1; R{1,1}(Tf, P{3}) = 1;
R{1,1}(P{2}, Tj) = 1; R{1,1}(P{3}, Tj) = 1; R{1,1}(Tj, P{4}) = 1;
R{1,1}(P{4}, Tb) = 1; R{1,1}(Tb, P{1}) = 1;
model.link(R);
for i = 1:4, P{i}.setState((i == 1) * njobs); end
end
