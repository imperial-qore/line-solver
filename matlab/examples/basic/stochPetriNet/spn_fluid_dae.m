% Fluid (mean-field) analysis of a stochastic Petri net, with SolverFLD.
%
% A GSPN is a density-dependent Markov population process: the marking is the
% population, a transition mode is a reaction, and the firing rate
% lambda*min(enabling degree, servers) is the same min() non-linearity the
% min-normal closure of SolverFLD exists to smooth. The 'dae' method is the one
% that can carry it, because a Petri net needs three things stated as EQUATIONS
% rather than integrated:
%
%   * the P-invariants, which hold to solver tolerance instead of integrator
%     tolerance -- and supply the rank the drift Jacobian is missing;
%   * the firing FLOW of an immediate transition, an algebraic unknown pinned by
%     the constraint that its input place holds no mass;
%   * a bounded place, a linear inequality on the marking.
%
% SolverFLD(model) resolves to 'dae' on any model holding a Transition node, so
% no method has to be named. Unlike every other solver of a Petri net in LINE it
% also returns a SECOND MOMENT: the marking covariance of the linear noise
% approximation, through getMoments().

clear model P T R jc solver

%% -- a closed net whose fluid answer is EXACT -------------------------------
% Every mode is infinite-server with one input arc, so min(m/w, Inf) = m and the
% drift is LINEAR: the fluid mean is then the exact mean, and the covariance the
% exact covariance (a binomial marking).
exact = Network('spn_fluid_exact');
P1 = Place(exact, 'P1');
P2 = Place(exact, 'P2');
T1 = Transition(exact, 'T1');
T2 = Transition(exact, 'T2');
jc = ClosedClass(exact, 'Class1', 4, P1, 0);
m1 = T1.addMode('Mode1'); T1.setNumberOfServers(m1, Inf); T1.setDistribution(m1, Exp(2));
T1.setEnablingConditions(m1, jc, P1, 1); T1.setFiringOutcome(m1, jc, P2, 1);
m2 = T2.addMode('Mode2'); T2.setNumberOfServers(m2, Inf); T2.setDistribution(m2, Exp(3));
T2.setEnablingConditions(m2, jc, P2, 1); T2.setFiringOutcome(m2, jc, P1, 1);
R = exact.initRoutingMatrix();
R.set(jc, jc, P1, T1, 1.0); R.set(jc, jc, T1, P2, 1.0);
R.set(jc, jc, P2, T2, 1.0); R.set(jc, jc, T2, P1, 1.0);
exact.link(R);
P1.setState(jc.population); P2.setState(0);

fld = SolverFLD(exact);
AvgTableFLD = fld.getAvgTable()
AvgTableCTMC = SolverCTMC(exact, 'cutoff', 6).getAvgTable()

% The exact marking is Binomial(4, 3/5), so the variance is 4*0.6*0.4 = 0.96.
mom = fld.getMoments();
fprintf('marking variance: %.6f %.6f  (exact 0.96)\n', mom.petri.markingVar(P1.index,1), ...
    mom.petri.markingVar(P2.index,1));
fprintf('invariant "%s" = %g, error %.2e\n', mom.petri.invariantLabel{1}, ...
    mom.petri.invariantValue(1), mom.petri.invariantError(1));

%% -- an immediate transition, as an algebraic flow --------------------------
% P1 -T1-> P2 -(immediate)-> P3 -T2-> P1. The vanishing place P2 holds exactly
% zero mass and the net answers as the reduced two-place net does, which is what
% the algebraic flow buys: an approximation of the immediate transition by a
% large finite rate would only approach it.
imm = Network('spn_fluid_immediate');
Q1 = Place(imm, 'P1'); Q2 = Place(imm, 'P2'); Q3 = Place(imm, 'P3');
U1 = Transition(imm, 'T1'); Ui = Transition(imm, 'Ti'); U3 = Transition(imm, 'T3');
jq = ClosedClass(imm, 'Class1', 4, Q1, 0);
a1 = U1.addMode('M1'); U1.setDistribution(a1, Exp(2));
U1.setEnablingConditions(a1, jq, Q1, 1); U1.setFiringOutcome(a1, jq, Q2, 1);
ai = Ui.addMode('Mi'); Ui.setDistribution(ai, Immediate());
Ui.setTimingStrategy(ai, TimingStrategy.IMMEDIATE);
Ui.setEnablingConditions(ai, jq, Q2, 1); Ui.setFiringOutcome(ai, jq, Q3, 1);
a3 = U3.addMode('M3'); U3.setDistribution(a3, Exp(3));
U3.setEnablingConditions(a3, jq, Q3, 1); U3.setFiringOutcome(a3, jq, Q1, 1);
Ri = imm.initRoutingMatrix();
Ri.set(jq, jq, Q1, U1, 1.0); Ri.set(jq, jq, U1, Q2, 1.0);
Ri.set(jq, jq, Q2, Ui, 1.0); Ri.set(jq, jq, Ui, Q3, 1.0);
Ri.set(jq, jq, Q3, U3, 1.0); Ri.set(jq, jq, U3, Q1, 1.0);
imm.link(Ri);
Q1.setState(jq.population); Q2.setState(0); Q3.setState(0);

AvgTableImmediate = SolverFLD(imm).getAvgTable()

% The same net with the immediate transition eliminated by hand.
red = Network('spn_fluid_reduced');
W1 = Place(red, 'P1'); W3 = Place(red, 'P3');
V1 = Transition(red, 'T1'); V3 = Transition(red, 'T3');
jr = ClosedClass(red, 'Class1', 4, W1, 0);
b1 = V1.addMode('M1'); V1.setDistribution(b1, Exp(2));
V1.setEnablingConditions(b1, jr, W1, 1); V1.setFiringOutcome(b1, jr, W3, 1);
b3 = V3.addMode('M3'); V3.setDistribution(b3, Exp(3));
V3.setEnablingConditions(b3, jr, W3, 1); V3.setFiringOutcome(b3, jr, W1, 1);
Rr = red.initRoutingMatrix();
Rr.set(jr, jr, W1, V1, 1.0); Rr.set(jr, jr, V1, W3, 1.0);
Rr.set(jr, jr, W3, V3, 1.0); Rr.set(jr, jr, V3, W1, 1.0);
red.link(Rr);
W1.setState(jr.population); W3.setState(0);

AvgTableReduced = SolverFLD(red).getAvgTable()
