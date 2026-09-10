% SPN_LPBOUNDS  Bound a stochastic Petri net by linear programming.
%
% SolverBA's 'spnlp' family is the first BOUNDING route LINE offers for a Petri
% net. SolverCTMC builds the explicit generator, SolverSSA and SolverLDES
% simulate, SolverFLD fluidises, and SolverNC 'rec' needs a product form; this
% one needs none of that. It relaxes the stationary chain to a MOMENT POLYTOPE
% -- the uniformized evolution equation written for E[X_p], E[X_p^2] and
% E[X_p1 X_p2], plus behavioural and probabilistic inequalities -- and then
% minimises and maximises each reported measure over it. Every stationary point
% of the true chain satisfies every row, so the two optima bracket the exact
% value.
%
%   SPN_LPBND    assembles the polytope and solves the LPs
%   SPN_SINVARIANTS  supplies the weighted place invariants, which are both the
%                    conservation equalities and the a priori per-place caps
%
% Reference: Z. Liu, "Performance Analysis of Stochastic Timed Petri Nets Using
% Linear Programming Approach", IEEE Trans. Software Engineering 24(11), 1998,
% 1014-1030.

%% The reference's own Fig. 2b, and its Table 2
%
% Four servers in a line, blocking before service: server i cannot start until
% the downstream buffer has a free slot. The buffers hold 3, 2 and 4, and each
% is a conserved pair of places -- (p5,p2), (p4,p1), (p3,p0) -- so the net is a
% strongly connected marked graph and all four transitions carry the same
% throughput.
%
% Table 2 of the paper reports four bound columns on five rate vectors. All
% four are reproduced below to the three decimals it prints. The one column not
% reproduced is its u.b.1, which is the upper side further tightened by the
% subnet-throughput theorems (its Thms 1 and 2); those are not implemented, so
% u.b.2 is the column to compare against.
mus = [1 1.25 2 0.5; 1 1.25 2 2.5; 1 1.25 1.25 2.5; 1 1.25 1.25 1; ...
    1.111 1.111 1.111 1.111];
pub = [1.165 1.951 2.000 0.930 2.000
       1.829 2.978 3.529 1.481 4.000
       1.581 2.873 3.333 1.333 4.000
       1.359 2.757 3.333 1.111 4.000
       1.350 2.667 2.963 1.111 4.444];

line_printf('\nLiu (1998) Table 2: total throughput of the production line\n');
line_printf('%-5s %-32s %-19s %-19s\n', 'case', 'Markovian LP', 'published', 'operational LP');
line_printf('%-5s %9s %9s %9s   %9s %9s   %9s %9s\n', ...
    '', 'lower', 'simul', 'upper', 'l.b.', 'u.b.2', 'o.l.b.', 'o.u.b.');
for c = 1:size(mus,1)
    sn = i_prodline(mus(c,:)).getStruct();
    % The liveness rows of the reference's Table 1 are OPT-IN, because they hold
    % only on a live net and SPN_LPBND cannot certify liveness. This one is
    % live: a strongly connected marked graph with a token on every cycle. They
    % are the whole of the lower side, so the published l.b. needs them.
    bLo = spn_lpbnd(sn, struct('markovian', true, 'assumelive', true));
    bUp = spn_lpbnd(sn, struct('markovian', true));
    bOp = spn_lpbnd(sn, struct('markovian', false, 'assumelive', true));
    line_printf('%-5d %9.4f %9.4f %9.4f   %9.3f %9.3f   %9.4f %9.4f\n', c, ...
        sum(bLo.modeTput(1,:)), pub(c,2), sum(bUp.modeTput(2,:)), ...
        pub(c,1), pub(c,3), sum(bOp.modeTput(1,:)), sum(bOp.modeTput(2,:)));
end

%% Through SolverBA, on a net with an inhibitor arc
%
% The four method names are spnlp.upper, spnlp.lower and their spnlp.op.*
% counterparts, which drop the second-moment, covariance and Little's-law
% families and so need only a mean firing time rather than an exponential one.
% They are the only family SolverBA offers on a Petri net, and the only one it
% withholds off a Petri net: every other family is parameterized by demands and
% a population, which a marking is not.
model = i_inhibiting(4);
line_printf('\nmethods offered on this net: %s\n', ...
    strjoin(SolverBA(model).listValidMethods(), ', '));

lo = SolverBA(model, 'spnlp.lower').getAvgTable();
up = SolverBA(model, 'spnlp.upper').getAvgTable();
ex = SolverCTMC(i_inhibiting(4)).getAvgTable();
line_printf('\nMean tokens per place, exact between the two sides:\n');
line_printf('%-8s %10s %10s %10s\n', 'place', 'lower', 'exact', 'upper');
for i = 1:height(ex)
    line_printf('%-8s %10.5f %10.5f %10.5f\n', string(ex.Station(i)), ...
        lo.QLen(i), ex.QLen(i), up.QLen(i));
end

% U = Q at a Place, which LINE models as an INF station -- the same convention
% SolverCTMC and SolverNC report. The paper's place utilization 1 - P(m = 0) is
% a different quantity and is not this column.
line_printf('\nsolver.citations():\n');
SolverBA(i_inhibiting(4), 'spnlp.upper').citations();

%% ------------------------------------------------------------------------
function model = i_prodline(mu)
% Fig. 2b: t1 -> p5 -> t2 -> p4 -> t3 -> p3 -> t4, with p2, p1 and p0 holding
% the free slots of the three finite buffers.
model = Network('liu98');
p5 = Place(model,'p5'); p4 = Place(model,'p4'); p3 = Place(model,'p3');
p2 = Place(model,'p2'); p1 = Place(model,'p1'); p0 = Place(model,'p0');
t1 = Transition(model,'t1'); t2 = Transition(model,'t2');
t3 = Transition(model,'t3'); t4 = Transition(model,'t4');
jc = ClosedClass(model,'Class1', 9, p2, 0);
m = t1.addMode('m1'); t1.setDistribution(m, Exp(mu(1)));
t1.setEnablingConditions(m, jc, p2, 1); t1.setFiringOutcome(m, jc, p5, 1);
m = t2.addMode('m2'); t2.setDistribution(m, Exp(mu(2)));
t2.setEnablingConditions(m, jc, p5, 1); t2.setEnablingConditions(m, jc, p1, 1);
t2.setFiringOutcome(m, jc, p4, 1); t2.setFiringOutcome(m, jc, p2, 1);
m = t3.addMode('m3'); t3.setDistribution(m, Exp(mu(3)));
t3.setEnablingConditions(m, jc, p4, 1); t3.setEnablingConditions(m, jc, p0, 1);
t3.setFiringOutcome(m, jc, p3, 1); t3.setFiringOutcome(m, jc, p1, 1);
m = t4.addMode('m4'); t4.setDistribution(m, Exp(mu(4)));
t4.setEnablingConditions(m, jc, p3, 1); t4.setFiringOutcome(m, jc, p0, 1);
R = model.initRoutingMatrix();
R.set(jc,jc,p2,t1,1.0); R.set(jc,jc,t1,p5,1.0);
R.set(jc,jc,p5,t2,1.0); R.set(jc,jc,p1,t2,1.0);
R.set(jc,jc,t2,p4,1.0); R.set(jc,jc,t2,p2,1.0);
R.set(jc,jc,p4,t3,1.0); R.set(jc,jc,p0,t3,1.0);
R.set(jc,jc,t3,p3,1.0); R.set(jc,jc,t3,p1,1.0);
R.set(jc,jc,p3,t4,1.0); R.set(jc,jc,t4,p0,1.0);
model.link(R);
p5.setState(0); p4.setState(0); p3.setState(0);
p2.setState(3); p1.setState(2); p0.setState(4);
end

function model = i_inhibiting(n)
% Three places, four modes, one inhibitor arc; the token count is conserved.
model = Network('spn');
P1 = Place(model,'P1'); P2 = Place(model,'P2'); P3 = Place(model,'P3');
T1 = Transition(model,'T1'); T2 = Transition(model,'T2'); T3 = Transition(model,'T3');
jc = ClosedClass(model,'Class1', n, P1, 0);
m = T1.addMode('Mode1'); T1.setDistribution(m, Exp(2));
T1.setEnablingConditions(m, jc, P1, 2); T1.setFiringOutcome(m, jc, P2, 2);
m = T1.addMode('Mode2'); T1.setDistribution(m, Exp(1));
T1.setEnablingConditions(m, jc, P1, 1); T1.setFiringOutcome(m, jc, P3, 1);
m = T2.addMode('Mode3'); T2.setDistribution(m, Exp(4));
T2.setEnablingConditions(m, jc, P2, 1); T2.setFiringOutcome(m, jc, P1, 1);
m = T3.addMode('Mode4'); T3.setDistribution(m, Exp(1));
T3.setEnablingConditions(m, jc, P3, 3); T3.setInhibitingConditions(m, jc, P2, 1);
T3.setFiringOutcome(m, jc, P1, 3);
R = model.initRoutingMatrix();
R.set(jc,jc,P1,T1,1.0); R.set(jc,jc,P2,T2,1.0);
R.set(jc,jc,P2,T3,1.0); R.set(jc,jc,P3,T3,1.0);
R.set(jc,jc,T1,P2,1.0); R.set(jc,jc,T1,P3,1.0);
R.set(jc,jc,T2,P1,1.0); R.set(jc,jc,T3,P1,1.0);
model.link(R);
P1.setState(n); P2.setState(0); P3.setState(0);
end
