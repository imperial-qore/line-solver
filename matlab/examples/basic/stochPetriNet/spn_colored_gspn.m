% Colored generalized stochastic Petri net (CGSPN).
%
% A generalized stochastic Petri net mixes IMMEDIATE transitions, which fire
% as soon as they are enabled and consume no time, with TIMED ones. It becomes
% COLORED when the tokens carry a type: in LINE a token color is a job class,
% so a place holds a marking per color and a transition declares one MODE per
% color it serves.
%
% Here two colors, Gold and Silver, circulate between a Buffer place and a
% Server place. The immediate transition admit moves one token into the empty
% server; its two modes are held apart by inhibiting arcs on both colors, so
% the server is a mutual-exclusion resource and the firing weights arbitrate
% between the colors whenever both are waiting. The timed transition serve
% returns the token to the buffer at a color-dependent rate.
%
% Admission is immediate, so the server is never idle and the model is exact
% by hand: the weights split the completions 2:1 in favour of Gold, hence
% X_Gold = 2*X_Silver, U_c = X_c/mu_c and U_Gold + U_Silver = 1.

model = Network('ColoredGSPN');

% Two places and two transitions. Declare every node before parameterizing
% any of it: the enabling, inhibiting and firing matrices are indexed by
% (node, class) over the nodes that exist when they are first set.
buf   = Place(model, 'Buffer');
srv   = Place(model, 'Server');
admit = Transition(model, 'admit');
serve = Transition(model, 'serve');

% Two token colors, both starting in the buffer.
gold   = ClosedClass(model, 'Gold',   2, buf, 0);
silver = ClosedClass(model, 'Silver', 2, buf, 0);

% Immediate transition: one mode per color, each admitting a token only when
% the server holds no token of either color. The weights decide the color mix.
mg = admit.addMode('gold');
admit.setTimingStrategy(mg, TimingStrategy.IMMEDIATE);
admit.setEnablingConditions(mg, gold, buf, 1);
admit.setInhibitingConditions(mg, gold, srv, 1);
admit.setInhibitingConditions(mg, silver, srv, 1);
admit.setFiringOutcome(mg, gold, srv, 1);
admit.setFiringWeights(mg, 2.0);

ms = admit.addMode('silver');
admit.setTimingStrategy(ms, TimingStrategy.IMMEDIATE);
admit.setEnablingConditions(ms, silver, buf, 1);
admit.setInhibitingConditions(ms, gold, srv, 1);
admit.setInhibitingConditions(ms, silver, srv, 1);
admit.setFiringOutcome(ms, silver, srv, 1);
admit.setFiringWeights(ms, 1.0);

% Timed transition: one mode per color, with color-dependent service rates.
sg = serve.addMode('gold');
serve.setDistribution(sg, Exp(3.0));
serve.setEnablingConditions(sg, gold, srv, 1);
serve.setFiringOutcome(sg, gold, buf, 1);

ss = serve.addMode('silver');
serve.setDistribution(ss, Exp(1.5));
serve.setEnablingConditions(ss, silver, srv, 1);
serve.setFiringOutcome(ss, silver, buf, 1);

% Topology: both colors follow Buffer -> admit -> Server -> serve -> Buffer.
P = model.initRoutingMatrix();
for c = {gold, silver}
    P.set(c{1}, c{1}, buf,   admit, 1.0);
    P.set(c{1}, c{1}, admit, srv,   1.0);
    P.set(c{1}, c{1}, srv,   serve, 1.0);
    P.set(c{1}, c{1}, serve, buf,   1.0);
end
model.link(P);

% Initial marking: two tokens of each color in the buffer, server empty.
buf.setState([2, 2]);
srv.setState([0, 0]);

% Exact solution of the underlying CTMC, and a simulation cross-check.
AvgTableCTMC = SolverCTMC(model, 'cutoff', 4).avgTable()
AvgTableSSA = SolverSSA(model, 'seed', 23000, 'samples', 2e5).avgTable()
