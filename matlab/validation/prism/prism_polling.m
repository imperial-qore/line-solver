function model = prism_polling()
% PRISM_POLLING  Two-station cyclic server polling system of Ibe and Trivedi.
%
% MODEL = PRISM_POLLING() builds LINE's encoding of the case study distributed
% with PRISM as prism-examples/ctmcs/polling/poll2.sm, whose source is
%
%   O. Ibe and K. Trivedi, "Stochastic Petri Net Models of Polling Systems",
%   IEEE Journal on Selected Areas in Communications, 8(9):1649-1657, 1990.
%
% The server walks between two stations, polling at rate gamma and serving at
% rate mu; each station holds at most one job and receives Poisson arrivals at
% rate lambda = mu/N. PRISM carries the server position s and its activity a as
% integer variables; here each is one-hot over a pair of places, which is what
% a Petri net encoding of a finite control variable looks like.
%
% Arrivals are modelled as a closed net rather than with a Source. PRISM enables
% each station's arrival command exactly while that station is empty, so a
% two-token pool feeding the two stations is an exact equivalent: a station
% holds at most one job, so the pool can only run dry when both are full, and
% then no arrival was possible anyway. A Source would have raised the question
% of what LINE does with the routing probability of a blocked destination,
% which is not a question this model should be asking.
%
% Arc convention: an enabling condition is a PRE arc and consumes its
% multiplicity, so a place that a mode only tests is given back with a firing
% outcome of +1.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

N = 2;
mu = 1;
gamma = 200;
lambda = mu / N;

model = Network('polling');

pool = Place(model, 'pool');   % jobs outside the two stations
S1 = Place(model, 'S1');   % server at station 1
S2 = Place(model, 'S2');   % server at station 2
A0 = Place(model, 'A0');   % polling
A1 = Place(model, 'A1');   % serving
Q1 = Place(model, 'Q1');   % station 1 full
Q2 = Place(model, 'Q2');   % station 2 full

loop1a = Transition(model, 'loop1a');
loop1b = Transition(model, 'loop1b');
serve1 = Transition(model, 'serve1');
loop2a = Transition(model, 'loop2a');
loop2b = Transition(model, 'loop2b');
serve2 = Transition(model, 'serve2');

arr1 = Transition(model, 'arr1');
arr2 = Transition(model, 'arr2');

jobs = ClosedClass(model, 'jobs', 4, pool);

for p = {S1, S2, A0, A1, Q1, Q2}
    p{1}.setClassCapacity(jobs, 1);
end
pool.setClassCapacity(jobs, 2);
pool.setMarking(2);
S1.setMarking(1);   % the server starts polling station 1
A0.setMarking(1);

% [] (s1=0) -> lambda : (s1'=1), and its station-2 twin
m = arr1.addMode('m');
arr1.setDistribution(m, Exp(lambda));
arr1.setEnablingConditions(m, jobs, pool, 1);
arr1.setInhibitingConditions(m, jobs, Q1, 1);
arr1.setFiringOutcome(m, jobs, Q1, 1);

m = arr2.addMode('m');
arr2.setDistribution(m, Exp(lambda));
arr2.setEnablingConditions(m, jobs, pool, 1);
arr2.setInhibitingConditions(m, jobs, Q2, 1);
arr2.setFiringOutcome(m, jobs, Q2, 1);

% [loop1a] (s=1)&(a=0)&(s1=0) -> gamma : (s'=2)
m = loop1a.addMode('m');
loop1a.setDistribution(m, Exp(gamma));
loop1a.setEnablingConditions(m, jobs, S1, 1);
loop1a.setEnablingConditions(m, jobs, A0, 1);
loop1a.setInhibitingConditions(m, jobs, Q1, 1);
loop1a.setFiringOutcome(m, jobs, S2, 1);
loop1a.setFiringOutcome(m, jobs, A0, 1);

% [loop1b] (s=1)&(a=0)&(s1=1) -> gamma : (a'=1)
m = loop1b.addMode('m');
loop1b.setDistribution(m, Exp(gamma));
loop1b.setEnablingConditions(m, jobs, S1, 1);
loop1b.setEnablingConditions(m, jobs, A0, 1);
loop1b.setEnablingConditions(m, jobs, Q1, 1);
loop1b.setFiringOutcome(m, jobs, S1, 1);
loop1b.setFiringOutcome(m, jobs, A1, 1);
loop1b.setFiringOutcome(m, jobs, Q1, 1);

% [serve1] (s=1)&(a=1)&(s1=1) -> mu : (s'=2)&(a'=0)&(s1'=0)
m = serve1.addMode('m');
serve1.setDistribution(m, Exp(mu));
serve1.setEnablingConditions(m, jobs, S1, 1);
serve1.setEnablingConditions(m, jobs, A1, 1);
serve1.setEnablingConditions(m, jobs, Q1, 1);
serve1.setFiringOutcome(m, jobs, S2, 1);
serve1.setFiringOutcome(m, jobs, A0, 1);
serve1.setFiringOutcome(m, jobs, pool, 1);   % the served job returns to the pool

% [loop2a] (s=2)&(a=0)&(s2=0) -> gamma : (s'=1)
m = loop2a.addMode('m');
loop2a.setDistribution(m, Exp(gamma));
loop2a.setEnablingConditions(m, jobs, S2, 1);
loop2a.setEnablingConditions(m, jobs, A0, 1);
loop2a.setInhibitingConditions(m, jobs, Q2, 1);
loop2a.setFiringOutcome(m, jobs, S1, 1);
loop2a.setFiringOutcome(m, jobs, A0, 1);

% [loop2b] (s=2)&(a=0)&(s2=1) -> gamma : (a'=1)
m = loop2b.addMode('m');
loop2b.setDistribution(m, Exp(gamma));
loop2b.setEnablingConditions(m, jobs, S2, 1);
loop2b.setEnablingConditions(m, jobs, A0, 1);
loop2b.setEnablingConditions(m, jobs, Q2, 1);
loop2b.setFiringOutcome(m, jobs, S2, 1);
loop2b.setFiringOutcome(m, jobs, A1, 1);
loop2b.setFiringOutcome(m, jobs, Q2, 1);

% [serve2] (s=2)&(a=1)&(s2=1) -> mu : (s'=1)&(a'=0)&(s2'=0)
m = serve2.addMode('m');
serve2.setDistribution(m, Exp(mu));
serve2.setEnablingConditions(m, jobs, S2, 1);
serve2.setEnablingConditions(m, jobs, A1, 1);
serve2.setEnablingConditions(m, jobs, Q2, 1);
serve2.setFiringOutcome(m, jobs, S1, 1);
serve2.setFiringOutcome(m, jobs, A0, 1);
serve2.setFiringOutcome(m, jobs, pool, 1);

R = model.initRoutingMatrix();
arcs = {pool, arr1; arr1, Q1; pool, arr2; arr2, Q2; ...
        serve1, pool; serve2, pool; ...
S1, loop1a; A0, loop1a; loop1a, S2; loop1a, A0; ...
        S1, loop1b; A0, loop1b; Q1, loop1b; loop1b, S1; loop1b, A1; loop1b, Q1; ...
        S1, serve1; A1, serve1; Q1, serve1; serve1, S2; serve1, A0; ...
        S2, loop2a; A0, loop2a; loop2a, S1; loop2a, A0; ...
        S2, loop2b; A0, loop2b; Q2, loop2b; loop2b, S2; loop2b, A1; loop2b, Q2; ...
        S2, serve2; A1, serve2; Q2, serve2; serve2, S1; serve2, A0};
for i = 1:size(arcs, 1)
    R.set(jobs, jobs, arcs{i,1}, arcs{i,2}, 1.0);
end
model.link(R);
end
