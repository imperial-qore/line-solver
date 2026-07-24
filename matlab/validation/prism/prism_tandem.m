function model = prism_tandem(c)
% PRISM_TANDEM  Tandem queueing network of Hermanns, Meyer-Kayser and Siegle.
%
% MODEL = PRISM_TANDEM(C) builds LINE's encoding of the case study distributed
% with PRISM as prism-examples/ctmcs/tandem/tandem.sm, whose source is
%
%   H. Hermanns, J. Meyer-Kayser and M. Siegle, "Multi-Terminal Binary Decision
%   Diagrams to Represent and Analyse Continuous Time Markov Chains", in Proc.
%   3rd International Workshop on the Numerical Solution of Markov Chains,
%   pp. 188-207, 1999.
%
% An M/Cox2/1/c station feeds an M/M/1/c station. A service completion at the
% first station is blocked while the second is full: PRISM disables the
% synchronised route action, which freezes the Coxian phase rather than letting
% the job finish and wait. That is not LINE's blocking-after-service, so the
% model is built as a Petri net, where an inhibitor arc reproduces the freeze
% exactly, and the Coxian phase is carried by two one-hot places rather than by
% a phase-type firing time (a phase-type mode would restart on unblocking).
%
% Arc convention: an enabling condition is a PRE arc and consumes its
% multiplicity, so a place that a mode only tests is given back with a firing
% outcome of +1.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 1
    c = 2;
end

lambda = 4*c;
mu1a = 0.1*2;
mu1b = 0.9*2;
mu2 = 2;
kappa = 4;

model = Network('tandem');

source = Source(model, 'source');
sink = Sink(model, 'sink');
sc = Place(model, 'sc');
sm = Place(model, 'sm');
ph1 = Place(model, 'ph1');
ph2 = Place(model, 'ph2');
tPhase = Transition(model, 'phase');    % ph 1 -> 2, no job movement
tRoute1 = Transition(model, 'route1');  % completion from phase 1
tRoute2 = Transition(model, 'route2');  % completion from phase 2
tServe = Transition(model, 'serve');    % departure from the second station

jobs = OpenClass(model, 'jobs');
source.setArrival(jobs, Exp(lambda));

sc.setClassCapacity(jobs, c);
sm.setClassCapacity(jobs, c);
ph1.setClassCapacity(jobs, 1);
ph2.setClassCapacity(jobs, 1);
ph1.setMarking(1);   % the Coxian starts in phase 1

% [] (sc>0) & (ph=1) -> mu1a : (ph'=2);
m = tPhase.addMode('a');
tPhase.setDistribution(m, Exp(mu1a));
tPhase.setEnablingConditions(m, jobs, sc, 1);
tPhase.setEnablingConditions(m, jobs, ph1, 1);
tPhase.setFiringOutcome(m, jobs, sc, 1);    % tested, not consumed
tPhase.setFiringOutcome(m, jobs, ph2, 1);

% [route] (sc>0) & (ph=1) -> mu1b : (sc'=sc-1);
m = tRoute1.addMode('b');
tRoute1.setDistribution(m, Exp(mu1b));
tRoute1.setEnablingConditions(m, jobs, sc, 1);
tRoute1.setEnablingConditions(m, jobs, ph1, 1);
tRoute1.setInhibitingConditions(m, jobs, sm, c);   % blocked while sm is full
tRoute1.setFiringOutcome(m, jobs, ph1, 1);         % phase unchanged
tRoute1.setFiringOutcome(m, jobs, sm, 1);

% [route] (sc>0) & (ph=2) -> mu2 : (ph'=1) & (sc'=sc-1);
m = tRoute2.addMode('c');
tRoute2.setDistribution(m, Exp(mu2));
tRoute2.setEnablingConditions(m, jobs, sc, 1);
tRoute2.setEnablingConditions(m, jobs, ph2, 1);
tRoute2.setInhibitingConditions(m, jobs, sm, c);
tRoute2.setFiringOutcome(m, jobs, ph1, 1);
tRoute2.setFiringOutcome(m, jobs, sm, 1);

% [] (sm>0) -> kappa : (sm'=sm-1);
m = tServe.addMode('d');
tServe.setDistribution(m, Exp(kappa));
tServe.setEnablingConditions(m, jobs, sm, 1);

R = model.initRoutingMatrix();
R.set(jobs, jobs, source, sc, 1.0);
R.set(jobs, jobs, sc, tPhase, 1.0);
R.set(jobs, jobs, tPhase, sc, 1.0);
R.set(jobs, jobs, ph1, tPhase, 1.0);
R.set(jobs, jobs, tPhase, ph2, 1.0);
R.set(jobs, jobs, sc, tRoute1, 1.0);
R.set(jobs, jobs, ph1, tRoute1, 1.0);
R.set(jobs, jobs, tRoute1, ph1, 1.0);
R.set(jobs, jobs, tRoute1, sm, 1.0);
R.set(jobs, jobs, sc, tRoute2, 1.0);
R.set(jobs, jobs, ph2, tRoute2, 1.0);
R.set(jobs, jobs, tRoute2, ph1, 1.0);
R.set(jobs, jobs, tRoute2, sm, 1.0);
R.set(jobs, jobs, sm, tServe, 1.0);
R.set(jobs, jobs, tServe, sink, 1.0);
model.link(R);
end
