% Joint probability of the per-station TOTAL queue lengths, all classes summed
% out, from SolverNC.getProbSysMarg.
%
% Compare with statepr_sys_aggr.m, which fixes the PER-CLASS population of
% every station: each probability here is the sum of that one over every
% per-class table with these row sums. The fibre grows combinatorially, so the
% quantity is evaluated as a permanent of the demand matrix replicated once per
% job rather than by enumerating it.
clear node jobclass solver

model = Network('model');

node{1} = Delay(model, 'Delay');
node{2} = Queue(model, 'Queue1', SchedStrategy.PS);
node{3} = Queue(model, 'Queue2', SchedStrategy.PS);

N = [2,1];
jobclass{1} = ClosedClass(model, 'Class1', N(1), node{1}, 0);
jobclass{2} = ClosedClass(model, 'Class2', N(2), node{1}, 0);

node{1}.setService(jobclass{1}, Exp(1/1.5));
node{1}.setService(jobclass{2}, Exp(1/2.0));
node{2}.setService(jobclass{1}, Exp(1/0.7));
node{2}.setService(jobclass{2}, Exp(1/0.4));
node{3}.setService(jobclass{1}, Exp(1/0.3));
node{3}.setService(jobclass{2}, Exp(1/0.9));

P = model.initRoutingMatrix;
P{jobclass{1}} = Network.serialRouting(node{1},node{2},node{3});
P{jobclass{2}} = Network.serialRouting(node{1},node{2},node{3});
model.link(P);

solver = SolverNC(model);

% Every way of splitting the closed population across the stations
states = multichoose(model.getNumberOfStations(), sum(N));

line_printf('\n  n(Delay) n(Queue1) n(Queue2)        P(n)\n');
p = zeros(size(states,1),1);
for j = 1:size(states,1)
    p(j) = solver.getProbSysMarg(states(j,:));
    line_printf('  %8d %9d %9d  %10.6f\n', states(j,1), states(j,2), states(j,3), p(j));
end
line_printf('  %s\n  sum %39.6f\n', repmat('-',1,42), sum(p));

% The law is exact, so its first moments are the queue lengths
line_printf('\n  E[n] from the joint law : %s\n', mat2str(p'*states, 8));
line_printf('  QLen from SolverCTMC    : %s\n', mat2str(sum(SolverCTMC(model).getAvgQLen(),2)', 8));

% The approximate permanent engines trade accuracy for cost on models whose
% class count makes the exact expansion dear. They need a demand matrix with
% full support and refuse a structural zero rather than flooring it.
pb = zeros(size(states,1),1);
for j = 1:size(states,1)
    pb(j) = solver.getProbSysMarg(states(j,:), 'bethe');
end
pb = pb / sum(pb);
line_printf('  Bethe engine, mean relative error : %.2f%%\n', 100*mean(abs(pb-p)./p));

solver.citations()
