% Mean busy period of order n for a subnetwork of a closed queueing network.
%
% The busy period of order n for a set of stations is the time from the
% instant a job entering the set finds n-1 jobs in it up to the next instant
% when fewer than n remain. SolverNC evaluates it exactly from the
% normalizing constants of the subnetwork and of its complement (H. Daduna,
% "Busy Periods for Subnetworks in Stochastic Networks: Mean Value Analysis",
% J. ACM 35(3), 1988); SolverLDES measures the same quantity along a
% simulated sample path.

clear node jobclass;

N = 5;
model = Network('busyPeriodModel');

node{1} = Queue(model, 'Queue1', SchedStrategy.FCFS);
node{2} = Queue(model, 'Queue2', SchedStrategy.FCFS);
node{3} = Queue(model, 'Queue3', SchedStrategy.FCFS);

jobclass{1} = ClosedClass(model, 'Class1', N, node{1}, 0);

node{1}.setService(jobclass{1}, Exp(1.5));
node{2}.setService(jobclass{1}, Exp(0.9));
node{3}.setService(jobclass{1}, Exp(2.0));

P = model.initRoutingMatrix();
P{jobclass{1}, jobclass{1}} = [0 0.6 0.4; 0.7 0 0.3; 0.5 0.5 0];
model.link(P);

%% Exact mean value analysis of the busy period
ncSolver = SolverNC(model);
fprintf(1, '\nMean busy period of order n (exact, SolverNC)\n');
fprintf(1, '%-18s %10s %10s %10s\n', 'subnetwork', 'n=1', 'n=3', 'n=5');
subnets = {1, 2, 3, [1 2]};
for s = 1:numel(subnets)
    I = subnets{s};
    b = ncSolver.getAvgBusyPeriod(I, [1 3 5]);
    fprintf(1, '%-18s %10.4f %10.4f %10.4f\n', mat2str(I), b(1), b(2), b(3));
end

%% The same quantity measured by simulation
ldesSolver = SolverLDES(model, 'samples', 2e6, 'seed', 23000);
fprintf(1, '\nMean busy period of order n (measured, SolverLDES)\n');
fprintf(1, '%-18s %10s %10s %10s\n', 'subnetwork', 'n=1', 'n=3', 'n=5');
for s = 1:numel(subnets)
    I = subnets{s};
    b = ldesSolver.getAvgBusyPeriod(I, -1, [1 3 5]);
    fprintf(1, '%-18s %10.4f %10.4f %10.4f\n', mat2str(I), b(1), b(2), b(3));
end
