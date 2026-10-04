% Open queue with a MAP_t SERVICE process, simulated by LDES.
%
% A MAP_t service is the only process LDES simulates that is both non-renewal
% and non-stationary. Segment k covers [breakpoints(k), breakpoints(k+1)) and
% carries the pair (D0{k}, D1{k}); D1 fires a service COMPLETION and D0 only
% moves the modulating phase. The engine walks that phase process forward from
% the instant the job ENTERS SERVICE -- not from the instant it arrived, which
% may fall under a different segment -- and resumes from the phase the previous
% completion left behind, so successive service times are correlated exactly as
% for an ordinary MAP.
%
% The schedule below is an MMPP whose two environment states hold the same
% sojourn rates throughout (Q_env) while their arrival intensities change from
% segment to segment, so the server alternates between a fast and a slow mode
% AND the pair of modes itself changes with the clock.
%
% Because the walk is a wall-clock one, it is exact only where service runs
% continuously at unit rate once started. LDES therefore refuses a MAP_t (or
% PH_t) service under processor sharing, under a preemptive discipline, with
% load dependence and with heterogeneous server pools, rather than returning a
% number from a sample path that does not model any of those. Use INF or a
% non-preemptive FCFS/LCFS-family discipline, as here.

model = Network('model');

source = Source(model,'Source');
queue = Queue(model, 'Queue', SchedStrategy.FCFS);
sink = Sink(model,'Sink');

jobclass = OpenClass(model, 'OpenClass', 0);

% Three segments held for 1, 1, 2 time units, repeating with period 4.
breakpoints = [0, 1, 2, 4];
Q_env = [-1, 1; 2, -2];
gamma = -diag(Q_env);
mu_values = [2, 8, 4; ...
             4, 2, 8];
D0 = cell(1, numel(breakpoints) - 1);
D1 = cell(1, numel(breakpoints) - 1);
for segment = 1:numel(D0)
    D0{segment} = Q_env;
    D0{segment}(1, 1) = -gamma(1) - mu_values(1, segment);
    D0{segment}(2, 2) = -gamma(2) - mu_values(2, segment);
    D1{segment} = diag(mu_values(:, segment));
end

source.setArrival(jobclass, Exp(10));
queue.setService(jobclass, MAPt(breakpoints, D0, D1, true));
queue.setCapacity(5);
queue.setNumberOfServers(1);

model.link(Network.serialRouting(source,queue,sink));

% The offered load exceeds what the server can clear, so the finite buffer
% drops: Tput falls short of the arrival rate and QLen sits just under the
% capacity. Both are signatures that the service process is being simulated --
% a service that sampled as zero would report QLen = Util = 0 with Tput = 10.
AvgTable = LDES(model,'seed',1234,'samples',1e6).getAvgTable;
AvgTable

% A MAP_t whose segments are identical carries no time dependence, so it must
% reproduce the ordinary MAP with those matrices. This is the degeneracy that
% pins the schedule machinery: it exercises every boundary crossing and still
% has to land on the time-homogeneous answer.
D0bar = {Q_env - diag([5 10]), Q_env - diag([5 10])};
D1bar = {diag([5 10]), diag([5 10])};
flat = Network('flat');
src2 = Source(flat,'Source');
q2 = Queue(flat,'Queue',SchedStrategy.FCFS);
snk2 = Sink(flat,'Sink');
cls2 = OpenClass(flat,'OpenClass',0);
src2.setArrival(cls2, Exp(1));
q2.setService(cls2, MAPt([0 0.25 0.5], D0bar, D1bar, true));
flat.link(Network.serialRouting(src2,q2,snk2));

homog = Network('homog');
src3 = Source(homog,'Source');
q3 = Queue(homog,'Queue',SchedStrategy.FCFS);
snk3 = Sink(homog,'Sink');
cls3 = OpenClass(homog,'OpenClass',0);
src3.setArrival(cls3, Exp(1));
q3.setService(cls3, MAP(D0bar{1}, D1bar{1}));
homog.link(Network.serialRouting(src3,q3,snk3));

tFlat = LDES(flat,'seed',1234,'samples',4e5).getAvgTable;
tHomog = LDES(homog,'seed',1234,'samples',4e5).getAvgTable;
fprintf('\nconstant schedule vs ordinary MAP (must agree):\n');
fprintf('  MAPt  QLen=%8.5f Util=%8.5f RespT=%8.5f\n', ...
    tFlat.QLen(2), tFlat.Util(2), tFlat.RespT(2));
fprintf('  MAP   QLen=%8.5f Util=%8.5f RespT=%8.5f\n', ...
    tHomog.QLen(2), tHomog.Util(2), tHomog.RespT(2));
