% renv_container_terminal - Daily-cycle container terminal in a random environment
%
% Advanced example for ENV's STATE-VECTOR analyzer (options.method =
% 'statevec'), which carries the full per-stage state distribution across
% environment switches instead of collapsing it to marginal mean queue lengths.
%
% Case study (after the fyp26 SOQN blending port): the Rotterdam container
% terminal of Dhingra et al. Container handling demand varies over the 24 hours
% of a day; each hour is one environment "stage" with its own demand intensity
% and (random, exponential) duration. The environment visits the 24 hourly
% stages in a fixed daily cycle 1->2->...->24->1.
%
% Here the semi-open SOQN of the original study is rendered as the equivalent
% finite-token CLOSED network required by ENV: N straddle carriers cycle
% between a yard staging Delay and a multi-server quay-crane Queue. The hourly
% demand modulates the yard staging rate (the closed-network analog of the
% time-varying external arrival intensity): busier hours stage containers to the
% cranes faster, so the cranes congest during the daily peaks.
%
% The script compares, for the day-averaged metrics:
%   (a) statevec  - full state-vector blending (this feature),
%   (b) meanfield - the default mean-field (marginal mean-queue-length) coupling,
%   (c) exact     - the stationary solution of the full joint
%                   (hour x network-state) random-environment CTMC.
% The state-vector blend tracks the exact joint solution far more closely than
% the mean-field collapse.
%
% A second section repeats the comparison with the internal terminal handling
% (quay cranes -> stacking cranes) collapsed by Norton's theorem into a single
% closed, load-dependent Flow-Equivalent Server (FES). This shows that a closed
% FES (load-dependent station, sn.lldscaling) is fully supported by the
% state-vector analyzer's CTMC backend.
%
% A third section analyses the OPEN counterpart: the hourly demand becomes an
% external Markov-modulated Poisson arrival stream (MMPP(24)) into a multi-server
% queue with an unbounded buffer, solved exactly by MAM as a QBD. This is
% the open / infinite-buffer regime that the closed-CTMC state-vector path cannot
% represent, and corresponds to the joint-MMPP baseline of the fyp26 study.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

clear; warning off;

%% Daily demand schedule (Dhingra Rotterdam terminal, 24 hours)
% Hourly container demand (moves/hr) and the mean duration of each hour-stage.
dailyRateHr     = [  6,  30,  40,  62,  76,  79, 119, 164, 152, 130,  79,  70, ...
                    57,  57, 113, 130, 162, 202, 148, 118,  92,  62,  36,   8];
dailyDurationHr = [0.51,0.84,0.76,0.98,0.67,2.33,1.80,0.62,0.36,1.30,1.02,0.98, ...
                   0.86,1.56,1.30,1.57,0.34,0.22,0.47,1.33,1.56,0.93,0.71,1.11];
E = numel(dailyRateHr);           % 24 hourly environment stages

%% Terminal physical parameters (kept small so the CTMC is exact and fast)
N        = 6;     % straddle carriers / containers circulating (closed tokens)
nCranes  = 2;     % quay cranes (multi-server FCFS queue)
craneRate = 16;   % moves/hr served per crane
% Yard staging rate per hour: scale demand so peak hours saturate the cranes.
stageRate = dailyRateHr / N;      % per-token yard completion rate in each hour

%% Build the random environment: one stage per hour, cyclic daily transitions
env = Environment('RotterdamDailyCycle');
for h = 1:E
    env.addStage(sprintf('Hour%02d', h-1), 'operational', ...
        terminalModel(stageRate(h), craneRate, nCranes, N));
end
% Deterministic-on-average daily cycle h -> h+1 with exponential hour durations.
for h = 1:E
    hNext = mod(h, E) + 1;
    env.addTransition(sprintf('Hour%02d', h-1), sprintf('Hour%02d', hNext-1), ...
        Exp(1/dailyDurationHr(h)));
end
env.init();

fprintf('Rotterdam container terminal: %d hourly stages, %d carriers, %d cranes.\n', E, N, nCranes);

%% Inner solver: CTMC with a finite transient horizon (required by statevec)
T = 12;                            % hours; covers the hour-duration CDF tails (>99%)
ctmcFactory = @(m) CTMC(m, 'timespan', [0,T], 'verbose', false);

baseOpt = Solver.defaultOptions;
baseOpt.iter_max = 100;
baseOpt.iter_tol = 1e-5;
baseOpt.verbose  = false;

%% (a) State-vector analyzer
optStatevec = baseOpt; optStatevec.method = 'statevec';
solverStatevec = ENV(env, ctmcFactory, optStatevec);
[Qsv, Usv, ~, Tsv] = solverStatevec.getAvg();

%% (b) Mean-field analyzer (default coupling)
optMeanfield = baseOpt; optMeanfield.method = 'meanfield';
solverMeanfield = ENV(env, ctmcFactory, optMeanfield);
[Qmf, Umf, ~, Tmf] = solverMeanfield.getAvg();

%% (c) Exact joint random-environment CTMC (ground truth)
[Qex, Uex, Tex] = exactJointMetrics(env, ctmcFactory);

%% Report day-averaged quay-crane metrics
crane = 2;   % station index of the QuayCranes queue
fprintf('\n=== Day-averaged quay-crane metrics (container terminal) ===\n');
fprintf('%-10s %12s %12s %12s\n', 'analyzer', 'QLen', 'Util', 'Tput');
fprintf('%-10s %12.5f %12.5f %12.5f\n', 'exact',    Qex(crane), Uex(crane), Tex(crane));
fprintf('%-10s %12.5f %12.5f %12.5f\n', 'statevec', Qsv(crane), Usv(crane), Tsv(crane));
fprintf('%-10s %12.5f %12.5f %12.5f\n', 'meanfield', Qmf(crane), Umf(crane), Tmf(crane));

fprintf('\nQLen error vs exact:  statevec = %.3e , meanfield = %.3e\n', ...
    abs(Qsv(crane)-Qex(crane)), abs(Qmf(crane)-Qex(crane)));
fprintf('Util error vs exact:  statevec = %.3e , meanfield = %.3e\n', ...
    abs(Usv(crane)-Uex(crane)), abs(Umf(crane)-Uex(crane)));

fprintf('\nDay-averaged crane-queue table (state-vector analyzer):\n');
solverStatevec.getAvgTable()

%% ===== Closed flow-equivalent-server (FES) variant =========================
% The internal terminal handling (quay cranes -> stacking cranes) is collapsed
% by Norton's theorem into a single load-dependent FES whose rate mu(n) is the
% throughput of that subnetwork at n containers in process. The load-dependent
% station populates sn.lldscaling, which CTMC honours, so the state-vector
% analyzer applies unchanged to this closed FES network. Each hour-stage is now
% Yard(stageRate_h) -> FES; the physical terminal is fixed across the day, so the
% FES rate curve is computed once and reused by every stage.
fprintf('\n=== Closed FES variant (Norton-aggregated terminal internals) ===\n');
muQuay = 16; muStack = 20;                       % internal quay / stacking rates
fesRate = fesRateCurve(muQuay, muStack, N);      % load-dependent mu(n), n=1..N
fprintf('Norton FES rate curve mu(n) = %s\n', mat2str(fesRate, 5));

envFES = Environment('RotterdamDailyCycleFES');
for h = 1:E
    envFES.addStage(sprintf('Hour%02d', h-1), 'operational', ...
        terminalFESModel(stageRate(h), fesRate));
end
for h = 1:E
    hNext = mod(h, E) + 1;
    envFES.addTransition(sprintf('Hour%02d', h-1), sprintf('Hour%02d', hNext-1), ...
        Exp(1/dailyDurationHr(h)));
end
envFES.init();

[Qsf, Usf, ~, Tsf] = ENV(envFES, ctmcFactory, optStatevec).getAvg();
[Qff, Uff, ~, Tff] = ENV(envFES, ctmcFactory, optMeanfield).getAvg();
[Qxf, Uxf, Txf]    = exactJointMetrics(envFES, ctmcFactory);

fes = 2;   % station index of the load-dependent FES
fprintf('%-10s %12s %12s %12s\n', 'analyzer', 'QLen', 'Util', 'Tput');
fprintf('%-10s %12.5f %12.5f %12.5f\n', 'exact',    Qxf(fes), Uxf(fes), Txf(fes));
fprintf('%-10s %12.5f %12.5f %12.5f\n', 'statevec', Qsf(fes), Usf(fes), Tsf(fes));
fprintf('%-10s %12.5f %12.5f %12.5f\n', 'meanfield', Qff(fes), Uff(fes), Tff(fes));
fprintf('\nQLen error vs exact:  statevec = %.3e , meanfield = %.3e\n', ...
    abs(Qsf(fes)-Qxf(fes)), abs(Qff(fes)-Qxf(fes)));
fprintf('Util error vs exact:  statevec = %.3e , meanfield = %.3e\n', ...
    abs(Usf(fes)-Uxf(fes)), abs(Uff(fes)-Uxf(fes)));

%% ===== Open system via MAM (MMPP(24)/M/c, infinite buffer) =================
% The semi-open SOQN's open counterpart: rather than a finite carrier pool, the
% 24-hour demand is an external arrival stream whose rate is modulated by the
% hour of day. This is a Markov-modulated Poisson process (MMPP): the 24-state
% ring is the environment generator and the per-state arrival rate is the hourly
% demand. The resulting MMPP(24)/M/c queue has an UNBOUNDED buffer, so it is
% outside the closed-CTMC state-vector path; MAM solves it exactly as a QBD
% (this is the joint-MMPP baseline of the fyp26 study). Because the open stream
% is not throttled by a finite token pool, more crane capacity is needed for
% stability than in the closed models.
fprintf('\n=== Open system via MAM (MMPP(24)/M/c, infinite buffer) ===\n');
Qenv = zeros(E);
for h = 1:E
    hNext = mod(h, E) + 1;
    Qenv(h, hNext) = 1/dailyDurationHr(h);   % cyclic hour transition
    Qenv(h, h)     = -1/dailyDurationHr(h);
end
mmpp  = MAP({Qenv - diag(dailyRateHr), diag(dailyRateHr)});  % rate-modulated arrivals
cOpen = 8;                                                   % cranes (open regime)

openModel = Network('OpenTerminal');
ships    = Source(openModel, 'Ships');
quayOpen = Queue(openModel, 'QuayCranes', SchedStrategy.FCFS);
done     = Sink(openModel, 'Departures');
oc = OpenClass(openModel, 'Containers');
ships.setArrival(oc, mmpp);
quayOpen.setService(oc, Exp(craneRate));
quayOpen.setNumberOfServers(cOpen);
openModel.link(Network.serialRouting(ships, quayOpen, done));

meanLambda = sum(dailyDurationHr .* dailyRateHr) / sum(dailyDurationHr);
fprintf('MMPP mean arrival = %.2f/hr, %d cranes @ %.0f/hr, rho = %.3f\n', ...
    meanLambda, cOpen, craneRate, meanLambda/(cOpen*craneRate));
MAM(openModel).getAvgTable

%% ----- local functions -----------------------------------------------------
function qn = terminalModel(stageRate, craneRate, nCranes, N)
% One hour-stage network: N carriers cycling yard-Delay -> quay-crane Queue.
qn = Network('Terminal');
yard   = Delay(qn, 'Yard');
cranes = Queue(qn, 'QuayCranes', SchedStrategy.FCFS);
cranes.setNumberOfServers(nCranes);
containers = ClosedClass(qn, 'Containers', N, yard);
yard.setService(containers,   Exp(stageRate));
cranes.setService(containers, Exp(craneRate));
qn.link(Network.serialRouting(yard, cranes));
end

function qn = terminalFESModel(stageRate, fesRate)
% One hour-stage with the internal terminal handling as a closed load-dependent
% FES: N containers cycling yard-Delay -> FES, where the FES serves at rate
% mu(n) = fesRate(n) when n containers are inside it (set via setLoadDependence).
N = numel(fesRate);
qn = Network('TerminalFES');
yard = Delay(qn, 'Yard');
fesq = Queue(qn, 'TerminalFES', SchedStrategy.PS);
containers = ClosedClass(qn, 'Containers', N, yard);
yard.setService(containers, Exp(stageRate));
fesq.setService(containers, Exp(fesRate(1)));     % base rate mu(1)
fesq.setLoadDependence(fesRate / fesRate(1));      % scaling alpha(n) = mu(n)/mu(1)
qn.link(Network.serialRouting(yard, fesq));
end

function mu = fesRateCurve(muQuay, muStack, N)
% Norton flow-equivalent rates: throughput of the isolated internal subnetwork
% (quay cranes -> stacking cranes, processor sharing) at n = 1..N containers,
% with the rest of the terminal short-circuited by a near-instantaneous delay.
mu = zeros(1, N);
mvaOpt = MVA.defaultOptions; mvaOpt.method = 'exact'; mvaOpt.verbose = false;
for n = 1:N
    sub   = Network('TerminalInternals');
    ref   = Delay(sub, 'ShortCircuit');
    quay  = Queue(sub, 'QuayCranes', SchedStrategy.PS);
    stack = Queue(sub, 'StackingCranes', SchedStrategy.PS);
    cls   = ClosedClass(sub, 'Containers', n, ref);
    ref.setService(cls,   Exp(1e6));               % ~instantaneous short-circuit
    quay.setService(cls,  Exp(muQuay));
    stack.setService(cls, Exp(muStack));
    sub.link(Network.serialRouting(ref, quay, stack));
    T = MVA(sub, mvaOpt).getAvgTable;
    mu(n) = T.Tput(1);                             % subnetwork throughput at n
end
end

function [QN, UN, TN] = exactJointMetrics(env, ctmcFactory)
% Day-averaged metrics from the exact joint (hour x network-state) CTMC.
solverX = ENV(env, ctmcFactory, Solver.defaultOptions);
renvQ   = solverX.getGenerator();          % flattened exact generator
piJoint = ctmc_solve(renvQ); piJoint = piJoint(:)';

mdl  = env.getEnsemble;
E    = numel(mdl);
M    = mdl{1}.getNumberOfStations;
K    = mdl{1}.getNumberOfClasses;
QN = zeros(M,K); UN = zeros(M,K); TN = zeros(M,K); probEnv = zeros(1,E);
off  = 0;
for e = 1:E
    opts = ctmcFactory(mdl{e}).getOptions;
    [Qe, SSe, SSae, ~, arve, depe, sne] = solver_ctmc(mdl{e}.getStruct, opts);
    ns  = size(Qe,1);
    blk = piJoint(off+1:off+ns); off = off + ns;
    probEnv(e) = sum(blk);
    if sum(blk) > 0, blk = blk/sum(blk); end
    [QNe, UNe, ~, TNe] = solver_ctmc_avg_from_pi(sne, blk, SSe, SSae, arve, depe);
    QN = QN + probEnv(e)*QNe;
    UN = UN + probEnv(e)*UNe;
    TN = TN + probEnv(e)*TNe;
end
end
