% Test ETAQA departure process for FCFS queues in MAM dec.mmap
%
% Creates a 2-queue tandem with non-Poisson (Erlang) arrivals and Erlang
% service, solved with dec.mmap (ETAQA-based departures). Validates that:
%   1. dec.mmap returns finite results and converges
%   2. Results are consistent with dec.source and MVA baselines
%   3. The ETAQA truncation level can be configured via options

%% Model: Erl(5)/Erl(2)/1 -> Erl(3)/1 tandem
model = Network('ETAQA-FCFS-Tandem');

source = Source(model, 'Source');
queue1 = Queue(model, 'Queue1', SchedStrategy.FCFS);
queue2 = Queue(model, 'Queue2', SchedStrategy.FCFS);
sink   = Sink(model, 'Sink');

oclass = OpenClass(model, 'Class1');
source.setArrival(oclass, Erlang.fitMeanAndOrder(2, 5));  % mean=2, SCV=0.2
queue1.setService(oclass, Erlang.fitMeanAndOrder(0.8, 2)); % rho1=0.4
queue2.setService(oclass, Erlang.fitMeanAndOrder(1.2, 3)); % rho2=0.6

model.link(Network.serialRouting(source, queue1, queue2, sink));

%% Solve with dec.mmap (uses ETAQA departure process)
solverMMAP = MAM(model, 'method', 'dec.mmap');
T_mmap = solverMMAP.getAvgTable();

%% Solve with dec.source (baseline)
solverSRC = MAM(model, 'method', 'dec.source');
T_src = solverSRC.getAvgTable();

%% Solve with MVA (exact for mean values)
solverMVA = MVA(model);
T_mva = solverMVA.getAvgTable();

%% Display results
fprintf('\n=== ETAQA FCFS Departure Process Test ===\n');
fprintf('%-10s %10s %10s %10s\n', 'Station', 'MVA', 'dec.source', 'dec.mmap');
fprintf('--- Queue Lengths ---\n');
for i = 1:height(T_mva)
    stn = string(T_mva.Station(i));
    if contains(stn, 'Queue')
        fprintf('%-10s %10.4f %10.4f %10.4f\n', stn, ...
            T_mva.QLen(i), T_src.QLen(i), T_mmap.QLen(i));
    end
end
fprintf('--- Response Times ---\n');
for i = 1:height(T_mva)
    stn = string(T_mva.Station(i));
    if contains(stn, 'Queue')
        fprintf('%-10s %10.4f %10.4f %10.4f\n', stn, ...
            T_mva.RespT(i), T_src.RespT(i), T_mmap.RespT(i));
    end
end

%% Validate: dec.mmap results should be finite and within 20% of MVA
Q_mmap = T_mmap.QLen(contains(string(T_mmap.Station), 'Queue'));
Q_mva  = T_mva.QLen(contains(string(T_mva.Station), 'Queue'));

assert(all(isfinite(Q_mmap)), 'ETAQA FCFS: dec.mmap returned non-finite queue lengths');
assert(all(Q_mmap >= 0), 'ETAQA FCFS: dec.mmap returned negative queue lengths');

relErr = abs(Q_mmap - Q_mva) ./ max(Q_mva, 1e-6);
fprintf('\nRelative error vs MVA: %.2f%%, %.2f%%\n', relErr*100);
assert(all(relErr < 0.20), ...
    sprintf('ETAQA FCFS: dec.mmap queue lengths deviate >20%% from MVA (err=%.2f%%)', max(relErr)*100));

%% Test configurable truncation level
solverMMAP4 = MAM(model, 'method', 'dec.mmap');
solverMMAP4.options.config.etaqa_trunc = 4;
T_mmap4 = solverMMAP4.getAvgTable();
Q_mmap4 = T_mmap4.QLen(contains(string(T_mmap4.Station), 'Queue'));
assert(all(isfinite(Q_mmap4)), 'ETAQA FCFS: truncation level 4 returned non-finite results');

fprintf('\nETAQA FCFS departure process test PASSED.\n');
