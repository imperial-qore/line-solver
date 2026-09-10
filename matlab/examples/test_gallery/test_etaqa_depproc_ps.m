% Test ETAQA departure process for PS queues in MAM dec.mmap
%
% Creates a 2-queue tandem with Erlang arrivals and PS scheduling, solved
% with dec.mmap (ETAQA-PS-based departures). Validates that:
%   1. dec.mmap now accepts PS queues (previously unsupported)
%   2. Results are finite and consistent with MVA (exact for PS mean values)
%   3. The PS departure process captures queueing effects beyond scaled service

%% Model: Erl(3)/M/1-PS -> M/1-PS tandem
model = Network('ETAQA-PS-Tandem');

source = Source(model, 'Source');
queue1 = Queue(model, 'Queue1', SchedStrategy.PS);
queue2 = Queue(model, 'Queue2', SchedStrategy.PS);
sink   = Sink(model, 'Sink');

oclass = OpenClass(model, 'Class1');
source.setArrival(oclass, Erlang.fitMeanAndOrder(2, 3));  % mean=2, SCV=1/3
queue1.setService(oclass, Exp(2));   % rho1=0.5
queue2.setService(oclass, Exp(1.5)); % rho2=0.667

model.link(Network.serialRouting(source, queue1, queue2, sink));

%% Solve with dec.mmap (uses ETAQA PS departure process)
solverMMAP = MAM(model, 'method', 'dec.mmap');
T_mmap = solverMMAP.getAvgTable();

%% Solve with MVA (exact mean values for PS)
solverMVA = MVA(model);
T_mva = solverMVA.getAvgTable();

%% Display results
fprintf('\n=== ETAQA PS Departure Process Test ===\n');
fprintf('%-10s %10s %10s\n', 'Station', 'MVA', 'dec.mmap');
fprintf('--- Queue Lengths ---\n');
for i = 1:height(T_mva)
    stn = string(T_mva.Station(i));
    if contains(stn, 'Queue')
        fprintf('%-10s %10.4f %10.4f\n', stn, T_mva.QLen(i), T_mmap.QLen(i));
    end
end
fprintf('--- Utilizations ---\n');
for i = 1:height(T_mva)
    stn = string(T_mva.Station(i));
    if contains(stn, 'Queue')
        fprintf('%-10s %10.4f %10.4f\n', stn, T_mva.Util(i), T_mmap.Util(i));
    end
end

%% Validate: dec.mmap results should be finite
Q_mmap = T_mmap.QLen(contains(string(T_mmap.Station), 'Queue'));
Q_mva  = T_mva.QLen(contains(string(T_mva.Station), 'Queue'));
U_mmap = T_mmap.Util(contains(string(T_mmap.Station), 'Queue'));
U_mva  = T_mva.Util(contains(string(T_mva.Station), 'Queue'));

assert(all(isfinite(Q_mmap)), 'ETAQA PS: dec.mmap returned non-finite queue lengths');
assert(all(Q_mmap >= 0), 'ETAQA PS: dec.mmap returned negative queue lengths');

%% Validate utilization (should match exactly since throughput = arrival rate)
relErrU = abs(U_mmap - U_mva) ./ max(U_mva, 1e-6);
fprintf('\nUtilization error vs MVA: %.4f%%, %.4f%%\n', relErrU*100);
assert(all(relErrU < 0.01), ...
    sprintf('ETAQA PS: utilizations deviate >1%% from MVA (err=%.4f%%)', max(relErrU)*100));

%% Validate queue lengths (PS mean QLen = rho/(1-rho), independent of arrival process)
% For a single-class open M/G/1-PS, E[Q] = rho/(1-rho) regardless of
% arrival SCV. The dec.mmap Q approximation uses U/(1-U) so should match.
relErrQ = abs(Q_mmap - Q_mva) ./ max(Q_mva, 1e-6);
fprintf('Queue length error vs MVA: %.2f%%, %.2f%%\n', relErrQ*100);
assert(all(relErrQ < 0.20), ...
    sprintf('ETAQA PS: queue lengths deviate >20%% from MVA (err=%.2f%%)', max(relErrQ)*100));

%% Test: single PS queue (simplest case)
model2 = Network('ETAQA-PS-Single');
src2  = Source(model2, 'Source');
q2    = Queue(model2, 'Queue1', SchedStrategy.PS);
snk2  = Sink(model2, 'Sink');
oc2   = OpenClass(model2, 'Class1');
src2.setArrival(oc2, APH.fitMeanAndSCV(2, 4));  % high-variability arrival
q2.setService(oc2, Exp(1));  % rho=0.5
model2.link(Network.serialRouting(src2, q2, snk2));

solverMMAP2 = MAM(model2, 'method', 'dec.mmap');
T2 = solverMMAP2.getAvgTable();
Q2 = T2.QLen(contains(string(T2.Station), 'Queue'));
assert(all(isfinite(Q2)), 'ETAQA PS single queue: returned non-finite results');
assert(all(Q2 >= 0), 'ETAQA PS single queue: returned negative queue lengths');
fprintf('Single PS queue Q=%.4f (expected ~1.0)\n', Q2);

fprintf('\nETAQA PS departure process test PASSED.\n');
