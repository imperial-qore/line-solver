%% QRF Method Comparison on a Simple Closed Queueing Network
% Single-class closed network: Delay + Queue(FCFS) with 2-phase MAP service.
% Compares exact CTMC solution against all available QRF methods.
%
% QRF approximation methods (qrf.mmi, qrf.mem, qrf.bethe, qrf.mmi.ld,
% qrf.mmi.linear)
% approximate performance metrics via nonlinear programming.
% QRF bounds methods (qrf.bas, qrf.rsrd) compute utilization bounds via LP.

clear;

%% Build model
model = Network('QRF Comparison');

node{1} = Delay(model, 'Delay');
node{2} = Queue(model, 'Queue1', SchedStrategy.FCFS);

N = 5;
jobclass{1} = ClosedClass(model, 'Class1', N, node{1}, 0);

% Delay: exponential think time (mean = 1)
node{1}.setService(jobclass{1}, Exp.fitMean(1));

% Queue: 2-phase MAP service (mean ~ 0.5, moderate variability)
D0 = [-4, 2; 1, -3];
D1 = [1, 1; 0.5, 1.5];
node{2}.setService(jobclass{1}, MAP(D0, D1));

% Serial routing: Delay -> Queue -> Delay
P = model.initRoutingMatrix();
P{1} = Network.serialRouting(node);
model.link(P);

%% Exact CTMC reference
fprintf('=== Exact CTMC ===\n');
solver_exact = SolverCTMC(model);
AvgTable_exact = solver_exact.getAvgTable();
disp(AvgTable_exact);

%% QRF methods (no-blocking variants)
qrf_methods = {'qrf.mmi', 'qrf.mem', 'qrf.bethe', 'qrf.mmi.ld', 'qrf.mmi.linear'};

results = struct();
results.method = {};
results.QN_delay = [];
results.QN_queue = [];
results.UN_delay = [];
results.UN_queue = [];
results.TN_queue = [];
results.RN_queue = [];

% Store exact results
results.method{end+1} = 'exact (CTMC)';
exact_Q = AvgTable_exact.QLen;
exact_U = AvgTable_exact.Util;
exact_T = AvgTable_exact.Tput;
exact_R = AvgTable_exact.RespT;
results.QN_delay(end+1) = exact_Q(1);
results.QN_queue(end+1) = exact_Q(2);
results.UN_delay(end+1) = exact_U(1);
results.UN_queue(end+1) = exact_U(2);
results.TN_queue(end+1) = exact_T(2);
results.RN_queue(end+1) = exact_R(2);

for m = 1:length(qrf_methods)
    method = qrf_methods{m};
    fprintf('\n=== %s ===\n', method);
    try
        solver = SolverCTMC(model, 'method', method);
        AvgTable = solver.getAvgTable();
        disp(AvgTable);

        results.method{end+1} = method;
        results.QN_delay(end+1) = AvgTable.QLen(1);
        results.QN_queue(end+1) = AvgTable.QLen(2);
        results.UN_delay(end+1) = AvgTable.Util(1);
        results.UN_queue(end+1) = AvgTable.Util(2);
        results.TN_queue(end+1) = AvgTable.Tput(2);
        results.RN_queue(end+1) = AvgTable.RespT(2);
    catch err
        fprintf('  FAILED: %s\n', err.message);
        results.method{end+1} = [method, ' (FAILED)'];
        results.QN_delay(end+1) = NaN;
        results.QN_queue(end+1) = NaN;
        results.UN_delay(end+1) = NaN;
        results.UN_queue(end+1) = NaN;
        results.TN_queue(end+1) = NaN;
        results.RN_queue(end+1) = NaN;
    end
end

%% QRF bounds methods (qrf.bas, qrf.rsrd)
% These use no-blocking defaults for an unbounded network
bounds_methods = {'qrf.bas', 'qrf.rsrd'};
for m = 1:length(bounds_methods)
    method = bounds_methods{m};
    fprintf('\n=== %s ===\n', method);
    try
        solver = SolverCTMC(model, 'method', method);
        AvgTable = solver.getAvgTable();
        disp(AvgTable);

        results.method{end+1} = method;
        results.QN_delay(end+1) = AvgTable.QLen(1);
        results.QN_queue(end+1) = AvgTable.QLen(2);
        results.UN_delay(end+1) = AvgTable.Util(1);
        results.UN_queue(end+1) = AvgTable.Util(2);
        results.TN_queue(end+1) = AvgTable.Tput(2);
        results.RN_queue(end+1) = AvgTable.RespT(2);
    catch err
        fprintf('  FAILED: %s\n', err.message);
        results.method{end+1} = [method, ' (FAILED)'];
        results.QN_delay(end+1) = NaN;
        results.QN_queue(end+1) = NaN;
        results.UN_delay(end+1) = NaN;
        results.UN_queue(end+1) = NaN;
        results.TN_queue(end+1) = NaN;
        results.RN_queue(end+1) = NaN;
    end
end

%% Summary comparison table
fprintf('\n\n========================================\n');
fprintf('     QRF Method Comparison Summary\n');
fprintf('========================================\n');
fprintf('%-20s %8s %8s %8s %8s %8s\n', ...
    'Method', 'QN(D)', 'QN(Q)', 'UN(Q)', 'TN(Q)', 'RN(Q)');
fprintf('%-20s %8s %8s %8s %8s %8s\n', ...
    '--------------------', '--------', '--------', '--------', '--------', '--------');
for i = 1:length(results.method)
    fprintf('%-20s %8.4f %8.4f %8.4f %8.4f %8.4f\n', ...
        results.method{i}, ...
        results.QN_delay(i), results.QN_queue(i), ...
        results.UN_queue(i), results.TN_queue(i), results.RN_queue(i));
end

% Relative errors vs exact
fprintf('\n%-20s %8s %8s %8s %8s %8s\n', ...
    'Relative Error (%)', 'QN(D)', 'QN(Q)', 'UN(Q)', 'TN(Q)', 'RN(Q)');
fprintf('%-20s %8s %8s %8s %8s %8s\n', ...
    '--------------------', '--------', '--------', '--------', '--------', '--------');
for i = 2:length(results.method)
    errQD = abs(results.QN_delay(i) - results.QN_delay(1)) / max(abs(results.QN_delay(1)), 1e-10) * 100;
    errQQ = abs(results.QN_queue(i) - results.QN_queue(1)) / max(abs(results.QN_queue(1)), 1e-10) * 100;
    errUQ = abs(results.UN_queue(i) - results.UN_queue(1)) / max(abs(results.UN_queue(1)), 1e-10) * 100;
    errTQ = abs(results.TN_queue(i) - results.TN_queue(1)) / max(abs(results.TN_queue(1)), 1e-10) * 100;
    errRQ = abs(results.RN_queue(i) - results.RN_queue(1)) / max(abs(results.RN_queue(1)), 1e-10) * 100;
    fprintf('%-20s %7.2f%% %7.2f%% %7.2f%% %7.2f%% %7.2f%%\n', ...
        results.method{i}, errQD, errQQ, errUQ, errTQ, errRQ);
end
