%% Test BMAP/PH/N queue simulation
% Simple test of BMAP arrivals with PH service using LINE

clear; clc;

fprintf('=== BMAP/PH/N Queue Test ===\n\n');

%% Create network
model = Network('BMAP_PH_N_Test');

%% Define nodes
source = Source(model, 'Source');
queue = Queue(model, 'Queue', SchedStrategy.FCFS);
queue.setNumberOfServers(5);  % N=5 servers
sink = Sink(model, 'Sink');

%% Define job class
jobclass = OpenClass(model, 'Customer');

%% Define BMAP arrival process (simplified from paper)
% Using a 2-phase MAP approximation for testing
D0 = [-2, 0.5; 0.3, -1];
D1 = [1.5, 0; 0, 0.7];

try
    % Try MAP first (simpler)
    map_arrival = MAP(D0, D1);
    source.setArrival(jobclass, map_arrival);
    fprintf('Using MAP arrivals\n');
catch
    % Fall back to exponential
    source.setArrival(jobclass, Exp(2));
    fprintf('Using Exp arrivals\n');
end

%% Define PH service (from paper)
beta = [0.5, 0.5];
S = [-2, 0; 0, -2/3];

try
    ph_service = PH(beta, S);
    queue.setService(jobclass, ph_service);
    fprintf('Using PH service with mean = %.4f\n', ph_service.getMean());
catch
    queue.setService(jobclass, Exp(1));
    fprintf('Using Exp service\n');
end

%% Link nodes
model.link(Network.serialRouting(source, queue, sink));

fprintf('\nNetwork created successfully.\n');

%% Solve with different methods
fprintf('\n--- Solving with JMT (Simulation) ---\n');

try
    options = Solver.defaultOptions();
    options.samples = 50000;
    options.seed = 23000;

    solver = SolverJMT(model, options);
    [QN, UN, RN, TN] = solver.getAvg();

    fprintf('Queue length: %.4f\n', QN(2));
    fprintf('Utilization: %.4f\n', UN(2));
    fprintf('Response time: %.4f\n', RN(2));
    fprintf('Throughput: %.4f\n', TN(2));
catch ME
    fprintf('JMT Error: %s\n', ME.message);
end

fprintf('\n--- Solving with MAM (Analytical) ---\n');

try
    solver_mam = SolverMAM(model);
    [QN_mam, UN_mam, RN_mam, TN_mam] = solver_mam.getAvg();

    fprintf('Queue length: %.4f\n', QN_mam(2));
    fprintf('Utilization: %.4f\n', UN_mam(2));
    fprintf('Response time: %.4f\n', RN_mam(2));
    fprintf('Throughput: %.4f\n', TN_mam(2));
catch ME
    fprintf('MAM Error: %s\n', ME.message);
end

fprintf('\nTest complete.\n');
