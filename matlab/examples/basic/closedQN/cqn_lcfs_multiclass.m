clear node jobclass solver AvgTable;
% Example: Multiclass closed queueing network with LCFS and LCFS-PR scheduling
%
% This example demonstrates a 2-station closed queueing network where:
%   - Station 1: LCFS (Last-Come-First-Served, non-preemptive)
%   - Station 2: LCFS-PR (LCFS with Preemption-Resume)
%
% Reference:
%   G. Casale, "A family of multiclass LCFS queueing networks with
%   order-dependent product-form solutions", QUESTA 2026.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

%% Define service rates
% mu(i,r) = service rate at station i for class r
mu = 1./[1, 3, 5; 2, 4, 6];  % 2 stations x 3 classes
R = size(mu, 2);  % number of classes

%% Create the network model
model = Network('LCFS Multiclass Model');

% Create nodes
node{1} = Queue(model, 'Queue1', SchedStrategy.LCFS);
node{2} = Queue(model, 'Queue2', SchedStrategy.LCFSPR);

% Create job classes (one job per class)
for r = 1:R
    jobclass{r} = ClosedClass(model, sprintf('Class%d', r), 1, node{1}, 0);
end

% Set service times (exponential distributions)
for r = 1:R
    node{1}.setService(jobclass{r}, Exp(mu(1, r)));
    node{2}.setService(jobclass{r}, Exp(mu(2, r)));
end

% Set up routing: jobs alternate between the two queues
P = model.initRoutingMatrix;
for r = 1:R
    P{r, r} = [0, 1; 1, 0];  % Queue1 -> Queue2 -> Queue1
end
model.link(P);

%% Solve with MVA
fprintf('Solving LCFS+LCFS-PR network with MVA...\n');
solver_mva = MVA(model, 'method', 'exact');
AvgTable_MVA = solver_mva.getAvgTable();
disp('MVA Results:');
disp(AvgTable_MVA);

%% Solve with CTMC for validation
fprintf('\nSolving with CTMC for validation...\n');
solver_ctmc = CTMC(model);
AvgTable_CTMC = solver_ctmc.getAvgTable();
disp('CTMC Results:');
disp(AvgTable_CTMC);

%% Compare results
fprintf('\n=== Comparison ===\n');
fprintf('Metric\t\t\tMVA\t\tCTMC\t\tDiff\n');
fprintf('------\t\t\t---\t\t----\t\t----\n');

% Compare queue lengths
for r = 1:R
    for k = 1:2
        idx = (r-1)*2 + k;
        q_mva = AvgTable_MVA.QLen(idx);
        q_ctmc = AvgTable_CTMC.QLen(idx);
        diff = abs(q_mva - q_ctmc);
        fprintf('Q(%d,%d)\t\t\t%.6f\t%.6f\t%.2e\n', k, r, q_mva, q_ctmc, diff);
    end
end

fprintf('\n');

% Compare utilizations
for r = 1:R
    for k = 1:2
        idx = (r-1)*2 + k;
        u_mva = AvgTable_MVA.Util(idx);
        u_ctmc = AvgTable_CTMC.Util(idx);
        diff = abs(u_mva - u_ctmc);
        fprintf('U(%d,%d)\t\t\t%.6f\t%.6f\t%.2e\n', k, r, u_mva, u_ctmc, diff);
    end
end

fprintf('\n');

% Compare throughputs
for r = 1:R
    for k = 1:2
        idx = (r-1)*2 + k;
        t_mva = AvgTable_MVA.Tput(idx);
        t_ctmc = AvgTable_CTMC.Tput(idx);
        diff = abs(t_mva - t_ctmc);
        fprintf('T(%d,%d)\t\t\t%.6f\t%.6f\t%.2e\n', k, r, t_mva, t_ctmc, diff);
    end
end

%% Verify results match within tolerance
tol = 1e-4;
qlen_match = all(abs(AvgTable_MVA.QLen - AvgTable_CTMC.QLen) < tol);
util_match = all(abs(AvgTable_MVA.Util - AvgTable_CTMC.Util) < tol);
tput_match = all(abs(AvgTable_MVA.Tput - AvgTable_CTMC.Tput) < tol);

fprintf('\n=== Validation ===\n');
if qlen_match && util_match && tput_match
    fprintf('SUCCESS: MVA and CTMC results match within tolerance %.0e\n', tol);
else
    fprintf('WARNING: Results differ beyond tolerance %.0e\n', tol);
    if ~qlen_match
        fprintf('  - Queue lengths do not match\n');
    end
    if ~util_match
        fprintf('  - Utilizations do not match\n');
    end
    if ~tput_match
        fprintf('  - Throughputs do not match\n');
    end
end
