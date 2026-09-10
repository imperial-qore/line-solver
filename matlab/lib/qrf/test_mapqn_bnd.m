%% test_mapqn_bnd - Test suite for mapqn_bnd_* functions
% Tests the 6 newly created MATLAB ports against known reference values.

clear;
clc;
fprintf('=== Testing mapqn_bnd_* functions ===\n\n');

nPass = 0;
nFail = 0;

%% Test 1: mapqn_bnd_lr_pf - Symmetric 2-queue tandem
fprintf('--- Test 1: mapqn_bnd_lr_pf (symmetric 2-queue tandem) ---\n');
try
    params = struct();
    params.M = 2;
    params.N = 3;
    params.mu = [1.0; 1.0];
    params.r = [0, 1; 1, 0];
    params.verbose = false;

    % MVA reference: U = 0.75 for symmetric tandem N=3
    [result_min] = mapqn_bnd_lr_pf(params, 1, 'min');
    [result_max] = mapqn_bnd_lr_pf(params, 1, 'max');

    fprintf('  U1 bounds: [%.6f, %.6f] (MVA ref: 0.7500)\n', result_min.objective, result_max.objective);
    if result_min.exitflag > 0 && result_max.exitflag > 0 && ...
       abs(result_min.objective - 0.75) < 0.01 && abs(result_max.objective - 0.75) < 0.01
        fprintf('  PASS (tight bounds match MVA)\n');
        nPass = nPass + 1;
    else
        fprintf('  FAIL (exitflags: %d, %d)\n', result_min.exitflag, result_max.exitflag);
        nFail = nFail + 1;
    end
catch e
    fprintf('  FAIL (error: %s)\n', e.message);
    nFail = nFail + 1;
end

%% Test 2: mapqn_bnd_lr_pf - 3-queue asymmetric
fprintf('\n--- Test 2: mapqn_bnd_lr_pf (3-queue asymmetric) ---\n');
try
    params = struct();
    params.M = 3;
    params.N = 3;
    params.mu = [1.0; 2.0; 0.5];
    params.r = [0.1, 0.6, 0.3; 1.0, 0, 0; 1.0, 0, 0];
    params.verbose = false;

    [r1] = mapqn_bnd_lr_pf(params, 1, 'min');
    [r2] = mapqn_bnd_lr_pf(params, 2, 'min');
    [r3] = mapqn_bnd_lr_pf(params, 3, 'min');

    fprintf('  U1 min: %.6f, U2 min: %.6f, U3 min: %.6f\n', r1.objective, r2.objective, r3.objective);
    if r1.exitflag > 0 && r2.exitflag > 0 && r3.exitflag > 0
        fprintf('  PASS (feasible)\n');
        nPass = nPass + 1;
    else
        fprintf('  FAIL\n');
        nFail = nFail + 1;
    end
catch e
    fprintf('  FAIL (error: %s)\n', e.message);
    nFail = nFail + 1;
end

%% Test 3: mapqn_bnd_lr - 2-phase MAP queue
fprintf('\n--- Test 3: mapqn_bnd_lr (2-phase MAP queue) ---\n');
try
    params = struct();
    params.M = 2;
    params.N = 3;
    params.K = [1; 2];
    params.mu = {1.0; [0.5, 0.3; 0.2, 0.8]};
    params.v = {0.0; [0, 0.1; 0.05, 0]};
    params.r = [0, 1; 1, 0];
    params.verbose = false;

    [result_max] = mapqn_bnd_lr(params, 1, 1, 'max');
    [result_min] = mapqn_bnd_lr(params, 1, 1, 'min');

    fprintf('  U1(phase 1) bounds: [%.6f, %.6f]\n', result_min.objective, result_max.objective);
    if result_max.exitflag > 0 && result_min.exitflag > 0 && ...
       result_max.objective > 0 && result_max.objective <= 1
        fprintf('  PASS (feasible, valid bounds)\n');
        nPass = nPass + 1;
    else
        fprintf('  FAIL (exitflags: %d, %d)\n', result_min.exitflag, result_max.exitflag);
        nFail = nFail + 1;
    end
catch e
    fprintf('  FAIL (error: %s)\n', e.message);
    nFail = nFail + 1;
end

%% Test 4: mapqn_bnd_lr_mva - MVA version
fprintf('\n--- Test 4: mapqn_bnd_lr_mva (MVA 3-queue) ---\n');
try
    params = struct();
    params.M = 3;
    params.N = 5;
    params.K = 2;
    params.muM = [2.0; 1.5];
    params.muMAP = [0.5, 0.1; 0.2, 0.8];
    params.r = [0, 0.5, 0.5; 1, 0, 0; 1, 0, 0];
    params.v = [0, 0.05; 0.03, 0];
    params.verbose = false;

    [result] = mapqn_bnd_lr_mva(params, 1, 1, 'max');

    fprintf('  UN(1,1) max: %.6f\n', result.objective);
    if result.exitflag > 0 && result.objective > 0 && result.objective <= 1
        fprintf('  PASS (feasible, valid bound)\n');
        nPass = nPass + 1;
    else
        fprintf('  FAIL (exitflag: %d)\n', result.exitflag);
        nFail = nFail + 1;
    end
catch e
    fprintf('  FAIL (error: %s)\n', e.message);
    nFail = nFail + 1;
end

%% Test 5: mapqn_bnd_qr - Quadratic reduction with 2-phase
fprintf('\n--- Test 5: mapqn_bnd_qr (2-queue 2-phase) ---\n');
try
    params = struct();
    params.M = 2;
    params.N = 2;  % small N to keep problem size manageable
    params.K = [2; 2];
    params.mu = {[0.8, 0.2; 0.1, 0.6]; [0.5, 0.1; 0.2, 0.7]};
    params.v = {[0, 0.1; 0.05, 0]; [0, 0.05; 0.1, 0]};
    params.r = [0, 1; 1, 0];
    params.verbose = false;

    [result_max] = mapqn_bnd_qr(params, 1, 1, 'max');
    [result_min] = mapqn_bnd_qr(params, 1, 1, 'min');

    fprintf('  U(1,1) bounds: [%.6f, %.6f]\n', result_min.objective, result_max.objective);
    if result_max.exitflag > 0 && result_min.exitflag > 0 && ...
       result_min.objective >= 0 && result_max.objective <= 1 && ...
       result_min.objective <= result_max.objective
        fprintf('  PASS (feasible, valid bounds)\n');
        nPass = nPass + 1;
    else
        fprintf('  FAIL (exitflags: %d, %d)\n', result_min.exitflag, result_max.exitflag);
        nFail = nFail + 1;
    end
catch e
    fprintf('  FAIL (error: %s)\n', e.message);
    nFail = nFail + 1;
end

%% Test 6: mapqn_bnd_qr_ld - Load-dependent with 2-phase
fprintf('\n--- Test 6: mapqn_bnd_qr_ld (2-queue 2-phase LD) ---\n');
try
    params = struct();
    params.M = 2;
    params.N = 2;
    params.K = [2; 2];
    params.mu = {[0.8, 0.2; 0.1, 0.6]; [0.5, 0.1; 0.2, 0.7]};
    params.v = {[0, 0.1; 0.05, 0]; [0, 0.05; 0.1, 0]};
    params.alpha = ones(2, 2);  % no load-dependence
    params.r = [0, 1; 1, 0];
    params.verbose = false;

    [result] = mapqn_bnd_qr_ld(params, 1, 1, 1, 'max');

    fprintf('  p2(1,1,1,1,1,1) max: %.6f\n', result.objective);
    if result.exitflag > 0 && result.objective >= 0 && result.objective <= 1
        fprintf('  PASS (feasible, valid bound)\n');
        nPass = nPass + 1;
    else
        fprintf('  FAIL (exitflag: %d)\n', result.exitflag);
        nFail = nFail + 1;
    end
catch e
    fprintf('  FAIL (error: %s)\n', e.message);
    nFail = nFail + 1;
end

%% Test 7: mapqn_bnd_qr_delay - Delay system with IS alpha
% Note: delay model requires alpha(M,n) = n for IS (infinite server) behavior
fprintf('\n--- Test 7: mapqn_bnd_qr_delay (2-queue delay with IS alpha) ---\n');
try
    params = struct();
    params.M = 2;
    params.N = 3;
    params.K = [1; 1];
    params.Z = 2.0;   % Z = 1/mu_delay = 1/0.5 = 2.0
    params.D1 = 1.0;  % D1 = 1/mu_server = 1/1.0 = 1.0
    params.mu = {1.0; 0.5};
    params.v = {0; 0};
    params.alpha = [ones(1,3); 1:3];  % IS: alpha(M,n) = n
    params.r = [0, 1; 1, 0];
    params.verbose = false;

    [result_max] = mapqn_bnd_qr_delay(params, 1, 1, 1, 'max');
    [result_min] = mapqn_bnd_qr_delay(params, 1, 1, 1, 'min');

    fprintf('  p2(1,1,1,1,1,1) bounds: [%.6f, %.6f]\n', result_min.objective, result_max.objective);
    if result_max.exitflag > 0 && result_min.exitflag > 0 && ...
       result_min.objective >= 0 && result_max.objective <= 1
        fprintf('  PASS (feasible, valid bounds)\n');
        nPass = nPass + 1;
    else
        fprintf('  FAIL (exitflags: %d, %d)\n', result_min.exitflag, result_max.exitflag);
        nFail = nFail + 1;
    end
catch e
    fprintf('  FAIL (error: %s)\n', e.message);
    nFail = nFail + 1;
end

%% Test 8: Cross-validation - bnd_lr_pf min == bnd_lr_pf max for exact case
fprintf('\n--- Test 8: mapqn_bnd_lr_pf (exact bounds for N=2, 2-queue) ---\n');
try
    params = struct();
    params.M = 2;
    params.N = 2;
    params.mu = [1.0; 1.0];
    params.r = [0, 1; 1, 0];
    params.verbose = false;

    % MVA reference: N=2, U = 2/3 = 0.6667
    [result_min] = mapqn_bnd_lr_pf(params, 1, 'min');
    [result_max] = mapqn_bnd_lr_pf(params, 1, 'max');

    fprintf('  U1 bounds: [%.6f, %.6f] (MVA ref: 0.6667)\n', result_min.objective, result_max.objective);
    if result_min.exitflag > 0 && result_max.exitflag > 0 && ...
       abs(result_min.objective - 2/3) < 0.01 && abs(result_max.objective - 2/3) < 0.01
        fprintf('  PASS (tight bounds)\n');
        nPass = nPass + 1;
    else
        fprintf('  FAIL\n');
        nFail = nFail + 1;
    end
catch e
    fprintf('  FAIL (error: %s)\n', e.message);
    nFail = nFail + 1;
end

%% Summary
fprintf('\n====================================\n');
fprintf('RESULTS: %d PASS, %d FAIL out of %d tests\n', nPass, nFail, nPass + nFail);
fprintf('====================================\n');
