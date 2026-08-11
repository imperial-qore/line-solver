function [estVal] = estimator_ekf(self, nodes)
node = nodes{1};
% Extended Kalman Filter resource demand estimator
% This demand estimator is based on the method proposed in:
%
% Kumar, Dinesh and Tantawi, Asser and Zhang, Li
% Real-Time Performance Modeling for Adaptive Software Systems 2009
%
% Uses a LINE solver (default: SolverAuto) to compute the predicted
% measurement function h(x), replacing the M/GI/1-PS analytical formula.
% The Jacobian is computed numerically via finite differences.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
% This code is released under the 3-Clause BSD License.

% This estimator at present can handle estimation on a single resource.

%%
% rescale utilization to be mean number of busy servers
    sn = self.model.getStruct;

    if isfinite(node.getNumberOfServers())
        U = self.getAggrUtil(node);
        if ~isempty(U)
            avgU = U.data * node.getNumberOfServers();
        end
    end

    % obtain per class metrics
    for r=1:sn.nclasses
        avgArvR{r} = self.getArvR(node, self.model.classes{r});
        if isempty(avgArvR{r})
            error('Arrival rate data for node %d in class %d is missing.', self.model.getNodeIndex(node), r);
        else
            avgArvR{r} = avgArvR{r}.data;
        end
        avgRespT{r} = self.getRespT(node, self.model.classes{r});
        if isempty(avgRespT{r})
            error('Response time data for node %d in class %d is missing.', self.model.getNodeIndex(node), r);
        else
            avgRespT{r} = avgRespT{r}.data;
        end

    end

    try
        avgA = cell2mat(avgArvR);
        avgR = cell2mat(avgRespT);
    catch me
        switch me.identifier
            case 'MATLAB:catenate:dimensionMismatch'
                error('Sampled metrics have different number of samples, use interpolate() before starting this estimation algorithm.');
        end
    end

    % Get solver class from options (default: @SolverAuto)
    if isfield(self.options, 'solver') && ~isempty(self.options.solver)
        solverType = self.options.solver;
    else
        solverType = @SolverAuto;
    end

    % Support warm-start from a previous estimate
    if isfield(self.options, 'x0') && ~isempty(self.options.x0)
        x0 = self.options.x0;
    else
        x0 = [];
    end

    [estVal] = ekf_data(avgU, avgR, avgA, node.getNumberOfServers(), self.options.iter_max, self.model, node, solverType, x0);
    estVal = estVal(:)';
end

% ekf procedure based on the common data format
function [demEst] = ekf_data(cpuUtil, rAvgTimes, avgArvR, numServers, ITERMAX, model, node, solverType, x0)
    a = isnan(cpuUtil);
    if sum(a) > 0
        disp('NaN values found for CPU Utilization. Removing NaN values.');
        cpuUtil = cpuUtil(a == 0);
        rAvgTimes = rAvgTimes(a == 0,:);
        avgArvR = avgArvR(a == 0,:);
    end

    a = sum(avgArvR,2) == 0;
    if sum(a) > 0
        disp('Removing sampling intervals with zero throughput for all request types.');
        cpuUtil = cpuUtil(a == 0);
        rAvgTimes = rAvgTimes(a == 0,:);
        avgArvR = avgArvR(a == 0,:);
    end

    %% number of resources
    M = 1;
    %% number of classes
    R = size(rAvgTimes,2);

    %% initial state
    % x(r) is the mean service demand of class r (visits are assumed unitary)
    if nargin >= 9 && ~isempty(x0)
        x = x0(:);
    else
        x = rand(R,1).*max(rAvgTimes); % randomize service demand in [0,max(avgRTime)] for each class
    end
    p = diag(x.^2);

    mCovNoise = diag(repmat(0.01, 1, R + 1));
    pCovNoise = eye(size(x,1)).*0.001;

    %% Initialize solver: get sn struct and solver options once
    % Use sn_set_service for fast parameter updates (no model rebuild)
    sn = model.getStruct();
    stIdx = node.stationIndex;

    % Determine solver analyzer function from solver type
    solver = solverType(model);
    solverOpts = solver.getOptions();
    solverAnalyzer = getSolverAnalyzer(solver, solverOpts);

    %% optimization program
    N = size(cpuUtil,1); % number of experiments
    stepBound = 0.6;
    a_min = zeros(size(x));
    a_max = zeros(size(x)) + inf;

    for n=1:min(N, ITERMAX)
        % predict
        % xn = Fn xn-1
        x_n = x;

        % Pn = Fn Pn-1 F'n + Qn
        P_n = p + pCovNoise;

        % update
        stepUtil = cpuUtil(n);
        stepResponse = rAvgTimes(n, :);

        [z_n, sn] = getPredictedMeasurement(x_n, R, sn, stIdx, solverAnalyzer);
        H_n = getJacobian(x_n, R, sn, stIdx, solverAnalyzer, z_n);
        z = getMeasurement(R, stepUtil, stepResponse);
        y_n = z - z_n;

        % Sn =  Hn Pn H'n + Rn
        H_nT = H_n';
        S_n = H_n * P_n * H_nT + mCovNoise;

        % Kalman Gain Kn = Pn H'n S^-1n
        K_n = P_n * H_nT * inv(S_n);

        % xnn = xn  + Kn y~n
        x = x_n + K_n * y_n;
        xlower = (stepBound * a_min) + (1-stepBound)*x;
        xupper = (stepBound * a_max) + (1-stepBound)*x;
        x = min(xupper, max(xlower, x));
        if (sum(x) < 0)
            x = x * -1;
        end
        % Pnn = (I - Kn Hn) Pn
        p = (eye(size(x)) - (K_n * H_n))  * P_n;

    end
    [demEst] = x(1:M);
end

function solverAnalyzer = getSolverAnalyzer(solver, solverOpts)
% Return a function handle that takes sn and returns [Q,U,R,T].
% Selects the appropriate low-level analyzer based on solver type,
% calling the solver's analyzer function directly to avoid OOP overhead.
    if isa(solver, 'SolverMVA')
        solverAnalyzer = @(sn_arg) solver_mva_analyzer(sn_arg, solverOpts);
    elseif isa(solver, 'SolverAUTO') || isa(solver, 'SolverAuto')
        mvaOpts = SolverMVA.defaultOptions;
        solverAnalyzer = @(sn_arg) solver_mva_analyzer(sn_arg, mvaOpts);
    elseif isa(solver, 'SolverNC')
        solverAnalyzer = @(sn_arg) solver_nc_analyzer(sn_arg, solverOpts);
    elseif isa(solver, 'SolverCTMC')
        solverAnalyzer = @(sn_arg) solver_ctmc_analyzer(sn_arg, solverOpts);
    elseif isa(solver, 'SolverFluid')
        solverAnalyzer = @(sn_arg) solver_fluid_analyzer(sn_arg, solverOpts);
    elseif isa(solver, 'SolverMAM')
        solverAnalyzer = @(sn_arg) solver_mam_analyzer(sn_arg, solverOpts);
    else
        % Fallback: wrap via the solver object
        solverAnalyzer = @(sn_arg) solverAnalyzerWrapper(solver, sn_arg);
    end
end

function [Q,U,R,T] = solverAnalyzerWrapper(solver, sn)
% Generic wrapper: update the model's cached sn and re-solve.
    solver.model.sn = sn;
    solver.result = [];
    Q = solver.getAvgQLen();
    U = solver.getAvgUtil();
    R = solver.getAvgRespT();
    T = solver.getAvgTput();
end

function [Hx] = getJacobian(x, R, sn, stIdx, solverAnalyzer, h0)
% Compute Jacobian numerically via finite differences.
    delta = 1e-6;
    Hx = zeros(R+1, R);
    for c = 1:R
        x_pert = x;
        x_pert(c) = x_pert(c) + delta;
        % Update sn for perturbation
        sn_pert = sn;
        for cc = 1:R
            sn_pert = sn_set_service_coc(sn_pert, stIdx, cc, 1/x_pert(cc));
        end
        [~,U_pert,R_pert] = solverAnalyzer(sn_pert);
        h_pert = zeros(R+1, 1);
        for cc = 1:R
            h_pert(cc) = R_pert(stIdx, cc);
        end
        h_pert(R+1) = sum(U_pert(stIdx, :));
        Hx(:, c) = (h_pert - h0) / delta;
    end
end

function [h] = getMeasurement(R, util, responseTimes)
    h = zeros(R+1,1);
    for c=1:R
        h(c) = responseTimes(c);
    end
    h(R+1) = util;
end

function [h, sn] = getPredictedMeasurement(x, R, sn, stIdx, solverAnalyzer)
% Compute predicted measurement using a LINE solver analyzer.
% Updates service rates in the sn struct via sn_set_service and solves.
    for c = 1:R
        sn = sn_set_service_coc(sn, stIdx, c, 1/x(c));
    end

    [~, U, RN] = solverAnalyzer(sn);

    h = zeros(R+1, 1);
    for c = 1:R
        h(c) = RN(stIdx, c);
    end
    h(R+1) = sum(U(stIdx, :));
end
