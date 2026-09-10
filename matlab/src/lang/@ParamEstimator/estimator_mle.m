function [estVal,fObjFun] = estimator_mle(self, nodes)
node = nodes{1};
% MLE Maximum likelihood demand estimator
%
% Uses a LINE solver (default: SolverAuto) to compute predicted
% response times and utilizations, replacing the M/GI/1-PS formula.
% Updates service rates efficiently via sn_set_service.
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

[estVal, fObjFun] = mle_data(avgU, avgR, avgA, node.getNumberOfServers(), self.options.iter_max, self.model, node, solverType);
estVal = estVal(:)';
end

% MLE procedure using LINE solver for predicted measurements
function [demEst,fObjFun] = mle_data(cpuUtil, rAvgTimes, avgArvR, numServers, ITERMAX, model, node, solverType)
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

%% number of classes
R = size(rAvgTimes,2);

%% Initialize solver: get sn struct and solver options once
sn = model.getStruct();
stIdx = node.stationIndex;

solver = solverType(model);
solverOpts = solver.getOptions();
solverAnalyzer = getSolverAnalyzer(solver, solverOpts);

%% initial point
x0 = rand(1,R).*max(rAvgTimes);

%% options
options = optimset();
options.Display = 'off';
options.LargeScale = 'off';
options.MaxIter = ITERMAX;
options.MaxFunEvals = 1e10;
options.MaxSQPIter = 5000;
options.TolCon = 1e-8;
options.Algorithm = 'interior-point';

XLB = x0*0 + options.TolCon;
XUB = max(rAvgTimes);

%% optimization program
N = size(cpuUtil,1);
w = avgArvR./(sum(avgArvR,2)*ones(1,R));
[demEst, fObjFun]=fmincon(@objfun,x0,[],[],[],[],XLB,XUB,[],options);

    function f = objfun(x)
        % Update service rates via sn_set_service
        for c = 1:R
            sn = sn_set_service_coc(sn, stIdx, c, 1/x(c));
        end

        % Solve model to get predicted response times and utilization
        [~, U_pred, R_pred] = solverAnalyzer(sn);

        % Predicted response times and utilization at the target station
        predR = R_pred(stIdx, 1:R);
        predU = sum(U_pred(stIdx, 1:R));

        % Compute objective: weighted response time error + utilization error
        deltaj = repmat(predR, N, 1) - rAvgTimes;
        epsi = predU * ones(N, 1) - cpuUtil;
        f = 0;
        for i = 1:N
            f = f + w(i,:) .* deltaj(i,:).^2;
        end
        f = sum(f(:));
        for i = 1:N
            f = f + epsi(i).^2;
        end
    end

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
