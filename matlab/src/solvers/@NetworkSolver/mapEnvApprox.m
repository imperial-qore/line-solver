function out = mapEnvApprox(self, options)
% OUT = MAPENVAPPROX(OPTIONS)
%
% Solver-agnostic random-environment approximation of a network with
% MAP/MMPP/MMAP arrival or service processes, for solvers that cannot consume
% a non-renewal process natively.
%
% MAP2RENV turns each modulated process into a set of environment stages in
% which that process is exponential with the phase-conditional intensity, and
% SolverENV recombines the stages. The stage models are exponential and carry
% the base model's structure unchanged, so the stage solver is the caller
% itself: nothing here reads a solver internal, and any NetworkSolver whose
% analyzer honours the getAvg contract can be driven through it.
%
% Three stage recombinations are available, and which ones a given solver can
% reach is decided by its transient capability:
%
%   'meanfield' - each stage is integrated over its sojourn from an entry
%                 marginal that mixes its predecessors' exit marginals, so
%                 the queue state is CARRIED across a phase switch instead of
%                 being restarted. Requires getTranAvg on the stage solver,
%                 i.e. supportsTransientAnalysis (FLD, CTMC, LDES, JMT);
%   'dec'       - quasi-stationary (slow-environment) limit, exact as the phase
%                 process slows down relative to the queueing dynamics;
%   'avg'       - rate-averaged (fast-environment) limit, exact as the phase
%                 process speeds up relative to the queueing dynamics.
%
% OPTIONS.config.map_env_method selects one; 'auto' (the default) takes
% 'meanfield' whenever the solver supports transient analysis, since it is the
% only recombination that models the phase switch rather than a limit of it,
% and otherwise compares the mean stage holding time with the model relaxation
% time to pick the limit whose regime the model is in.
%
% OUT carries QN, UN, RN, TN, CN, XN, WN, runtime, method and actualmethod, and
% the averages are also stored on the solver through finalizeAvgResults.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

Tstart = tic;

[envModel, info] = map2renv(self.model, options);

envMethod = 'auto';
if isfield(options,'config') && isfield(options.config,'map_env_method') ...
        && ~isempty(options.config.map_env_method)
    envMethod = options.config.map_env_method;
end
if strcmpi(envMethod,'auto')
    if self.supportsTransientAnalysis()
        envMethod = 'meanfield';
    else
        envMethod = selectEnvLimit(self.model, envModel);
    end
end
if ~any(strcmpi(envMethod, {'dec','avg','meanfield'}))
    line_error(mfilename, sprintf(['options.config.map_env_method=''%s'' is not a supported environment ' ...
        'recombination. Use ''meanfield'', ''dec'', ''avg'' or ''auto''.'], envMethod));
end
if strcmpi(envMethod,'meanfield') && ~self.supportsTransientAnalysis()
    line_error(mfilename, sprintf(['The mean-field environment coupling integrates each stage over its ' ...
        'sojourn, so it needs transient averages from the stage solver, which %s does not produce. Use ' ...
        'options.config.map_env_method=''dec'' or ''avg''.'], class(self)));
end

% Stage solvers are instances of the calling solver. The recursion guard is
% redundant on a correct transformation (no stage model carries a MAP) but
% keeps a mis-detected process from re-entering this driver.
innerOptions = options;
innerOptions.config.map_env = 'off';
if strcmpi(envMethod,'meanfield')
    % The stage transients are weighted by the sojourn density over the
    % integration grid, so the horizon must cover the sojourn distribution;
    % beyond it the weights vanish and the extra span is inert.
    innerOptions.timespan = [0, 20 * maxHoldTime(envModel)];
end
solverClass = class(self);
solverFactory = @(stageModel) feval(solverClass, stageModel, innerOptions);

envOptions = SolverENV.defaultOptions;
envOptions.method = lower(envMethod);
if strcmpi(envMethod,'meanfield')
    envOptions.method = 'default'; % SolverENV selects the mean-field analyzer
end
envOptions.verbose = options.verbose;
envOptions.iter_max = options.iter_max;
envOptions.iter_tol = options.iter_tol;

line_debug(options, 'MAP/MMPP random-environment approximation: %d stages, limit ''%s'', %s image.', ...
    info.nstages, lower(envMethod), mapImageKind(info));

envSolver = SolverENV(envModel, solverFactory, envOptions);
[QN, UN, ~, TN] = envSolver.getAvg();

sn = self.model.getStruct();
M = sn.nstations;
K = sn.nclasses;
RN = zeros(M,K);
RN(TN > 0) = QN(TN > 0) ./ TN(TN > 0);
XN = zeros(1,K);
CN = zeros(1,K);
for k = 1:K
    XN(k) = TN(sn.refstat(k), k);
    if XN(k) > 0
        CN(k) = sum(QN(:,k)) / XN(k);
    end
end
WN = sn_get_residt_from_respt(sn, RN, self.getAvgResidTHandles());

% The system metrics read the reference station, which for an open class is
% the Source. A stage solver whose transient does not report Source
% throughput (SolverFluid) leaves it at zero under the mean-field coupling,
% and the zero propagates into XN and CN. Report that rather than
% substituting the arrival rate, which would hide whose metric is missing.
openZero = find(isinf(sn.njobs(:)') & XN == 0 & any(QN > 0, 1));
if ~isempty(openZero)
    line_warning(mfilename, ['The %s environment coupling returned no reference-station throughput for open ' ...
        'class(es) %s, so their system throughput and system response time are reported as zero. This stage ' ...
        'solver does not measure Source throughput in transient mode; use options.config.map_env_method=''dec'' ' ...
        'or ''avg'' for system-level metrics.\n'], lower(envMethod), mat2str(openZero));
end

runtime = toc(Tstart);
actualmethod = ['env.', lower(envMethod)];
self.finalizeAvgResults(QN, UN, RN, TN, CN, XN, runtime, options.method, NaN, actualmethod, [], WN);

line_warning(mfilename, ['This solver has no native support for the non-renewal (MAP/MMPP) processes of this model; ' ...
    'the reported averages come from its %s random-environment approximation (%d stages, %s image). ' ...
    'Set options.config.map_env=''off'' to reject the model instead.\n'], ...
    lower(envMethod), info.nstages, mapImageKind(info));

out = struct('QN',QN,'UN',UN,'RN',RN,'TN',TN,'CN',CN,'XN',XN,'WN',WN, ...
    'runtime',runtime,'method',options.method,'actualmethod',actualmethod,'info',info);
end

function kind = mapImageKind(info)
% An MMPP image preserves the modulating chain exactly; a general MAP image
% aggregates the event-epoch phase jumps into the phase generator.
if info.isMMPP
    kind = 'exact-modulation';
else
    kind = 'intensity-matched';
end
end

function T = maxHoldTime(envModel)
% Longest mean stage sojourn of the environment.
E = numel(envModel.holdTime);
T = 0;
for e = 1:E
    T = max(T, map_mean(envModel.holdTime{e}));
end
if ~isfinite(T) || T <= 0
    line_error(mfilename, 'The environment image has no finite stage sojourn, so no integration horizon can be set.');
end
end

function limit = selectEnvLimit(model, envModel)
% Timescale test: compare the mean stage holding time of the environment with
% the relaxation time of the model, taken as the time the slowest station needs
% to clear the jobs it can hold. A stage that outlives the relaxation time lets
% each stage reach its own steady state, which is the quasi-stationary regime
% ('dec'); a stage that expires first leaves the model responding to the mean
% rate only, which is the rate-averaged regime ('avg'). The closed population
% enters the relaxation time because a closed queue drains in N services, so a
% populous model relaxes far more slowly than one service time.
E = height(envModel.envGraph.Nodes);
exitRate = zeros(1,E);
for e = 1:E
    for h = 1:E
        if e ~= h && ~isa(envModel.env{e,h}, 'Disabled')
            exitRate(e) = exitRate(e) + envModel.env{e,h}.getRate();
        end
    end
end
exitRate = exitRate(exitRate > 0);
if isempty(exitRate)
    limit = 'dec'; % absorbing environment: every stage is its own steady state
    return
end
tauEnv = mean(1 ./ exitRate);

sn = model.getStruct();
rates = sn.rates(:);
rates = rates(isfinite(rates) & rates > 0);
if isempty(rates)
    limit = 'dec';
    return
end
njobs = sn.njobs(:);
njobs = njobs(isfinite(njobs));
tauSys = (1 + sum(njobs)) / min(rates);

if tauEnv >= tauSys
    limit = 'dec';
else
    limit = 'avg';
end
end
