function [runtime, tranSysState, tranSync] = runAnalyzer(self, options)
% [RUNTIME, TRANSYSSTATE] = RUNANALYZER()
%
% Run the LDES solver on the queueing network model.
%
% Communication with the LDES engine is fully JSON-mediated: the model is
% serialized to model.json (linemodel_save) and the engine is run as a
% subprocess ("solve model.json -o result.json"); the result JSON is parsed
% back here. No in-process Java object marshalling (JPype/JLINE) is used for
% regular Network models. LayeredNetwork (LQN) models are still simulated by
% the Java LDES ensemble backend (self.obj), which is a separate path.

T0 = tic;
if nargin < 2
    options = self.getOptions;
end

% options.events (DES event budget) overrides options.samples when set;
% samples remains accepted as a deprecated alias for the event budget.
if isfield(options,'events') && ~isempty(options.events) && isfinite(options.events) && options.events > 0
    options.samples = round(options.events);
end

tranSysState = [];
tranSync = [];

if isa(self.model, 'LayeredNetwork')
    % The gate's own sentence before the backend's: ldes_ln_refusal is what
    % supports() asks, so a refusal reads the same here as in model.help.
    [okln, whyln] = ldes_ln_refusal(self.model);
    if ~okln
        line_error(mfilename, whyln);
    end
    self.obj.getAvg(); % runs the LN LDES analyzer (Java ensemble backend)
    runtime = toc(T0);
    return
end

verboseGuard = self.runAnalyzerChecks(options); %#ok<NASGU> restores the caller verbosity on return
Solver.resetRandomGeneratorSeed(options.seed);

line_debug('LDES solver starting: samples=%d, seed=%d', options.samples, options.seed);

sn = self.model.getStruct;
M = sn.nstations;
R = sn.nclasses;
[confintEnabled, ~] = Solver.parseConfInt(options.confint);

isTransient = isfield(options, 'timespan') && numel(options.timespan) >= 2 && ...
              isfinite(options.timespan(2));

extraFlags = {};
if isTransient
    extraFlags = {'--timespan', ...
        sprintf('%.10g,%.10g', options.timespan(1), options.timespan(2)), ...
        '--trajectory'};
    % Ensemble transient: a single realization is not the transient mean
    % E[N](t) (no time-ergodicity at fixed t). When replications > 1 is
    % requested, drive the engine's parallel analyzer via --replications so
    % the per-bucket QNt/UNt/TNt are averaged across independent replications.
    % Read from options.replications or options.config.replications; default
    % (absent/<=1) keeps the single-path run.
    reps = 0;
    if isfield(options, 'replications') && ~isempty(options.replications)
        reps = options.replications;
    elseif isfield(options, 'config') && isfield(options.config, 'replications') ...
            && ~isempty(options.config.replications)
        reps = options.config.replications;
    end
    if reps > 1
        extraFlags{end+1} = '--replications';
        extraFlags{end+1} = sprintf('%d', round(reps));
        if isfield(options, 'config') && isfield(options.config, 'numthreads') ...
                && ~isempty(options.config.numthreads) && options.config.numthreads > 0
            extraFlags{end+1} = '--numthreads';
            extraFlags{end+1} = sprintf('%d', round(options.config.numthreads));
        end
    end
end

LineConsole.step('serializing the model to JSON and running the LDES engine');
[data, engine] = self.solveCli(options, extraFlags);
LineConsole.substep('engine that answered: %s', engine);
% NAME THE ENGINE THAT ANSWERED. The constructor pins options.lang='java' to
% keep this solver off the JLINE/CPPLINE dispatch (getAvg.m and getAvgTable.m
% exclude SolverLDES by name for the same reason), but the banner prints
% options.lang, so every solve claimed 'java' while common/ldes -- the NATIVE
% C++ engine since 2026-08-01 -- was what actually ran. Set after the solve, so
% the dispatch decisions upstream are untouched.
self.options.lang = engine;
if ~isstruct(data) || ~isfield(data, 'metrics')
    line_error(mfilename, 'LDES engine returned no metrics.');
end

LineConsole.step('parsing the LDES result document');
if isTransient
    parseTransient(self, data, M, R, options);
else
    parseSteady(self, data, sn, M, R, confintEnabled, options);
end

runtime = toc(T0);
end

% =========================================================================

function parseSteady(self, data, sn, M, R, confintEnabled, options)
% PARSESTEADY Marshal steady-state result JSON onto the solver.

QN = ldesJson2mat(data.metrics.QN, M, R);
UN = ldesJson2mat(data.metrics.UN, M, R);
RN = ldesJson2mat(data.metrics.RN, M, R);
TN = ldesJson2mat(data.metrics.TN, M, R);
AN = ldesJson2mat(data.metrics.AN, M, R);
WN = ldesJson2mat(data.metrics.WN, M, R);
CN = ldesJson2mat(data.metrics.CN, 1, R);
XN = ldesJson2mat(data.metrics.XN, 1, R);
if isempty(QN)
    line_error('runAnalyzer', 'LDES result metrics could not be parsed.');
end

% A closed chain must conserve its population in the reported queue lengths.
% Engine builds that drop blocked closed jobs at a finite-capacity station
% (pre-guard ldes.jar) leak population silently, and the reported means then
% describe a different model. Surface it here, at the wrapper, so a stale
% engine cannot fail silently. FCR models are excluded: a WAITQ waiter is
% parked outside the station rows and would trip the check spuriously.
if sn.nregions == 0
    for nc = 1:sn.nchains
        inchain = find(sn.chains(nc,:));
        njobs_chain = sum(sn.njobs(inchain));
        if isfinite(njobs_chain) && njobs_chain > 0
            qtot = sum(sum(QN(:,inchain), 'omitnan'));
            if abs(qtot - njobs_chain) > max(0.05*njobs_chain, 0.05)
                line_warning(mfilename, sprintf(['LDES returned %.3f jobs for closed chain %d, whose population is %g: ' ...
                    'the engine dropped blocked jobs at a finite-capacity station, so these averages describe a different ' ...
                    'model. Use SolverCTMC/SolverSSA for capped closed networks, or rebuild common/ldes.jar.\n'], ...
                    qtot, nc, njobs_chain));
            end
        end
    end
end

% Fork-Join quorum sibling-drop rate (station-indexed), folded into
% getAvgLossTable at the Join rows. Absent unless the model has a quorum Join.
DropRateJoin = ldesJson2mat(ldesGetField(data.metrics, 'DropRateJoin', []), M, R);
if ~isempty(DropRateJoin)
    self.result.DropRateJoin = DropRateJoin;
end

if isfield(data, 'runtime') && isscalar(data.runtime) && isnumeric(data.runtime)
    runtime = data.runtime;
else
    runtime = 0;
end

% Finite capacity region (FCR) metrics: append region rows after the stations.
F = sn.nregions;
WeightNfcr = [];
MemOccNfcr = [];
if F > 0 && isfield(data, 'fcr')
    fcr = data.fcr;
    QNfcr = ldesJson2mat(ldesGetField(fcr, 'QNfcr', []), F, R);
    UNfcr = ldesJson2mat(ldesGetField(fcr, 'UNfcr', []), F, R);
    RNfcr = ldesJson2mat(ldesGetField(fcr, 'RNfcr', []), F, R);
    TNfcr = ldesJson2mat(ldesGetField(fcr, 'TNfcr', []), F, R);
    WNfcr = ldesJson2mat(ldesGetField(fcr, 'WNfcr', []), F, R);
    if ~isempty(QNfcr)
        ANfcr = NaN(F, R);  % arrival rate not applicable to FCR rows
        QN = [QN; QNfcr];
        UN = [UN; UNfcr];
        RN = [RN; RNfcr];
        TN = [TN; TNfcr];
        WN = [WN; WNfcr];
        AN = [AN; ANfcr];
        WeightNfcr = ldesJson2mat(ldesGetField(fcr, 'WeightNfcr', []), F, R);
        MemOccNfcr = ldesJson2mat(ldesGetField(fcr, 'MemOccNfcr', []), F, R);
        if isempty(WeightNfcr), WeightNfcr = zeros(F, R); end
        if isempty(MemOccNfcr), MemOccNfcr = zeros(F, R); end
        % Region carried throughput and drop rate for getAvgRegionLossTable.
        DropRateNfcr = ldesJson2mat(ldesGetField(fcr, 'DropRateNfcr', []), F, R);
        if isempty(DropRateNfcr), DropRateNfcr = zeros(F, R); end
        self.result.FCR.TNfcr = TNfcr;
        self.result.FCR.DropRateNfcr = DropRateNfcr;
    end
end

if options.verbose
    if LineConsole.isActive()
        LineConsole.substep('sample path of %d samples drawn', options.samples);
    else
        line_printf('LDES samples: %8d\n', options.samples);
    end
end

self.setAvgResults(QN, UN, RN, TN, AN, WN, CN, XN, runtime, options.method, options.samples);

% FCR-only weighted-occupation (Total Weight) and memory-occupation metrics:
% NaN at station rows, populated at the FCR rows (M+1 .. M+F).
self.result.Avg.Weight = NaN(M + F, R);
self.result.Avg.MemOcc = NaN(M + F, R);
if F > 0 && ~isempty(WeightNfcr)
    self.result.Avg.Weight(M+1:M+F, :) = WeightNfcr;
    self.result.Avg.MemOcc(M+1:M+F, :) = MemOccNfcr;
end

% Marshal cache hit/miss/latency from the result JSON onto the Cache nodes.
% Clear any stale split first, then set each field independently so a single
% missing field does not discard the others (order-independence).
haveCacheMetrics = isfield(data, 'cacheMetrics');
for ind = 1:sn.nnodes
    if sn.nodetype(ind) == NodeType.Cache
        mcache = self.model.nodes{ind};
        mcache.setResultHitProb(sparse([]));
        mcache.setResultMissProb(sparse([]));
        mcache.setResultDelayedHitProb(sparse([]));
        mcache.setResultResidT(sparse([]));
        mcache.setResultListCost([]);
        if haveCacheMetrics
            cname = matlab.lang.makeValidName(char(mcache.getName()));
            if isfield(data.cacheMetrics, cname)
                cm = data.cacheMetrics.(cname);
                hp  = ldesJson2mat(ldesGetField(cm, 'hit', []), [], []);
                mp  = ldesJson2mat(ldesGetField(cm, 'miss', []), [], []);
                dhp = ldesJson2mat(ldesGetField(cm, 'delayed', []), [], []);
                lp  = ldesJson2mat(ldesGetField(cm, 'latency', []), [], []);
                hpl = ldesJson2mat(ldesGetField(cm, 'hitList', []), [], []);
                ip  = ldesJson2mat(ldesGetField(cm, 'itemProb', []), [], []);
                lc  = ldesJson2mat(ldesGetField(cm, 'listCost', []), [], []);
                if ~isempty(hp),  mcache.setResultHitProb(hp);         end
                if ~isempty(mp),  mcache.setResultMissProb(mp);        end
                if ~isempty(dhp), mcache.setResultDelayedHitProb(dhp); end
                if ~isempty(lp),  mcache.setResultResidT(lp);          end
                if ~isempty(hpl), mcache.setResultHitProbList(hpl);    end
                if ~isempty(ip),  mcache.setResultItemProb(ip);        end
                if ~isempty(lc),  mcache.setResultListCost(lc);        end
            end
        end
    end
end

% Confidence intervals (the result JSON always carries them when requested).
if confintEnabled && isfield(data, 'confidenceIntervals')
    ci = data.confidenceIntervals;
    QNCI = ldesJson2mat(ldesGetField(ci, 'QNCI', []), M, R);
    UNCI = ldesJson2mat(ldesGetField(ci, 'UNCI', []), M, R);
    RNCI = ldesJson2mat(ldesGetField(ci, 'RNCI', []), M, R);
    TNCI = ldesJson2mat(ldesGetField(ci, 'TNCI', []), M, R);
    ANCI = ldesJson2mat(ldesGetField(ci, 'ANCI', []), M, R);
    WNCI = ldesJson2mat(ldesGetField(ci, 'WNCI', []), M, R);
    self.setAvgResultsCI(QNCI, UNCI, RNCI, TNCI, ANCI, WNCI, [], []);
    % HOW LONG THE RUN SHOULD HAVE BEEN, when the caller asked for it. The
    % engine's half-width at the configured confidence over the events it ran
    % pins the ASYMPTOTIC variance, which is the quantity a run length is
    % planned from -- not the stationary variance, which on M/M/1 differs from
    % it by a factor blowing up like (1-rho)^-2.
    if isfield(options,'config') && isstruct(options.config) ...
            && isfield(options.config,'runLengthPlan') && ~isempty(options.config.runLengthPlan)
        spec = options.config.runLengthPlan;
        relprecision = 0.05;
        [~, confidence] = Solver.parseConfInt(options.confint);
        if isstruct(spec)
            if isfield(spec,'relprecision') && ~isempty(spec.relprecision)
                relprecision = spec.relprecision;
            end
            if isfield(spec,'confidence') && ~isempty(spec.confidence)
                confidence = spec.confidence;
            end
        elseif isnumeric(spec) && isscalar(spec) && spec > 0
            relprecision = spec;
        end
        % The ACTUAL number of events simulated where the engine reports it,
        % not the budget: LDES stops early on convergence, and planning from a
        % budget it never spent would overstate N and so overstate sigma^2.
        used = 0;
        if isfield(data,'totalSimulatedEvents') && data.totalSimulatedEvents > 0
            used = data.totalSimulatedEvents;
        elseif isfield(data,'events') && data.events > 0
            used = data.events;
        elseif isfield(options,'samples')
            used = options.samples;
        end
        if used > 0
            self.result.runLengthPlan = sim_runlength_plan(self.result.Avg.Q, QNCI, used, ...
                'relprecision', relprecision, 'confidence', confidence);
        end
    end
end
end

% =========================================================================

function parseTransient(self, data, M, R, options)
% PARSETRANSIENT Marshal the transient trajectory JSON onto the solver.
% data.transient carries t (time vector) and QNt/UNt/TNt as nested
% [station][class] arrays, each an (nTimePoints x 2) [value, time] matrix.

runtime = 0;
if isfield(data, 'runtime') && isscalar(data.runtime) && isnumeric(data.runtime)
    runtime = data.runtime;
end

if ~isfield(data, 'transient') || isempty(data.transient) || ~isfield(data.transient, 't')
    if options.verbose
        line_printf('LDES transient analysis produced no trajectory.\n');
    end
    return;
end
tran = data.transient;

QNt = ldesTrajCell(ldesGetField(tran, 'QNt', []), M, R);
UNt = ldesTrajCell(ldesGetField(tran, 'UNt', []), M, R);
TNt = ldesTrajCell(ldesGetField(tran, 'TNt', []), M, R);
RNt = cell(M, R);  % response-time transient not computed by LDES
CNt = cell(1, R);
XNt = cell(1, R);
% Empty cells become scalar NaN so downstream indexing (getTranAvg) is uniform.
for ist = 1:M
    for r = 1:R
        if isempty(QNt{ist, r}), QNt{ist, r} = NaN; end
        if isempty(UNt{ist, r}), UNt{ist, r} = NaN; end
        if isempty(TNt{ist, r}), TNt{ist, r} = NaN; end
        RNt{ist, r} = NaN;
    end
end
for r = 1:R
    CNt{1, r} = NaN;
    XNt{1, r} = NaN;
end

self.setTranAvgResults(QNt, UNt, RNt, TNt, CNt, XNt, runtime);

% Steady-state estimates from the final transient values.
QN = NaN(M, R); UN = NaN(M, R); RN = NaN(M, R); TN = NaN(M, R);
WN = NaN(M, R); AN = NaN(M, R); CN = NaN(1, R); XN = NaN(1, R);
for ist = 1:M
    for r = 1:R
        if ~isscalar(QNt{ist, r}) && ~isnan(QNt{ist, r}(end, 1))
            QN(ist, r) = QNt{ist, r}(end, 1);
        end
        if ~isscalar(UNt{ist, r}) && ~isnan(UNt{ist, r}(end, 1))
            UN(ist, r) = UNt{ist, r}(end, 1);
        end
        if ~isscalar(TNt{ist, r}) && ~isnan(TNt{ist, r}(end, 1))
            TN(ist, r) = TNt{ist, r}(end, 1);
        end
    end
end
self.setAvgResults(QN, UN, RN, TN, AN, WN, CN, XN, runtime, options.method, 0);

if options.verbose
    line_printf('LDES transient analysis complete, timespan = [%g, %g]\n', ...
        options.timespan(1), options.timespan(2));
end
end

% =========================================================================
