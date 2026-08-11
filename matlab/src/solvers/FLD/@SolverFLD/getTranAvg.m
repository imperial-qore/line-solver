function [QNclass_t, UNclass_t, TNclass_t] = getTranAvg(self,Qt,Ut,Tt)
% [QNCLASS_T, UNCLASS_T, TNCLASS_T] = GETTRANAVG(SELF,QT,UT,TT)
% Returns transient mean performance metrics over time
%
% @brief Computes transient queue length, utilization, and throughput time series
%
% This method returns transient performance metrics (time-dependent averages)
% for each station and job class. The transient analysis tracks how metrics
% evolve from the initial state toward steady state over the specified time span.
%
% The Fluid solver uses ODE-based mean-field approximations to compute transient
% behavior efficiently. The method returns different data structures depending on
% whether optional handles are provided.
%
% @param self SolverFluid instance (must have timespan configured, e.g.,
%             SolverFluid(model, 'timespan', [0, 50]))
% @param Qt (optional) Queue length handle from previous computation. If provided,
%           returns cached results. Otherwise computes new QN(t).
% @param Ut (optional) Utilization handle. If provided, uses cached results.
%           Otherwise computes new UN(t).
% @param Tt (optional) Throughput handle. If provided, uses cached results.
%           Otherwise computes new TN(t).
%
% @return QNclass_t Nested cell array of queue length time series
%         - Structure: {station_1, station_2, ...} where each station contains
%           {class_1_timeseries, class_2_timeseries, ...}
%         - Each element is a vector of queue lengths at each time point
%         - Shape: [num_stations][num_classes] with each element a time-series array
%
% @return UNclass_t Nested cell array of utilization time series
%         - Same structure as QNclass_t but containing utilization values
%
% @return TNclass_t Nested cell array of throughput time series
%         - Same structure as QNclass_t but containing throughput values
%
% @note Transient analysis is only available for certain solvers (CTMC, Fluid, JMT).
%       Must configure timespan when creating solver instance:
%       solver = SolverFluid(model, 'timespan', [0, 100]);
%
% @warning The Fluid approximation is mean-field based and provides estimates
%          rather than exact transient probabilities. Use CTMC for small models
%          requiring exact transient analysis.
%
% @see getTranHandles - Get result handles for transient metrics
% @see getAvg - Get steady-state average metrics
%
% Example:
% @code
% model = Network('example');
% % ... model construction ...
% solver = SolverFluid(model, 'timespan', [0, 50]);
% [QN_t, UN_t, TN_t] = solver.getTranAvg();
%
% % Access queue length time series for station 0, class 1
% qlen_time_series = QN_t{1}{2};
% plot(qlen_time_series);
% @endcode
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

% temporarily switch to closing method
if nargin == 1
    [Qt,Ut,Tt] = self.getTranHandles;
end

options = self.options;
% Force MATLAB path for transient analysis (Java path only produces steady-state)
self.options.lang = 'matlab';

% see _kb/06-solver-catalog.md for rationale
sn = self.getStruct;
hasCache = any(sn.nodetype == NodeType.Cache);

% see _kb/06-solver-catalog.md for rationale
self.options.config.nhpp_sched = local_detect_nhpp(self.model, sn);

switch options.method
    case {'default', 'matrix', 'closing', 'rmf', 'fluid.rmf'}
        % These methods can switch to closing silently for the queueing
        % transient (RMF caches keep their transient via cacheqn_tran below).
        self.options.method = 'closing';
    case {'tbi','fluid.tbi'}
        % TBI integrates the closing ODEs by cell decomposition and
        % produces a genuine transient; keep the method as is.
    otherwise
        line_warning(mfilename,'getTranAvg is not offered by the specified method. Setting the solution method to ''''closing''''.\n');
        self.options.method = 'closing';
        self.reset();
end

if hasCache
    % The cache carries its own transient through the RMF drift.
    tranopts = self.options;
    tranopts.method = 'rmf';
    tranopts.timespan = options.timespan;
    [tcache, hitprob_t, missprob_t, caches, arate] = solver_fld_cacheqn_tran(sn, tranopts);
    self.result.CacheTran = struct('t', tcache, 'hitprob', hitprob_t, ...
        'missprob', missprob_t, 'nodes', caches, 'arate', arate);
    % see _kb/06-solver-catalog.md for rationale
    try
        [QNclass_t, UNclass_t, TNclass_t] = getTranAvg@NetworkSolver(self,Qt,Ut,Tt);
    catch
        [QNclass_t, UNclass_t, TNclass_t] = local_const_tran(self, options.timespan);
    end
else
    [QNclass_t, UNclass_t, TNclass_t] = getTranAvg@NetworkSolver(self,Qt,Ut,Tt);
end

self.options = options;
end

function sched = local_detect_nhpp(model, sn)
% Build the nhpp_sched struct array (station index in sn station space, class,
% process handle) for every EXT/source station carrying a non-homogeneous
% arrival process (identified by the getRateSchedule method, as in
% refreshProcessRepresentations). Empty when there is none.
sched = [];
for i = 1:sn.nstations
    if sn.sched(i) ~= SchedStrategy.EXT
        continue
    end
    node = model.nodes{sn.stationToNode(i)};
    if ~isa(node, 'Source')
        continue
    end
    for c = 1:sn.nclasses
        if numel(node.input.sourceClasses) < c || isempty(node.input.sourceClasses{c})
            continue
        end
        proc = node.input.sourceClasses{c}{end};
        if ismethod(proc, 'getRateSchedule')
            entry = struct('station', i, 'class', c, 'nhpp', proc);
            if isempty(sched)
                sched = entry;
            else
                sched(end+1) = entry; %#ok<AGROW>
            end
        end
    end
end
end

function [QNt, UNt, TNt] = local_const_tran(self, timespan)
% Build constant transient handles from the steady-state solution, for cache
% networks whose queueing part has no fluid dynamics to integrate.
solver = FLD(self.model, 'method', 'rmf');
[QN, UN, ~, TN] = solver.getAvg();
M = size(QN, 1); K = size(QN, 2);
t0 = timespan(1); if isinf(t0), t0 = 0; end
tcol = [t0; timespan(2)];
QNt = cell(M, K); UNt = cell(M, K); TNt = cell(M, K);
for i = 1:M
    for r = 1:K
        QNt{i,r} = [QN(i,r); QN(i,r)]; QNt{i,r}(:,2) = tcol;
        UNt{i,r} = [UN(i,r); UN(i,r)]; UNt{i,r}(:,2) = tcol;
        TNt{i,r} = [TN(i,r); TN(i,r)]; TNt{i,r}(:,2) = tcol;
    end
end
end
