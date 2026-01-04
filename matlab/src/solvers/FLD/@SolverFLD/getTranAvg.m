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
switch options.method
    case {'default', 'matrix', 'closing'}
        % These methods can switch to closing silently for transient analysis
        self.options.method = 'closing';
    otherwise
        line_warning(mfilename,'getTranAvg is not offered by the specified method. Setting the solution method to ''''closing''''.\n');
        self.options.method = 'closing';
        self.reset();
end

[QNclass_t, UNclass_t, TNclass_t] = getTranAvg@NetworkSolver(self,Qt,Ut,Tt);
self.options = options;
end
