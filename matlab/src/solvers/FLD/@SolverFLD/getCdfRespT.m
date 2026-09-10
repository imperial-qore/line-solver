function RD = getCdfRespT(self, R)
% RD = GETCDFRESPT(R) Returns cumulative distribution function of response times
%
% @brief Computes steady-state response time distribution for each station and class
%
% This method returns the cumulative distribution function (CDF) of response times
% at each station for each job class. The response time represents the time from
% when a job arrives at a station until it completes service and departs.
%
% The Fluid solver uses ODE-based mean-field approximations to estimate response
% time distributions. The distributions are fitted based on the computed steady-state
% moments and assume phase-type characteristics.
%
% @param self SolverFluid instance
% @param R (optional) Response time computation handle. If provided, uses cached results.
%         Otherwise computes new CDFs by solving the passage time ODEs.
%
% @return RD Nested cell array of CDF data
%         - RD is a {nstations x nclasses} cell array
%         - RD{i,k} contains CDF data for station i, job class k
%         - Each non-empty element is [n x 2] matrix:
%           - Column 1: CDF values (cumulative probability from 0 to 1)
%           - Column 2: Response time values (increasing sequence)
%         - Empty cells indicate metric is disabled or not applicable
%
% @note Response time CDF depends on the queue length distribution and service
%       times at each station. Requires getAvg() to have been called first to
%       establish steady-state baseline metrics.
%
% @warning The Fluid approximation is mean-field based. For exact distributions,
%          use CTMC solver on small models or SSA for general models.
%          Results may not be accurate for systems with high variability or
%          strong state dependencies.
%
% @see getAvg - Computes steady-state average metrics (prerequisite)
% @see getCdfPassT - Returns first passage time distributions
%
% Example:
% @code
% model = Network('example');
% % ... model construction ...
% solver = SolverFluid(model);
% solver.getAvg();  % Compute steady-state first
%
% % Get response time CDFs
% cdf_respt = solver.getCdfRespT();
%
% % Access CDF for station 0, class 1
% if ~isempty(cdf_respt{1,2})
%     cdf_matrix = cdf_respt{1,2};
%     probs = cdf_matrix(:,1);
%     times = cdf_matrix(:,2);
%
%     % Plot response time CDF
%     plot(times, probs, 'LineWidth', 2);
%     xlabel('Response Time');
%     ylabel('Cumulative Probability');
%     grid on;
% end
% @endcode


% lang='cpp' takes the fluid passage-time law from line-cli (-s fluid -a cdf),
% which fits an APH to the passage time LSODA integrates. It is the same
% quantity this function builds below, so recomputing it here would report a
% MATLAB integration as a C++ one.
if isfield(self.options,'lang') && strcmp(self.options.lang,'cpp')
    RD = CPPLINE.cdfRespT(self.name, self.model, self.options);
    return
end

sn = self.getStruct;
if GlobalConstants.DummyMode
    RD = cell(sn.nstations,sn.nclasses);    
    return
end

T0 = tic;
if nargin<2 %~exist('R','var')
    R = self.getAvgRespTHandles;
    % to do: check if some R are disabled
end
% Force MATLAB path for passage time computation (needs odeStateVec)
origLang = self.options.lang;
self.options.lang = 'matlab';
self.result = []; % Clear cached results to force MATLAB re-computation
self.getAvg; % get steady-state solution with odeStateVec populated
self.options.lang = origLang;
options = self.getOptions;
options.init_sol = self.result.solverSpecific.odeStateVec;
% The passage-time ODE is a SECOND solve on the fixed point the mean solve
% reached, so it must be driven by the same drift. Under a moment closure that
% drift closes min(n_i,c_i) at the fixed-point variance; the first-order min()
% drains a station the mean solve holds below capacity at full rate, and the
% distribution then contradicts the mean the same solver reports. Only the
% per-STATION variance transfers: the transient class only relabels a station
% population, whereas the coordinate covariance is indexed by a state layout
% that class changes, so the class share stays the plug-in ratio.
options.config.moment_sigma2 = [];
options.config.moment_cov = {};
if isfield(self.result.solverSpecific,'moments') && ...
        isstruct(self.result.solverSpecific.moments) && ...
        isfield(self.result.solverSpecific.moments,'sigma2Drift')
    options.config.moment_sigma2 = self.result.solverSpecific.moments.sigma2Drift;
end
% we need to pass the modified sn as the number of phases may have changed
% during the fluid iterations, affecting the size of odeStateVec
RD = solver_fluid_passage_time(self.result.solverSpecific.sn, options);
runtime = toc(T0);
self.setDistribResults(RD, runtime);
end