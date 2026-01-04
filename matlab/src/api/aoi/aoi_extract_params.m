function [aoiParams, aoiInfo] = aoi_extract_params(sn, aoiInfo, options)
%AOI_EXTRACT_PARAMS Extract parameters from LINE network for aoi-fluid solvers
%
% [aoiParams, aoiInfo] = AOI_EXTRACT_PARAMS(sn, aoiInfo, options)
%
% Maps LINE model representation to aoi-fluid function inputs:
% - For bufferless (capacity=1): extract (tau, T, sigma, S, p)
% - For single-buffer (capacity=2): extract (lambda, sigma, S, r)
%
% Parameters:
%   sn (struct): Network structure from Network.getStruct()
%   aoiInfo (struct): Topology information from aoi_is_aoi()
%   options (struct): Solver options, may contain:
%       - config.aoi_preemption: Override preemption probability (0-1)
%
% Returns:
%   aoiParams (struct): Parameters for aoi-fluid solver:
%       For bufferless:
%           .tau: Arrival initial probability vector
%           .T: Arrival sub-generator matrix
%           .sigma: Service initial probability vector
%           .S: Service sub-generator matrix
%           .p: Preemption probability (0=FCFS, 1=preemptive)
%       For single-buffer:
%           .lambda: Arrival rate (Poisson)
%           .sigma: Service initial probability vector
%           .S: Service sub-generator matrix
%           .r: Replacement probability (0=FCFS, 1=replacement)
%   aoiInfo (struct): Updated with extracted parameters

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 3
    options = struct();
end

aoiParams = struct();

% Get station indices
sourceStation = aoiInfo.sourceStation;
queueStation = aoiInfo.queueStation;

% Find the open class
openClasses = find(isinf(sn.njobs));
classIdx = openClasses(1);

% Get capacity and system type
capacity = aoiInfo.capacity;
schedStrategy = aoiInfo.schedStrategy;

% =========================================================================
% EXTRACT SERVICE PROCESS (same for both bufferless and single-buffer)
% =========================================================================

serviceProc = sn.proc{queueStation}{classIdx};

if isempty(serviceProc) || (iscell(serviceProc) && isnan(serviceProc{1}(1)))
    % Simple exponential service
    mu = sn.rates(queueStation, classIdx);
    aoiParams.sigma = 1;
    aoiParams.S = -mu;
else
    % Phase-type service
    [sigma, S] = aoi_dist2ph(serviceProc);
    aoiParams.sigma = sigma;
    aoiParams.S = S;
end

% =========================================================================
% EXTRACT ARRIVAL PROCESS AND PREEMPTION PARAMETER
% =========================================================================

if capacity == 1
    % BUFFERLESS: PH/PH/1/1 or PH/PH/1/1*
    % Need: tau, T, sigma, S, p

    arrivalProc = sn.proc{sourceStation}{classIdx};

    if isempty(arrivalProc) || (iscell(arrivalProc) && isnan(arrivalProc{1}(1)))
        % Simple exponential arrival (Poisson)
        lambda = sn.rates(sourceStation, classIdx);
        aoiParams.tau = 1;
        aoiParams.T = -lambda;
    else
        % Phase-type arrival
        [tau, T] = aoi_dist2ph(arrivalProc);
        aoiParams.tau = tau;
        aoiParams.T = T;
    end

    % Determine preemption probability p
    if isfield(options, 'config') && isfield(options.config, 'aoi_preemption')
        % User-specified preemption probability
        aoiParams.p = options.config.aoi_preemption;
    else
        % Automatic based on scheduling strategy
        switch schedStrategy
            case SchedStrategy.FCFS
                aoiParams.p = 0;  % No preemption
            case SchedStrategy.LCFS
                aoiParams.p = 0;  % Non-preemptive LCFS
            case SchedStrategy.LCFSPR
                aoiParams.p = 1;  % Preemptive LCFS
            otherwise
                aoiParams.p = 0;  % Default to FCFS behavior
        end
    end

    aoiInfo.arrivalType = 'PH';

else
    % SINGLE-BUFFER: M/PH/1/2 or M/PH/1/2*
    % Need: lambda, sigma, S, r

    % Arrival must be Poisson (exponential interarrival times)
    aoiParams.lambda = sn.rates(sourceStation, classIdx);

    % Determine replacement probability r
    if isfield(options, 'config') && isfield(options.config, 'aoi_preemption')
        % User-specified replacement probability
        aoiParams.r = options.config.aoi_preemption;
    else
        % Automatic based on scheduling strategy
        switch schedStrategy
            case SchedStrategy.FCFS
                aoiParams.r = 0;  % No replacement (FCFS)
            case SchedStrategy.LCFS
                aoiParams.r = 1;  % Replacement (LCFS)
            case SchedStrategy.LCFSPR
                aoiParams.r = 1;  % Preemptive also maps to replacement
            otherwise
                aoiParams.r = 0;  % Default to FCFS behavior
        end
    end

    aoiInfo.arrivalType = 'M';
end

% Store extracted info
aoiInfo.aoiParams = aoiParams;

end
