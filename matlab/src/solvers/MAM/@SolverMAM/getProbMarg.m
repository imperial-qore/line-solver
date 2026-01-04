function [Pmarg, logPmarg] = getProbMarg(self, ist, jobclass, state_m)
% [PMARG, LOGPMARG] = GETPROBMARG(IST, JOBCLASS, STATE_M)
%
% Probability distribution for queue length of a SINGLE class at a station.
% Returns P(n jobs of class r) for n=0,1,...,N(r).
%
% Compare with getProbAggr: returns probability of a specific per-class
% distribution, e.g., P(2 class-1, 1 class-2) as a scalar.
%
% Input:
%   ist      - Station index
%   jobclass - Class index
%   state_m  - (optional) Specific states to query, or [] for all
%
% Output:
%   Pmarg    - Vector where Pmarg(n+1) = P(n jobs of this class)
%   logPmarg - Log probabilities for numerical stability

if nargin < 3
    line_error(mfilename,'getProbMarg requires station and job class parameters.');
end
if nargin < 4
    state_m = [];
end

sn = self.getStruct;

% Check if station index is valid
if ist > sn.nstations
    line_error(mfilename,'Station number exceeds the number of stations in the model.');
end
if jobclass > sn.nclasses
    line_error(mfilename,'Job class index exceeds the number of classes in the model.');
end

% Check if this is a network (more than one queue station)
queueStations = 0;
for i=1:sn.nstations
    if sn.nodetype(sn.stationToNode(i)) == NodeType.Queue
        queueStations = queueStations + 1;
    end
end

if queueStations > 1
    line_error(mfilename,'getProbMarg is currently only supported for single queue models in SolverMAM.\nFor networks, this method is not yet implemented.');
end

% Check if the model has only one queue
if queueStations == 0
    line_error(mfilename,'Model does not contain any queue stations.');
end

% Ensure results are available
if isempty(self.result)
    self.run;
end

% Get model parameters needed for QBD solution
K = sn.nclasses;
N = sn.njobs';

% Check if this is a closed model
if all(isfinite(N))
    % Closed model - compute probability distribution using QBD
    maxLevel = sum(N(isfinite(N))) + 1;

    % Build the arrival and service processes
    PH = sn.proc;
    pie = cell(1,K);
    D0 = cell(1,K);

    % Extract service process parameters
    for k=1:K
        PH{ist}{k} = map_scale(PH{ist}{k}, 1./sn.rates(ist,k)/sn.nservers(ist));
        pie{k} = map_pie(PH{ist}{k});
        D0{k} = PH{ist}{k}{1};
        if any(isnan(D0{k}))
            D0{k} = -GlobalConstants.Immediate;
            pie{k} = 1;
        end
    end

    % Build aggregate arrival process (approximation using throughput)
    % This is a simplified approach - for closed networks we use the throughput
    % to approximate the arrival process
    T = self.result.Avg.T;
    lambda_total = sum(T(ist,:));

    if lambda_total < GlobalConstants.FineTol
        % No traffic at this station
        Pmarg = zeros(1, maxLevel);
        Pmarg(1) = 1.0; % All probability at state 0
        logPmarg = log(Pmarg);
        logPmarg(Pmarg == 0) = -Inf;

        if ~isempty(state_m)
            if max(state_m) > length(Pmarg) - 1
                line_error(mfilename,'Requested state exceeds maximum population.');
            end
            Pmarg = Pmarg(state_m + 1); % +1 for MATLAB indexing
            logPmarg = logPmarg(state_m + 1);
        end
        return;
    end

    % Build a simple MMAP approximation based on class throughputs
    D_approx = cell(1, K+1);
    D_approx{1} = -lambda_total * eye(1); % D0
    for k=1:K
        D_approx{k+1} = T(ist,k) * ones(1,1); % Dk - arrivals of class k
    end

    try
        % Compute the queue length distribution using MMAPPH1FCFS
        [pdistr] = MMAPPH1FCFS(D_approx, {pie{:}}, {D0{:}}, 'ncDistr', maxLevel);

        % Extract distribution for the specific class
        % This is an approximation: we extract the portion corresponding to this class
        if jobclass <= K && N(jobclass) > 0
            % Normalize and extract portion for this class
            pdistr = abs(pdistr);
            pdistr = pdistr / sum(pdistr);

            % For the specific class, we approximate the marginal
            % by assuming independence and using the class's share
            class_share = N(jobclass) / sum(N);
            max_class_level = N(jobclass) + 1;

            Pmarg = zeros(1, max_class_level);
            for n=0:N(jobclass)
                if n+1 <= length(pdistr)
                    Pmarg(n+1) = pdistr(n+1);
                end
            end

            % Renormalize
            Pmarg = Pmarg / sum(Pmarg);
            logPmarg = log(Pmarg);
            logPmarg(Pmarg == 0) = -Inf;
        else
            line_error(mfilename,'Invalid class index or zero population for class.');
        end

        % Filter based on state_m if provided
        if ~isempty(state_m)
            if max(state_m) > length(Pmarg) - 1
                line_error(mfilename,'Requested state exceeds maximum population for this class.');
            end
            Pmarg = Pmarg(state_m + 1); % +1 for MATLAB indexing
            logPmarg = logPmarg(state_m + 1);
        end

    catch ME
        line_error(mfilename,'Failed to compute marginalized probabilities: %s', ME.message);
    end

else
    % Open model - compute probability distribution using MAM (MMAPPH1FCFS)

    % Get model parameters
    PH = sn.proc;
    pie = cell(1, K);
    D0 = cell(1, K);

    % Extract service process parameters for the queue station
    for k = 1:K
        PH{ist}{k} = map_scale(PH{ist}{k}, 1./sn.rates(ist,k)/sn.nservers(ist));
        pie{k} = map_pie(PH{ist}{k});
        D0{k} = PH{ist}{k}{1};
        if any(isnan(D0{k}))
            D0{k} = -GlobalConstants.Immediate;
            pie{k} = 1;
        end
    end

    % Build the arrival process from source
    % Find the source station (reference station for open classes)
    nChains = sn.nchains;
    V = cellsum(sn.visits);

    % Get arrival process from source and build aggregate MMAP
    refstat = sn.refstat(1);  % Source station
    PH_src = sn.proc{refstat};

    % Build arrival process cell array: {D0, D1, D2, ...} for K classes
    D_arr = cell(1, K + 1);

    % Check if all arrivals are from a single source with simple MAPs
    totalLambda = 0;
    for k = 1:K
        if ~isnan(PH_src{k}{1})
            totalLambda = totalLambda + map_lambda(PH_src{k});
        end
    end

    if totalLambda < GlobalConstants.FineTol
        % No arrivals at this station
        Pmarg = zeros(1, 100);
        Pmarg(1) = 1.0;
        logPmarg = log(Pmarg);
        logPmarg(Pmarg == 0) = -Inf;
    else
        % Build the aggregate arrival MMAP
        % For simple exponential arrivals, construct directly
        arrMaps = cell(1, K);
        for k = 1:K
            if ~isnan(PH_src{k}{1})
                arrMaps{k} = PH_src{k};
            else
                arrMaps{k} = map_exponential(Inf);  % No arrivals
            end
        end

        % Superpose all arrival processes
        if K == 1
            % Single class - use the arrival process directly
            arrProcess = arrMaps{1};
            D_arr{1} = arrProcess{1};  % D0
            D_arr{2} = arrProcess{2};  % D1
        else
            % Multiple classes - superpose MAPs
            superMAP = arrMaps{1};
            for k = 2:K
                superMAP = map_super({superMAP, arrMaps{k}});
            end
            % Build MMAP with class marking based on arrival rates
            D_arr{1} = superMAP{1};  % D0
            for k = 1:K
                lambdaK = map_lambda(arrMaps{k});
                D_arr{k + 1} = (lambdaK / totalLambda) * superMAP{2};
            end
        end

        % Compute queue length distribution using MMAPPH1FCFS
        maxLevel = 100;  % Maximum queue length to compute
        if ~isempty(self.options.cutoff) && isfinite(self.options.cutoff) && self.options.cutoff > 0
            maxLevel = self.options.cutoff;
        end

        % Check stability
        aggrUtil = 0;
        for k = 1:K
            lambdaK = map_lambda(arrMaps{k});
            if ~isnan(sn.rates(ist, k)) && sn.rates(ist, k) > 0
                aggrUtil = aggrUtil + lambdaK / (sn.rates(ist, k) * sn.nservers(ist));
            end
        end

        if aggrUtil >= 1 - GlobalConstants.FineTol
            line_warning(mfilename, 'System utilization %.2f%% is too high, distribution may be inaccurate.', aggrUtil * 100);
        end

        try
            [pdistr] = MMAPPH1FCFS(D_arr, {pie{:}}, {D0{:}}, 'ncDistr', maxLevel);
            pdistr = abs(pdistr);
            pdistr = pdistr / sum(pdistr);

            Pmarg = pdistr(:)';
            logPmarg = log(Pmarg);
            logPmarg(Pmarg == 0) = -Inf;
        catch ME
            line_error(mfilename, sprintf('Failed to compute queue length distribution: %s', ME.message));
        end
    end

    % Filter based on state_m if provided
    if ~isempty(state_m)
        if max(state_m) > length(Pmarg) - 1
            line_error(mfilename, 'Requested state exceeds maximum queue length.');
        end
        Pmarg = Pmarg(state_m + 1);
        logPmarg = logPmarg(state_m + 1);
    end
end

end
