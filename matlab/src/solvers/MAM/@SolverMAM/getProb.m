function Pstate = getProb(self, node, state)
% PSTATE = GETPROB(NODE, STATE)
%
% Returns state probability for the specified node and state
% For QBD models, state is a 2-element vector [level, phase]
%
% Parameters:
%   node - node/station index
%   state - state vector [level, phase] or structure with .level and .phase fields
%           If state is omitted or empty, returns full probability matrix
%
% Returns:
%   Pstate - probability of the state, or full probability matrix if state not specified

if nargin < 2
    line_error(mfilename,'getProb requires a node parameter.');
end
if nargin < 3
    state = [];
end

sn = self.getStruct;

% Convert node to station index if needed
if node > sn.nnodes
    line_error(mfilename,'Node number exceeds the number of nodes in the model.');
end
ist = sn.nodeToStation(node);
if ist == 0
    line_error(mfilename,'Specified node is not a station.');
end

% Check if this is a network (more than one queue station)
queueStations = 0;
for i=1:sn.nstations
    if sn.nodetype(sn.stationToNode(i)) == NodeType.Queue
        queueStations = queueStations + 1;
    end
end

if queueStations > 1
    line_error(mfilename,'getProb is currently only supported for single queue models in SolverMAM.\nFor networks, this method is not yet implemented.');
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
    T = self.result.Avg.T;
    lambda_total = sum(T(ist,:));

    if lambda_total < GlobalConstants.FineTol
        % No traffic at this station
        if isempty(state)
            % Return full probability matrix - all probability at state (0,1)
            Pstate = zeros(maxLevel, 1);
            Pstate(1, 1) = 1.0;
        else
            % Parse state
            if isstruct(state)
                level = state.level;
                phase = state.phase;
            elseif length(state) >= 2
                level = state(1);
                phase = state(2);
            else
                line_error(mfilename,'State must be a 2-element vector [level, phase] or structure with .level and .phase fields.');
            end

            if level == 0 && phase == 1
                Pstate = 1.0;
            else
                Pstate = 0.0;
            end
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
        % For getProb, we need the joint distribution over (level, phase)
        % The MMAPPH1FCFS ncDistr option gives us marginal over levels
        % To get the full joint distribution, we would need to access
        % the internal QBD solution
        %
        % For now, we approximate using the marginal and assuming
        % uniform distribution over phases (or using initial distribution)

        [pdistr] = MMAPPH1FCFS(D_approx, {pie{:}}, {D0{:}}, 'ncDistr', maxLevel);
        pdistr = abs(pdistr);
        pdistr = pdistr / sum(pdistr);

        % Build phase distribution - for simplicity, use the steady-state
        % phase distribution from the service process
        % This is an approximation
        nPhases = size(D0{1}, 1);
        for k=2:K
            nPhases = max(nPhases, size(D0{k}, 1));
        end

        % Construct joint probability matrix: rows = levels, cols = phases
        % For now, approximate by assuming phases are independent of level
        % and use the service process steady-state distribution
        avgPie = zeros(1, nPhases);
        for k=1:K
            if size(pie{k}, 2) == nPhases
                avgPie = avgPie + pie{k} * T(ist,k) / lambda_total;
            end
        end

        if sum(avgPie) == 0 || any(isnan(avgPie))
            avgPie = ones(1, nPhases) / nPhases;
        else
            avgPie = avgPie / sum(avgPie);
        end

        % Joint distribution: P(level, phase) ≈ P(level) * P(phase)
        Pstate = zeros(min(maxLevel, length(pdistr)), nPhases);
        for level=1:min(maxLevel, length(pdistr))
            Pstate(level, :) = pdistr(level) * avgPie;
        end

        % If specific state requested, extract it
        if ~isempty(state)
            if isstruct(state)
                level = state.level;
                phase = state.phase;
            elseif length(state) >= 2
                level = state(1);
                phase = state(2);
            else
                line_error(mfilename,'State must be a 2-element vector [level, phase] or structure with .level and .phase fields.');
            end

            % Check bounds (MATLAB indexing: level+1, phase)
            if level+1 > size(Pstate, 1) || phase > size(Pstate, 2) || level < 0 || phase < 1
                Pstate = 0.0;
            else
                Pstate = Pstate(level+1, phase);
            end
        end

    catch ME
        line_error(mfilename,'Failed to compute state probabilities: %s', ME.message);
    end

else
    % Open model - compute joint probability distribution using MAM (MMAPPH1FCFS)

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
    refstat = sn.refstat(1);  % Source station
    PH_src = sn.proc{refstat};

    % Build arrival process cell array: {D0, D1, D2, ...} for K classes
    D_arr = cell(1, K + 1);

    % Check total arrival rate
    totalLambda = 0;
    for k = 1:K
        if ~isnan(PH_src{k}{1})
            totalLambda = totalLambda + map_lambda(PH_src{k});
        end
    end

    % Compute queue length distribution using MMAPPH1FCFS
    maxLevel = 100;  % Maximum queue length to compute
    if ~isempty(self.options.cutoff) && isfinite(self.options.cutoff) && self.options.cutoff > 0
        maxLevel = self.options.cutoff;
    end

    if totalLambda < GlobalConstants.FineTol
        % No arrivals at this station
        if isempty(state)
            Pstate = zeros(maxLevel, 1);
            Pstate(1, 1) = 1.0;
        else
            if isstruct(state)
                level = state.level;
                phase = state.phase;
            elseif length(state) >= 2
                level = state(1);
                phase = state(2);
            else
                line_error(mfilename, 'State must be a 2-element vector [level, phase] or structure.');
            end
            if level == 0 && phase == 1
                Pstate = 1.0;
            else
                Pstate = 0.0;
            end
        end
    else
        % Build the aggregate arrival MMAP
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

        try
            [pdistr] = MMAPPH1FCFS(D_arr, pie, D0, 'ncDistr', maxLevel);
            pdistr = abs(pdistr);
            pdistr = pdistr / sum(pdistr);

            % Build phase distribution - use the steady-state phase distribution
            nPhases = size(D0{1}, 1);
            for k = 2:K
                nPhases = max(nPhases, size(D0{k}, 1));
            end

            % Compute weighted average phase distribution
            V = cellsum(sn.visits);
            avgPie = zeros(1, nPhases);
            for k = 1:K
                if size(pie{k}, 2) <= nPhases
                    pieK = [pie{k}, zeros(1, nPhases - size(pie{k}, 2))];
                    lambdaK = map_lambda(arrMaps{k});
                    avgPie = avgPie + pieK * lambdaK / totalLambda;
                end
            end

            if sum(avgPie) == 0 || any(isnan(avgPie))
                avgPie = ones(1, nPhases) / nPhases;
            else
                avgPie = avgPie / sum(avgPie);
            end

            % Joint distribution: P(level, phase) ≈ P(level) * P(phase)
            Pstate = zeros(min(maxLevel, length(pdistr)), nPhases);
            for level = 1:min(maxLevel, length(pdistr))
                Pstate(level, :) = pdistr(level) * avgPie;
            end

            % If specific state requested, extract it
            if ~isempty(state)
                if isstruct(state)
                    level = state.level;
                    phase = state.phase;
                elseif length(state) >= 2
                    level = state(1);
                    phase = state(2);
                else
                    line_error(mfilename, 'State must be a 2-element vector [level, phase] or structure.');
                end

                if level + 1 > size(Pstate, 1) || phase > size(Pstate, 2) || level < 0 || phase < 1
                    Pstate = 0.0;
                else
                    Pstate = Pstate(level + 1, phase);
                end
            end
        catch ME
            line_error(mfilename, sprintf('Failed to compute state probabilities: %s', ME.message));
        end
    end
end

end
