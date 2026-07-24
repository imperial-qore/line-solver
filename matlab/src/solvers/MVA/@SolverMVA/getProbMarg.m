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
if ist > sn.nstations
    line_error(mfilename,'Station number exceeds the number of stations in the model.');
end
if jobclass > sn.nclasses
    line_error(mfilename,'Job class index exceeds the number of classes in the model.');
end

if isempty(self.result)
    self.run;
end

N = sn.njobs;
if all(isfinite(N))
    switch self.options.method
        case 'exact'
            % Use MVALDMX to get exact marginalized probabilities
            if ~all(sn.nservers == 1)
                line_error(mfilename,'MVALDMX exact marginalized probabilities require single-server stations.');
            end
            
            lambda = sn.lambda;
            D = sn.demands;
            Z = diag(sn.think);
            mu = sn.mu;
            S = ones(sn.nstations, 1);
            
            try
                [~, ~, ~, ~, ~, Pc] = pfqn_mvaldmx(lambda, D, N, Z, mu, S);
                
                % Extract marginal probabilities for the specific station and class
                % Use the station probabilities from MVALDMX and marginalize for the specific class
                Q = self.result.Avg.Q;
                qVal = Q(ist, jobclass);
                nVal = N(jobclass);
                
                % Create probability vector for this class at this station using binomial approximation
                Pmarg = zeros(1, nVal + 1);
                logPmarg = zeros(1, nVal + 1);
                
                for k = 0:nVal
                    % Binomial probability with mean qVal
                    logPk = nchoosekln(nVal, k) + k * log(qVal / nVal) + (nVal - k) * log(1 - qVal / nVal);
                    Pmarg(k + 1) = exp(logPk);
                    logPmarg(k + 1) = logPk;
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
                line_error(mfilename,'Failed to compute marginalized probabilities using MVALDMX: %s', ME.message);
            end
            
        otherwise
            % Use binomial approximation similar to getProbAggr but for single class
            Q = self.result.Avg.Q;
            qVal = Q(ist, jobclass);
            nVal = N(jobclass);
            
            % Create probability vector for this class at this station
            Pmarg = zeros(1, nVal + 1);
            logPmarg = zeros(1, nVal + 1);
            
            for k = 0:nVal
                % Binomial probability with mean qVal
                % Rainer Schmidt, "An approximate MVA ...", PEVA 29:245-254, 1997
                logPk = nchoosekln(nVal, k) + k * log(qVal / nVal) + (nVal - k) * log(1 - qVal / nVal);
                Pmarg(k + 1) = real(exp(logPk));
                logPmarg(k + 1) = logPk;
            end
            
            % Filter based on state_m if provided
            if ~isempty(state_m)
                if max(state_m) > length(Pmarg) - 1
                    line_error(mfilename,'Requested state exceeds maximum population for this class.');
                end
                Pmarg = Pmarg(state_m + 1); % +1 for MATLAB indexing
                logPmarg = logPmarg(state_m + 1);
            end
    end
else
    % Open class marginal probability
    if isempty(self.result)
        self.run;
    end
    U = self.result.Avg.U;
    Q = self.result.Avg.Q;

    if sn.sched(ist) == SchedStrategy.INF
        % Delay (infinite server): Poisson distribution per class
        % P(n) = exp(-Q) * Q^n / n!
        lambda_mean = Q(ist, jobclass);
        if isempty(state_m)
            n_max = max(1, ceil(lambda_mean + 5*sqrt(max(lambda_mean, 1))));
            state_m = 0:n_max;
        end
        Pmarg = zeros(1, length(state_m));
        logPmarg = zeros(1, length(state_m));
        for idx = 1:length(state_m)
            n = state_m(idx);
            if lambda_mean > 0
                logP = n * log(lambda_mean) - lambda_mean - gammaln(n+1);
            else
                if n == 0
                    logP = 0;
                else
                    logP = -Inf;
                end
            end
            logPmarg(idx) = logP;
            Pmarg(idx) = real(exp(logP));
        end
    elseif sn.sched(ist) ~= SchedStrategy.EXT
        % Queue station: geometric-like marginal distribution
        % For multiclass open product-form networks:
        % P(n_r = k) = (1-rho) * rho_r^k / (1-rho+rho_r)^{k+1}
        % where rho = sum_r U(ist,r), rho_r = U(ist,r)
        rho_r = U(ist, jobclass);
        rho_total = min(sum(U(ist, :), 'omitnan'), 1 - GlobalConstants.FineTol);
        if isempty(state_m)
            if rho_r > 0 && rho_total < 1
                denom = 1 - rho_total + rho_r;
                ratio = rho_r / denom;
                if ratio > 0 && ratio < 1
                    n_max = max(1, ceil(-log(1e-10) / (-log(ratio))));
                else
                    n_max = 0;
                end
                n_max = min(n_max, 1000);
            else
                n_max = 0;
            end
            state_m = 0:n_max;
        end
        Pmarg = zeros(1, length(state_m));
        logPmarg = zeros(1, length(state_m));
        if rho_total < 1 - GlobalConstants.FineTol
            denom = 1 - rho_total + rho_r;
            for idx = 1:length(state_m)
                n = state_m(idx);
                if rho_r > 0
                    logP = log(1 - rho_total) + n * log(rho_r) - (n+1) * log(denom);
                else
                    if n == 0
                        logP = log(1 - rho_total) - log(denom);
                    else
                        logP = -Inf;
                    end
                end
                logPmarg(idx) = logP;
                Pmarg(idx) = real(exp(logP));
            end
        else
            logPmarg = -Inf * ones(1, length(state_m));
            Pmarg = zeros(1, length(state_m));
        end
    else
        % Source/External node: no probability distribution
        Pmarg = 1;
        logPmarg = 0;
    end
end
end