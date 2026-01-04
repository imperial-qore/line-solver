%{ @file fj_dist2fj.m
 %  @brief Convert LINE PH/MAP distribution to FJ_codes format
 %
 %  @author LINE Development Team
%}

%{
 % @brief Convert LINE PH/MAP distribution to FJ_codes format
 %
 % @details
 % Converts LINE's MAP representation {D0, D1} to FJ_codes format.
 % For arrivals: returns (lambda, lambda0, lambda1, ma, Ia)
 % For service: returns (mu, ST, St, tau_st, SerChoice)
 %
 % Supported distributions:
 % - Exponential (Exp)
 % - 2-phase Hyperexponential (HyperExp with 2 phases)
 % - 2-phase Erlang (Erlang with 2 phases)
 % - 2-phase MAP (MAP with 2 states)
 %
 % @par Syntax:
 % @code
 % fjDist = fj_dist2fj(lineMAP, distType, sn, ist, r)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>lineMAP<td>LINE MAP representation {D0, D1}
 % <tr><td>distType<td>'arrival' or 'service'
 % <tr><td>sn<td>Network structure (for getting rate and process type)
 % <tr><td>ist<td>Station index
 % <tr><td>r<td>Class index
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>fjDist<td>FJ_codes compatible structure
 % </table>
 %
 % @par Reference:
 % Z. Qiu, J.F. Pérez, and P. Harrison, "Beyond the Mean in Fork-Join Queues:
 % Efficient Approximation for Response-Time Tails", IFIP Performance 2015.
 % Copyright 2015 Imperial College London
%}
function fjDist = fj_dist2fj(lineMAP, distType, sn, ist, r)

fjDist = struct();

% Get process type and mean rate
procType = sn.procid(ist, r);
mean_rate = sn.rates(ist, r);

% Extract D0 and D1 from LINE MAP format
D0 = lineMAP{1};
D1 = lineMAP{2};
nPhases = size(D0, 1);

% Check number of phases
if nPhases > 2
    line_error(mfilename, 'FJ_codes only supports distributions with at most 2 phases. Found %d phases.', nPhases);
end

% Compute mean arrival/service rate (lambda or mu)
lambda = map_lambda(lineMAP);
mean_time = 1 / lambda;

if strcmp(distType, 'arrival')
    % ===== ARRIVAL PROCESS =====
    fjDist.lambda = lambda;
    fjDist.lambda0 = D0;
    fjDist.lambda1 = D1;
    fjDist.ma = nPhases;
    fjDist.Ia = eye(nPhases);

    % Determine arrival choice for logging/debugging
    if procType == ProcessType.EXP
        fjDist.ArrChoice = 1; % Exp
    elseif procType == ProcessType.HYPEREXP && nPhases == 2
        fjDist.ArrChoice = 2; % HE2
    elseif procType == ProcessType.ERLANG && nPhases == 2
        fjDist.ArrChoice = 3; % ER2
    elseif procType == ProcessType.MAP && nPhases == 2
        fjDist.ArrChoice = 4; % MAP2
    else
        line_error(mfilename, 'Unsupported arrival distribution type: %s with %d phases.', ...
            ProcessType.toText(procType), nPhases);
    end

else
    % ===== SERVICE PROCESS =====
    % Convert MAP to PH representation for service
    % PH is represented as (tau_st, ST) where:
    % - tau_st is initial probability vector
    % - ST is the rate matrix
    % - St is the exit vector (computed as -ST * ones)

    % For service, we need to compute the PH representation from MAP
    % The initial distribution for service is the stationary distribution
    pie = map_pie(lineMAP);

    fjDist.mu = lambda; % Mean service rate
    fjDist.ST = D0;  % Rate matrix
    fjDist.St = -D0 * ones(nPhases, 1);  % Exit vector
    fjDist.tau_st = pie(:)';  % Initial probability (as row vector)

    % Determine service choice for logging/debugging
    if procType == ProcessType.EXP
        fjDist.SerChoice = 1; % Exp
    elseif procType == ProcessType.HYPEREXP && nPhases == 2
        fjDist.SerChoice = 2; % HE2
    elseif procType == ProcessType.ERLANG && nPhases == 2
        fjDist.SerChoice = 3; % ER2
    else
        line_error(mfilename, 'Unsupported service distribution type: %s with %d phases.', ...
            ProcessType.toText(procType), nPhases);
    end
end

end
