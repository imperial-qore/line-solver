%{ @file sn_refresh_process_fields.m
 %  @brief Refreshes process fields based on rate and SCV values
 %
 %  @author LINE Development Team
%}

%{
 % @brief Refreshes process fields for a station-class pair
 %
 % @details
 % Updates mu, phi, proc, pie, phases based on current rate and SCV values.
 % - SCV = 1.0: Exponential (1 phase)
 % - SCV < 1.0: Erlang approximation
 % - SCV > 1.0: Hyperexponential(2) approximation
 %
 % @par Syntax:
 % @code
 % sn = sn_refresh_process_fields(sn, stationIdx, classIdx)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>sn<td>Network structure
 % <tr><td>stationIdx<td>Station index (1-based)
 % <tr><td>classIdx<td>Class index (1-based)
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>sn<td>Modified network structure
 % </table>
%}
function sn = sn_refresh_process_fields(sn, stationIdx, classIdx)

rate = sn.rates(stationIdx, classIdx);
scv = sn.scv(stationIdx, classIdx);

% Skip if rate is invalid
if isnan(rate) || rate <= 0 || isinf(rate)
    return;
end

mean = 1.0 / rate;

% Determine process type and create MAP based on SCV
if isnan(scv) || abs(scv - 1.0) < 1e-10
    % Exponential
    MAP = map_exponential(mean);
    nPhases = 1;
    procType = ProcessType.EXP;
elseif scv < 1.0
    % Erlang: k = ceil(1/scv)
    k = max(1, ceil(1.0 / scv));
    MAP = map_erlang(mean, k);
    nPhases = k;
    procType = ProcessType.ERLANG;
else
    % Hyperexponential (scv > 1)
    MAP = map_hyperexp(mean, scv);
    if ~isempty(MAP)
        nPhases = 2;
        procType = ProcessType.HYPEREXP;
    else
        % Fallback to exponential
        MAP = map_exponential(mean);
        nPhases = 1;
        procType = ProcessType.EXP;
    end
end

% Update process fields
D0 = MAP{1};
D1 = MAP{2};

% Update proc
sn.proc{stationIdx, classIdx} = MAP;

% Update procid
sn.procid(stationIdx, classIdx) = procType;

% Update phases
sn.phases(stationIdx, classIdx) = nPhases;

% Update phasessz
sn.phasessz(stationIdx, classIdx) = max(nPhases, 1);

% Recompute phaseshift for this station
cumSum = 0;
sn.phaseshift(stationIdx, 1) = 0;
for c = 1:sn.nclasses
    cumSum = cumSum + sn.phasessz(stationIdx, c);
    if c + 1 <= size(sn.phaseshift, 2)
        sn.phaseshift(stationIdx, c + 1) = cumSum;
    end
end

% Update mu (rates from -diag(D0))
muVec = -diag(D0);
sn.mu{stationIdx, classIdx} = muVec;

% Update phi (completion probabilities)
phiVec = zeros(nPhases, 1);
for i = 1:nPhases
    d1RowSum = sum(D1(i, :));
    d0Diag = -D0(i, i);
    if d0Diag ~= 0
        phiVec(i) = d1RowSum / d0Diag;
    end
end
sn.phi{stationIdx, classIdx} = phiVec;

% Update pie (initial phase distribution)
sn.pie{stationIdx, classIdx} = map_pie(MAP);

end
