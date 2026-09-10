%{ @file sn_get_residt_from_respt.m
 %  @brief Computes residence times from response times
 %
 %  @author LINE Development Team
%}

%{
 % @brief Computes residence times from response times
 %
 % @details
 % This function converts response times to residence times by accounting
 % for visit ratios at each station.
 %
 % @par Syntax:
 % @code
 % WN = sn_get_residt_from_respt(sn, RN, WH)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>sn<td>Network structure
 % <tr><td>RN<td>Average response times
 % <tr><td>WH<td>Residence time handles
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>WN<td>Average residence times
 % </table>
%}
function WN=sn_get_residt_from_respt(sn, RN, WH)

M = sn.nstations;
K = sn.nclasses;
WN = zeros(M, K);

% sn.visits is cell array of stateful-indexed matrices (nstateful x nclasses)
% We need to convert to station-indexed matrix (M x K) using statefulToStation
Vstateful = cellsum(sn.visits);
V = zeros(M, K);
for sf = 1:sn.nstateful
    stationIdx = sn.statefulToStation(sf);
    if ~isnan(stationIdx) && stationIdx >= 1 && stationIdx <= M
        V(stationIdx, :) = V(stationIdx, :) + Vstateful(sf, :);
    end
end

for ist = 1:M
    for k = 1:K
        if isempty(WH) || WH{ist,k}.disabled
            WN(ist,k) = NaN;
        elseif ~isempty(RN) && RN(ist,k) > 0
            if RN(ist,k) < GlobalConstants.FineTol
                WN(ist,k) = RN(ist,k);
            else
                c = find(sn.chains(:, k));
                refclass = sn.refclass(c);
                if refclass > 0
                    WN(ist,k) = RN(ist,k) * V(ist,k) / sum(V(sn.refstat(k), refclass));
                else
                    WN(ist,k) = RN(ist,k) * V(ist,k) / sum(V(sn.refstat(k), sn.inchain{c}));
                end
            end
        end
    end
end
WN(isnan(WN)) = 0;
WN(WN < 10 * GlobalConstants.FineTol) = 0;
WN(WN < GlobalConstants.FineTol) = 0;
end