%{ @file sn_set_arrival.m
 %  @brief Sets arrival rate for a class at the Source station
 %
 %  @author LINE Development Team
%}

%{
 % @brief Sets arrival rate for a class at the Source station
 %
 % @details
 % Directly modifies the arrival rate in NetworkStruct. Finds the Source
 % station and updates its rate for the specified class.
 %
 % @par Syntax:
 % @code
 % sn = sn_set_arrival(sn, classIdx, rate)
 % sn = sn_set_arrival(sn, classIdx, rate, scv)
 % sn = sn_set_arrival(sn, classIdx, rate, scv, autoRefresh)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>sn<td>Network structure
 % <tr><td>classIdx<td>Class index (1-based)
 % <tr><td>rate<td>New arrival rate (lambda)
 % <tr><td>scv<td>Squared coefficient of variation (default 1.0)
 % <tr><td>autoRefresh<td>If true, refresh process fields (default false)
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>sn<td>Modified network structure
 % </table>
%}
function sn = sn_set_arrival(sn, classIdx, rate, scv, autoRefresh)

if nargin < 4 || isempty(scv)
    scv = 1.0;
end

if nargin < 5 || isempty(autoRefresh)
    autoRefresh = false;
end

% Find Source station index
sourceIdx = find(sn.nodetype == NodeType.Source);
if isempty(sourceIdx)
    error('sn_set_arrival: No Source station found in network');
end

% Convert node index to station index
stationIdx = sn.nodeToStation(sourceIdx);

% Delegate to sn_set_service
sn = sn_set_service(sn, stationIdx, classIdx, rate, scv, autoRefresh);

end
