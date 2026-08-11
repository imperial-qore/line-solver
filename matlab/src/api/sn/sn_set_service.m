%{ @file sn_set_service.m
 %  @brief Sets service rate at a specific station and class
 %
 %  @author LINE Development Team
%}

%{
 % @brief Sets service rate at a specific station and class
 %
 % @details
 % Directly modifies the service rate in NetworkStruct without rebuilding
 % the full Network object. Useful for fast parameter updates in optimization.
 %
 % @par Syntax:
 % @code
 % sn = sn_set_service(sn, stationIdx, classIdx, rate)
 % sn = sn_set_service(sn, stationIdx, classIdx, rate, scv)
 % sn = sn_set_service(sn, stationIdx, classIdx, rate, scv, autoRefresh)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>sn<td>Network structure
 % <tr><td>stationIdx<td>Station index (1-based)
 % <tr><td>classIdx<td>Class index (1-based)
 % <tr><td>rate<td>New service rate (must be positive)
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
function sn = sn_set_service(sn, stationIdx, classIdx, rate, scv, autoRefresh)

if nargin < 5 || isempty(scv)
    scv = 1.0;
end

if nargin < 6 || isempty(autoRefresh)
    autoRefresh = false;
end

% Update rates matrix
sn.rates(stationIdx, classIdx) = rate;

% Update scv matrix
sn.scv(stationIdx, classIdx) = scv;

% Auto-refresh process fields if requested
if autoRefresh
    sn = sn_refresh_process_fields(sn, stationIdx, classIdx);
end

end
