%{ @file sn_set_service_batch.m
 %  @brief Sets service rates for multiple station-class pairs
 %
 %  @author LINE Development Team
%}

%{
 % @brief Sets service rates for multiple station-class pairs
 %
 % @details
 % Batch update of service rates. NaN values are skipped (not updated).
 % More efficient than calling sn_set_service multiple times.
 %
 % @par Syntax:
 % @code
 % sn = sn_set_service_batch(sn, rates)
 % sn = sn_set_service_batch(sn, rates, scvs)
 % sn = sn_set_service_batch(sn, rates, scvs, autoRefresh)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>sn<td>Network structure
 % <tr><td>rates<td>Matrix of new rates (nstations x nclasses), NaN = skip
 % <tr><td>scvs<td>Matrix of new SCVs (optional)
 % <tr><td>autoRefresh<td>If true, refresh process fields (default false)
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>sn<td>Modified network structure
 % </table>
%}
function sn = sn_set_service_batch(sn, rates, scvs, autoRefresh)

if nargin < 3
    scvs = [];
end

if nargin < 4 || isempty(autoRefresh)
    autoRefresh = false;
end

M = sn.nstations;
K = sn.nclasses;

% Track updated pairs for auto-refresh
updatedPairs = [];

% Update rates
for i = 1:M
    for j = 1:K
        if ~isnan(rates(i, j))
            sn.rates(i, j) = rates(i, j);
            updatedPairs = [updatedPairs; i, j]; %#ok<AGROW>
        end
    end
end

% Update SCVs if provided
if ~isempty(scvs)
    for i = 1:M
        for j = 1:K
            if ~isnan(scvs(i, j))
                sn.scv(i, j) = scvs(i, j);
            end
        end
    end
end

% Auto-refresh if requested
if autoRefresh
    for idx = 1:size(updatedPairs, 1)
        sn = sn_refresh_process_fields(sn, updatedPairs(idx, 1), updatedPairs(idx, 2));
    end
end

end
