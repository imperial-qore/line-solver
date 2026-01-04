%{ @file sn_set_population.m
 %  @brief Sets the number of jobs for a closed class
 %
 %  @author LINE Development Team
%}

%{
 % @brief Sets the number of jobs for a closed class
 %
 % @details
 % Directly modifies the job population in NetworkStruct and recalculates
 % nclosedjobs.
 %
 % @par Syntax:
 % @code
 % sn = sn_set_population(sn, classIdx, nJobs)
 % sn = sn_set_population(sn, classIdx, nJobs, autoRefresh)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>sn<td>Network structure
 % <tr><td>classIdx<td>Class index (1-based)
 % <tr><td>nJobs<td>Number of jobs (non-negative)
 % <tr><td>autoRefresh<td>If true, refresh visit ratios (default false)
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>sn<td>Modified network structure
 % </table>
%}
function sn = sn_set_population(sn, classIdx, nJobs, autoRefresh)

if nargin < 4 || isempty(autoRefresh)
    autoRefresh = false;
end

% Update njobs
sn.njobs(classIdx) = nJobs;

% Recalculate nclosedjobs
sn.nclosedjobs = sum(sn.njobs(isfinite(sn.njobs)));

% Auto-refresh visit ratios if requested
if autoRefresh
    [~, ~, sn] = sn_refresh_visits(sn, sn.chains, sn.rt, sn.rtnodes);
end

end
