%{ @file sn_has_blocking.m
 %  @brief Checks if the network holds jobs back at a finite buffer or region
 %
 %  @author LINE Development Team
%}

%{
 % @brief Checks if the network has finite-buffer blocking or loss
 %
 % @details
 % Returns true when some station can refuse a job, either because its own
 % buffer BINDS (Kendall's K below the population that can reach it, whatever
 % the drop rule: WAITQ, DROP, BAS, BBS, RSRD) or because a finite capacity
 % region caps a set of stations jointly. Such a network is not product form:
 % the truncation couples the station occupancies, so no BCMP factorization
 % of the equilibrium distribution exists.
 %
 % Only a buffer that can actually BIND counts, which is what
 % SN_GET_BUFFER_SIZE decides: refreshCapacity derives a finite classcap (the
 % chain population) at every station of every closed model, so a plain
 % finiteness test would call every closed model blocking.
 %
 % Two shapes are exempt. A Cache builds its own capped retrieval queues
 % (classCap = 1), which the cache analyzers solve rather than treat as a
 % buffer constraint, the same exemption NetworkSolver.checkBindingCapacity
 % makes. And the single-station M/M/1/K loss system keeps the truncated
 % geometric distribution, a product form over its one station, which
 % qsys_mm1k_loss and qsys_mg1k_loss_mgs evaluate in closed form.
 %
 % @par Syntax:
 % @code
 % bool = sn_has_blocking(sn)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>sn<td>Network structure
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>bool<td>True if the network has binding finite buffers or regions
 % </table>
%}
function bool = sn_has_blocking(sn)
bool = false;

% a finite capacity region caps a SET of stations, which no per-station
% capacity can express and no product form survives
if isfield(sn,'nregions') && ~isempty(sn.nregions) && sn.nregions > 0
    bool = true;
    return
end

if isfield(sn,'nodetype') && any(sn.nodetype == NodeType.Cache)
    return
end

if sn_is_mm1k_loss(sn)
    return
end

for ist = 1:sn.nstations
    if isfinite(sn_get_buffer_size(sn, ist))
        bool = true;
        return
    end
end
end
