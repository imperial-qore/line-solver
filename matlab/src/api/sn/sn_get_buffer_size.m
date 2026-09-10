%{ @file sn_get_buffer_size.m
 %  @brief Physical buffer size of a station, in jobs
 %
 %  @author LINE Development Team
%}

%{
 % @brief Returns the number of jobs a station can hold, in service included
 %
 % @details
 % Kendall's K: the total occupancy bound of station IST, obtained as the
 % tighter of the station capacity sn.cap(ist) and the per-class capacities
 % sn.classcap(ist,:). Returns Inf when the station is unbounded. Both
 % fields are populated by refreshCapacity, which already folds setCapacity,
 % setClassCapacity, a finite orbit and the closed-chain population into
 % them, so this is the single place that decides whether a buffer BINDS.
 %
 % Only a buffer that can actually BIND is reported. refreshCapacity derives
 % a FINITE classcap (the chain population) for EVERY closed model, so a
 % plain finiteness test would report a buffer at every station of every
 % closed model; a capacity at least as large as the total population can
 % never refuse a job and is returned as Inf. sum(njobs) is Inf as soon as
 % one class is open, so any finite capacity reachable by an open class
 % binds.
 %
 % @par Syntax:
 % @code
 % N = sn_get_buffer_size(sn, ist)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>sn<td>Network structure
 % <tr><td>ist<td>Station index
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>N<td>Buffer size in jobs, Inf if unbounded
 % </table>
%}
function N = sn_get_buffer_size(sn, ist)
N = Inf;
if isfield(sn, 'cap') && ~isempty(sn.cap) && numel(sn.cap) >= ist
    if sn.cap(ist) >= 0
        N = min(N, sn.cap(ist));
    end
end
if isfield(sn, 'classcap') && ~isempty(sn.classcap) && size(sn.classcap, 1) >= ist
    ccap = sn.classcap(ist, :);
    ccap = ccap(ccap > 0); % a zero marks a class that is not served here
    if ~isempty(ccap)
        N = min(N, sum(ccap));
    end
end
if isfield(sn, 'njobs') && ~isempty(sn.njobs)
    % The population that can actually REACH this station, not the model's
    % total. A capacity is unreachable when it is at least as large as the
    % jobs able to arrive here, and only the classes served here can arrive.
    %
    % Using sum(njobs) over every class made a MIXED model read as finite
    % buffered at stations an open class never visits: the open class puts Inf
    % into the total, so the test below could not fire, while classcap at such
    % a station is the CLOSED chain's population, which refreshCapacity derives
    % and which can never refuse a job of that chain. On Source/Delay/Queue with
    % one closed class of 2 and one open class, the Delay reported a buffer of 2
    % and solver_nc_mem_supports then rejected the whole model with "MEM
    % supports finite station buffers only in open models", making the mixed MEM
    % branch unreachable.
    reach = true(1, numel(sn.njobs));
    if isfield(sn, 'classcap') && ~isempty(sn.classcap) && size(sn.classcap, 1) >= ist
        reach = sn.classcap(ist, :) > 0;
    end
    if any(reach)
        reachableJobs = sum(sn.njobs(reach));
    else
        reachableJobs = 0;
    end
    if N >= reachableJobs
        N = Inf; % declared but unreachable: the buffer can never refuse a job
    end
end
end
