%{ @file sn_to_qrf_capacity.m
 %  @brief Per-station occupancy bound F(i) for the QRF bounds
 %
 %  @author LINE Development Team
%}

%{
 % @brief Returns the QRF capacity vector F and which buffers actually bind
 %
 % @details
 % The QRF bounds index every marginal by 0..F(i), so F is an occupancy bound
 % rather than a declared capacity: it is the station's buffer where that
 % buffer BINDS, and the population N everywhere else, since no queue of a
 % closed model can hold more than N jobs.
 %
 % Binding is decided by SN_GET_BUFFER_SIZE, the single place in LINE that
 % makes that call: refreshCapacity derives a finite classcap (the chain
 % population) at every station of every closed model, so a plain finiteness
 % test on sn.cap would report a buffer at every station.
 %
 % Both QRF blocking bounds need this. 'qrf.bas' needs it beside the blocking
 % tables SN_TO_QRF_BLOCKING derives; 'qrf.rsrd' needs it ALONE, since its PBB
 % constraint reads only which queues can be full and it carries no blocking
 % tables at all.
 %
 % @par Syntax:
 % @code
 % [F, binding, msg] = sn_to_qrf_capacity(sn)
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
 % <tr><td>F<td>(nstations x 1) occupancy bound of each station, in jobs
 % <tr><td>binding<td>(nstations x 1) logical, true where the buffer can refuse a job
 % <tr><td>msg<td>Empty on success, otherwise why F is not defined for this model
 % </table>
%}
function [F, binding, msg] = sn_to_qrf_capacity(sn)
msg = '';
M = sn.nstations;
F = zeros(M, 1);
binding = false(M, 1);

N = sum(sn.njobs);
if ~isfinite(N) || N < 1
    msg = 'the QRF bounds need a closed model with a finite population.';
    return
end

for i = 1:M
    b = sn_get_buffer_size(sn, i);
    binding(i) = isfinite(b);
    if ~isfinite(b) || b > N
        F(i) = N;
    else
        F(i) = b;
    end
    if F(i) < 1
        msg = sprintf(['station %d has capacity %d: the QRF bounds need every queue to be able ' ...
            'to hold at least one job.'], i, F(i));
        return
    end
end
end
