%{ @file ba_blocking_default.m
 %  @brief What 'default'/'auto' must mean on a model with finite buffers
 %
 %  @author LINE Development Team
%}

%{
 % @brief Selects the QRF BAS bound as the default of a blocked model
 %
 % @details
 % 'default' resolves to the geometric upper bound, which is parameterized by
 % demands and a population alone and therefore bounds a blocked model as if
 % its buffers were unbounded. BA_IGNORES_BLOCKING refuses that, which is
 % right; but refusing is not the whole answer, because 'qrf.bas' DOES model
 % the finite buffer and the model itself says everything that bound needs.
 % So a blocked model that fits the QRF shape gets 'qrf.bas' as its default,
 % exactly as SolverMVA routes a BAS model to 'sqd' through SN_IS_BAS_MODEL.
 %
 % The shape is the one `listValidMethods` already calls "reducible" and
 % `solver_ba_qrf_analyzer` gates on: single class, closed, no delay station,
 % no multiserver station. On top of that the blocking tables must actually
 % derive, which is asked of SN_TO_QRF_BLOCKING rather than re-tested here --
 % it owns the single-finite-buffer rule and the size guard, and a second copy
 % of either is how the two drift apart.
 %
 % WHY is the derivation's own reason when the routing does not apply, so the
 % caller is told what about THIS model rules the blocking bound out (two
 % binding buffers, an oversized enumeration) instead of a generic refusal.
 % It is empty when the model is simply not blocked.
 %
 % Only the UPPER side is routed. `qrf.bas` is solved in the 'max' direction
 % alone by the analyzer, so 'auto.lower' has no blocking counterpart and
 % keeps refusing rather than being answered with the wrong side.
 %
 % @par Syntax:
 % @code
 % [method, why] = ba_blocking_default(sn)
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
 % <tr><td>method<td>'qrf.bas' when the blocked model admits it, otherwise empty
 % <tr><td>why<td>Reason the blocking bound does not apply, empty when it does or when the model is unblocked
 % </table>
%}
function [method, why] = ba_blocking_default(sn)
method = '';
why = '';

if ~sn_has_blocking(sn)
    return
end

% Same premise listValidMethods calls "reducible": the QRF reduction bounds
% are derived for a single-class closed network of single servers.
if any(isinf(sn.njobs)) || sn.nclasses ~= 1
    why = ['the QRF blocking bounds are derived for a single-class closed network, ' ...
        'which this model is not.'];
    return
end
if any(sn.sched == SchedStrategy.INF)
    why = ['the QRF blocking bounds model every station as a single server and have no ' ...
        'infinite-server notion, so a delay station rules them out.'];
    return
end
if any(sn.nservers(sn.sched ~= SchedStrategy.INF) > 1)
    why = ['the QRF blocking bounds model every station as a single server, so a ' ...
        'multiserver station rules them out.'];
    return
end

[~, blkMsg] = sn_to_qrf_blocking(sn);
if ~isempty(blkMsg)
    why = blkMsg;
    return
end

method = 'qrf.bas';
end
