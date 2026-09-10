%{ @file sn_is_mm1k_loss.m
 %  @brief Checks if the network is a single-station M/M/1/K queue with tail drop
 %
 %  @author LINE Development Team
%}

function bool = sn_is_mm1k_loss(sn)
% BOOL = SN_IS_MM1K_LOSS(SN)
% Returns true for a single-class open Source-Queue-Sink system whose queue is
% a single-server exponential M/M/1/K with tail drop (DropStrategy.DROP). This
% is the exact regime shared by the closed-form loss scripts qsys_mm1k_loss
% (probability-based, used by SolverNC) and qsys_mg1k_loss_mgs (moment-based,
% used by SolverMVA), so both solvers gate their finite-capacity loss branch on
% this predicate.
bool = false;
if sn.nclasses ~= 1 || sn.nclosedjobs ~= 0 || numel(sn.nodetype) ~= 3
    return
end
if ~all(sort(sn.nodetype(:))' == sort([NodeType.Source,NodeType.Queue,NodeType.Sink]))
    return
end
qist = sn.nodeToStation(sn.nodetype == NodeType.Queue);
sist = sn.nodeToStation(sn.nodetype == NodeType.Source);
if sn.nservers(qist) ~= 1
    return
end
if isempty(sn.droprule) || sn.droprule(qist,1) ~= DropStrategy.DROP
    return
end
if ~isfinite(sn.cap(qist)) || sn.cap(qist) <= 0
    return
end
if abs(sn.scv(sist,1)-1) > 1e-6 || abs(sn.scv(qist,1)-1) > 1e-6
    return
end
bool = true;
end
