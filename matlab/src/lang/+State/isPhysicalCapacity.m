function tf = isPhysicalCapacity(sn, ist, class)
% TF = ISPHYSICALCAPACITY(SN, IST, CLASS)
%
% True when the capacity bound at station IST for CLASS is a PHYSICAL finite
% capacity (setCapacity/setClassCapacity), as opposed to a state-space CUTOFF
% imposed on an open class only to bound enumeration.
%
% This distinction matters because the producer's capacity/classcap arguments
% have the open-class cutoff folded in: solver_ssa overwrites sn.cap/sn.classcap
% with min(cutoff, physical), so at the cutoff boundary they are finite even
% when there is no physical cap. Treating a cutoff boundary as a physical one
% would turn a state-space truncation into a self-loop loss (wrong ArvR and a
% perturbed sample path). Only a physical cap should trigger the loss/block
% refusal semantics; a cutoff-only refusal must fall back to the pre-change
% truncation (place the job, let the en_o capacity filter delete the row).
%
% The reliable in-producer signal is the DROP RULE. refreshCapacity sets a
% finite-capacity drop rule (DROP, or a blocking/retrial rule) exactly when the
% station has a physical finite capacity for the class; an open class bounded
% only by the cutoff keeps the WAITQ default. So a non-WAITQ, non-unset drop
% rule marks a physical capacity. (A user who explicitly sets WAITQ on a
% physical cap is the one ambiguous case; it is already ill-defined -- CTMC
% drops, JMT blocks -- and is treated here as a cutoff, i.e. truncated.)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if isempty(sn.droprule) || size(sn.droprule,1) < ist || size(sn.droprule,2) < class
    tf = false;
    return
end
dr = sn.droprule(ist,class);
tf = dr ~= DropStrategy.WAITQ && dr ~= 0;
end
