function [isReneging, info] = mam_reneging_applicable(sn)
% [ISRENEGING, INFO] = MAM_RENEGING_APPLICABLE(SN)
%
% Is this model the MAP/M/s+G shape the reneging analyzer of solver_mam_retrial
% solves (Gursoy, Mehr, Akar, "The MAP/M/s + G Call Center Model with
% Generally Distributed Patience Times")? INFO carries the indices the
% analyzer reads (sourceIdx, queueIdx, classIdx, nServers, serviceRate) and,
% on a refusal, errorMsg naming the requirement that failed.
%
% Requirements:
% - Open model, single class
% - Single Queue station with reneging patience configured
% - MAP/BMAP arrival at the Source
% - Exponential service at the Queue (single-phase)
% - FCFS scheduling
%
% ONE PREDICATE, THREE CALLERS. solver_mam_retrial asks it to decide which of
% its two analyzers runs; MAM_RETRIAL_APPLICABLE asks it so that 'retrial' is
% offered only on a shape the analyzer answers; SolverMAM.supportsModelMethod
% asks it for 'default' on a reneging model, which the default routes here.
% Until 2026-09-05 it was a local function of solver_mam_retrial and the gate
% tested PRESENCE (MAM_HAS_RENEGING_PATIENCE) alone, so a two-queue reneging
% model was offered and then refused by the analyzer.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

info = struct();
info.sourceIdx = [];
info.queueIdx = [];
info.classIdx = [];
info.nServers = [];
info.serviceRate = [];
info.errorMsg = '';
isReneging = false;

% Check open model
if ~sn_is_open_model(sn)
    info.errorMsg = 'MAPMsG requires open queueing model.';
    return;
end

% Check single class (current limitation)
if sn.nclasses > 1
    info.errorMsg = 'MAPMsG currently supports single class only.';
    return;
end
info.classIdx = 1;

% Find source and queue stations
sourceIdx = [];
queueIdx = [];
for ist = 1:sn.nstations
    nodeIdx = sn.stationToNode(ist);
    if sn.nodetype(nodeIdx) == NodeType.Source
        sourceIdx = ist;
    elseif sn.nodetype(nodeIdx) == NodeType.Queue
        if isempty(queueIdx)
            queueIdx = ist;
        else
            % Multiple queues - not supported
            info.errorMsg = 'MAPMsG requires single queue station.';
            return;
        end
    end
end

if isempty(sourceIdx)
    info.errorMsg = 'No Source node found.';
    return;
end
if isempty(queueIdx)
    info.errorMsg = 'No Queue node found.';
    return;
end

info.sourceIdx = sourceIdx;
info.queueIdx = queueIdx;

% Check for reneging patience configuration
if ~isfield(sn, 'impatienceClass') || isempty(sn.impatienceClass)
    info.errorMsg = 'No patience/impatience configuration found.';
    return;
end

if sn.impatienceClass(queueIdx, info.classIdx) ~= ImpatienceType.RENEGING
    info.errorMsg = 'Queue does not have reneging configured.';
    return;
end

% Check patience distribution exists
if ~isfield(sn, 'patienceProc') || isempty(sn.patienceProc)
    info.errorMsg = 'No patience distribution found.';
    return;
end
if isempty(sn.patienceProc{queueIdx, info.classIdx})
    info.errorMsg = 'No patience distribution for this class.';
    return;
end

% Check FCFS scheduling
if sn.sched(queueIdx) ~= SchedStrategy.FCFS
    info.errorMsg = 'MAPMsG requires FCFS scheduling.';
    return;
end

% Check exponential service (single-phase)
serviceProc = sn.proc{queueIdx}{info.classIdx};
if isempty(serviceProc) || ~iscell(serviceProc) || length(serviceProc) < 2
    info.errorMsg = 'Invalid service process.';
    return;
end
% For exponential, the service process should be 1x1 matrices
if size(serviceProc{1}, 1) ~= 1
    info.errorMsg = 'MAPMsG requires exponential service (single-phase).';
    return;
end

% Extract service rate
info.serviceRate = -serviceProc{1}(1,1);
info.nServers = sn.nservers(queueIdx);

% Check MAP arrival process
arrivalProc = sn.proc{sourceIdx}{info.classIdx};
if isempty(arrivalProc) || ~iscell(arrivalProc) || length(arrivalProc) < 2
    info.errorMsg = 'Invalid arrival process.';
    return;
end

% All checks passed
isReneging = true;

end
