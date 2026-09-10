function tf = mam_transient_qbd_applicable(sn)
% MAM_TRANSIENT_QBD_APPLICABLE True when transient analysis of the model should
% use the Laplace-domain transient QBD solver (solver_mam_transient_qbd) rather
% than the libQBD/expm fast path (solver_mam_ldqbd_transient).
%
% The Laplace solver is selected for single-server open queues whose arrival is
% non-Poisson (correlated/PH-renewal MAP with >1 phase) or whose service is a
% correlated MAP (non-renewal), which the fast path cannot represent exactly.
% Poisson arrival with PH/exp service (M/M/1, M/PH/1) and M/M/c stay on the fast
% path.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
tf = false;

if sn.nclasses ~= 1 || ~isinf(sn.njobs')
    return;   % Laplace solver handles single-class open models only
end

sourceIdx = find(sn.sched == SchedStrategy.EXT);
queueIdx  = find(sn.sched == SchedStrategy.FCFS);
if numel(sourceIdx) ~= 1 || numel(queueIdx) ~= 1
    return;
end
if sn.nservers(queueIdx) ~= 1
    return;   % Laplace solver is single-server; multiserver stays on fast path
end

arrProc = sn.proc{sourceIdx}{1};
svcProc = sn.proc{queueIdx}{1};
Da0 = arrProc{1}; Da1 = arrProc{2};
Ds0 = svcProc{1}; Ds1 = svcProc{2};

arrivalIsPoisson  = (size(Da0, 1) == 1);
serviceIsRenewal  = is_renewal_map(Ds0, Ds1);

tf = ~arrivalIsPoisson || ~serviceIsRenewal;
end

function tf = is_renewal_map(D0, D1)
% A PH/renewal service, embedded as a MAP, has D1 = t * alpha (rank one): the
% post-completion phase does not depend on the completing phase.
ns = size(D0, 1);
if ns == 1
    tf = true;
    return;
end
tExit = D1 * ones(ns, 1);
alpha = map_pie({D0, D1});
tf = norm(D1 - tExit * alpha, 'fro') < 1e-9 * max(1, norm(D1, 'fro'));
end
