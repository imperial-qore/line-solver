function [bool, reason] = mam_dt_supports(sn)
% [BOOL, REASON] = MAM_DT_SUPPORTS(SN)
%
% Can the discrete-time (slotted) path of SolverMAM represent this model?
% REASON is '' when it can and otherwise names the construct it cannot.
%
% SOLVER_MAM_ANALYZER routes EVERY method to SOLVER_MAM_DT when
% SN_IS_DISCRETE_TIME says the laws live on a slot lattice, before the method
% name is read, so these rules bind every MAM method on such a model:
%   - open models only (a closed slotted model needs a level-dependent
%     discrete chain the Q-MAM discrete-time catalogue does not cover);
%   - one class only (independent per-class lattice sources fire in the same
%     slot with positive probability, and a batch of simultaneous arrivals of
%     different classes is not an MMAP[K]);
%   - FCFS single-server stations and the Source only (a slotted multiserver
%     queue needs the level-dependent boundary of Geo/Geo/c).
%
% ONE PREDICATE, TWO CALLERS: solver_mam_dt raises with it and
% SolverMAM.supportsModelMethod answers with it. It was a local function of
% solver_mam_dt until 2026-09-05, so the report offered every MAM method on a
% closed or multiclass slotted model and the run refused it.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

bool = false;
if ~sn_is_open_model(sn)
    reason = ['The discrete-time path supports open models only. A closed slotted ' ...
        'model needs a level-dependent discrete chain, which the Q-MAM discrete-time catalogue ' ...
        'does not cover.'];
    return
end
if sn.nclasses > 1
    reason = ['The discrete-time path supports one class only. Independent per-class ' ...
        'lattice sources fire in the same slot with positive probability, and a batch of ' ...
        'simultaneous arrivals of different classes is not an MMAP[K], which is what ' ...
        'Q_DT_MMAPK_PHK_1 consumes. Supply a single class, or a marked discrete arrival process ' ...
        'once DMMAP support lands.'];
    return
end
for ist = 1:sn.nstations
    switch sn.sched(ist)
        case SchedStrategy.EXT
            % source
        case SchedStrategy.FCFS
            if sn.nservers(ist) > 1
                reason = sprintf(['Station %d has %d servers. The discrete-time path models ' ...
                    'one server per station: a slotted multiserver queue needs the level-dependent ' ...
                    'boundary of Geo/Geo/c, which is not implemented.'], ist, sn.nservers(ist));
                return
            end
        otherwise
            reason = sprintf(['Station %d uses scheduling %s. The discrete-time path supports ' ...
                'FCFS single-server stations and the Source only.'], ist, SchedStrategy.toText(sn.sched(ist)));
            return
    end
end
bool = true;
reason = '';
end
