function reason = mam_buffer_refusal(sn)
% REASON = MAM_BUFFER_REFUSAL(SN)
%
% Why no MAM method can answer a model whose finite buffer a CLOSED class can
% fill, or '' when every binding buffer is reached by open classes only.
%
% A MAM analyzer represents a finite buffer as a LOSS buffer: solver_mam_basic
% (default, dec.source, dec.poisson) solves an M/M/c/K or an MMAP[K]/G/1/K, and
% solver_mam (dec.mmap), solver_mam_basic_mmap (dec.source.mmap) and
% solver_mna_open (mna) truncate and renormalize the same way. That is the
% right model for an OPEN class, whose refused arrival is lost. A closed job
% that finds no room BLOCKS instead -- LINE disables the upstream departure
% and holds the job where it is -- and no MAM analyzer blocks: the loss
% formulas answer a different system, and the closed routes of 'default'
% (solver_mam_ldqbd, solver_mam_bgchain, solver_mna_closed) read no sn.cap or
% sn.classcap at all. So the pair is refused rather than answered.
%
% ONE PREDICATE, TWO CALLERS: SolverMAM.supportsModelMethod reports it and
% solver_mam_analyzer raises it. Only a buffer that can BIND counts, which is
% what SN_GET_BUFFER_SIZE decides: refreshCapacity derives a finite classcap
% (the chain population) at every station of every closed model, so a plain
% finiteness test would refuse every closed model. A Cache builds its own
% capped retrieval queues and is exempt, as in sn_has_blocking.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

reason = '';
if isfield(sn, 'nodetype') && any(sn.nodetype == NodeType.Cache)
    return
end
if ~isfield(sn, 'classcap') || isempty(sn.classcap)
    return
end
closedClass = isfinite(sn.njobs(:)');
for ist = 1:sn.nstations
    if ~isfinite(sn_get_buffer_size(sn, ist))
        continue
    end
    served = sn.classcap(ist, :) > 0;
    r = find(closedClass & served, 1);
    if isempty(r)
        continue
    end
    reason = sprintf(['Station %s carries a finite capacity that binds for the closed class %s. ' ...
        'SolverMAM represents a finite buffer as a LOSS buffer (M/M/c/K, MMAP[K]/G/1/K, ' ...
        'truncate-and-renormalize), which is the open-class model: a closed job that finds ' ...
        'no room blocks instead, and no MAM analyzer blocks, while the closed routes of the ' ...
        'default method (ldqbd, bgchain, mna) read no capacity at all. Use SolverCTMC, ' ...
        'SolverSSA or SolverLDES, or SolverMVA with method ''sqd'' for blocking after service.'], ...
        sn.nodenames{sn.stationToNode(ist)}, sn.classnames{r});
    return
end
end
