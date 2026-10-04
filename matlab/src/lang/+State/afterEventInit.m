function ctx = afterEventInit(sn)
% CTX = AFTEREVENTINIT(SN)
% Precomputes the loop-invariant setup of State.afterEvent so that hot
% callers (e.g., the solver_ssa Gillespie loop, which re-evaluates every
% synchronization at every step) avoid re-deriving it on each call.
%
% IMPORTANT: ctx caches sn.nservers, sn.cap and sn.classcap. It must be
% built AFTER any caller-side rewrite of these fields (solver_ssa rewrites
% them in its preamble: Inf servers at delay nodes and open-class cutoffs)
% and AFTER sn_nonmarkov_toph, and the SAME sn must be passed to afterEvent
% alongside ctx. Building ctx from a stale sn silently drops blocking and
% server-count semantics.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

M = sn.nstations;
R = sn.nclasses;

ctx.M = M;
ctx.R = R;
ctx.S = sn.nservers;
ctx.phasessz = sn.phasessz;
ctx.phaseshift = sn.phaseshift;
ctx.pie = sn.pie;
ctx.mu = sn.mu;
ctx.phi = sn.phi;
ctx.proc = sn.proc;
ctx.capacity = sn.cap;
ctx.classcap = sn.classcap;

lldscaling = sn.lldscaling;
if isempty(lldscaling)
    lldlimit = max(sum(sn.nclosedjobs),1);
    lldscaling = ones(M,lldlimit);
else
    lldlimit = size(lldscaling,2);
end
ctx.lldscaling = lldscaling;
ctx.lldlimit = lldlimit;

% sn.cdscaling holds a class-dependence handle only for the stations that
% declare one; the others are left empty (see getLimitedClassDependence), since
% a constant 1 is not a valid class-dependent RATE. The state machinery indexes
% every station unconditionally, so fill the gaps with the neutral scaling here.
cdscaling = sn.cdscaling;
if isempty(cdscaling)
    cdscaling = cell(M,1);
end
neutral = @(ni) 1;
for i = 1:M
    if i > numel(cdscaling) || isempty(cdscaling{i})
        cdscaling{i} = neutral;
    end
end
% Joint-dependence handles (sn.jdscaling, non-product-form eta_i) are applied
% by the same state machinery. cd and jd are evaluated identically, so fold
% them into a single effective per-station handle eta_i(n)*beta_i(n) (scalar
% .* vector broadcasts). For stations with only one of the two, the other is
% the neutral 1, so the product reproduces the single-mechanism case exactly.
jdscaling = sn.jdscaling;
if ~isempty(jdscaling)
    for i = 1:M
        if i <= numel(jdscaling) && ~isempty(jdscaling{i})
            cdh = cdscaling{i};
            jdh = jdscaling{i};
            cdscaling{i} = @(ni) cdh(ni) .* jdh(ni);
        end
    end
end
ctx.cdscaling = cdscaling;

ctx.ismkvmod = false(sn.nnodes,1);
ctx.ismkvmodclass = cell(sn.nnodes,1);
for ind=1:sn.nnodes
    if sn.isstation(ind)
        ist = sn.nodeToStation(ind);
        ctx.ismkvmod(ind) = any(sn.procid(ist,:)==ProcessType.MAP | sn.procid(ist,:)==ProcessType.MMPP2);
        ismkvmodclass = zeros(R,1);
        for r=1:R
            ismkvmodclass(r) = any(sn.procid(ist,r)==ProcessType.MAP | sn.procid(ist,r)==ProcessType.MMPP2);
        end
        ctx.ismkvmodclass{ind} = ismkvmodclass;
    end
end
end
