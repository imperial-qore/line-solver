function tf = mam_gk1_applicable(sn, ist, K)
% MAM_GK1_APPLICABLE  True when station IST should be answered by MMAP[K]/G[K]/1.
%
% The generic MMAPPH1FCFS path reads the service law out of SN.PROC, which holds
% its PHASE-TYPE FIT: for a Uniform, a Gamma, a Pareto, a Weibull, a Lognormal or
% a Det that fit matches the mean and, once the SCV exceeds one, nothing else. He
% (2001) needs only the TRANSFORM of the original law, which SN.LST carries, so
% wherever a class declares one of those laws the exact analysis is available and
% the fit is not needed.
%
% Two structural conditions. The result is a /1, so a multiserver station is out.
% And every class must carry a usable transform handle, since the analysis takes
% all K of them together rather than one at a time.
%
% A matrix-exponential service qualifies too, its transform being rational and
% SN.LST evaluating it, which is why the ME warning about falling back to a
% phase-type approximation is suppressed when this returns true. A RAP does NOT
% qualify. He's analysis assumes INDEPENDENT service times, so reading a
% correlated service through its marginal transform would discard exactly the
% autocorrelation the RAP was declared to carry.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

tf = false;
if sn.nservers(ist) ~= 1
    return;
end
if isempty(sn.lst) || numel(sn.lst) < ist || isempty(sn.lst{ist})
    return;
end
nonPh = [ProcessType.DET, ProcessType.UNIFORM, ProcessType.GAMMA, ...
    ProcessType.PARETO, ProcessType.WEIBULL, ProcessType.LOGNORMAL, ...
    ProcessType.ME];
% The DECLARED tags, snapshotted by solver_mam_analyzer before sn_nonmarkov_toph:
% that conversion retags procid APH/ME/MAP, so reading the live procid would see
% the surrogate and leave this path unreachable for every law but Det.
if isfield(sn, 'procidDeclared') && ~isempty(sn.procidDeclared)
    procid = sn.procidDeclared;
else
    procid = sn.procid;
end
if ~any(ismember(procid(ist, 1:K), nonPh))
    return;
end
for k = 1:K
    if numel(sn.lst{ist}) < k || isempty(sn.lst{ist}{k})
        return;
    end
end
tf = true;
end
