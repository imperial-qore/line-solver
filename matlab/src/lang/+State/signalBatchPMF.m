function [kvals, kprobs] = signalBatchPMF(sn, class, ntot)
% [KVALS, KPROBS] = SIGNALBATCHPMF(SN, CLASS, NTOT)
%
% Batch-size distribution of the jobs removed by signal class CLASS when
% NTOT eligible jobs are present. Without a removal distribution a signal
% removes exactly one job. With one, the pmf is clipped at NTOT: a batch
% larger than the eligible population empties it instead of driving the
% queue negative, so the tail P(B >= NTOT) lumps onto "remove NTOT". This
% matches the min(B, n) clipping in SolverLDES and the tail term that
% SolverMAM puts on the empty state.

kvals = 1;
kprobs = 1;
if ~isfield(sn, 'signalremdist') || numel(sn.signalremdist) < class || isempty(sn.signalremdist{class})
    return
end
dist = sn.signalremdist{class};

support = 0:ntot;
pmf = zeros(1, numel(support));
for i = 1:numel(support)
    pmf(i) = dist.evalPMF(support(i));
end
% Everything at or beyond ntot removes the whole eligible population.
head = pmf(1:end-1);           % B = 0 .. ntot-1
tail = max(0, 1 - sum(head));  % P(B >= ntot)
kvals = [support(1:end-1), ntot];
kprobs = [head, tail];

keep = kprobs > 0;
kvals = kvals(keep);
kprobs = kprobs(keep);
if isempty(kvals)
    kvals = 1;
    kprobs = 1;
    return
end
kprobs = kprobs / sum(kprobs);
end
