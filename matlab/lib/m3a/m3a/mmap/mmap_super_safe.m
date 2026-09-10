function SUP = mmap_super_safe(MMAPs, maxorder, method)
% SUP = mmap_super_safe(MMAPs, maxorder)
% Superposition of marked MAPS

if nargin == 2
    method = 'default';
end
empty = cellfun(@isempty, MMAPs);
MMAPs(empty)=[];
% A component with an all-zero aggregate arrival matrix (D1 == 0) has zero
% arrival rate: it contributes nothing to the superposition. When such a
% component carries more than one phase (e.g. the transient slow phase of a
% high-SCV APH fit that reaches a station the flow never visits), its phase
% generator D0+D1 is absorbing, so map_scv/map_pie/map_prob -> ctmc_solve fail
% with "no recurrent state". Canonicalize it to the equivalent order-1 null
% (matching the marking count), whose moments are well defined and whose
% superposition is an identity. Detection uses the arrival-matrix norm, not
% mmap_lambda, because mmap_lambda itself calls map_prob -> ctmc_solve.
for iz = 1:numel(MMAPs)
    if size(MMAPs{iz}{1},1) > 1 && norm(full(MMAPs{iz}{2}),1) < 1e-13
        MMAPs{iz} = mmap_exponential(zeros(1,numel(MMAPs{iz})-2), 1);
    end
end
scv_unmarked = cellfun(@map_scv, MMAPs);
[~,Iset]=sort(scv_unmarked); % sort flows with small scv first
% Mark provenance, so the SCV sort cannot permute the marking. mmap_super
% concatenates the marks of its operands in FOLD order, and every caller reads
% mark k as its own k-th class, so folding low-SCV first silently renamed the
% classes whenever the components' SCVs did not happen to be in caller order.
% A zero-rate component has an infinite mean, hence SCV NaN, which sort() puts
% last, so a chain that never visits a station moved its marks ahead of a chain
% that does. The sort stays, as the numerical heuristic it was meant to be, and
% the marks are permuted back at the end.
markcounts = cellfun(@(x) numel(x)-2, MMAPs);
markbase = cumsum([0, markcounts(:)']);
outorder = [];
SUP = {};
for i=Iset(:)' % low-SCV first
    outorder = [outorder, markbase(i)+(1:markcounts(i))]; %#ok<AGROW>
    % Bound the order of each individual flow to maxorder. A single flow
    % whose order already exceeds maxorder (e.g. a high-order Erlang from a
    % near-deterministic APH fit) would otherwise pass through uncapped as
    % the superposition base and blow up downstream matrix-analytic solves.
    if length(MMAPs{i}{1}) > maxorder
        if maxorder >= 2
            MMAPs{i} = mamap2m_fit_gamma_fb_mmap(MMAPs{i});
        else
            MMAPs{i} = mmap_exponential(mmap_lambda(MMAPs{i}));
        end
    end
    if isempty(SUP) % is this is the first MMAP
        SUP = MMAPs{i};
        if maxorder == 1
            %  then treat it as a Poisson process if the limit is 1
            SUP = mmap_exponential(mmap_lambda(MMAPs{i}));
        end
    else
        if length(SUP{1}) * length(MMAPs{i}{1}) > maxorder
            %  otherwise treat it as a marked AMAP(2) process if the limit is >1
            if length(SUP{1}) * 2 <= maxorder
                SUP = mmap_super(SUP, mamap2m_fit_gamma_fb_mmap(MMAPs{i}), method);
            else
                SUP = mmap_super(SUP, mmap_exponential(mmap_lambda(MMAPs{i})), method);
            end
        else
            SUP = mmap_super(SUP,MMAPs{i}, method);
        end
    end

end
% Restore the caller's mark order; 'match' keeps one mark per class and fails
% the count test, so it is left alone.
if numel(SUP)-2 == numel(outorder) && ~issorted(outorder)
    [~,perm] = sort(outorder);
    SUP = {SUP{1}, SUP{2}, SUP{2+perm}};
end
end