function nstates = mam_bgchain_states(sn, options)
% NSTATES = MAM_BGCHAIN_STATES(SN, OPTIONS)
%
% Number of states of the background-chain CTMC that SOLVER_MAM_BGCHAIN would
% build on this model, WITHOUT building it. The size is what decides whether
% bgchain is affordable, and MAM_BGCHAIN_CTMC only discovers it after the
% partition is fixed, so the default-method chooser needs it up front.
%
% The count mirrors the partition SOLVER_MAM_BGCHAIN uses: a pass carries the
% tagged closed chain as background class 1 and the demand-similar groups of
% the other closed chains as classes 2..1+G, each class enumerating the
% compositions of its population over the stations its members visit. Merging
% two chains onto the UNION of their supports can raise the count as easily as
% lower it, so the passes are enumerated rather than bounded, and the largest
% is returned: that is the one MAM_BGCHAIN_CTMC would refuse.
%
% The per-class count below is the ROW COUNT OF STATE.SPACECLOSEDSINGLE(m, Nb),
% the primitive MAM_BGCHAIN_CTMC enumerates each block with, evaluated in closed
% form as nchoosek(Nb+m-1, m-1) so that the size can be compared without paying
% for the enumeration.
%
% Returns 0 when bgchain does not apply to the model at all (no closed chain,
% or no station visited by one), which leaves the caller to its own guard.
%
% See also SOLVER_MAM_BGCHAIN, MAM_BGCHAIN_CTMC, MAM_BGCHAIN_GROUPS.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

nstates = 0;

if nargin < 2 || ~isstruct(options) || ~isfield(options, 'config') || isempty(options.config)
    bgaggr = 1;
else
    if isfield(options.config, 'bgaggr') && ~isempty(options.config.bgaggr)
        bgaggr = options.config.bgaggr;
    else
        bgaggr = 1;
    end
end

C = sn.nchains;
N = sn.njobs';
[Lchain, ~, Vchain, ~, Nchain] = sn_get_demands_chain(sn);

isopenchain = false(1, C);
for c = 1:C
    isopenchain(c) = any(isinf(N(sn.inchain{c})));
end
closedChains = find(~isopenchain & Nchain > 0);
R = numel(closedChains);
if R == 0
    return;
end

% Support of the background chain: the stations the closed chains visit
inbg = false(1, sn.nstations);
for c = closedChains
    inbg = inbg | (Vchain(:, c)' > GlobalConstants.Zero);
end
cst = find(inbg);
if isempty(cst)
    return;
end

naggr = min(max(round(bgaggr), 1), max(R - 1, 1));
noAggr = (R == 1) || (naggr >= R - 1);

if noAggr
    passes = {num2cell(closedChains)};
else
    passes = cell(1, R);
    for ridx = 1:R
        r = closedChains(ridx);
        others = setdiff(closedChains, r);
        grp = mam_bgchain_groups(Lchain(cst, others), naggr);
        members = cell(1, 1 + naggr);
        members{1} = r;
        for g = 1:naggr
            members{1 + g} = others(grp == g);
        end
        passes{ridx} = members;
    end
end

for pidx = 1:numel(passes)
    members = passes{pidx};
    n = 1;
    for b = 1:numel(members)
        mem = members{b};
        Nb = sum(Nchain(mem));
        sb = false(1, numel(cst));
        for oi = 1:numel(mem)
            sb = sb | (Vchain(cst, mem(oi))' > GlobalConstants.Zero);
        end
        m = sum(sb);
        if m == 0
            m = 1;   % an empty class still needs one slot to be indexed by
        end
        n = n * bgchain_binomial(round(Nb) + m - 1, m - 1);
        if ~isfinite(n)
            nstates = Inf;
            return;
        end
    end
    nstates = max(nstates, n);
end

end

function c = bgchain_binomial(n, k)
% nchoosek in floating point, so a chain far above any usable size still
% compares instead of raising the exact-arithmetic warning of nchoosek.
if k < 0 || k > n
    c = 0;
    return;
end
kk = min(k, n - k);
c = 1;
for i = 1:kk
    c = c * (n - kk + i) / i;
end
end
