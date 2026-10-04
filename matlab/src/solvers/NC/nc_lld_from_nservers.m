function lldscaling = nc_lld_from_nservers(sn, latticeMax)
% LLDSCALING = NC_LLD_FROM_NSERVERS(SN, LATTICEMAX)
%
% Represents every finite multiserver station of a closed model as the exact
% load-dependent rate lattice mu(n)=min(n,c), the form SOLVER_NCLD consumes and
% PFQN_NCLD solves exactly. Returns [] when the conversion does not apply, in
% which case the caller keeps Seidmann's approximation:
%
%   - no finite station has more than one server (forcing the lattice on an
%     all-single-server model mishandles a self-looping closed chain)
%   - the total population is not finite (an open or mixed model)
%   - the per-chain population lattice exceeds LATTICEMAX
%
% LATTICEMAX prices the exact enumeration SOLVER_NCLD performs, whose cost is
% unbounded in N while the multiserver test is topological. Callers pass the
% same 6000-state budget SolverCTMC uses for exact enumeration.
%
% The server count c is left on sn.nservers: SOLVER_NCLD keeps reading it to
% normalize utilization, and LLD_ENCODES_MULTISERVER recognises the pairing.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

lldscaling = [];

if ~isempty(sn.lldscaling)
    return % an explicit limited load-dependence already occupies the slot
end
if ~any(sn.nservers(isfinite(sn.nservers)) > 1)
    return
end
if ~all(isfinite(sn.njobs))
    return
end

Nt = sum(sn.njobs);
if ~isfinite(Nt) || Nt < 1
    return
end

if nargin >= 2 && ~isempty(latticeMax)
    Nchain = sn.njobs(:)';
    if ~isempty(sn.chains)
        Nchain = zeros(1, size(sn.chains,1));
        for c = 1:size(sn.chains,1)
            inc = find(sn.chains(c,:) > 0);
            Nchain(c) = sum(sn.njobs(inc(isfinite(sn.njobs(inc)))));
        end
    end
    if prod(1 + Nchain(isfinite(Nchain))) > latticeMax
        return
    end
end

lldscaling = ones(sn.nstations, Nt);
for i = 1:sn.nstations
    if sn.nservers(i) > 1 && isfinite(sn.nservers(i))
        lldscaling(i,:) = min(1:Nt, sn.nservers(i));
    end
end
end
