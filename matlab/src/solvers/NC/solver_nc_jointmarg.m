function [Pr,lPr,lG,runtime] = solver_nc_jointmarg(sn, options, nvec, engine, lG)
% [PR,LPR,LG,RUNTIME] = SOLVER_NC_JOINTMARG(SN, OPTIONS, NVEC, ENGINE, LG)
%
% Joint probability that station I holds NVEC(I) jobs IN TOTAL, all classes
% summed out. This is NOT solver_nc_jointaggr, which fixes the per-class
% population of every station: each state here is the sum of jointaggr over
% the whole fibre of per-class tables with these row sums, and that fibre
% grows combinatorially. The permanent evaluates the sum in closed form
% (pfqn_jointmarg).
%
% LG, when supplied, is the log normalizing constant of the model, the same
% one that @SolverNC/getProbAggr caches as result.Prob.logNormConstAggr.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

Tstart = tic;

if nargin < 4 || isempty(engine)
    engine = 'exact';
end

M = sn.nstations;
nvec = nvec(:)';
if numel(nvec) ~= M
    line_error(mfilename,'The occupancy vector has %d entries but the model has %d stations.', numel(nvec), M);
end

[Lchain,~,~,~,Nchain] = sn_get_demands_chain(sn);
Lchain(~isfinite(Lchain)) = 0;

solver_nc_jointmarg_supports(sn, Nchain);

infset = find(isinf(sn.nservers));
infset = infset(:)';

if nargin < 5 || isempty(lG) || ~isfinite(lG)
    lG = [];
end

[Pr,lPr] = pfqn_jointmarg(nvec, Lchain, Nchain, infset, lG, engine);

if isempty(lG)
    % Recover the constant the API computed, so the caller can cache it and
    % the sweep over a whole lattice pays for it once.
    isinfrow = false(1,M);
    isinfrow(infset) = true;
    if any(isinfrow)
        Z = sum(Lchain(isinfrow,:),1);
    else
        Z = zeros(1,size(Lchain,2));
    end
    [~,lG] = pfqn_ca(Lchain(~isinfrow,:), Nchain, Z);
end

runtime = toc(Tstart);
end

function solver_nc_jointmarg_supports(sn, Nchain)
% The permanent identity supplies one n_i! per queueing station and none per
% infinite server. A multiserver or load-dependent station has neither, so it
% is refused by name rather than approximated.
if any(~isfinite(Nchain)) || any(isinf(sn.njobs))
    line_error(mfilename,'getProbSysMarg requires a closed model: the joint law of the total queue lengths is not defined when a class has an infinite population.');
end
if ~isempty(sn.lldscaling)
    line_error(mfilename,'getProbSysMarg does not support load-dependent stations (sn.lldscaling is set): the permanent identity supplies exactly one n_i! per queueing station.');
end
if ~isempty(sn.cdscaling)
    line_error(mfilename,'getProbSysMarg does not support class-dependent scaling (sn.cdscaling is set).');
end
ms = find(isfinite(sn.nservers) & sn.nservers > 1);
if ~isempty(ms)
    line_error(mfilename,'getProbSysMarg does not support the multiserver station %d (%d servers): the permanent identity supplies exactly one n_i! per queueing station.', ms(1), sn.nservers(ms(1)));
end
end
