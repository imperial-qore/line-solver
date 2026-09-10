function fcn = fluid_conservation_guard(phases, njobs, chains, tol)
% FCN = FLUID_CONSERVATION_GUARD(PHASES, NJOBS, CHAINS, TOL)
%
% An ODE OutputFcn that halts an integration whose state has left the model.
%
% WHY IT EXISTS. The moment-closure drift can leave the simplex: on a station
% where min(n,c) is not the identity the Gaussian correction to the per-class
% share can drive a coordinate negative, the drift is conservative so another
% grows to match, and the trajectory runs away. `odeset` carries `NonNegative`
% over every coordinate, so the excursion is CLAMPED rather than reported --
% which INJECTS mass, collapses the step size, and leaves the window never
% returning. That is not a slow solve, it is a solve that does not terminate:
% one MATLAB suite run sat in `test_CQN_Cox_CS_9` for 3h16m, and the 2026-08-27
% run was killed after `test11_interlock_lqnx` had held the suite for 100
% minutes, taking every block after it down with it. See
% _kb/06-solver-catalog.md.
%
% THE TEST IS AN EXACT INVARIANT, not a heuristic bound on time or magnitude.
% The drift conserves the population of every CLOSED CHAIN exactly, so any
% deviation is the clamp injecting mass and nothing else. TOL is therefore a
% generous fraction of that population rather than a numerical tolerance:
% the integrator's own error is ~1e-4 relative, while the documented excursion
% reaches 5.2e4 against a true population of 0.05. A closed model whose
% population has moved by TOL is no longer solving the model, whatever it is
% converging to.
%
% THE CHAIN IS THE CONSERVED UNIT, NOT THE CLASS, and the difference is the
% whole correctness of this check. SN.NJOBS(k) is the population class k STARTS
% with; class switching then moves jobs between the classes of one chain, so
% only the chain total is invariant. Watching classes instead condemns every
% class-switching model out of hand -- measured on `cqn_twoclass_hyperl` (313
% of 447 accepted states), on `init_state_ps` (286 of 310) and on every one of
% the 162 fluid layers an LQN builds under the `srvn.cs` encoding, where the
% chain sum never moved at all. A cache model is the same story with the
% hit/miss classes.
%
% A wall-clock budget would have caught the same hang and was rejected: it makes
% the answer depend on how busy the host is, so the same model would fall back
% on one machine and not on another. This invariant is deterministic.
%
% Parameters:
%   phases - (nstations,nclasses) phases per (station,class), as SOLVER_FLUID
%            builds it from Mu; fixes the coordinate blocks of the state vector
%   njobs  - (1,nclasses) class populations; a chain summing to Inf is open,
%            its population is not conserved, and it is therefore not watched
%   chains - (nchains,nclasses) membership, i.e. SN.CHAINS. Empty degrades to
%            one chain per class, which is the safe reading rather than a
%            guess: with no chain map there is no class switching to merge
%            classes, so each class IS its own conserved unit
%   tol    - relative deviation that trips the guard (0.1 = 10%)
%
% Returns:
%   fcn - handle for odeset('OutputFcn', ...). Returns status 1 to halt, which
%         the caller reads as a short final time and turns into
%         'LINE:FluidNonHyperbolic' so the existing fallback ladder answers.
%
% See also SOLVER_FLUID_ITERATION, SOLVER_FLUID_MOMENTS, FLUID_LYAPUNOV.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

[M, K] = size(phases);
if nargin < 3 || isempty(chains)
    chains = eye(K);
end
nchains = size(chains, 1);
% Coordinate blocks, in the layout SOLVER_FLUID reads QN back out of:
% station-major, then class, then that pair's phases. One block per CHAIN,
% pooling the columns of every class the chain holds.
idx = cell(1, nchains);
target = zeros(1, nchains);
watched = false(1, nchains);
for ch = 1:nchains
    members = find(chains(ch,:));
    members = members(members <= numel(njobs));
    total = sum(njobs(members));
    if isempty(members) || ~isfinite(total) || total <= 0
        continue % open, or absent: no conserved population to check
    end
    cols = [];
    for k = members
        for ist = 1:M
            if phases(ist,k) <= 0
                continue
            end
            shift = sum(sum(phases(1:ist-1,:))) + sum(phases(ist,1:k-1)) + 1;
            cols = [cols, shift:(shift + phases(ist,k) - 1)]; %#ok<AGROW>
        end
    end
    if isempty(cols)
        continue
    end
    idx{ch} = sort(cols);
    target(ch) = total;
    watched(ch) = true;
end
watchList = find(watched);

fcn = @guard;

    function status = guard(~, y, flag)
        status = 0;
        if strcmp(flag, 'done') || isempty(y) || isempty(watchList)
            return
        end
        % 'init' hands the initial state, '' a column per accepted step.
        for c = 1:size(y, 2)
            yc = y(:, c);
            for kk = watchList
                if idx{kk}(end) > numel(yc)
                    continue % a caller integrating a different state vector
                end
                mass = sum(yc(idx{kk}));
                if abs(mass - target(kk)) > tol * max(1, target(kk))
                    status = 1;
                    return
                end
            end
        end
    end
end
