function placement = warmStartPlacement(initSolver, sn)
% PLACEMENT = WARMSTARTPLACEMENT(INITSOLVER, SN)
%
% Integer job placement (nstations x nclasses) decided by the steady-state
% solution of an auxiliary solver.
%
% If the auxiliary solver is a SolverCTMC, the exact stationary distribution
% over the aggregate state space is computed and the placement is the mode of
% that distribution (the most probable aggregate state). For any other network
% solver, the steady-state mean queue lengths are used instead and rounded to
% an integer placement that conserves each closed-class population.
%
% Only service stations (Queue/Delay) receive an initial population.

M = sn.nstations;
K = sn.nclasses;
if isa(initSolver, 'SolverCTMC')
    placement = placementFromCtmcSteadyState(initSolver, M, K);
else
    placement = placementFromMeanQLen(initSolver, sn, M, K);
end
end

function placement = placementFromCtmcSteadyState(ctmcSolver, M, K)
% Placement from the exact CTMC stationary distribution: solve the CTMC,
% aggregate the stationary probabilities over the aggregate (per-station,
% per-class job count) state space, and return the aggregate state of maximum
% stationary probability.
snc = ctmcSolver.getStruct();
[Q,~,SSq,~,~,~,snc] = solver_ctmc(snc, ctmcSolver.getOptions());
pi = ctmc_solve_reducible(Q);
pi = pi(:);
pi(pi < GlobalConstants.Zero) = 0;

% Aggregate the stationary probability over identical aggregate states and
% locate the mode of the aggregate distribution.
[uRows, ~, ic] = unique(SSq, 'rows');
aggrProb = accumarray(ic, pi);
[~, bidx] = max(aggrProb);
modeState = uRows(bidx, :);

% Map the aggregate-state columns (nclasses per stateful station, in station
% order) onto the placement matrix.
placement = zeros(M, K);
col = 1;
for i = 1:snc.nstations
    ind = snc.stationToNode(i);
    if ~snc.isstateful(ind)
        continue;
    end
    isService = snc.nodetype(ind) == NodeType.Queue || snc.nodetype(ind) == NodeType.Delay;
    for r = 1:K
        v = modeState(col);
        col = col + 1;
        if isService
            placement(i, r) = v;
        end
    end
end
end

function placement = placementFromMeanQLen(initSolver, sn, M, K)
% Placement from the steady-state mean queue lengths of a generic network
% solver: floor the per-station means and distribute the residual closed-class
% jobs by largest remainder so each closed population is conserved.
QN = initSolver.getAvgQLen();
placement = zeros(M, K);
for r = 1:K
    njobs = sn.njobs(r);
    isClosed = isfinite(njobs);
    floors = zeros(1, M);
    fracs = -ones(1, M);
    eligible = false(1, M);
    for i = 1:M
        ind = sn.stationToNode(i);
        eligible(i) = sn.nodetype(ind) == NodeType.Queue || sn.nodetype(ind) == NodeType.Delay;
        if ~eligible(i)
            continue;
        end
        m = max(0, QN(i, r));
        if isClosed
            floors(i) = floor(m);
            fracs(i) = m - floors(i);
        else
            floors(i) = round(m);
        end
        placement(i, r) = floors(i);
    end
    if isClosed
        % Largest-remainder apportionment of the residual jobs.
        residual = round(njobs) - sum(floors);
        while residual > 0
            [bf, bi] = max(fracs);
            if bf < 0
                % All remainders consumed: place the leftover jobs at the
                % reference station so the population is conserved.
                refStation = sn.refstat(r);
                if ~eligible(refStation)
                    refStation = find(eligible, 1);
                end
                placement(refStation, r) = placement(refStation, r) + residual;
                residual = 0;
                break;
            end
            placement(bi, r) = placement(bi, r) + 1;
            fracs(bi) = -1;
            residual = residual - 1;
        end
    end
end
end
