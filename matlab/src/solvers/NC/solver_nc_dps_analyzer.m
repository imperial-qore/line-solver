function [Q,U,R,T,C,X,lG,runtime,iter,method] = solver_nc_dps_analyzer(sn, options)
% [Q,U,R,T,C,X,LG,RUNTIME,ITER,METHOD] = SOLVER_NC_DPS_ANALYZER(SN, OPTIONS)
%
% Heavy-usage asymptotic analysis of the closed two-station network with one
% think (infinite-server) station and one discriminatory processor-sharing
% station, by the generating-function expansion of
%
%   J.A. Morrison, "Asymptotic analysis of a large closed queueing network
%   with discriminatory processor sharing", Queueing Systems 9 (1991) 191-214.
%
% Admitted only on the exact shape NC_IS_DPS_MODEL tests for. The kernel is
% NPFQN_DPS_MORRISON; this function maps the model struct onto it and lifts the
% per-class DPS results into the station-by-class arrays the NC analyzers
% return.
%
% THERE IS NO NORMALIZING CONSTANT HERE. A DPS station is not product-form --
% that is the premise of the paper -- so lG is returned as NaN, as it is on the
% maximum-entropy route. NC hosts this method because NC is where LINE keeps
% the asymptotic expansions of generating functions and normalizing-constant
% integrals (pana, mmint2, le, ble, gleint, rayint), which is the family
% Morrison's expansion belongs to, not because a constant is being computed.
%
% Response times come from Little's law on the queue-length result rather than
% from the expanded RESULT 2 (eq. 4.17), so that Q = R*T holds exactly in the
% returned table; the two agree to the order of the approximation, since
% Morrison derives (4.17) as the ratio (4.11)/(4.15).
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
Tstart = tic;
iter = 1;
method = 'morrison';
lG = NaN;                                  % no product form: no normalizing constant

% The shape is re-checked here, not assumed from the caller. This analyzer is
% reachable from runAnalyzer, solver_nc_analyzer, ncDispatch and directly from
% user code, and every quantity below -- the think time, the DPS service time,
% the weights, the visit ratios -- is meaningless off the shape the expansion
% was derived for. A wrong number is worse than a refusal.
if ~nc_is_dps_model(sn)
    line_error(mfilename, ['solver_nc_dps_analyzer applies only to a CLOSED network of exactly ' ...
        'two stations, one infinite-server (think) station and one single-server DPS station with ' ...
        'exponential service and one visit each per cycle (see nc_is_dps_model).']);
end

M = sn.nstations;
K = sn.nclasses;
iInf = find(sn.sched == SchedStrategy.INF);
iDps = find(sn.sched == SchedStrategy.DPS);

% ---- per-class visits, normalized at the reference station ----------------
V = zeros(M, K);
for r = 1:K
    c = find(sn.chains(:, r));
    vis = sn.visits{c};
    for ist = 1:M
        V(ist, r) = vis(sn.stationToStateful(ist), r);
    end
    vref = V(sn.refstat(r), r);
    if vref > 0
        V(:, r) = V(:, r) / vref;
    end
end

% ---- Morrison's inputs: populations, think times, DPS service times, weights
Npop = sn.njobs(:).';
Z = 1 ./ sn.rates(iInf, :);
St = 1 ./ sn.rates(iDps, :);
w = sn.schedparam(iDps, :);

line_debug(options, 'NC DPS analyzer: N=%s, Z=%s, S=%s, w=%s', ...
    mat2str(Npop), mat2str(Z, 4), mat2str(St, 4), mat2str(w, 4));

[Qdps, Rmorr, Xr, aux] = npfqn_dps_morrison(Npop, Z, St, w);

if any(Qdps < 0) || any(Qdps > Npop) || any(~isfinite(Qdps))
    line_warning(mfilename, ['The Morrison expansion returned queue lengths outside [0,N] ' ...
        '(usage rho=%g, a=%g). The model is outside the moderately-heavy regime the expansion ' ...
        'assumes; treat the result as unreliable and cross-check with SolverMVA or SolverCTMC.'], ...
        aux.rho, aux.a);
    Qdps = min(max(Qdps, 0), Npop);
    Xr = (Npop - Qdps) ./ Z;
end

line_debug(options, 'NC DPS analyzer: rho=%g, a=%g, leading E[n]=%s, two-term E[n]=%s, RESULT 2 E[T]=%s', ...
    aux.rho, aux.a, mat2str(aux.Qlead, 4), mat2str(Qdps, 4), mat2str(Rmorr, 4));

% ---- lift into the station-by-class arrays --------------------------------
Q = zeros(M, K); U = zeros(M, K); R = zeros(M, K); T = zeros(M, K);
X = Xr;

Q(iDps, :) = Qdps;
Q(iInf, :) = Npop - Qdps;                  % population conservation (exact, closed)

for ist = 1:M
    T(ist, :) = X .* V(ist, :);
end

U(iInf, :) = Q(iInf, :);                   % INF utilization convention
c = sn.nservers(iDps);
if ~isfinite(c) || c <= 0, c = 1; end
U(iDps, :) = X .* (V(iDps, :) .* St) / c;

for ist = 1:M
    for r = 1:K
        if T(ist, r) > 0
            R(ist, r) = Q(ist, r) / T(ist, r);
        end
    end
end

C = zeros(1, K);
for r = 1:K
    if X(r) > 0
        C(r) = Npop(r) / X(r);
    end
end

runtime = toc(Tstart);
end
