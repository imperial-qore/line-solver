function rates = ode_rates_closing_factors(x, M, K, enabled, q_indices, Kic, nservers, w, sched_id, sigma2, lld, covblk)
% RATES = ODE_RATES_CLOSING_FACTORS(...) per-coordinate service share, before
% event indexing and before the constant rate factors are applied. Kept
% separate so the moment-closure methods can read the same service shares the
% ODE integrated (see SOLVER_FLUID_MOMENTS).
%
% SIGMA2 is the per-station closure variance (0 or omitted: first-order
% closure). LLD is `sn.lldscaling` (empty: no load dependence). COVBLK is a
% per-station cell of coordinate covariance blocks, which closes the DPS
% capacity-share ratio at second order (see FLUID_SHARE_CLOSURE). All three
% leave the legacy code path bit-identical when absent, so the untouched
% methods are unaffected.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 10 || isempty(sigma2)
    sigma2 = zeros(M,1);
end
if nargin < 11
    lld = [];
end
if nargin < 12
    covblk = {};
end
gaussian = any(sigma2 > 0);
hasLD = ~isempty(lld);

rates = x; % basic vector valid for INF and PS case min(ni,nservers(i))=ni
for i = 1:M
    if hasLD
        lldrow = lld(i,:);
        if all(lldrow == 1)
            lldrow = []; % station without load dependence keeps the plain branch
        end
    else
        lldrow = [];
    end
    switch sched_id(i) % source
        case SchedStrategy.INF
            % without load dependence each job is served at its own rate and
            % the share is the identity; alpha(n_i) scales the whole station
            if ~isempty(lldrow)
                idxIni = q_indices(i,1);
                idxEnd = q_indices(i,K) + Kic(i,K) - 1;
                ni = sum( x(idxIni:idxEnd) );
                if ni > 0
                    h = fluid_capacity_closure(ni, nservers(i), sigma2(i), lldrow, true);
                    rates(idxIni:idxEnd) = x(idxIni:idxEnd)/ni * h;
                end
            end
        case SchedStrategy.EXT  %EXT
            % this is treated by a delay except that we require mass
            % conservation in the local population
            for k=1:K
                idxIni = q_indices(i,k);
                idxEnd = q_indices(i,k) + Kic(i,k) - 1;
                if enabled(i,k)
                    rates(idxIni) = 1-sum(x(idxIni+1:idxEnd)); % keep total mass 1 into the source for all classes at all times, not needed for idxIni+1:idxEnd as rates is initialized equal to x
                end
            end
        case {SchedStrategy.PS, SchedStrategy.FCFS}
            idxIni = q_indices(i,1);
            idxEnd = q_indices(i,K) + Kic(i,K) - 1;
            blk = idxIni:idxEnd;
            ni = sum( x(blk) );
            if (gaussian || ~isempty(lldrow)) && ni > 0
                [h, dh] = fluid_capacity_closure(ni, nservers(i), sigma2(i), lldrow, false);
                Ci = [];
                if numel(covblk) >= i
                    Ci = covblk{i};
                end
                if isempty(Ci)
                    rates(blk) = x(blk)/ni * h;
                else
                    % THE SHARE AND THE CAPACITY ARE CLOSED JOINTLY. What the
                    % station clears is S_j*psi(N), and both factors move with
                    % N, so the product needs Cov(S_j,N)*psi'(n) on top of the
                    % two separate closures; see FLUID_SHARE_CLOSURE. With unit
                    % weights this is the DPS branch below.
                    [s, ~, cn] = fluid_share_closure(x(blk), ones(numel(blk),1), Ci);
                    rates(blk) = local_project_rate(s * h + dh * cn, x(blk), lldrow, h);
                end
            elseif ni > nservers(i) % case  min = ni handled by rates = x
                rates(blk) = x(blk)/ni * nservers(i);
            end
        case SchedStrategy.DPS %DPS
            % DPS is PS with a weighted share: the class-k coordinates get
            % w_k*x/ni of the station capacity psi(xi) instead of x/xi of it.
            % The denominator used to carry an ADDITIVE mean(w) term as a
            % divide-by-zero guard. That term never cancels, so the shares
            % summed to 1 - mean(w)/ni instead of 1 and utilization was
            % depressed by that factor (measured 24.9% on a two-class model
            % with weights [1 4]); the capacity was also taken as the full c
            % rather than psi(xi) = min(xi,c)*alpha(xi), so an underloaded
            % station was served at full rate. Both are now handled exactly as
            % in the PS/FCFS branch above, guarding with ni > 0 the way that
            % branch guards with xi > 0. Equal weights now reduce DPS to PS
            % identically, which they did not before.
            % The share itself is a RATIO, so evaluating it at the mean is a
            % separate closure from the min(): FLUID_SHARE_CLOSURE corrects it
            % at second order when a covariance block is supplied, and the
            % product of the two is closed jointly as in the PS branch.
            w(i,:) = w(i,:)/sum(w(i,:));

            lo = q_indices(i,1);
            hi = q_indices(i,K) + Kic(i,K) - 1;
            blk = lo:hi;
            wv = zeros(numel(blk),1);
            for k=1:K
                if enabled(i,k)
                    idx = (q_indices(i,k):(q_indices(i,k)+Kic(i,k)-1)) - lo + 1;
                    wv(idx) = w(i,k);
                end
            end
            xi = sum(x(blk));
            if xi > 0 && wv' * x(blk) > 0
                [psi, dpsi] = fluid_capacity_closure(xi, nservers(i), sigma2(i), lldrow, false);
                Ci = [];
                if numel(covblk) >= i
                    Ci = covblk{i};
                end
                [s, ~, cn] = fluid_share_closure(x(blk), wv, Ci);
                rates(blk) = local_project_rate(s * psi + dpsi * cn, x(blk), lldrow, psi);
            end
        case SchedStrategy.GPS
            % GPS splits the server by WEIGHT among the BACKLOGGED classes,
            % then equally among that class's own jobs. The share is a function
            % of the backlog indicator, so FLUID_GPS_SHARE closes it over the
            % 2^K patterns using P(X_k >= 1). No capacity term multiplies it:
            % GPS is single-server and the indicator already carries the idle
            % server, so the shares sum to 1 - P(station empty) by design.
            if nservers(i) > 1
                line_error(mfilename,'Multi-server GPS stations are not supported yet.');
            end
            [xk, vk, blk_k] = local_class_moments(x, i, K, enabled, q_indices, Kic, covblk);
            wk = w(i,:)';
            sk = fluid_gps_share(xk, wk, vk);
            a = 1;
            if ~isempty(lldrow)
                a = fluid_lld_scaling(lldrow, sum(x(q_indices(i,1):(q_indices(i,K)+Kic(i,K)-1))));
            end
            for k = 1:K
                if enabled(i,k) && xk(k) > 0
                    rates(blk_k{k}) = x(blk_k{k})/xk(k) * sk(k) * a;
                end
            end
    end
end
end

function [xk, vk, blk_k] = local_class_moments(x, i, K, enabled, q_indices, Kic, covblk)
% per-class population and variance at station i, plus the coordinate blocks
xk = zeros(K,1);
vk = zeros(K,1);
blk_k = cell(K,1);
lo = q_indices(i,1);
Ci = [];
if numel(covblk) >= i
    Ci = covblk{i};
end
for k = 1:K
    if ~enabled(i,k)
        continue
    end
    blk_k{k} = q_indices(i,k):(q_indices(i,k)+Kic(i,k)-1);
    xk(k) = sum(x(blk_k{k}));
    if ~isempty(Ci)
        loc = blk_k{k} - lo + 1;
        vk(k) = max(0, sum(sum(Ci(loc,loc))));
    end
end
end

function r = local_project_rate(r, xb, lldrow, tot)
% R = LOCAL_PROJECT_RATE(R, XB, LLDROW, TOT) project a jointly-closed
% per-coordinate service share onto the set it has to live in: R >= 0, R <= XB
% where that bound applies, and sum(R) = TOT.
%
% THE JOINT CLOSURE IS AN EXPANSION AND CAN LEAVE THAT SET. R = S*psi +
% psi'*Cov(S,N) adds a term that sums to ZERO over the coordinates, so it moves
% mass between them and its entries can push one past either bound; the
% first-order share X_j/n_i*psi cannot, being X_j scaled by psi/n_i <= 1.
% Either breach ends the same way, because SOLVER_FLUID_ITERATION integrates
% with NonNegative set on every coordinate: R_j > X_j drains coordinate j
% faster than it holds, the state goes negative, and the integrator CLAMPS it
% to zero -- which INJECTS mass. On cqn_twoclass_hyperl (Delay + PS, closed
% population 4) the total queue length came back 4.6719, with Queue1/Class2 at
% 3.4161 against the exact 2.7857, and on cqn_scheduling_dps (DPS) 3.1273
% against a population of 3. 'closing' and 'matrix', which carry no covariance
% and so no Cov(S,N) term, conserve exactly on both.
%
% THE UPPER BOUND HOLDS ONLY WITHOUT LOAD DEPENDENCE, which is what LLDROW
% selects. R is an expected NUMBER in service, so R_j <= X_j; but psi(n) =
% min(n,c)*alpha(n) folds the load-dependent scaling into the same variable, and
% with alpha > 1 the first-order share itself exceeds X_j. Capping there would
% move the load-dependent models that are already correct.
%
% Clip, then move the residual onto the coordinates that still have slack, in
% proportion to it, so sum(R) = TOT survives and the station still clears what
% its capacity closure says it clears. This is the rule FLUID_SHARE_CLOSURE
% applies to the shares one level up -- clip and renormalise the survivors,
% falling back to the first-order share when nothing survives -- and it is a
% NO-OP whenever the expansion stayed inside the set, which is why the methods
% and models already inside it are bit-identical.
%
% See also FLUID_SHARE_CLOSURE, FLUID_CAPACITY_CLOSURE, SOLVER_FLUID_ITERATION.

% Copyright (c) 2012-2026, QORE Lab, Imperial College London
% All rights reserved.

r = r(:);
xb = xb(:);
zt = GlobalConstants.Zero;
capped = isempty(lldrow);
if all(r >= -zt) && (~capped || all(r <= xb + zt))
    return
end
r = max(r, 0);
if capped
    r = min(r, xb);
end
for it = 1:(numel(r)+1)
    d = tot - sum(r);
    if abs(d) <= zt
        break
    end
    if d > 0
        if capped
            slack = xb - r;   % room to grow
        else
            slack = ones(numel(r),1);
        end
    else
        slack = r;            % room to shrink
    end
    tsl = sum(slack);
    if tsl <= zt
        break
    end
    r = max(r + d*slack/tsl, 0);
    if capped
        r = min(r, xb);
    end
end
if abs(sum(r) - tot) > zt
    % nothing had slack: fall back to the first-order share, which satisfies
    % every bound by construction
    sx = sum(xb);
    if sx > zt
        r = xb * (tot / sx);
    end
end
end
