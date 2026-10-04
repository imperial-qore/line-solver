function [A, G] = fluid_drift_jacobian(x, M, K, enabled, q_indices, Kic, nservers, w, sched_id, rateBase, eventIdx, D, sigma2, lld, covblk)
% [A, G] = FLUID_DRIFT_JACOBIAN(X, M, K, ENABLED, Q_INDICES, KIC, NSERVERS, W, SCHED_ID, RATEBASE, EVENTIDX, D, SIGMA2, LLD, COVBLK)
%
% Analytic Jacobian of the fluid drift F(x) = D*(rateBase .* g(x)(eventIdx)),
% where g is the rate-factor vector built by ODE_RATES_CLOSING. The Jacobian
% drives the covariance (Lyapunov) equation of the linear noise approximation
% and the 1/N refinement, so it must mirror ODE_RATES_CLOSING branch by
% branch: any policy without a case there keeps g = x and contributes the
% identity here.
%
% With sigma2 = 0 the derivative of the occupancy factor is the indicator of
% the unsaturated region, i.e. the a.e. derivative of the first-order
% closure. With sigma2 > 0 it is the smooth derivative of the Gaussian
% closure returned by FLUID_MIN_CLOSURE.
%
% Parameters:
%   x         - state vector (phase-resolved populations)
%   D         - (n x nevents) jump matrix from ODE_JUMPS_NEW
%   sigma2    - (M x 1) station population variances, 0 for the first-order closure
%   remaining arguments as in ODE_RATES_CLOSING
%
% Returns:
%   A - (n x n) Jacobian dF/dx
%   G - (n x n) Jacobian dg/dx of the rate factors
%
% See also ODE_RATES_CLOSING, FLUID_MIN_CLOSURE, FLUID_LYAPUNOV.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 13 || isempty(sigma2)
    sigma2 = zeros(M,1);
end
if nargin < 14
    lld = [];
end
if nargin < 15
    covblk = {};
end
gaussian = any(sigma2 > 0);

n = numel(x);
G = eye(n); % INF, EXT phases 2..end, and every policy without a case below

for i = 1:M
    if ~isempty(lld)
        lldrow = lld(i,:);
        if all(lldrow == 1)
            lldrow = [];
        end
    else
        lldrow = [];
    end
    switch sched_id(i)
        case SchedStrategy.INF
            % g = x, identity already in place unless a load dependence
            % scales the whole station
            if ~isempty(lldrow)
                blk = q_indices(i,1):(q_indices(i,K) + Kic(i,K) - 1);
                ni = sum(x(blk));
                if ni > 0
                    [h, dh] = fluid_capacity_closure(ni, nservers(i), sigma2(i), lldrow, true);
                    f = h/ni;
                    fp = (ni*dh - h)/ni^2;
                    G(blk,:) = 0;
                    G(blk,blk) = f*eye(numel(blk)) + x(blk)*fp*ones(1,numel(blk));
                end
            end
        case SchedStrategy.EXT
            for k = 1:K
                if enabled(i,k)
                    idxIni = q_indices(i,k);
                    idxEnd = q_indices(i,k) + Kic(i,k) - 1;
                    G(idxIni,:) = 0;
                    G(idxIni, idxIni+1:idxEnd) = -1;
                end
            end
        case {SchedStrategy.PS, SchedStrategy.FCFS}
            idxIni = q_indices(i,1);
            idxEnd = q_indices(i,K) + Kic(i,K) - 1;
            blk = idxIni:idxEnd;
            ni = sum(x(blk));
            if ni <= 0
                continue % g = x on an empty station
            end
            if gaussian || ~isempty(lldrow)
                [h, dh] = fluid_capacity_closure(ni, nservers(i), sigma2(i), lldrow, false);
                % the share closure is reached only from this branch, exactly
                % as in ODE_RATES_CLOSING_FACTORS, so the two never disagree
                Ci = [];
                if numel(covblk) >= i
                    Ci = covblk{i};
                end
                if ~isempty(Ci)
                    % g = s(x_blk)*h(ni) + h'(ni)*cn(x_blk), the joint closure
                    % of the share and the capacity; differentiating it with C
                    % held fixed adds h'*dcn and h''*cn to the product rule
                    [~, ~, d2h] = fluid_capacity_closure(ni, nservers(i), sigma2(i), lldrow, false);
                    [s, ds, cn, dcn] = fluid_share_closure(x(blk), ones(numel(blk),1), Ci);
                    e = ones(1,numel(blk));
                    G(blk,:) = 0;
                    G(blk,blk) = ds*h + s*(dh*e) + dh*dcn + cn*(d2h*e);
                    continue
                end
            elseif ni > nservers(i) - GlobalConstants.FineTol*max(1, ni)
                % THE SATURATION TEST CARRIES A BAND, and it is a cross-codebase
                % requirement: a saturated fixed point sits exactly at ni = c, and
                % each engine's ODE stops on its own residual (MATLAB 1.0004, the
                % C++ port 1 - 1.8e-13 on the same model). A strict ni > c reads
                % saturated in one and unsaturated in the other, which flips this
                % whole station block between a zero row and the identity, and with
                % it the hyperbolicity verdict FLUID_LYAPUNOV returns and the method
                % SolverFLD ends up answering with. See FLUID_MIN_CLOSURE, whose
                % degenerate branch carries the same band.
                h = nservers(i); dh = 0;
            else
                continue % g = x, identity already in place
            end
            % g_j = x_j*h(ni)/ni  ->  dg_j/dx_m = delta_jm*f + x_j*f',
            % f = h/ni, f' = (ni*dh - h)/ni^2
            f = h/ni;
            fp = (ni*dh - h)/ni^2;
            G(blk,:) = 0;
            G(blk,blk) = f*eye(numel(blk)) + x(blk)*fp*ones(1,numel(blk));
        case SchedStrategy.DPS
            % g = s(x_blk)*psi(xi) + psi'(xi)*cn(x_blk), the joint closure of
            % the capacity share of FLUID_SHARE_CLOSURE with the capacity,
            % mirroring the DPS branch of ODE_RATES_CLOSING_FACTORS:
            %   dg_j/dx_m = ds_j/dx_m*psi + s_j*dpsi + dpsi*dcn_j/dx_m
            %               + cn_j*d2psi  (m inside the station)
            wi = w(i,:)/sum(w(i,:));
            stblk = q_indices(i,1):(q_indices(i,K) + Kic(i,K) - 1);
            wv = zeros(numel(stblk),1); % per-coordinate DPS weight
            for k = 1:K
                if enabled(i,k)
                    idx = (q_indices(i,k):(q_indices(i,k)+Kic(i,k)-1)) - stblk(1) + 1;
                    wv(idx) = wi(k);
                end
            end
            xi = sum(x(stblk));
            if xi <= 0 || wv' * x(stblk) <= 0
                continue % g = x on an empty station
            end
            [psi, dpsi, d2psi] = fluid_capacity_closure(xi, nservers(i), sigma2(i), lldrow, false);
            Ci = [];
            if numel(covblk) >= i
                Ci = covblk{i};
            end
            [s, ds, cn, dcn] = fluid_share_closure(x(stblk), wv, Ci);
            e = ones(1,numel(stblk));
            G(stblk,:) = 0;
            G(stblk,stblk) = ds*psi + s*(dpsi*e) + dpsi*dcn + cn*(d2psi*e);
        case SchedStrategy.GPS
            % g_j = (x_j/x_k)*s_k(x_1..x_K)*a  for coordinate j of class k, so
            %   dg_j/dx_l = [delta_jl/x_k - x_j/x_k^2]*s_k*a        (l in class k)
            %             + (x_j/x_k)*ds_k/dx_m*a                   (l in class m)
            %             + (x_j/x_k)*s_k*da/dxi                    (l in station)
            if nservers(i) > 1
                line_error(mfilename,'Multi-server GPS stations are not supported yet.');
            end
            stblk = q_indices(i,1):(q_indices(i,K) + Kic(i,K) - 1);
            xk = zeros(K,1); vk = zeros(K,1); blk_k = cell(K,1);
            Ci = [];
            if numel(covblk) >= i
                Ci = covblk{i};
            end
            for k = 1:K
                if enabled(i,k)
                    blk_k{k} = q_indices(i,k):(q_indices(i,k)+Kic(i,k)-1);
                    xk(k) = sum(x(blk_k{k}));
                    if ~isempty(Ci)
                        loc = blk_k{k} - stblk(1) + 1;
                        vk(k) = max(0, sum(sum(Ci(loc,loc))));
                    end
                end
            end
            [sk, dsk] = fluid_gps_share(xk, w(i,:)', vk);
            xi = sum(x(stblk));
            a = 1; da = 0;
            if ~isempty(lldrow)
                [a, da] = fluid_lld_scaling(lldrow, xi);
            end
            G(stblk,:) = 0;
            for k = 1:K
                if ~enabled(i,k) || xk(k) <= 0
                    continue
                end
                bk = blk_k{k};
                G(bk,bk) = G(bk,bk) + (sk(k)*a/xk(k))*eye(numel(bk)) ...
                    - (sk(k)*a/xk(k)^2)*x(bk)*ones(1,numel(bk));
                for m = 1:K
                    if enabled(i,m)
                        G(bk,blk_k{m}) = G(bk,blk_k{m}) ...
                            + (a*dsk(k,m)/xk(k))*x(bk)*ones(1,numel(blk_k{m}));
                    end
                end
                if da ~= 0
                    G(bk,stblk) = G(bk,stblk) + (sk(k)*da/xk(k))*x(bk)*ones(1,numel(stblk));
                end
            end
    end
end

A = D * (rateBase .* G(eventIdx,:));
end
