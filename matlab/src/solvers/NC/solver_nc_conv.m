function [Q,U,R,T,C,X,lG,runtime,iter,method] = solver_nc_conv(sn, options)
% [Q,U,R,T,C,X,LG,RUNTIME,ITER,METHOD] = SOLVER_NC_CONV(SN, OPTIONS)
%
% Exact normalizing constant solver for closed networks with Limited
% class-dependent (cdscaling) service rates, using the multichain
% convolution algorithm of Sauer (1983), Section 5.2.
%
% This solver handles models where some stations have class-dependent scaling
% (e.g., Flow-Equivalent Servers from aggregateFES).

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

Tstart = tic;
method = 'conv';
iter = 1;

M = sn.nstations;
K = sn.nclasses;
NK = sn.njobs';
nservers = sn.nservers;

V = cellsum(sn.visits);
ST = 1 ./ sn.rates;
ST(isnan(ST)) = 0;

% Demands: L(ist,k) = V(ist,k) * ST(ist,k)
Ldemand = V .* ST;

% Separate delay and queue stations
isDelay = isinf(nservers);
delayIdx = find(isDelay);
queueIdx = find(~isDelay);
nQueues = length(queueIdx);

% Build pfqn_conv inputs
% Z: total delay demand per class
Z_conv = zeros(1, K);
for ist = delayIdx(:)'
    Z_conv = Z_conv + Ldemand(ist, :);
end

% L_conv: demands for queue stations only
L_conv = Ldemand(queueIdx, :);

% Class-dependence handles beta_{i,r}(n) for the queue stations. Each handle
% takes the per-class population vector at its station and returns either a
% scalar (chain-independent) or a length-R vector of per-class rates.
cdscaling_conv = cell(nQueues, 1);
if ~isempty(sn.cdscaling)
    for qi = 1:nQueues
        ist = queueIdx(qi);
        if ist <= length(sn.cdscaling) && ~isempty(sn.cdscaling{ist})
            cdscaling_conv{qi} = sn.cdscaling{ist};
        end
    end
end
% Fold joint-dependence handles (sn.jdscaling, non-product-form eta_i) into the
% same per-station handle used by the convolution recursion. cd and jd are
% evaluated identically; the product reproduces the single-mechanism case when
% only one is present (the other is treated as absent below).
if ~isempty(sn.jdscaling)
    for qi = 1:nQueues
        ist = queueIdx(qi);
        if ist <= length(sn.jdscaling) && ~isempty(sn.jdscaling{ist})
            jdh = sn.jdscaling{ist};
            if isempty(cdscaling_conv{qi})
                cdscaling_conv{qi} = jdh;
            else
                cdh = cdscaling_conv{qi};
                cdscaling_conv{qi} = @(ni) cdh(ni) .* jdh(ni);
            end
        end
    end
end

%% Compute G(N)
[G_N, lG] = pfqn_conv(L_conv, NK, Z_conv, cdscaling_conv);

%% Compute G(N - e_k) for each class -> throughput
XN = zeros(1, K);
for k = 1:K
    if NK(k) > 0
        NK_minus = NK;
        NK_minus(k) = NK_minus(k) - 1;
        [G_Nk, ~] = pfqn_conv(L_conv, NK_minus, Z_conv, cdscaling_conv);
        XN(k) = G_Nk / G_N;
    end
end

%% Compute per-station throughput
TN = V .* repmat(XN, M, 1);

%% Compute queue lengths
QN = zeros(M, K);

% Delay stations: Q = L * X
for ist = delayIdx(:)'
    for k = 1:K
        QN(ist, k) = Ldemand(ist, k) * XN(k);
    end
end

% Queue stations: use marginal distribution
% P_m(n|N) = X_m(n) * G_{-m}(N-n) / G(N)
% Q_m_k = sum_{n: n_k>=1} n_k * P_m(n|N)
%
% G_{-m}(N-n) is computed by pfqn_conv on all stations except m
stateSpaceSize = prod(NK + 1);

for qi = 1:nQueues
    ist = queueIdx(qi);

    % Build X_m(n) for this station
    Xm = zeros(stateSpaceSize, 1);
    Xm(1) = 1; % X_m(0) = 1
    isCdStation = ~isempty(cdscaling_conv{qi});

    n = pprod_init(NK);
    while n(1) >= 0
        idx = hashpop(n, NK);
        if sum(n) > 0
            if isCdStation
                % class-dependent: X_m(n) = (L/mu_km(n)) * X_m(n-e_k) via eq. (40)
                % Pick any k with n_k > 0 (result is path-independent)
                for r = 1:K
                    if n(r) > 0
                        % beta_{qi,r}(n): DIMENSIONLESS scaling of the demand,
                        % so the effective demand is L/beta. Handle returns a
                        % scalar (shared by all classes) or a length-R vector.
                        bval = cdscaling_conv{qi}(n);
                        if numel(bval) > 1
                            beta = bval(r);
                        else
                            beta = bval;
                        end
                        % X_m(n) = (|n|/n_r) * (L/beta) * X_m(n-e_r); at beta=1
                        % this is the load-independent multinomial recurrence.
                        tot = sum(n);
                        nr = n(r);
                        n(r) = n(r) - 1;
                        idx_prev = hashpop(n, NK);
                        n(r) = n(r) + 1;
                        if beta > 0
                            Xm(idx) = (tot / nr) * (L_conv(qi, r) / beta) * Xm(idx_prev);
                        end
                        break
                    end
                end
            else
                % LI: X_m(n) = Σ_r L(m,r) * X_m(n-e_r) (multinomial form)
                for r = 1:K
                    if n(r) > 0
                        n(r) = n(r) - 1;
                        idx_prev = hashpop(n, NK);
                        n(r) = n(r) + 1;
                        Xm(idx) = Xm(idx) + L_conv(qi, r) * Xm(idx_prev);
                    end
                end
            end
        end
        n = pprod_next(n, NK);
    end

    % Build complement: all stations except qi
    L_comp = L_conv; L_comp(qi, :) = [];
    cd_comp = cdscaling_conv; cd_comp(qi) = [];

    % Compute Q_m_k using marginal
    n = pprod_init(NK);
    while n(1) >= 0
        if any(n > 0)
            idx = hashpop(n, NK);
            nmi = NK - n;
            if all(nmi >= 0)
                % G_{-m}(N-n) with delay
                [G_comp, ~] = pfqn_conv(L_comp, nmi, Z_conv, cd_comp);
                prob = Xm(idx) * G_comp / G_N;
                for k = 1:K
                    QN(ist, k) = QN(ist, k) + n(k) * prob;
                end
            end
        end
        n = pprod_next(n, NK);
    end
end

%% Compute remaining metrics
RN = QN ./ TN;
RN(TN == 0) = 0;
UN = TN .* ST;

% Utilization at a class-dependent station is normalized by the peak service
% capacity (sn.cdscalingpeak), not T*ST, so U<=1 by construction; see
% _kb/06-solver-catalog.md (NC section, convolution/beta scaling)
for qi = 1:nQueues
    ist = queueIdx(qi);
    if isempty(cdscaling_conv{qi})
        continue
    end
    for r = 1:K
        % Effective peak = product of the class- and joint-dependence peaks
        % declared at the station (a missing one contributes 1).
        bmax = 1;
        haspeak = false;
        if ~isempty(sn.cdscaling) && ist <= length(sn.cdscaling) && ~isempty(sn.cdscaling{ist})
            bmax = bmax * sn.cdscalingpeak(ist,r); haspeak = true;
        end
        if ~isempty(sn.jdscaling) && ist <= length(sn.jdscaling) && ~isempty(sn.jdscaling{ist})
            bmax = bmax * sn.jdscalingpeak(ist,r); haspeak = true;
        end
        if haspeak && bmax > 0
            UN(ist,r) = UN(ist,r) / bmax;
        end
    end
end
CN = NK ./ XN;
CN(XN == 0) = 0;
CN = CN - sum(Z_conv .* (repmat(1, M, 1) .* isDelay(:)), 1); % subtract delay

% Output
Q = QN;
U = UN;
R = RN;
T = TN;
C = CN;
X = XN;
runtime = toc(Tstart);
end

%% --- Local helper functions ---


function idx = hashpop(n, N)
idx = 1;
R = length(N);
for r = 1:R
    idx = idx + prod(N(1:r-1) + 1) * n(r);
end
end

function n = pprod_init(N)
n = zeros(size(N));
end

function n = pprod_next(n, N)
R = length(N);
if all(n == N)
    n = -1 * ones(1, R);
    return
end
s = R;
while s > 0 && n(s) == N(s)
    n(s) = 0;
    s = s - 1;
end
if s > 0
    n(s) = n(s) + 1;
end
end
