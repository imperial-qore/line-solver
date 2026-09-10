function [Qtable, Utable] = fes_compute_metrics(L, mi, isDelay, cutoffs, options)
% [QTABLE, UTABLE] = FES_COMPUTE_METRICS(L, MI, ISDELAY, CUTOFFS, OPTIONS)
%
% Per-station metrics of the ISOLATED subnetwork at every population state,
% the companion of FES_COMPUTE_THROUGHPUTS.
%
% FES_COMPUTE_THROUGHPUTS returns the aggregate throughput X(n) that becomes
% the flow-equivalent server's rate. That is all the REDUCED model needs, but
% it is not enough to report the collapsed stations' own metrics: those are
% recovered by conditioning on the FES population,
%
%   E[Q_i] = sum_n P(N_fes = n) * Q_i(n),
%
% which is the Chandy-Herzog-Woo hierarchical decomposition and is EXACT when
% the subnetwork is product-form. This function supplies the Q_i(n) and
% U_i(n) that the sum is taken over.
%
% Input:
%   L        - (M_sub x K) service demands of the isolated subnetwork
%   mi       - (1 x M_sub) servers per station, Inf at a Delay
%   isDelay  - (1 x M_sub) logical, true at a Delay
%   cutoffs  - (1 x K) per-class population cutoffs
%   options  - struct with optional field .verbose
%
% Output (each a cell array indexed by LJD_LINEARIZE, as
% FES_COMPUTE_THROUGHPUTS indexes its scaling table):
%   Qtable{idx} - (M_sub x K) queue lengths at that population
%   Utable{idx} - (M_sub x K) utilizations
%
% Throughput needs no table: flow through a station is fixed by the ROUTING,
% so the caller derives it from the original model's visit ratios and the
% chain throughput the reduced solve already reports.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 5
    options = struct();
end
if ~isfield(options, 'verbose')
    options.verbose = false;
end

[M_sub, K] = size(L);
queueIdx = find(~isDelay);
delayIdx = find(isDelay);
M_queue = length(queueIdx);

if M_queue > 0
    L_queue = L(queueIdx, :);
    mi_queue = mi(queueIdx);
else
    L_queue = zeros(0, K);
    mi_queue = [];
end
if ~isempty(delayIdx)
    Z = sum(L(delayIdx, :), 1);
else
    Z = zeros(1, K);
end

tableSize = prod(cutoffs + 1);
Qtable = cell(1, tableSize);
Utable = cell(1, tableSize);

for idx = 1:tableSize
    Qtable{idx} = zeros(M_sub, K);
    Utable{idx} = zeros(M_sub, K);

    nvec = ljd_delinearize(idx, cutoffs);
    if sum(nvec) == 0
        continue % an empty subnetwork holds nothing and serves nothing
    end

    % MI is the additive C=L*(mi+Qarv) term, not a server count: multiservers
    % go through PFQN_MVAMS, as in FES_COMPUTE_THROUGHPUTS
    [XN, QNq] = pfqn_mvams(zeros(1,K), L_queue, nvec, Z, ...
        ones(M_queue,1), mi_queue(:));
    Q = zeros(M_sub, K); U = zeros(M_sub, K);
    for a = 1:M_queue
        i = queueIdx(a);
        Q(i,:) = QNq(a,:);
        % PFQN_MVAMS reports UN per STATION on its closed multiserver branch and
        % per station-class elsewhere, so it is not read here: utilization is
        % recomputed analytically as U=X*L/S, the [0,1] convention LINE uses at
        % every queueing station whatever its multiplicity (SOLVER_MVA does the
        % same with the same reason).
        U(i,:) = XN(:)' .* L_queue(a,:) / max(1, mi_queue(a));
    end
    % A Delay holds X_k * Z_i(k) jobs by Little's law on a station with no
    % queueing, and its INF "utilization" is that same population.
    for a = 1:length(delayIdx)
        i = delayIdx(a);
        Q(i,:) = XN(:)' .* L(i,:);
        U(i,:) = Q(i,:);
    end
    Qtable{idx} = Q; Utable{idx} = U;
end

if options.verbose
    line_printf(sprintf('Computed FES conditional metrics for %d population states\n', tableSize));
end
end
