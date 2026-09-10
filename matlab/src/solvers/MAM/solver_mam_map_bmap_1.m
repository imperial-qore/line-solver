function [QN, UN, RN, TN, piAgg, R] = solver_mam_map_bmap_1(mapArr, bmapSvc, varargin)
% SOLVER_MAM_MAP_BMAP_1 Solve MAP/BMAP/1 queue using GI/M/1-type analysis.
%
% Solves a MAP/BMAP/1 queue where:
%   - Arrivals follow a Markovian Arrival Process (MAP)
%   - Service follows a Batch Markovian Arrival Process (BMAP) for
%     batch service completions
%   - Single server
%
% The queue is modeled as a GI/M/1-type Markov chain because:
%   - MAP arrivals increase level by exactly 1
%   - BMAP service can decrease level by 1, 2, 3, ... (batch sizes)
%
% USAGE:
%   [QN, UN, RN, TN] = solver_mam_map_bmap_1(mapArr, bmapSvc)
%   [QN, UN, RN, TN, piAgg, R] = solver_mam_map_bmap_1(mapArr, bmapSvc)
%
% INPUT:
%   mapArr  - MAP object or cell array {C0, C1} for arrivals
%   bmapSvc - BMAP object or cell array {D0, D1, D2, ...} for batch service
%
% OUTPUT:
%   QN    - Mean queue length E[N]
%   UN    - Server utilization rho
%   RN    - Mean response time E[R]
%   TN    - Throughput (arrival rate)
%   piAgg - Aggregated stationary probabilities [pi0, pi1, piStar]
%   R     - R matrix from GI/M/1-type solution
%
% The GI/M/1-type structure for MAP/BMAP/1 is:
%       B1  A0  0   0   ...
%       B2  A1  A0  0   ...
%   Q = B3  A2  A1  A0  ...
%       ...
%
% Where:
%   A0 = C1 ⊗ I_ms           (MAP arrival, level +1)
%   A1 = C0 ⊗ I_ms + I_ma ⊗ D0  (phase changes, level 0)
%   A_{k+1} = I_ma ⊗ D_k     (batch size k service, level -k)
%
% REFERENCES:
%   Riska, A., & Smirni, E. (2003). ETAQA: An Efficient Technique for the
%   Analysis of QBD-Processes by Aggregation. Performance Evaluation.
%
% See also: solver_mam_bmap_m_1, MAP, BMAP

%% Parse inputs
if isa(mapArr, 'MAP')
    C0 = mapArr.getD0();
    C1 = mapArr.getD1();
elseif iscell(mapArr)
    C0 = mapArr{1};
    C1 = mapArr{2};
else
    error('mapArr must be a MAP object or cell array {C0, C1}');
end

if isa(bmapSvc, 'BMAP')
    K = bmapSvc.getNumberOfPhases() - 1;  % Max batch size
    D = cell(1, K + 1);
    D{1} = bmapSvc.getD0();
    for k = 1:K
        D{k+1} = bmapSvc.D(1, k);  % D_k for batch size k
    end
elseif iscell(bmapSvc)
    D = bmapSvc;
    K = length(D) - 1;  % Max batch size
else
    error('bmapSvc must be a BMAP object or cell array {D0, D1, D2, ...}');
end

%% Get dimensions
ma = size(C0, 1);  % Number of MAP phases
ms = size(D{1}, 1);  % Number of BMAP phases
m = ma * ms;  % Combined phases per level

%% Validate inputs
% Check MAP is valid generator
rowSumsC = sum(C0 + C1, 2);
if max(abs(rowSumsC)) > 1e-10
    error('MAP matrices C0 + C1 must have zero row sums');
end

% Check BMAP is valid generator
DTotal = D{1};
for k = 1:K
    DTotal = DTotal + D{k+1};
end
rowSumsD = sum(DTotal, 2);
if max(abs(rowSumsD)) > 1e-10
    error('BMAP matrices D0 + D1 + ... + DK must have zero row sums');
end

%% Compute arrival and service rates
% Arrival rate from MAP
piC = map_prob({C0, C1});
lambdaArr = piC * C1 * ones(ma, 1);

% Service rate from BMAP (total customers served per unit time)
piD = map_prob({D{1}, DTotal - D{1}});  % D0, sum of D1...DK
muTotal = 0;
for k = 1:K
    muTotal = muTotal + k * (piD * D{k+1} * ones(ms, 1));
end

% Utilization
rho = lambdaArr / muTotal;

if rho >= 1.0
    warning('System is unstable (rho = %.4f >= 1). Results may be invalid.', rho);
end

%% Construct A matrix for GI/M/1-type (repeating part)
% A = [A0; A1; A2; ...; A_{K+1}] with (K+2) blocks of size m x m
%
% A0 = C1 ⊗ I_ms           (level +1: arrival)
% A1 = C0 ⊗ I_ms + I_ma ⊗ D0  (level 0: phase changes)
% A_{k+1} = I_ma ⊗ D_k     (level -k: batch service of k)

A = zeros(m * (K + 2), m);

% A0: MAP arrival (level +1)
A0 = kron(C1, eye(ms));
A(1:m, :) = A0;

% A1: Phase changes only (level 0)
A1 = kron(C0, eye(ms)) + kron(eye(ma), D{1});
A(m+1:2*m, :) = A1;

% A_{k+1}: Batch service of size k (level -k)
for k = 1:K
    Ak = kron(eye(ma), D{k+1});
    A((k+1)*m+1:(k+2)*m, :) = Ak;
end

%% Construct B matrix for boundary levels
% For levels 0, 1, ..., K-1, we need special handling because
% batch service of size k is only possible when queue length >= k
%
% B = [B1; B2; ...; B_{K+1}] where B_i is m x m
%
% The boundary structure accounts for:
% - Level 0: no service possible (empty queue)
% - Level 1: only batch size 1 service possible
% - Level k: batch sizes 1, 2, ..., k possible

B = zeros(m * (K + 1) + m, m);  % B has mb=m boundary rows, then K+1 blocks

% B1 (level 0 -> level 0): Phase changes only, no service from empty
% Arrivals from level 0 are handled by B0 = A0
B1 = kron(C0, eye(ms)) + kron(eye(ma), D{1});
% Add service transitions that would go to negative levels back to level 0
for k = 1:K
    B1 = B1 + kron(eye(ma), D{k+1});  % Service from empty stays at 0
end
B(1:m, :) = B1;

% B2, B3, ..., B_{K+1}: Transitions to level 0 from levels 1, 2, ..., K
for j = 1:K
    % From level j, we can service batches of size 1, 2, ..., j going to level 0
    Bj = zeros(m, m);
    % Only batch sizes <= j are valid; larger batches get truncated
    for k = j:K
        % Batch size k from level j: if k > j, customer stays at 0
        % But for GI/M/1 B matrix, we accumulate transitions to level 0
        % B_{j+1} handles transition from level j to level 0
        Bj = Bj + kron(eye(ma), D{k+1});
    end
    B(m + (j-1)*m + 1 : m + j*m, :) = Bj;
end

%% Compute R matrix using GI/M/1 ETAQA
R = GIM1_R_ETAQA(A);

%% Compute stationary probabilities using GI/M/1 ETAQA
piAgg = GIM1_pi_ETAQA(B, A, R, 'Boundary', A0);

%% Compute performance metrics
% Extract probability components
pi0 = piAgg(1:m);
pi1 = piAgg(m+1:2*m);
piStar = piAgg(2*m+1:3*m);

% Mean queue length using ETAQA moments
QN = GIM1_qlen_ETAQA(B, A, R, piAgg, 1, 'Boundary', A0);

% Throughput equals arrival rate (in steady state)
TN = lambdaArr;

% Utilization
UN = rho;

% Mean response time via Little's Law
RN = QN / TN;

end
