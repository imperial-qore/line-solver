function result = qsys_tandem_lindley(A, S, W0)
% QSYS_TANDEM_LINDLEY Tandem network Lindley recursion on a sample path.
%
% RESULT = QSYS_TANDEM_LINDLEY(A, S) propagates the waiting times of a series of
% K single-server FCFS stations in tandem, driven by the primitives of the
% sample path. A is the length-N vector of interarrival times at the first
% station, A(n) separating customers n and n+1, and S is the N-by-K matrix of
% service times, S(n,k) being the service time of customer n at station k. All
% waiting times start from zero.
%
% RESULT = QSYS_TANDEM_LINDLEY(A, S, W0) starts customer 1 from the given
% length-K vector of waiting times instead of from an empty network.
%
% At the first station this is Lindley's recursion,
%   W(n+1,1) = max(W(n,1) + S(n,1) - A(n), 0).
% Downstream the interarrival time is not a primitive: the arrival epoch of
% customer n at station k is its departure epoch from station k-1, so the
% interarrival time at station k is the interdeparture time upstream. Writing
% G(n,k) for the interarrival time at station k between customers n and n+1,
% with G(n,1) = A(n), the exact interdeparture identity is
%   G(n,k+1) = G(n,k) + W(n+1,k) - W(n,k) + S(n+1,k) - S(n,k),
% equivalently and more transparently
%   G(n,k+1) = max(G(n,k) - W(n,k) - S(n,k), 0) + S(n+1,k),
% an idle period at station k followed by the next customer's service there.
% The recursion at station k is then
%   W(n+1,k) = max(W(n,k) + S(n,k) - G(n,k), 0).
%
% Nothing here is distributional, so the recursion is exact for arbitrary
% interarrival and service times, dependent or not, and is the reference a
% simulated tandem sample path can be checked against directly. It reproduces a
% direct event-driven tandem simulation to 1e-12 over four stations.
%
% Note that proposition 1 of the reference states this identity without the
% S(n+1,k) - S(n,k) term, which makes it wrong as a sample-path identity: the
% omitted difference has mean zero, so the mean interdeparture time survives,
% but individual waiting times do not. Implementing it as published gives
% station-1 waiting times that are correct and downstream ones that are not,
% by up to several mean service times. The form above is used instead.
%
% Returns a struct with fields:
%   W        - N-by-K waiting times, W(n,k) for customer n at station k
%   G        - N-by-K interarrival times, G(n,k) between customers n and n+1 at
%              station k, so G(:,1) is A; the last row is NaN, there being no
%              customer N+1 to separate from
%   T        - N-by-K sojourn times, W + S
%   departure- N-by-K departure epochs of each customer from each station
%   analyzer - Identifier string
%
% Examples:
%   A = exprnd(1/0.8, 1000, 1);
%   S = exprnd(1, 1000, 2);
%   r = qsys_tandem_lindley(A, S);
%   mean(r.W(:,1))    % upstream mean wait
%   mean(r.W(:,2))    % downstream mean wait
%
% Reference: S. Palomo, J. Pender, "Learning the Tandem Network Lindley
% Recursion", Proc. Winter Simulation Conference, 2021, equations 2 and 3 and
% proposition 1, the last corrected as described above; D. V. Lindley, "The
% Theory of Queues with a Single Server", Proc. Camb. Phil. Soc. 48, 1952.
%
% See also QSYS_MM1_TANDEM_LINDLEY, QSYS_MM1_LINDLEY, QSYS_HH1_LINDLEY
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

A = A(:);
N = numel(A);
if N < 1
    line_error(mfilename, 'A must hold at least one interarrival time');
end
if size(S, 1) ~= N
    line_error(mfilename, ...
        'S must have one row per customer, got %d rows for %d interarrivals', ...
        size(S, 1), N);
end
K = size(S, 2);
if K < 1
    line_error(mfilename, 'S must have at least one column, one per station');
end
if ~isreal(A) || any(A < 0) || any(~isfinite(A))
    line_error(mfilename, 'A must hold finite nonnegative interarrival times');
end
if ~isreal(S) || any(S(:) < 0) || any(~isfinite(S(:)))
    line_error(mfilename, 'S must hold finite nonnegative service times');
end

if nargin < 3 || isempty(W0)
    W0 = zeros(1, K);
end
W0 = W0(:)';
if numel(W0) ~= K
    line_error(mfilename, 'W0 must hold one waiting time per station, %d of them', K);
end
if ~isreal(W0) || any(W0 < 0) || any(~isfinite(W0))
    line_error(mfilename, 'W0 must hold finite nonnegative waiting times');
end

W = zeros(N, K);
G = nan(N, K);
W(1, :) = W0;

for n = 1:N - 1
    gap = A(n);
    for k = 1:K
        G(n, k) = gap;
        W(n + 1, k) = max(W(n, k) + S(n, k) - gap, 0);
        % the interdeparture time here is the interarrival time one station on
        gap = max(gap - W(n, k) - S(n, k), 0) + S(n + 1, k);
    end
end

T = W + S;
departure = zeros(N, K);
arrivalEpoch = cumsum([0; A(1:N - 1)]);
departure(:, 1) = arrivalEpoch + T(:, 1);
for k = 2:K
    departure(:, k) = departure(:, k - 1) + T(:, k);
end

result = struct('W', W, 'G', G, 'T', T, 'departure', departure, ...
    'analyzer', 'qsys_tandem_lindley');
end
