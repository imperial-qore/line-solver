%{ @file ldqbd_pi.m
 %  @brief Computes stationary distribution for level-dependent QBD processes
 %
 %  @author QORE Lab, Imperial College London (Algorithm by Phung-Duc et al.)
%}

%{
 % @brief Computes the stationary distribution for a level-dependent QBD
 %
 % @details
 % This function implements Algorithm 3 from the Phung-Duc et al. paper for
 % computing the stationary distribution of level-dependent QBD processes.
 %
 % The algorithm uses:
 % 1. Boundary equation: pi_0 * (Q1^(0) + R^(1)*Q2^(1)) = 0
 % 2. Forward recursion: pi_n = pi_{n-1} * R^(n)
 % 3. Normalize: sum(pi) = 1
 %
 % @par Syntax:
 % @code
 % pi = ldqbd_pi(R, Q0, Q1, Q2)
 % pi = ldqbd_pi(R, Q0, Q1, Q2, options)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>R<td>Cell array {R^(1), R^(2), ..., R^(N)} - rate matrices from ldqbd_R
 % <tr><td>Q0<td>Cell array {Q0^(0), Q0^(1), ..., Q0^(N-1)} - upward transitions
 % <tr><td>Q1<td>Cell array {Q1^(0), Q1^(1), ..., Q1^(N)} - local transitions
 % <tr><td>Q2<td>Cell array {Q2^(1), Q2^(2), ..., Q2^(N)} - downward transitions
 % <tr><td>options<td>Optional struct with fields: verbose (default false)
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>pi<td>Row vector [pi_0, pi_1, ..., pi_N] - stationary distribution
 %              (aggregated over phases at each level)
 % <tr><td>pi_cell<td>(Optional) Cell array with full phase-level distributions
 % </table>
 %
 % @par Reference:
 % T. Phung-Duc, H. Masuyama, S. Kasahara, Y. Takahashi, "A Simple Algorithm
 % for the Rate Matrices of Level-Dependent QBD Processes", QTNA 2010.
%}
function [pi, pi_cell] = ldqbd_pi(R, Q0, Q1, Q2, varargin)

% Parse options
if nargin >= 5 && ~isempty(varargin{1})
    options = varargin{1};
else
    options = struct();
end

if ~isfield(options, 'verbose'), options.verbose = false; end

% Determine number of levels
N = length(Q1) - 1;  % Maximum level

if options.verbose
    fprintf('LD-QBD Stationary Distribution: N=%d levels\n', N);
end

% Check if all levels have same dimension (homogeneous)
isScalar = all(cellfun(@(x) numel(x) == 1, Q1));
dims = cellfun(@(x) size(x, 1), Q1);
isHomogeneous = (length(unique(dims)) == 1);

if isScalar && isHomogeneous
    % Scalar case: direct computation
    pi = zeros(1, N+1);
    pi_cell = cell(N+1, 1);

    % Start with pi_0 = 1 (will normalize later)
    pi_cell{1} = 1;

    % Forward recursion: pi_n = pi_{n-1} * R^(n)
    for n = 1:N
        if ~isempty(R{n})
            pi_cell{n+1} = pi_cell{n} * R{n};
        else
            pi_cell{n+1} = 0;
        end
    end

    % Convert to vector and normalize
    for n = 0:N
        pi(n+1) = pi_cell{n+1};
    end
    pi = pi / sum(pi);

    % Update pi_cell with normalized values
    for n = 0:N
        pi_cell{n+1} = pi(n+1);
    end
else
    % Heterogeneous or matrix case: use cell array for pi with varying dimensions
    pi_cell = cell(N+1, 1);

    % Initialize pi_0 based on its dimension
    if numel(Q1{1}) == 1
        % Scalar level 0: start with pi_0 = 1
        pi_cell{1} = 1;
    else
        % Matrix level 0: solve boundary equation
        Q1_0 = Q1{1};
        Q2_1 = Q2{1};
        A = Q1_0 + R{1} * Q2_1;

        % Find null space of A' (left null space of A)
        [V, D] = eig(A');
        [~, idx] = min(abs(diag(D)));
        pi0 = real(V(:, idx))';

        % Ensure positive and normalize
        if sum(pi0) < 0
            pi0 = -pi0;
        end
        pi0 = abs(pi0);
        pi0 = pi0 / sum(pi0);
        pi_cell{1} = pi0;
    end

    % Forward recursion: pi_n = pi_{n-1} * R^(n)
    for n = 1:N
        pi_cell{n+1} = pi_cell{n} * R{n};
    end

    % Normalize and convert to vector (sum over phases at each level)
    total = sum(cellfun(@sum, pi_cell));
    pi = zeros(1, N+1);
    for n = 0:N
        pi_cell{n+1} = pi_cell{n+1} / total;
        pi(n+1) = sum(pi_cell{n+1});  % Sum over phases if vector
    end
end

if options.verbose
    fprintf('  Stationary distribution computed, sum = %.6f\n', sum(pi));
end

end
