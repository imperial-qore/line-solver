%{ @file ldqbd_R.m
 %  @brief Computes rate matrices for level-dependent QBD processes
 %
 %  @author QORE Lab, Imperial College London (Algorithm by Phung-Duc et al.)
%}

%{
 % @brief Computes all rate matrices R^(n) for a level-dependent QBD
 %
 % @details
 % This function implements the backward recursion algorithm for computing
 % rate matrices of level-dependent QBD processes with heterogeneous dimensions.
 %
 % For a level-dependent QBD with levels 0, 1, ..., N:
 %   R^(N) = Q0^(N-1) * (-Q1^(N))^{-1}
 %   R^(n) = Q0^(n-1) * (-Q1^(n) - R^(n+1) * Q2^(n+1))^{-1}  for n = N-1,...,1
 %
 % Dimensions:
 %   R^(n) is (states at level n-1) x (states at level n)
 %   Q0^(n-1) is (states at level n-1) x (states at level n)
 %   Q1^(n) is (states at level n) x (states at level n)
 %   Q2^(n) is (states at level n) x (states at level n-1)
 %
 % @par Syntax:
 % @code
 % R = ldqbd_R(Q0, Q1, Q2)
 % R = ldqbd_R(Q0, Q1, Q2, options)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>Q0<td>Cell array {Q0^(0), Q0^(1), ..., Q0^(N-1)} - upward transitions
 % <tr><td>Q1<td>Cell array {Q1^(0), Q1^(1), ..., Q1^(N)} - local transitions
 % <tr><td>Q2<td>Cell array {Q2^(1), Q2^(2), ..., Q2^(N)} - downward transitions
 % <tr><td>options<td>Optional struct with fields: verbose (default false)
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>R<td>Cell array {R^(1), R^(2), ..., R^(N)} - rate matrices
 % </table>
 %
 % @par Reference:
 % T. Phung-Duc, H. Masuyama, S. Kasahara, Y. Takahashi, "A Simple Algorithm
 % for the Rate Matrices of Level-Dependent QBD Processes", QTNA 2010.
%}
function R = ldqbd_R(Q0, Q1, Q2, varargin)

% Parse options
if nargin >= 4 && ~isempty(varargin{1})
    options = varargin{1};
else
    options = struct();
end

if ~isfield(options, 'verbose'), options.verbose = false; end

% Determine number of levels
N = length(Q1) - 1;  % Maximum level

if options.verbose
    fprintf('LD-QBD Rate Matrix Computation: N=%d levels\n', N);
end

% Initialize R cell array
R = cell(N, 1);

% Start at level N (boundary): R^(N) = Q0^(N-1) * (-Q1^(N))^{-1}
Q0_Nminus1 = Q0{N};      % Q0^(N-1)
Q1_N = Q1{N+1};          % Q1^(N)

U = -Q1_N;
if numel(U) == 1
    if abs(U) > 1e-14
        R{N} = Q0_Nminus1 / U;
    else
        R{N} = Q0_Nminus1 * 0;
    end
else
    if abs(det(U)) > 1e-14
        R{N} = Q0_Nminus1 / U;
    else
        R{N} = Q0_Nminus1 * pinv(U);
    end
end

% Backward recursion for n = N-1 down to 1
for n = N-1:-1:1
    % R^(n) = Q0^(n-1) * (-Q1^(n) - R^(n+1) * Q2^(n+1))^{-1}
    Q0_nminus1 = Q0{n};      % Q0^(n-1): (states at n-1) x (states at n)
    Q1_n = Q1{n+1};          % Q1^(n): (states at n) x (states at n)
    Q2_nplus1 = Q2{n+1};     % Q2^(n+1): (states at n+1) x (states at n)
    R_nplus1 = R{n+1};       % R^(n+1): (states at n) x (states at n+1)

    % Compute R^(n+1) * Q2^(n+1): (states at n) x (states at n)
    RQ2 = R_nplus1 * Q2_nplus1;

    % Compute U = -Q1^(n) - R^(n+1) * Q2^(n+1)
    U = -Q1_n - RQ2;

    % Compute R^(n) = Q0^(n-1) * U^{-1}
    if numel(U) == 1
        if abs(U) > 1e-14
            R{n} = Q0_nminus1 / U;
        else
            R{n} = Q0_nminus1 * 0;
        end
    else
        if abs(det(U)) > 1e-14
            R{n} = Q0_nminus1 / U;
        else
            R{n} = Q0_nminus1 * pinv(U);
        end
    end
end

if options.verbose
    for n = 1:N
        fprintf('  R^(%d) computed (%dx%d)\n', n, size(R{n}, 1), size(R{n}, 2));
    end
end

end
