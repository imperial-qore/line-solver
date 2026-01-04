%{ @file ldqbd.m
 %  @brief Solves level-dependent QBD processes using matrix continued fractions
 %
 %  @author QORE Lab, Imperial College London (Algorithm by Phung-Duc et al.)
%}

%{
 % @brief Solves level-dependent QBD processes using matrix continued fractions
 %
 % @details
 % This function implements Algorithm 1 from "A Simple Algorithm for the Rate
 % Matrices of Level-Dependent QBD Processes" by Phung-Duc, Masuyama, Kasahara,
 % and Takahashi (2010), QTNA Conference.
 %
 % It computes the rate matrices R^(n) for an ergodic level-dependent QBD
 % process with a finite number of levels, and optionally the stationary
 % distribution.
 %
 % For a level-dependent QBD with levels 0, 1, ..., N, the infinitesimal
 % generator has block-tridiagonal structure:
 %
 %   Q^(0)_1  Q^(0)_0   O        O      ...   O
 %   Q^(1)_2  Q^(1)_1  Q^(1)_0   O      ...   O
 %    O       Q^(2)_2  Q^(2)_1  Q^(2)_0 ...   O
 %   ...      ...      ...      ...    ...  ...
 %    O        O        O       ...   Q^(N)_2  Q^(N)_1
 %
 % where Q_0^(n) are upward transitions (level n to n+1), Q_1^(n) are local
 % transitions (within level n), and Q_2^(n) are downward transitions
 % (level n to n-1).
 %
 % @par Syntax:
 % @code
 % [R, pi] = ldqbd(Q0, Q1, Q2)
 % [R, pi] = ldqbd(Q0, Q1, Q2, options)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>Q0<td>Cell array {Q0^(0), Q0^(1), ..., Q0^(N-1)} - upward transitions
 % <tr><td>Q1<td>Cell array {Q1^(0), Q1^(1), ..., Q1^(N)} - local transitions
 % <tr><td>Q2<td>Cell array {Q2^(1), Q2^(2), ..., Q2^(N)} - downward transitions
 % <tr><td>options<td>Optional struct with fields: epsilon (default 1e-10),
 %                    maxIter (default 1000), verbose (default false)
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>R<td>Cell array {R^(1), R^(2), ..., R^(N)} - rate matrices
 % <tr><td>pi<td>Row vector [pi_0, pi_1, ..., pi_N] - stationary distribution
 % </table>
 %
 % @par Reference:
 % T. Phung-Duc, H. Masuyama, S. Kasahara, Y. Takahashi, "A Simple Algorithm
 % for the Rate Matrices of Level-Dependent QBD Processes", QTNA 2010.
 %
 % @see ldqbd_R, ldqbd_pi
%}
function [R, pi] = ldqbd(Q0, Q1, Q2, varargin)

% Parse options
if nargin >= 4 && ~isempty(varargin{1})
    options = varargin{1};
else
    options = struct();
end

if ~isfield(options, 'epsilon'), options.epsilon = 1e-10; end
if ~isfield(options, 'maxIter'), options.maxIter = 1000; end
if ~isfield(options, 'verbose'), options.verbose = false; end

% Determine number of levels
N = length(Q1) - 1;  % Maximum level

if options.verbose
    fprintf('LD-QBD Solver: N=%d levels, epsilon=%.2e\n', N, options.epsilon);
end

% Compute all R matrices using ldqbd_R
R = ldqbd_R(Q0, Q1, Q2, options);

% Compute stationary distribution if requested using ldqbd_pi
if nargout >= 2
    pi = ldqbd_pi(R, Q0, Q1, Q2, options);
end

end
