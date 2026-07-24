%{ @file fj_char_max.m
 %  @brief Characteristic maximum M_K for order statistics
 %
 %  @author LINE Development Team
%}

%{
 % @brief Characteristic maximum M_K for order statistics
 %
 % @details
 % Computes the characteristic maximum M_K, which provides a bound for
 % the expected value of the maximum of K i.i.d. random variables.
 %
 % Let m_K be the greatest lower bound such that P(X > m_K) <= 1/K.
 % For continuous distributions with a density: P(X > m_K) = 1/K.
 %
 % The characteristic maximum is:
 %   M_K = m_K + K * integral_{m_K}^{inf} P(X > x) dx
 %
 % For specific distributions:
 %   - Exponential: m_K = ln(K)/μ, M_K = H_K/μ
 %   - Erlang-k: m_K solves exp(-μ*m_K) * sum_{i=0}^{k-1} (μ*m_K)^i/i! = 1/K
 %               M_K = (k/μ) * [1 + K*exp(-μ*m_K)*(μ*m_K)^k/k!]
 %
 % @par Syntax:
 % @code
 % [MK, mK] = fj_char_max(K, mu)                     % Exponential
 % [MK, mK] = fj_char_max(K, mu, 'exp')              % Exponential (explicit)
 % [MK, mK] = fj_char_max(K, [k, mu], 'erlang')      % Erlang-k
 % [MK, mK] = fj_char_max(K, S_X, 'general')         % General (S_X is survival function)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>K<td>Number of random variables (positive integer)
 % <tr><td>param<td>Distribution parameters (rate mu, or [k, mu] for Erlang)
 % <tr><td>dist_type<td>Optional: 'exp' (default), 'erlang', or 'general'
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>MK<td>Characteristic maximum M_K
 % <tr><td>mK<td>Threshold m_K where P(X > m_K) = 1/K
 % </table>
 %
 % @par Reference:
 % A. Thomasian, "Analysis of Fork/Join and Related Queueing Systems",
 % ACM Computing Surveys, Vol. 47, No. 2, Article 17, July 2014.
 % Eq. (47) on page 17:24.
 %
 % A. Gravey, "A Simple Construction of an Upper Bound for the Mean of
 % the Maximum of n Identically Distributed Random Variables", 1985.
%}
function [MK, mK] = fj_char_max(K, param, dist_type)

if nargin < 3
    dist_type = 'exp';
end

if K < 1
    line_error(mfilename, 'K must be a positive integer. Got K=%d.', K);
end

switch lower(dist_type)
    case 'exp'
        % Exponential distribution with rate mu
        mu = param;
        if mu <= 0
            line_error(mfilename, 'Rate mu must be positive.');
        end
        % m_K = ln(K) / mu
        mK = log(K) / mu;
        % M_K = H_K / mu (exact for exponential)
        MK = fj_harmonic(K) / mu;

    case 'erlang'
        % Erlang-k distribution
        if numel(param) ~= 2
            line_error(mfilename, 'For Erlang, param must be [k, mu].');
        end
        k = param(1);
        mu = param(2);

        if k < 1
            line_error(mfilename, 'Erlang stages k must be positive.');
        end
        if mu <= 0
            line_error(mfilename, 'Rate mu must be positive.');
        end

        % Find m_K by solving: exp(-mu*m_K) * sum_{i=0}^{k-1} (mu*m_K)^i/i! = 1/K
        mK = find_mK_erlang(K, k, mu);

        % M_K = (k/mu) * [1 + K * exp(-mu*m_K) * (mu*m_K)^k / k!]
        MK = (k / mu) * (1 + K * exp(-mu * mK) * (mu * mK)^k / factorial(k));

    case 'general'
        % General distribution: param is survival function S_X(x) = P(X > x)
        S_X = param;
        if ~isa(S_X, 'function_handle')
            line_error(mfilename, 'For general distribution, param must be a survival function handle.');
        end

        % Find m_K such that S_X(m_K) = 1/K
        mK = find_mK_general(K, S_X);

        % M_K = m_K + K * integral_{m_K}^{inf} S_X(x) dx
        integrand = @(x) S_X(x);
        upper_limit = find_upper_limit(S_X, 1e-10);
        tail_integral = integral(integrand, mK, upper_limit, 'RelTol', 1e-8, 'AbsTol', 1e-10);
        MK = mK + K * tail_integral;

    otherwise
        line_error(mfilename, 'Unknown distribution type: %s. Valid: exp, erlang, general.', dist_type);
end

end

function mK = find_mK_erlang(K, k, mu)
% Find m_K for Erlang-k by solving:
% exp(-mu*m_K) * sum_{i=0}^{k-1} (mu*m_K)^i/i! = 1/K

% Erlang survival function
S_erlang = @(x) erlang_survival(x, k, mu);

% Use fzero to find the root
objective = @(x) S_erlang(x) - 1/K;

% Initial guess: for large K, m_K is approximately the mean plus adjustment
initial_guess = k / mu + log(K) / mu;

options = optimset('TolX', 1e-10, 'Display', 'off');
mK = fzero(objective, initial_guess, options);
end

function S = erlang_survival(x, k, mu)
% Erlang-k survival function: P(X > x) = exp(-mu*x) * sum_{i=0}^{k-1} (mu*x)^i/i!
if x <= 0
    S = 1;
else
    S = 0;
    for i = 0:(k-1)
        S = S + (mu * x)^i / factorial(i);
    end
    S = exp(-mu * x) * S;
end
end

function mK = find_mK_general(K, S_X)
% Find m_K such that S_X(m_K) = 1/K
objective = @(x) S_X(x) - 1/K;

% Find bounds
lower = 0;
upper = 1;
while S_X(upper) > 1/K
    upper = upper * 2;
end

options = optimset('TolX', 1e-10, 'Display', 'off');
mK = fzero(objective, [lower, upper], options);
end

function u = find_upper_limit(S_X, target)
% Find u such that S_X(u) <= target
u = 1;
max_iter = 100;
for iter = 1:max_iter
    if S_X(u) <= target
        break;
    end
    u = u * 2;
end
end
