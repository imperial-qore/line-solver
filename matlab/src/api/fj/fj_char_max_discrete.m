%{ @file fj_char_max_discrete.m
 %  @brief Characteristic maximum of a discrete random variable
 %
 %  @author LINE Development Team
%}

%{
 % @brief Characteristic maximum of a discrete random variable
 %
 % @details
 % Gravey's characteristic maximum for lattice distributions. Let m_K be the
 % smallest integer with P(X > m_K) <= 1/K; then
 %
 %   M_K = m_K + K * sum_{k >= m_K} P(X > k),
 %
 % which upper bounds the expected maximum of K i.i.d. copies of X and costs
 % O(1) instead of the alternating binomial sum. Two lattice laws admit a
 % closed form for the tail sum:
 %
 %   geometric, P(X = k) = (1-p) p^k:
 %       m_K = ceil( -ln(K) / ln(p) ),  M_K = m_K + K * p^(m_K+1) / (1-p),
 %       and the exact maximum is
 %       E[Y_K] = sum_{k=1..K} binom(K,k) (-1)^(k+1) p^k / (1 - p^k);
 %
 %   Poisson, P(X = k) = exp(-theta) theta^k / k!:
 %       M_K = m_K * (1 - K*P(X > m_K)) + K*theta*P(X > m_K - 1),
 %   which is the same tail sum rewritten through E[(X-m)^+] = theta*P(X>m-1)
 %   - m*P(X>m).
 %
 % The exact expected maximum is returned alongside: in closed form for the
 % geometric law, and by summing 1 - F(k)^K over the lattice for the Poisson.
 %
 % @par Syntax:
 % @code
 % [MK, mK] = fj_char_max_discrete(K, dist_type, param)
 % [MK, mK, exact] = fj_char_max_discrete(K, dist_type, param)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>K<td>Number of i.i.d. copies (positive integer)
 % <tr><td>dist_type<td>'geometric' or 'poisson'
 % <tr><td>param<td>Success-complement p in (0,1) for the geometric, mean theta > 0 for the Poisson
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>MK<td>Characteristic maximum, an upper bound on the expected maximum
 % <tr><td>mK<td>Smallest integer with P(X > mK) <= 1/K
 % <tr><td>exact<td>Exact expected maximum of the K copies
 % </table>
 %
 % @par Reference:
 % A. Thomasian, "Analysis of Fork/Join and Related Queueing Systems",
 % ACM Computing Surveys, Vol. 47, No. 2, Article 17, July 2014, Eq. (48) and
 % the lattice cases on page 17:25.
 %
 % Original: A. Gravey, "A Simple Construction of an Upper Bound for the Mean
 % of the Maximum of N Identically Distributed Random Variables",
 % J. Applied Probability 22(4), 1985.
%}
function [MK, mK, exact] = fj_char_max_discrete(K, dist_type, param)

if K < 1 || K ~= round(K)
    line_error(mfilename, 'K must be a positive integer. Got K=%g.', K);
end

switch lower(dist_type)
    case {'geometric', 'geom'}
        p = param;
        if p <= 0 || p >= 1
            line_error(mfilename, 'The geometric parameter p must lie in (0,1). Got p=%g.', p);
        end
        % Smallest integer k with p^k <= 1/K
        mK = ceil(-log(K) / log(p));
        mK = max(mK, 0);
        MK = mK + K * p^(mK + 1) / (1 - p);
        % Exact maximum by inclusion-exclusion on the geometric tail
        exact = 0;
        for k = 1:K
            exact = exact + nchoosek(K, k) * (-1)^(k + 1) * p^k / (1 - p^k);
        end

    case {'poisson', 'pois'}
        theta = param;
        if theta <= 0
            line_error(mfilename, 'The Poisson mean theta must be positive. Got theta=%g.', theta);
        end
        % Walk the lattice until the tail drops to 1/K, then to numerical zero
        kmax = ceil(theta + 12 * sqrt(theta) + 40);
        k = (0:kmax);
        pmf = exp(-theta + k * log(theta) - gammaln(k + 1));
        cdf = cumsum(pmf);
        cdf = min(cdf, 1);
        tail = 1 - cdf;
        mK = find(tail <= 1 / K, 1, 'first');
        if isempty(mK)
            line_error(mfilename, 'The Poisson lattice truncation at %d did not reach a tail of 1/K.', kmax);
        end
        mK = mK - 1;
        if mK == 0
            tail_prev = 1;
        else
            tail_prev = tail(mK);
        end
        MK = mK * (1 - K * tail(mK + 1)) + K * theta * tail_prev;
        % Exact maximum as the sum over the lattice of 1 - F(k)^K
        exact = sum(1 - cdf.^K);

    otherwise
        line_error(mfilename, 'Unsupported discrete distribution "%s". Use geometric or poisson.', dist_type);
end

end
