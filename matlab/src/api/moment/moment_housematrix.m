function T = moment_housematrix(edge, n)
% T = moment_housematrix(edge, n)
%
% Conversion matrix of one edge of the house of moments.
%
% The edge is returned as a linear map on the moment subspace {m_0 = 1}. Four
% edges (the Lah pair and the shifted-binomial pair) pin their zeroth output to
% 1 rather than propagating element 0, so as maps of the whole space they are
% affine. Here the offset is folded into column 0, which is empty for those
% edges, making every edge a genuine matrix. On a moment vector, whose element
% 0 is 1 by definition, the two agree. This is also what makes those edges
% usable dimension by dimension in the joint conversions.
%
% Input:
%   edge: one of 'factorial_from_raw', 'raw_from_factorial',
%         'upfactorial_from_raw', 'raw_from_upfactorial',
%         'binomial_from_factorial', 'factorial_from_binomial',
%         'negbinomial_from_upfactorial', 'upfactorial_from_negbinomial',
%         'factorial_from_upfactorial', 'upfactorial_from_factorial',
%         'negbinomial_from_binomial', 'binomial_from_negbinomial',
%         'binomial_from_tail', 'tail_from_binomial'
%   n: maximum order of the mode (n >= 0)
%
% Output:
%   T: (n+1)x(n+1) conversion matrix
%
% Example:
%   T = moment_housematrix('factorial_from_raw', 4);
%
% Reference:
% A. Heindl and A. van de Liefvoort. Moment conversions for discrete
% distributions. PMCCS, 2003.

if ~isscalar(n) || n < 0 || n ~= round(n)
    line_error(mfilename,'The maximum order n must be a nonnegative integer.');
end
switch edge
    case 'factorial_from_raw'
        T = moment_stirling1(n);
    case 'raw_from_factorial'
        T = moment_stirling2(n);
    case 'upfactorial_from_raw'
        T = moment_stirlingcycle(n);
    case 'raw_from_upfactorial'
        S = moment_stirling2(n);
        T = zeros(n+1,n+1);
        for i = 0:n
            for j = 0:i
                T(i+1,j+1) = (-1)^(i-j) * S(i+1,j+1);
            end
        end
    case {'binomial_from_factorial','negbinomial_from_upfactorial'}
        T = diag(1 ./ factorial(0:n));
    case {'factorial_from_binomial','upfactorial_from_negbinomial'}
        T = diag(factorial(0:n));
    case {'factorial_from_upfactorial','upfactorial_from_factorial'}
        L = moment_lah(n);
        T = zeros(n+1,n+1);
        T(1,1) = 1;
        for i = 1:n
            for k = 1:i
                if strcmp(edge,'upfactorial_from_factorial')
                    T(i+1,k+1) = L(i+1,k+1);
                else
                    T(i+1,k+1) = (-1)^(i-k) * L(i+1,k+1);
                end
            end
        end
    case {'negbinomial_from_binomial','binomial_from_negbinomial'}
        T = zeros(n+1,n+1);
        T(1,1) = 1;
        for i = 1:n
            for k = 1:i
                c = nchoosek(i-1,k-1);
                if strcmp(edge,'negbinomial_from_binomial')
                    T(i+1,k+1) = c;
                else
                    T(i+1,k+1) = (-1)^(i-k) * c;
                end
            end
        end
    case {'binomial_from_tail','tail_from_binomial'}
        T = zeros(n+1,n+1);
        T(1,1) = 1;
        for i = 1:n
            for k = i:n
                c = nchoosek(k-1,i-1);
                if strcmp(edge,'binomial_from_tail')
                    T(i+1,k+1) = c;
                else
                    T(i+1,k+1) = (-1)^(k-i) * c;
                end
            end
        end
    otherwise
        line_error(mfilename,sprintf('Unknown edge %s.',edge));
end
end
