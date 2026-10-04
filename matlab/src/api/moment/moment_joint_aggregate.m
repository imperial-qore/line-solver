function f = moment_joint_aggregate(F)
% f = moment_joint_aggregate(F)
%
% Factorial moments of a total count from the joint factorial moments of its
% parts. For N = N_1+...+N_d the Vandermonde convolution of falling factorials
% gives
%
%   f_n = sum_(|a|=n) (n! / prod_j a_j!) * F_a
%
% which holds for ANY joint law of the parts, marked or not, and is the inverse
% direction of moment_joint_marking whenever the marking is multinomial. The
% order reached is limited by the smallest per-class order in F, since the term
% a = n*e_j must be available for every j.
%
% Input:
%   F: array of size (n_1+1)x...x(n_d+1) holding the joint factorial moments of
%      the parts
%
% Output:
%   f: column vector of length min_j(n_j)+1 holding f_0,...,f_min_j(n_j), the
%      factorial moments of the total
%
% Example:
%   f = moment_joint_aggregate(moment_joint_marking([1, 2, 4], [0.3, 0.7], [1, 1]));
%
% Reference:
% A. Heindl and A. van de Liefvoort. Moment conversions for discrete
% distributions. PMCCS, 2003.

sz = moment_tensorsize(F);
d = numel(sz);
nel = prod(sz);
nmax = min(sz) - 1;
Fv = reshape(F, [], 1);
f = zeros(nmax+1,1);
a = ones(1,d);
for ia = 1:nel
    ord = a - 1;
    n = sum(ord);
    if n <= nmax
        f(n+1) = f(n+1) + factorial(n) / prod(factorial(ord)) * Fv(ia);
    end
    for l = 1:d
        a(l) = a(l) + 1;
        if a(l) <= sz(l)
            break
        end
        a(l) = 1;
    end
end
end
