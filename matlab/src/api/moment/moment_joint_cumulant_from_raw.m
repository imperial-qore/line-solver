function kappa = moment_joint_cumulant_from_raw(m)
% kappa = moment_joint_cumulant_from_raw(m)
%
% Converts the joint power (raw) moments of a random vector (N_1,...,N_d) into
% its joint cumulants, the coefficients of the joint cumulant generating
% function
%
%   log E[exp(s_1 N_1 + ... + s_d N_d)] = sum_(a ~= 0) kappa_a prod_j
%   s_j^(a_j) / a_j!
%
% They obey the multivariate exponential formula, equivalently the
% Leonov-Shiryaev partition formula. With j the first dimension in which the
% multi-index a is nonzero,
%
%   m_a = sum_(0<b<=a) prod_l nchoosek(a_l-[l=j], b_l-[l=j]) kappa_b m_(a-b)
%
% which isolates kappa_a because the b = a term has unit coefficient and
% m_0 = 1. Unlike every other conversion in the house, this one does not factor
% into a product of univariate transforms: the cumulant of multi-order (1,1) is
% the covariance, which mixes the dimensions.
%
% Input:
%   m: array of size (n_1+1)x...x(n_d+1) holding the joint power moments, with
%      element 1 equal to 1
%
% Output:
%   kappa: array of the same size holding the joint cumulants, element 1 being
%          kappa_0 = 0
%
% Example:
%   kappa = moment_joint_cumulant_from_raw(m); % kappa(2,2) is the covariance
%
% Reference:
% V. P. Leonov and A. N. Shiryaev. On a method of calculation of
% semi-invariants. Theory of Probability and its Applications,
% 4(3):319-329, 1959.

sz = moment_tensorsize(m);
d = numel(sz);
nel = prod(sz);
stride = cumprod([1, sz(1:end-1)]);
mv = reshape(m, [], 1);
kv = zeros(nel,1);
a = ones(1,d);
for ia = 1:nel
    if any(a > 1)
        j = find(a > 1, 1);
        acc = 0;
        b = ones(1,d);
        for ib = 1:prod(a)
            if any(b > 1) && ~isequal(b,a) && b(j) > 1
                c = 1;
                for l = 1:d
                    alpha = a(l) - 1 - (l == j);
                    beta = b(l) - 1 - (l == j);
                    c = c * nchoosek(alpha, beta);
                end
                if c ~= 0
                    acc = acc + c * kv(1 + sum((b-1).*stride)) * mv(1 + sum((a-b).*stride));
                end
            end
            for l = 1:d
                b(l) = b(l) + 1;
                if b(l) <= a(l)
                    break
                end
                b(l) = 1;
            end
        end
        kv(ia) = mv(ia) - acc;
    end
    for l = 1:d
        a(l) = a(l) + 1;
        if a(l) <= sz(l)
            break
        end
        a(l) = 1;
    end
end
kappa = reshape(kv, size(m));
end
