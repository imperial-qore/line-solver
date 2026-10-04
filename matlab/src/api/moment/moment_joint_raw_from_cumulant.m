function m = moment_joint_raw_from_cumulant(kappa)
% m = moment_joint_raw_from_cumulant(kappa)
%
% Converts the joint cumulants of a random vector into its joint power (raw)
% moments, by running the multivariate exponential-formula recursion forward,
%
%   m_a = sum_(0<b<=a) prod_l nchoosek(a_l-[l=j], b_l-[l=j]) kappa_b m_(a-b)
%
% with m_0 = 1 and j the first dimension in which a is nonzero. Inverse of
% moment_joint_cumulant_from_raw.
%
% Input:
%   kappa: array of size (n_1+1)x...x(n_d+1) holding the joint cumulants;
%          element 1 is ignored
%
% Output:
%   m: array of the same size holding the joint power moments, element 1 being
%      1
%
% Example:
%   m = moment_joint_raw_from_cumulant(moment_joint_cumulant_from_raw(m0));
%
% Reference:
% V. P. Leonov and A. N. Shiryaev. On a method of calculation of
% semi-invariants. Theory of Probability and its Applications,
% 4(3):319-329, 1959.

sz = moment_tensorsize(kappa);
d = numel(sz);
nel = prod(sz);
stride = cumprod([1, sz(1:end-1)]);
kv = reshape(kappa, [], 1);
mv = zeros(nel,1);
mv(1) = 1;
a = ones(1,d);
for ia = 1:nel
    if any(a > 1)
        j = find(a > 1, 1);
        acc = 0;
        b = ones(1,d);
        for ib = 1:prod(a)
            if any(b > 1) && b(j) > 1
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
        mv(ia) = acc;
    end
    for l = 1:d
        a(l) = a(l) + 1;
        if a(l) <= sz(l)
            break
        end
        a(l) = 1;
    end
end
m = reshape(mv, size(kappa));
end
