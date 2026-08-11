function F = moment_joint_marking(f, p, dims)
% F = moment_joint_marking(f, p, dims)
%
% Joint factorial moments of the per-class counts under multinomial marking.
%
% If a count N is marked independently, every event receiving class j with
% probability p_j, then the per-class counts (N_1,...,N_d) have joint factorial
% moments
%
%   E[prod_j (N_j)_(a_j)] = (prod_j p_j^(a_j)) * f_(|a|)
%
% where f is the factorial moment sequence of the aggregate count N and
% |a| = a_1+...+a_d. This is the counting-process counterpart of the marking
% (class-splitting) formulas of the M3A fitters, and it is exact for the
% per-class counts of a MAP marked in this i.i.d. way, in particular for an
% MMAP whose marking probabilities do not depend on the phase.
%
% Input:
%   f: vector of length n+1 holding the factorial moments f_0,...,f_n of the
%      aggregate count
%   p: vector of length d holding the marking probabilities
%   dims: vector of length d holding the maximum order per class. Their sum
%         must not exceed n, since an entry of multi-order a consumes the
%         aggregate moment of order |a|
%
% Output:
%   F: array of size (dims(1)+1)x...x(dims(d)+1) holding the joint factorial
%      moments of the per-class counts
%
% Example:
%   F = moment_joint_marking([1, 2, 4, 8], [0.3, 0.7], [1, 1]);
%
% Reference:
% A. Heindl and A. van de Liefvoort. Moment conversions for discrete
% distributions. PMCCS, 2003.

fcol = f(:);
pv = p(:).';
dv = dims(:).';
d = numel(pv);
if numel(dv) ~= d
    line_error(mfilename,'The vectors p and dims must have the same length.');
end
if any(dv < 0) || any(dv ~= round(dv))
    line_error(mfilename,'The maximum orders must be nonnegative integers.');
end
if sum(dv) > length(fcol)-1
    line_error(mfilename,'The aggregate factorial moments must reach order sum(dims).');
end
sz = dv + 1;
nel = prod(sz);
Fv = zeros(nel,1);
a = ones(1,d);
for ia = 1:nel
    ord = a - 1;
    Fv(ia) = prod(pv .^ ord) * fcol(sum(ord)+1);
    for l = 1:d
        a(l) = a(l) + 1;
        if a(l) <= sz(l)
            break
        end
        a(l) = 1;
    end
end
F = reshape(Fv, [sz 1]);
end
