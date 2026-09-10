function bmax = cd_peak_scaling(beta, NK, K) %#ok<INUSD>
% CD_PEAK_SCALING Peak of a class-dependence handle over the population lattice
%
% bmax = CD_PEAK_SCALING(beta, NK, K)
%
% Peak of the class-dependence handle over the reachable population lattice
% 0 <= n(r) <= NK(r). The handle returns either a scalar (shared by every
% class) or a length-K vector, so the peak is taken over both the states and
% the classes: utilization is a per-station quantity, so the whole station
% shares one normalizer, as it does for max(lldscaling(ist,:)).
%
% This is the single normalizer used to report utilization at stations with
% limited class dependence, U = T*S/bmax, so that every solver follows the
% same convention as solver_ncld does for lldscaling (U/max(lldscaling)).

bmax = 0;
n = pprod_init(NK);
while n(1) >= 0
    if sum(n) > 0
        v = beta(n);
        v = v(isfinite(v));
        if ~isempty(v)
            bmax = max(bmax, max(v));
        end
    end
    n = pprod_next(n, NK);
end
end

function n = pprod_init(N)
n = zeros(size(N));
end

function n = pprod_next(n, N)
R = length(N);
if all(n == N)
    n = -1 * ones(1, R);
    return
end
s = R;
while s > 0 && n(s) == N(s)
    n(s) = 0;
    s = s - 1;
end
if s > 0
    n(s) = n(s) + 1;
end
end
