%{ @file cache_lrum_map_levelstats.m
 %  @brief Per-item level statistics for the LRU(m)-MAP TTL approximation
 %
 %  @author LINE Development Team
%}

%{
 % @brief Level statistics of one item's embedded (list, phase) chain
 %
 % @details
 % Evaluates eqs. (5)-(9) of Gast and Van Houdt (Performance Evaluation
 % 2017) for a single item with MAP request process (D0, D1) and
 % characteristic times T: the R-recursions of the embedded chain, the
 % level-0 boundary vector pi_0 (left Perron vector of R_1 expm(D0 T_1)),
 % the time-stationary level probabilities, and the request-weighted hit
 % fractions per list.
 %
 % @par Syntax:
 % @code
 % [prob, occ, hitfrac] = cache_lrum_map_levelstats(D0, D1, T)
 % @endcode
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>prob<td>(1,h+1) time-stationary probability of level 0..h
 % <tr><td>occ<td>(1,h) occupancy contribution of lists 1..h (= prob(2:end))
 % <tr><td>hitfrac<td>(1,h) fraction of the item's requests hitting in list l
 % </table>
%}
function [prob, occ, hitfrac] = cache_lrum_map_levelstats(D0, D1, T)
h = numel(T);
d = size(D0, 1);
iD0 = inv(-D0);

E = cell(1, h);   % expm(D0*T_l)
A = cell(1, h);   % A_l of the paper
Nh = cell(1, h);
for l = 1:h
    E{l} = expm(D0 * T(l));
    Nh{l} = (eye(d) - E{l}) * iD0; %#ok<MINV>
    A{l} = Nh{l} * D1;
end
A0 = iD0 * D1; %#ok<MINV>
N0 = iD0;

% R recursion, eqs. (6)-(7): R{l} is the paper's R_l over lists 1..h
R = cell(1, h);
for l = h:-1:1
    if l == h
        if h == 1
            Aprev = A0;
        else
            Aprev = A{h-1};
        end
        R{l} = Aprev / (eye(d) - A{h});
    elseif l == 1
        R{l} = A0 / (eye(d) - R{l+1} * E{l+1});
    else
        R{l} = A{l-1} / (eye(d) - R{l+1} * E{l+1});
    end
end

% pi_0: left Perron vector of R_1 expm(D0 T_1) (level-0 balance)
M = R{1} * E{1};
[V, D] = eig(M');
[~, idx] = max(real(diag(D)));
pi0 = real(V(:, idx))';
pi0 = pi0 / sum(pi0);

pih = cell(1, h);
pih{1} = pi0 * R{1};
for l = 2:h
    pih{l} = pih{l-1} * R{l};
end

e = ones(d, 1);
holding = zeros(1, h+1);
holding(1) = pi0 * N0 * e;
for l = 1:h
    holding(1+l) = pih{l} * Nh{l} * e;
end
denom = sum(holding);
prob = holding / denom;
occ = prob(2:end);

% request-weighted hit fractions: throughput of hits in list l over the
% item's stationary request rate
Q = D0 + D1;
piphase = ctmc_solve(Q);
lam = piphase * D1 * e;
hitfrac = zeros(1, h);
for l = 1:h
    hitfrac(l) = (pih{l} * Nh{l} * D1 * e) / denom;
end
if lam > 0
    hitfrac = hitfrac / lam;
else
    hitfrac = zeros(1, h);
end
end
