function [MMAP, EXACT] = mmap2k_fit(M1, M2, M3, GAMMA, P, F, B)
% [MMAP, EXACT] = mmap2k_fit(M1, M2, M3, GAMMA, P, F, B)
% Closed-form fit of an MMAP(2,K): a marked MAP of second order with K classes.
%
% The fit factorizes into two independent inverse problems:
%  1. the underlying MAP(2) from (M1, M2, M3, GAMMA), by the acyclic canonical
%     inverse of amap2_fit_gamma. Order two loses nothing: every MAP(2) is
%     equivalent to one of the two acyclic canonical forms;
%  2. the marking, which is LINEAR in the class characteristics. Writing the
%     marking as fractions of the aggregate D1,
%        form 1:  D1{c} = D1 .* [q1c 0; q2c q3c]
%        form 2:  D1{c} = D1 .* [0 q1c; q2c q3c]
%     the map (q1c,q2c,q3c) -> (p_c, p_c F_c, p_c B_c) is linear, with a 3x3
%     matrix depending only on (h1,h2,r1,r2) and IDENTICAL for every class. Its
%     inverse is therefore pre-computed once below and applied per class, so the
%     cost does not grow with K and no quadratic program is needed.
%
% The inverse was derived symbolically in SageMath, see
% io/sage/proofs/mmap2k_marking_inverse.py. It is singular exactly on the six
% degenerate loci r1 in {0,1}, r2 in {0,1}, h1-h2+h2*r1 = 0 and, per form,
% h1*r2-h2 = 0 or h1*r1*r2-h1*r1+h1-h2 = 0, which are the branches
% mamap2m_fit_fb_multiclass handles. There, and whenever the closed form leaves
% the unit box, this function falls back to that quadratic program.
%
% Input:
%  M1,M2,M3: raw moments of the inter-arrival times
%  GAMMA: autocorrelation decay rate
%  P: class probabilities (sum to one)
%  F: first-order forward moments (sum_c P(c) F(c) = M1)
%  B: first-order backward moments (sum_c P(c) B(c) = M1)
%
% Output:
%  MMAP: the fitted MMAP as {D0, D1, D11, ..., D1K}
%  EXACT: true when the closed form matched P, F and B exactly

degentol = 1e-8;
feastol = 1e-8;

P = P(:); F = F(:); B = B(:);
K = length(P);
if length(F) ~= K || length(B) ~= K
    error('mmap2k_fit: P, F and B must have the same length');
end

[~, AMAPS] = amap2_fit_gamma(M1, M2, M3, GAMMA);

BEST = []; BESTERR = Inf;
for j = 1:length(AMAPS)
    D0 = AMAPS{j}{1};
    D1 = AMAPS{j}{2};
    if size(D0,1) ~= 2 || D0(2,1) ~= 0
        continue
    end
    if D1(1,2) == 0
        form = 1;
    elseif D1(1,1) == 0
        form = 2;
    else
        continue
    end
    h1 = -1/D0(1,1);
    h2 = -1/D0(2,2);
    r1 = D0(1,2) * h1;
    r2 = D1(2,2) * h2;

    q = zeros(3,K);
    ok = true;
    for c = 1:K
        [q(1,c), q(2,c), q(3,c), good] = marking_inverse(form, h1, h2, r1, r2, P(c), F(c), B(c), degentol);
        if ~good
            ok = false;
            break
        end
    end
    if ~ok
        continue
    end
    % feasible when the fractions lie in the unit box and close per phase
    viol = max(0, max(max(-q))) + max(0, max(max(q-1))) + max(abs(sum(q,2) - 1));
    if viol < BESTERR
        BESTERR = viol;
        BEST = {D0, D1, form, q};
    end
end

if ~isempty(BEST) && BESTERR <= feastol
    D0 = BEST{1}; D1 = BEST{2}; form = BEST{3}; q = min(max(BEST{4},0),1);
    MMAP = cell(1,2+K);
    MMAP{1} = D0;
    MMAP{2} = D1;
    for c = 1:K
        if form == 1
            MMAP{2+c} = D1 .* [q(1,c) 0; q(2,c) q(3,c)];
        else
            MMAP{2+c} = D1 .* [0 q(1,c); q(2,c) q(3,c)];
        end
    end
    EXACT = true;
    return
end

% degenerate underlying form, or characteristics outside the feasible set
MMAP = mamap2m_fit_gamma_fb(M1, M2, M3, GAMMA, P, F, B);
EXACT = false;

end

function [q1, q2, q3, good] = marking_inverse(form, h1, h2, r1, r2, p, Fc, Bc, degentol)
% Pre-computed inverse of the per-class linear map, see the header.
q1 = 0; q2 = 0; q3 = 0;
good = false;
if form == 1
    d1 = (r2-1)*(r1-1)*(h1*r2-h2);
    d2 = (r2-1)*r1;
    d3 = r1*r2*(h1+h2*r1-h2);
    if min([abs(d1) abs(d2) abs(d3) abs(h1+h2*r1-h2) abs(h1*r2-h2)]) < degentol
        return
    end
    W = r1*r2 - r2 + 1;
    q1 = p * W * ((h1*r2 - h1 - h2) + Bc) / d1;
    q2 = p * W * ((h1^2*(r2-1) + h1*h2*r1*(r2-1) - h2^2*r1)/((h1+h2*r1-h2)*(h1*r2-h2)) ...
                  - Fc/(h1+h2*r1-h2) + Bc/(h1*r2-h2)) / d2;
    q3 = p * W * ((h1 + h2*r1) - Fc) / d3;
else
    U = h1*r1*r2 - h1*r1 + h1 - h2;
    d1 = (r2-1)*(r1-1)*U;
    d2 = (r2-1)*(h1+h2*r1-h2);
    if min([abs(d1) abs(d2) abs(r2) abs(U) abs(h1+h2*r1-h2)]) < degentol
        return
    end
    V = r1*r2 - r1 - r2 + 2;
    q1 = p * V * ((h1*r1*r2 - h1*r1 - h2) + Bc) / d1;
    q2 = p * V * (h2 - Fc) / d2;
    q3 = p * V * ((h1^2 + h1*h2*r1*r2 - h2^2)/((h1+h2*r1-h2)*U) ...
                  - Fc/(h1+h2*r1-h2) - Bc/U) / r2;
end
good = true;
end
