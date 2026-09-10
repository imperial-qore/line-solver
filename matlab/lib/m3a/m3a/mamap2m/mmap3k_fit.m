function [MMAP, EXACT] = mmap3k_fit(D0, D1, P, F, B, B2)
% [MMAP, EXACT] = mmap3k_fit(D0, D1, P, F, B, B2)
% Closed-form marking fit of an MMAP(3,K): a marked MAP of third order.
%
% The MMAP(2,K) argument does not depend on the order. Two facts carry over
% (verified symbolically, see io/sage/proofs/mmap3k_marking_inverse.py):
%  1. every per-class characteristic in which the class matrix appears exactly
%     once is LINEAR in the marking fractions, so z = nnz(D1) fractions per
%     class are determined by z characteristics through a square system;
%  2. that system is BLOCK DIAGONAL in the classes, so one z-by-z block is
%     built once and reused for every class: the cost does not grow with K.
%
% What changes with the order is WHICH characteristics are needed. At order two
% (p_c, F_c, B_c) suffice; at order three the independent set of lowest total
% order is
%    (a,b) = (1,0), (1,1), (2,0), (3,0)
%    i.e.    p_c,   F_c,   B_c,   B_c^(2)
% with a the backward and b the forward order of pie*A^a*D1c*A^b*1. Alternating
% forward and backward orders does NOT stay independent at higher orders.
%
% The block is assembled exactly by evaluating the linear map on unit markings:
% no finite differences and no computer algebra, since the map is linear with a
% zero offset. This is preferred to inlining the symbolic inverse, which at
% order three is about 5 KB of expressions and is tied to one sparsity pattern.
%
% The underlying MAP(3) is an input: unlike order two there is no canonical
% inverse here that turns moments and autocorrelation into an order-3 MAP.
%
% Input:
%  D0, D1: the underlying MAP (any order; validated at orders two and three)
%  P: class probabilities (sum to one)
%  F: first-order forward moments
%  B: first-order backward moments
%  B2: second-order backward moments (required from order three)
%
% Output:
%  MMAP: the marked MAP as {D0, D1, D11, ..., D1K}
%  EXACT: true when the closed-form marking is feasible

feastol = 1e-8;

P = P(:); F = F(:); B = B(:);
K = length(P);
n = size(D0,1);

[nzi, nzj] = find(D1);
z = length(nzi);
orders = marking_orders(n, z);
if nargin < 6 || isempty(B2)
    if z > 3
        error('mmap3k_fit: order %d needs the second-order backward moments B2', n);
    end
    B2 = zeros(K,1);
else
    B2 = B2(:);
end

A = inv(-D0);
P_emb = A * D1;
T = (P_emb' - eye(n));
T(n,:) = ones(1,n);
rhs = zeros(n,1); rhs(n) = 1;
pie = (T \ rhs)';

% z-by-z block: column j is the characteristic vector of a unit marking on the
% j-th nonzero of D1
M = zeros(z,z);
for jj = 1:z
    Dc = zeros(n,n);
    Dc(nzi(jj), nzj(jj)) = D1(nzi(jj), nzj(jj));
    for ii = 1:z
        a = orders(ii,1); b = orders(ii,2);
        M(ii,jj) = pie * (A^a) * Dc * (A^b) * ones(n,1);
    end
end

if abs(det(M)) < 1e-12 * max(1, max(abs(M(:)))^z)
    error('mmap3k_fit: the underlying MAP is on the degenerate locus of the marking system');
end

q = zeros(z,K);
for c = 1:K
    y = zeros(z,1);
    for ii = 1:z
        a = orders(ii,1); b = orders(ii,2);
        if a == 1 && b == 0
            y(ii) = P(c);
        elseif a == 1 && b == 1
            y(ii) = P(c) * F(c);
        elseif a == 2 && b == 0
            y(ii) = P(c) * B(c);
        elseif a == 3 && b == 0
            y(ii) = P(c) * B2(c);
        else
            error('mmap3k_fit: no target supplied for the characteristic (a=%d, b=%d)', a, b);
        end
    end
    q(:,c) = M \ y;
end

viol = max([0; -q(:); q(:)-1; abs(sum(q,2)-1)]);
EXACT = viol <= feastol;

qc = min(max(q,0),1);
MMAP = cell(1,2+K);
MMAP{1} = D0;
MMAP{2} = D1;
for c = 1:K
    Dc = zeros(n,n);
    for jj = 1:z
        Dc(nzi(jj), nzj(jj)) = D1(nzi(jj), nzj(jj)) * qc(jj,c);
    end
    MMAP{2+c} = Dc;
end

end

function orders = marking_orders(n, z)
% (backward, forward) orders of an independent characteristic set
if n == 2
    orders = [1 0; 1 1; 2 0];
elseif n == 3
    orders = [1 0; 1 1; 2 0; 3 0];
else
    orders = [1 0; 1 1];
    a = 2;
    while size(orders,1) < n+1
        orders = [orders; a 0]; %#ok<AGROW>
        a = a + 1;
    end
end
orders = orders(1:z,:);
end
