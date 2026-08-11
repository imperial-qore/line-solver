%{
%{
 % @file pfqn_cbh.m
 % @brief Convolutional Bound Hierarchy (Dowdy, Eager, Gordon, Saxton 1984)
 %        for single-class closed product-form networks.
%}
%}

function [Xlo,Xhi] = pfqn_cbh(L,N,Z,level)
%{
%{
 % @brief Level-`level` convolutional bound hierarchy on throughput. Fills
 %        column c = M-level of Buzen's g-array with BJB-derived estimates
 %        (e_0=1, e_i=e_{i-1}/BOUND(i)), then convolves the remaining `level`
 %        servers exactly (g(n,m)=g(n,m-1)+L_m g(n-1,m)). Bounds tighten
 %        monotonically with level and equal the exact solution at level M
 %        (Dowdy et al. 1984). BJB upper fill yields the upper bound; BJB
 %        lower fill the lower bound.
 % @fn pfqn_cbh(L, N, Z, level)
 % @param L Service demand vector (M x 1).
 % @param N Population (scalar).
 % @param Z Think time (scalar, default 0); folded as an extra IS column.
 % @param level Number of exactly-convolved servers, 1..M (default 2).
 % @return Xlo Lower throughput bound.
 % @return Xhi Upper throughput bound.
%}
%}
L = L(:);
if nargin < 3 || isempty(Z), Z = 0; end
if nargin < 4 || isempty(level), level = 2; end
M = numel(L);
level = max(1, min(level, M));
c = max(1, M - level);   % c=1 (single-server column) is already exact
Xlo = cbh_hier(L,N,Z,c,'lower');
Xhi = cbh_hier(L,N,Z,c,'upper');
end

function X = cbh_hier(L,N,Z,c,side)
% BJB-filled column c, then exact convolution of the remaining queueing
% servers and (separately, always exactly) the IS think-time station.
M = numel(L);
Rc = sum(L(1:c)); Lbc = max(L(1:c)); Lac = mean(L(1:c));
% Column-c fill from a BJB estimate of the first c queueing servers. The
% think time is NOT folded here: the delay is an infinite-server station and
% is convolved exactly below, so folding Z into the BJB fill would corrupt
% both the bound and the exact-at-level-M property.
e = zeros(1, N+1); e(1) = 1;
for i = 1:N
    switch side
        case 'upper', B = i/(Rc + (i-1)*Lac);   % BJB upper fill -> upper bound
        case 'lower', B = i/(Rc + (i-1)*Lbc);   % BJB lower fill -> lower bound
    end
    e(i+1) = e(i)/B;
end
if c == 1
    e = L(1).^(0:N);                             % single-server column is exact
end
g = e;
for m = c+1:M                                    % exact convolution of servers c+1..M
    for n = 1:N
        g(n+1) = g(n+1) + L(m)*g(n);
    end
end
if Z > 0                                         % convolve the IS delay exactly: g_Z(j)=Z^j/j!
    gd = Z.^(0:N) ./ factorial(0:N);
    gfull = zeros(1, N+1);
    for n = 0:N
        acc = 0;
        for j = 0:n
            acc = acc + g(j+1)*gd(n-j+1);
        end
        gfull(n+1) = acc;
    end
    g = gfull;
end
X = g(N)/g(N+1);
end
