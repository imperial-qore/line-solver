%{
%{
 % @file pfqn_ssd.m
 % @brief Server-Station Disaggregation throughput bounds (Suri-Dallery 1986)
 %        for single-class closed networks with multiserver stations.
%}
%}

function [Xlo,Xhi] = pfqn_ssd(L,N,Z,nservers)
%{
%{
 % @brief SSD multiserver bounds (SIGMETRICS 1986, Theorem 5). Each C_k-server
 %        station of loading L_k is bracketed by disaggregations: C_k balanced
 %        single-server stations of loading L_k/C_k (lower) and one single
 %        server of loading L_k/C_k (upper). With R_l=sum L_k, Y_l=max L_k/C_k,
 %        R_u=sum L_k/C_k, Y_u=R_u/K:
 %          X_l = N/(R_l+(N-1)Y_l) <= X(N) <= N/(R_u+(N-1)Y_u) = X_u,
 %        the upper bound taken jointly with the ABA bound min(N/R_l, C_b/L_b).
 %        O(K) cost, same order as BJB on single-server networks.
 %        With Z>0 the queueing terms carry the terminal-workload correction of
 %        Lazowska et al. 1984, Table 5.2: (N-1)Y_l/(1+Z/(N R_l)) on the lower
 %        bound and (N-1)Y_u/(1+Z/R_u) on the upper. Adding Z without it does
 %        not yield a bound.
 % @fn pfqn_ssd(L, N, Z, nservers)
 % @param L Service demand vector (M x 1).
 % @param N Population (scalar).
 % @param Z Think time (scalar, default 0).
 % @param nservers Per-station server counts C_k (M x 1, default all 1).
 % @return Xlo Lower throughput bound (Theorem 5).
 % @return Xhi Upper throughput bound (Theorem 5, joint with ABA).
%}
%}
L = L(:);
K = numel(L);
if nargin < 3 || isempty(Z), Z = 0; end
if nargin < 4 || isempty(nservers), nservers = ones(K,1); end
C = nservers(:);
if isscalar(C), C = C*ones(K,1); end

Rl = sum(L);        Yl = max(L./C);
Ru = sum(L./C);     Yu = Ru/K;
[~, b] = max(L./C);

Xlo = N/(Rl + Z + (N-1)*Yl/(1 + Z/(N*Rl)));
Xhi = min([ N/(Ru + Z + (N-1)*Yu/(1 + Z/Ru)), ...     % Theorem 5 upper
            C(b)/L(b), ...                 % ABA capacity bound (eq. 3)
            N/(Rl + Z) ]);                 % ABA population bound
end
