%{
%{
 % @file pfqn_mcub.m
 % @brief Multiclass Composite Upper Bound (Kerola 1986) on per-class
 %        throughput for closed product-form networks.
%}
%}

function [Xub,Xlb] = pfqn_mcub(L,N,Z)
%{
%{
 % @brief Kerola's composite bound method (Perf. Eval. 6:1-9, eqs. 10-16).
 %        Given per-class multiclass BJB lower bounds X_s^-, the residual-
 %        utilization composite upper bound is
 %          X_r <= min_k [1 - sum_{s!=r} X_s^- L_ks] / L_kr,
 %        computed at O(KR). Much tighter than per-class ABA at moderate load.
 %        (Named pfqn_mcub because pfqn_cub is the unrelated cubature NC method.)
 % @fn pfqn_mcub(L, N, Z)
 % @param L Service demand matrix, station x class (M x R).
 % @param N Population vector (1 x R).
 % @param Z Think time vector (1 x R, default zeros).
 % @return Xub Composite upper throughput bound per class (1 x R).
 % @return Xlb Multiclass BJB lower throughput bound per class (1 x R, eq. 10).
%}
%}
[M,R] = size(L);
N = N(:)';
if nargin < 3 || isempty(Z), Z = zeros(1,R); end
Z = Z(:)';
Ntot = sum(N);
R0 = sum(L,1);                       % per-class total demand (1 x R)
Lb = max(L,[],1);                    % per-class bottleneck demand (1 x R)

% eq (10): multiclass Balanced Job Bounds lower throughput bound.
Xlb = N ./ (R0 + Z + (Ntot - 1).*Lb);

% eqs (13)-(16): composite upper bound per class.
Xub = zeros(1,R);
for r = 1:R
    Uoth = zeros(M,1);              % utilization at each device by other classes
    for s = 1:R
        if s ~= r
            Uoth = Uoth + Xlb(s) * L(:,s);
        end
    end
    Ucub = 1 - Uoth;               % residual utilization available to class r
    dev = inf(M,1);
    for k = 1:M
        if L(k,r) > 0
            dev(k) = Ucub(k) / L(k,r);
        end
    end
    Xub(r) = min(dev);
end
end
