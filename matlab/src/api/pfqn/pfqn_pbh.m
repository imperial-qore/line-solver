%{
%{
 % @file pfqn_pbh.m
 % @brief Performance Bound Hierarchy (Eager-Sevcik 1983) for single-class
 %        closed product-form networks.
%}
%}

function [Xlo,Xhi,Qlo,Qhi] = pfqn_pbh(L,N,Z,level)
%{
%{
 % @brief Level-`level` Performance Bound Hierarchy throughput/queue bounds.
 %        i MVA steps from an ABA-initialized residence give nested
 %        optimistic/pessimistic bounds converging to exact MVA as level->N
 %        (Eager-Sevcik 1983, ACM TOCS 1(2):99-115, eqs. 5-13). Level 1 with
 %        Z=0 equals the BJB optimistic bound; level 2 with Z=0 equals the
 %        closed form 1 + S(N-1), S = sum L_k^2.
 % @fn pfqn_pbh(L, N, Z, level)
 % @param L Service demand vector (M x 1).
 % @param N Population (scalar).
 % @param Z Think time (scalar, default 0).
 % @param level Hierarchy level >= 0 (default 1); clamped to N.
 % @return Xlo Lower throughput bound (pessimistic residence).
 % @return Xhi Upper throughput bound (optimistic residence, joint with the
 %             asymptotic bound and 1/max(L)).
 % @return Qlo Per-station lower queue-length bound (M x 1), Little bracket.
 % @return Qhi Per-station upper queue-length bound (M x 1), Little bracket.
%}
%}
L = L(:);
if nargin < 3 || isempty(Z), Z = 0; end
if nargin < 4 || isempty(level), level = 1; end
Lmax = max(L);

Ro = pbh_residence(L,N,Z,level,'opt');    % per-station optimistic residence
Rp = pbh_residence(L,N,Z,level,'pess');   % per-station pessimistic residence

% joint with the asymptotic residence lower bound: R(N) >= max(N*Lmax-Z, R(1)),
% R(1)=sum(L) (the Eager-Sevcik model is normalized to sum(L)=1).
RoC = max(sum(Ro), max(N*Lmax - Z, sum(L)));
Xhi = min(1/Lmax, N/(Z + RoC));
Xlo = N/(Z + sum(Rp));

% Little's law bracket: Q_k = X*R_k, with X in [Xlo,Xhi], R_k in [Ro_k,Rp_k].
Qlo = Xlo * Ro;
Qhi = Xhi * Rp;
end

function Rk = pbh_residence(L,N,Z,level,side)
% Per-station residence-time vector of the level-`level` PBH bound.
K = numel(L);
[~, b] = max(L);
level = min(level, N);
n0 = N - level;
if strcmp(side,'pess')
    % Pessimistic start, carried in QUEUE LENGTHS: all n0 customers at the
    % bottleneck (eq 13). Seeding a residence instead makes the assumed
    % population n0*Rb/(Z+Rtot) < n0 once Z > 0, so the pessimism is diluted
    % and the resulting Xlo stops being a bound (violated exact on 14% of
    % random delay models, worst 43%).
    Q = zeros(K,1);
    Q(b) = n0;
    Rk = L .* (1 + Q);
    for n = n0+1:N
        Rk = L .* (1 + Q);                               % eq (5), MVA step
        Q = (n / (Z + sum(Rk))) * Rk;
    end
    return
end
Rk = ones(K,1) * max(n0*L(b) - Z, sum(L))/K;              % eq (7), ABA opt
if n0 == 0
    Rk = zeros(K,1);                                      % exact empty base
end
for n = n0+1:N
    Rtot = sum(Rk);
    if n == 1 || Z + Rtot == 0
        Rk = L;                                          % empty-network base step
    else
        Rk = L .* (1 + (n-1) * Rk / (Z + Rtot));         % eq (5)
    end
end
end
