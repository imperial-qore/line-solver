%{
%{
 % @file pfqn_sib.m
 % @brief Srinivasan (1985) Successively Improving Bounds (SIB) on cycle time
 %        and throughput for single-class product-form closed networks.
%}
%}

function [Xlo,Xhi,Wlo,Whi] = pfqn_sib(L,N,Z,level)
%{
%{
 % @brief Successively Improving Bounds (Srinivasan, "Successively Improving
 %        Bounds on Performance Measures for Product Form Queueing Networks",
 %        IEEE ToC 1987 / TR 85-2). Closed-form hierarchy of upper/lower bounds
 %        on the cycle time W(N) and throughput X(N) of a single-class closed
 %        network of fixed-rate (and delay) stations, based on the MVA relation
 %        W(N)=L*(1+phi(N-1)) with phi(K)=sum_m rho_m Q_m(K). Level 1 is the
 %        closed form of Thm 2.1; higher levels use the S_i power sums
 %        (S_i=sum_m rho_m^i) via Thms 3.5 (upper) and 3.6 (lower), tightening
 %        monotonically toward exact. Bounds are always at least as tight as the
 %        Balanced Job Bounds. O(M) to compute; level n needs S_2..S_{n+2}.
 % @fn pfqn_sib(L, N, Z, level)
 % @param L Fixed-rate service demand vector (M x 1). Delay demand goes in Z.
 % @param N Total population (scalar, N>=2).
 % @param Z Think time (scalar; added to the cycle-time bracket).
 % @param level Bound level (>=1, default 3). Higher = tighter, more S_i terms.
 % @return Xlo Lower bound on throughput X(N).
 % @return Xhi Upper bound on throughput X(N).
 % @return Wlo Lower bound on cycle time (residence + think) W(N).
 % @return Whi Upper bound on cycle time W(N).
%}
%}

L = L(:);
if nargin < 3 || isempty(Z), Z = 0; end
Z = sum(Z(:));
if nargin < 4 || isempty(level), level = 3; end
level = max(1, round(level));

% see _kb/03-api-layer.md (pfqn/ family: scaling, log-domain switches, dispatch gates)
if Z > 0
    error('pfqn_sib:delayUnsupported', ...
        'pfqn_sib supports Z=0 only (delay needs the Section-3.2 demand substitution, not yet implemented).');
end

Lsum = sum(L);
rho = L/Lsum;
rho_u = max(rho);
% Power sums S_i, i=1..level+3 (S_1=1).
imax = level+3;
S = zeros(1,imax);
for i = 1:imax
    S(i) = sum(rho.^i);
end
S2 = S(2);

% alpha_i coefficients (eq. 3.5-3.6), 0-indexed: alpha(k+1) = alpha_k.
alpha = zeros(1,level+1);
alpha(1) = S2;                       % alpha_0
for i = 1:level
    acc = 0;
    for j = 0:i-1
        acc = acc + S(i+1-j) * alpha(j+1);
    end
    alpha(i+1) = S(i+2) - acc;       % alpha_i
end

    function v = phi_u1(K)
        % Level-1 upper bound on phi(K) (Thm 3.5 with n=1 / eq. 3.18).
        if K <= 0
            v = 0;
        elseif K == 1
            v = S2;                  % phi(1) = sum rho_m^2 = S_2 (exact)
        else
            etaK = (K-1)/K;   % was eta: this is a NESTED function, so assigning
                              % eta clobbered the parent's eta=(N-2)/(N-1) used
                              % in the Theorem-3.5 prefactor 0.5/eta below, and
                              % the SIB hierarchy went inert above level 1
                              % (levels 2..4 all returned the level-0 baseline,
                              % with level 2 LOOSER than level 1).
            T1 = (K-1)*rho_u - 1;
            v = 0.5/etaK * (T1 + sqrt(T1^2 + 4*(K-1)*S2));
        end
    end

    function s = sigma(NN,i)
        % eq. (3.22c): NN plays the role of (N-1); Dbar_{N-3}=1+phi_u1(N-3).
        s = 0;
        if i <= 0, return; end
        Dbar = 1 + phi_u1(NN-2);     % N-3 = (N-1)-2 = NN-2
        pnum = 1;
        for j = 1:i
            pnum = pnum * (NN-1-(j-1));          % prod_{m=0}^{j-1}(N-2-m)
            s = s + (rho_u*S(j+1) - S(j+2)) * pnum / Dbar^j;
        end
    end

    function b = betaL(NN,i)
        % eq. (3.23c): NN plays the role of (N-1).
        b = 0;
        if i <= 0, return; end
        % term1: sum_{j=1}^{i-1} alpha_j prod_{m=2}^{j}((N-1-m)/Dbar_{N-1-m})
        for j = 1:i-1
            p = 1;
            for m = 2:j
                p = p * (NN-m)/(1 + phi_u1(NN-m));
            end
            b = b + alpha(j+1)*p;
        end
        % term2: alpha_i prod_{m=2}^{i}(...)(1 + (alpha_i/alpha_{i-1})(N-i-2)/(1+(N-i-2)alpha_0))
        p = 1;
        for m = 2:i
            p = p * (NN-m)/(1 + phi_u1(NN-m));
        end
        Nim2 = NN-1-i-1;             % N-i-2 = (NN+1)-i-2 = NN-i-1
        corr = 1 + (alpha(i+1)/alpha(i)) * Nim2/(1 + Nim2*alpha(1));
        b = b + alpha(i+1)*p*corr;
    end

% phi(N-1): NN = N-1.
NN = N-1;
% Section-2 baseline bounds (always valid).
phi_lo = (N-1)*S2;
T1s2 = (N-1)*rho_u - 1;
phi_hi = 0.5*(T1s2 + sqrt(T1s2^2 + 4*(N-1)*S2));

if N >= 3
    % Section-3 level-n upper (Thm 3.5) and lower (Thm 3.6); take tightest valid.
    eta = (N-2)/(N-1);
    T1u = (N-2)*rho_u - 1;
    su = sigma(NN, level-1);
    phi_u_n = 0.5/eta * (T1u + sqrt(max(0, T1u^2 + 4*(N-2)*(S2 - su))));
    phi_hi = min(phi_hi, phi_u_n);

    T1l = (N-2)*S2 - 1;
    bl = betaL(NN, level-1);
    phi_l_n = (T1l + sqrt(max(0, T1l^2 + 4*(N-2)*(S2 + (N-2)*bl))))/(2*level);
    phi_lo = max(phi_lo, phi_l_n);
end

% Guard validity of the bracket.
phi_lo = max(0, phi_lo);
if phi_hi < phi_lo
    phi_hi = phi_lo;
end

% Cycle time and throughput. Delay Z adds to the cycle time.
Wlo = Lsum*(1 + phi_lo) + Z;
Whi = Lsum*(1 + phi_hi) + Z;
Xlo = N/Whi;      % larger cycle time -> lower throughput
Xhi = N/Wlo;
end
