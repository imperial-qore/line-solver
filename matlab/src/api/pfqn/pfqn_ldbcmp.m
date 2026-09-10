%{
%{
 % @file pfqn_ldbcmp.m
 % @brief Anselmi-Cremonesi (2008) lower throughput bound for closed
 %        single-class BCMP networks with load-dependent stations.
%}
%}

function [Xlo,Rhi,Qhat] = pfqn_ldbcmp(L,N,Z,c,varargin)
%{
%{
 % @brief Lower bound on system throughput (and upper bound on response time)
 %        for a closed, single-class BCMP network with load-dependent stations,
 %        via the asymptotic closed-open equivalence of Anselmi and Cremonesi,
 %        "Bounding the Performance of BCMP Networks with Load-Dependent
 %        Stations" (2008). The bound (their eq. 15) exploits the monotonicity
 %        of system throughput and the fact that a closed BCMP network is, in
 %        the limit N -> inf, equivalent to the open network obtained by
 %        removing the bottleneck and injecting arrivals at rate 1/D_max. It is
 %        applicable when N >= Qhat and is asymptotically exact; Algorithm 1
 %        refines it to a monotone fixed point.
 % @fn pfqn_ldbcmp(L, N, Z, c, tol)
 % @param L Fixed-rate (limiting) service demand vector (M x 1). For a Heffes
 %          LD station, L(i) is the limiting demand D_i = lim_n D_i(n).
 % @param N Total population (scalar).
 % @param Z Think time (scalar; modeled as a non-bottleneck delay station).
 % @param c Per-station Heffes load-dependence coefficient (M x 1, default 0).
 %          c(i)=0 marks a fixed-rate (LI) station with open queue
 %          rho_i/(1-rho_i); c(i)>0 a Heffes LD station with open queue
 %          (c(i)+1)*rho_i/(1-rho_i) (their eq. 22-23). The bottleneck is
 %          assumed fixed-rate (population transform (7) reduces to N'=N).
 % @param varargin Optional trailing arguments, in order: tol, the fixed-point
 %          tolerance for Algorithm 1 (default 1e-10).
 % @return Xlo Lower bound on system throughput X(N); NaN if N < Qhat.
 % @return Rhi Upper bound on system response+think time, N/Xlo (Little).
 % @return Qhat Sum of non-bottleneck limiting queue lengths (eq. 11).
%}
%}

L = L(:);
M = numel(L);
if nargin < 3 || isempty(Z)
    Z = 0;
end
Z = sum(Z(:));
if nargin < 4 || isempty(c)
    c = zeros(M,1);
end
c = c(:);
tol = 1e-10;
if numel(varargin) >= 1 && ~isempty(varargin{1})
    tol = varargin{1};
end

% Limiting effective demands and bottleneck (eqs. 3-5).
Dstar = L;
Dm = max(Dstar);
isbott = abs(Dstar - Dm) <= 1e-12*Dm;
bmax = sum(isbott);

% Qhat: sum of non-bottleneck queue lengths in the equivalent open BCMP
% network at arrival rate lambda = 1/Dm (eqs. 11, 17, 22-23).
lambda = 1/Dm;
Qhat = 0;
for i = 1:M
    if isbott(i)
        continue
    end
    rho_i = lambda*Dstar(i);
    if rho_i >= 1
        Xlo = NaN; Rhi = NaN; return
    end
    Qhat = Qhat + (c(i)+1)*rho_i/(1-rho_i);
end
% Delay station (infinite server): open queue lambda*Z.
Qhat = Qhat + lambda*Z;

% Applicability (eqs. 8-10 require N - Qhat >= 0).
if N - Qhat < 0
    Xlo = NaN;
    Rhi = NaN;
    return
end

% Algorithm 1: iterate the lower bound eq. (15) to a monotone fixed point.
a = N - Qhat;
Xprime = 0;
Xlo = 0;
for it = 1:10000
    Xprev = Xlo;
    denom = Dm*(bmax + N - Qhat) - bmax*(Dm*Xprime)^N*Dm;
    Xlo = a/denom;
    Xprime = Xlo;
    if Xprev > 0 && abs(Xprev - Xlo)/Xprev <= tol
        break
    end
end

Rhi = N/Xlo;   % response+think upper bound (Little, eq. 16)
end
