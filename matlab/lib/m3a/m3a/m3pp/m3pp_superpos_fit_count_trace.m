function [fit,m3pps] = m3pp_superpos_fit_count_trace(T, A, t, tinf)
% Superposes k M3PP to fit a multi-class trace with m classes.
% INPUT
% - T: inter-arrival time
% - A: labels
% - t: finite time scale
% - tinf: near-infinite time scale

% by default choose time scales as follows
if nargin < 3
    t = 10*mean(T);
    tinf = max(10*t, (sum(T)-T(1)) / 100);
end

% rate
a = 1/mean(T);

% per-class rates
L = unique(A);
m = length(L);
pv = zeros(m,1);
for j = 1:m
    pv(j) = sum(A==L(j))/length(A);
end
av = pv * a;

% compute counting process at resolution t1 and tinf
Nt = mtrace_iat2counts(T, A, t);
Ninf = mtrace_iat2counts(T, A, tinf);

% compute IDC(t) e IDC(inf)
btv = zeros(m,1);
binfv = zeros(m,1);
m3tv = zeros(m,1);
for i = 1:m
    btv(i) = var(Nt(:,i))/(av(i)*t);
    binfv(i) = var(Ninf(:,i))/(av(i)*tinf);
    mt = [mean(Nt(:,i)), mean(Nt(:,i).^2), mean(Nt(:,i).^3)];
    m3tv(i) = mt(3) - 3*mt(2)*mt(1) + 2*mt(1)^3;
end

% fit superposition
[fit, m3pps] = m3pp_superpos_fit_count(av, btv, binfv, m3tv, t, tinf);
end