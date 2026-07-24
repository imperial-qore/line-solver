function I = map_count_idc(map,t)
% I = map_count_idc(MAP,t) - Index of dispersion for counts (IDC) of a MAP
% evaluated at the time points t.
%
% The IDC of the counting process A(t) associated to the MAP is
%
%   I_a(t) = Var(A(t)) / E[A(t)],  t > 0,
%
% i.e. the scaled variance-time curve. It is a function of t that interpolates
% between I_a(0+) = SCV of the interarrival time (for a renewal MAP) and the
% asymptotic value I_a(Inf) = map_idc(MAP). See W. Whitt and W. You, "A Robust
% Queueing Network Analyzer Based on Indices of Dispersion", eq. (1).
%
% Input:
% - map: Markovian Arrival Process in the form {D0,D1}
% - t:   vector of time points (t>0)
% Output:
% - I:   column vector of IDC values, one per element of t

t = t(:);
I = zeros(numel(t),1);
m = map_count_mean(map,t);
v = map_count_var(map,t);
nz = m > 0;
I(nz) = v(nz) ./ m(nz);
% For an orderly point process the counting IDC tends to 1 as t->0 (locally
% Poisson: Var(N(t)) ~ E[N(t)]). This limit is only hit at t==0 exactly.
if any(~nz)
    I(~nz) = 1;
end
end
