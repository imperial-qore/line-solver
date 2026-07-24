function W=polling_qsys_decrementing(arvMAPs,svcMAPs,switchMAPs)
% W=polling_qsys_decrementing(arvMAPs,svcMAPs,switchMAPs)
%
% Exact mean waiting time solution of a symmetric polling system with open
% (Poisson) arrivals and decrementing (semiexhaustive) service. The server
% serves a queue until the number of jobs present drops to one less than the
% number found at the polling instant.
%
% Pittel (1973); Takagi (1984). See Takagi, ACM Computing Surveys, Vol. 20,
% No. 1, March 1988, eq (28). No exact closed form for the individual E[W_i]
% is known for asymmetric decrementing systems, so this analysis is restricted
% to the symmetric case.
%
% Example:
% W=polling_qsys_decrementing({map_exponential(1/0.4),map_exponential(1/0.4)},{map_exponential(1),map_exponential(1)},{map_exponential(0.1),map_exponential(0.1)})

n = length(arvMAPs); % number of classes
lambda = zeros(1,n); b = zeros(1,n); b2 = zeros(1,n);
r1 = zeros(1,n); delta2 = zeros(1,n);
for i=1:n
    lambda(i) = map_lambda(arvMAPs{i});
    b(i) = map_mean(svcMAPs{i});
    b2(i) = map_moment(svcMAPs{i},2);
    r1(i) = map_mean(switchMAPs{i});
    delta2(i) = map_var(switchMAPs{i});
end

% The exact result of Takagi (1984) requires a symmetric system.
tol = 1e-6;
for i=2:n
    if reldiff(lambda(i),lambda(1)) > tol || reldiff(b(i),b(1)) > tol || ...
            reldiff(b2(i),b2(1)) > tol || reldiff(r1(i),r1(1)) > tol || ...
            reldiff(delta2(i),delta2(1)) > tol
        line_error(mfilename,['MVA analysis for decrementing polling is only ' ...
            'available for symmetric systems (identical arrival, service and ' ...
            'switchover parameters across all queues).']);
    end
end

N = n;
lam = lambda(1);
b2s = b2(1);
r = r1(1);
d2 = delta2(1);
rho = N*lam*b(1);

denom = 2*(1 - rho - lam*r*(N - rho));
if denom <= 0
    line_error(mfilename,'Decrementing polling system is unstable: rho + lambda*r*(N-rho) >= 1.');
end

if r > 0
    residualSwitchover = d2/(2*r);
else
    residualSwitchover = 0;
end

Wval = residualSwitchover + ...
    (N*lam*b2s*(1 - lam*r) + (r + lam*d2)*(N - rho)) / denom;

W = Wval*ones(1,n);
end

function d=reldiff(a,b)
scale = max(abs(a),abs(b));
if scale < 1e-30
    d = 0;
else
    d = abs(a-b)/scale;
end
end
