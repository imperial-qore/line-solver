function result = qsys_tandem_ub_ciucu(x, lst, p, mu, dlst)
% QSYS_TANDEM_UB_CIUCU Tail bounds for a GI/Hn/1 -> ./Hn/1 tandem.
%
% RESULT = QSYS_TANDEM_UB_CIUCU(X, LST, P, MU) returns polynomial-exponential
% upper bounds on the tails of the end-to-end waiting time W and sojourn time S
% of a tandem of two single-server FCFS stations, fed by a renewal arrival
% process and serving the same hyperexponential law at both stations. X is the
% vector of thresholds at which the tails are bounded, LST is a function handle
% evaluating the interarrival Laplace-Stieltjes transform E[e^{-s X}] for s >= 0,
% and P, MU are the phase probabilities and rates of the service law
% Y, Z ~ sum_i P(i) Exp(MU(i)); a scalar P = 1 gives exponential service.
%
% RESULT = QSYS_TANDEM_UB_CIUCU(X, LST, P, MU, DLST) also takes a handle
% evaluating E[X e^{-s X}], i.e. minus the derivative of LST. Without it that
% derivative is obtained by a Richardson-extrapolated central difference, which
% costs about four extra LST evaluations and loses roughly four digits.
%
% The bounds are those of Ciucu and Mehri: with theta the positive root of
% E[e^{theta (Y-X)}] = 1 (Lundberg/Kingman exponent of the first station, which
% by identical service is also the exponent of the tandem) and
% alpha = E[X e^{-theta X}], the test function
%   gamma(u,v) = 1{0<=u<=v} [1 - A e^{-theta u} - (B + C u + D v) e^{-theta v}]
% is made to satisfy the integral inequality of their Theorem 1(b) by the five
% sufficient conditions of their Lemma 4, which fix
%   A = 1,  C = theta sum_i p_i/(mu_i-theta) / sum_i p_i mu_i/(mu_i-theta)^2,
%   D = (-C E[U e^{theta V}]/E[V e^{theta V}]) v 0,  U = Y-X, V = Z-X,
%   B = C (1/mu_1 - alpha E[e^{theta Z}])        if D = 0,
%     = (C+D)/(mu_1-theta) - theta/mu_1          if D > 0,
% with mu_1 the smallest service rate. Their Corollary 2 then gives
%   P(S > x) <= sum_i p_i { e^{-mu_i x}
%                 + mu_i/(mu_i-theta) (A+B) (e^{-theta x} - e^{-mu_i x})
%                 + mu_i/(mu_i-theta)^2 (C+D) (((mu_i-theta)x-1) e^{-theta x}
%                                              + e^{-mu_i x}) }
% and the corresponding closed form for W when the service is exponential.
% E[V e^{theta V}] is positive at any stable load, so D is always well defined:
% h(s) = E[e^{s(Z-X)}] is convex with h(0) = h(theta) = 1, hence h'(theta) > 0.
%
% The two exponentials mix a polynomial of degree one in x, which is what lets
% the bound follow the concave bend of the tail on a linear-log scale where a
% purely exponential bound cannot. In the M/M/1 -> ./M/1 case the five
% inequalities hold as equalities, so gamma is the exact joint distribution and
% both bounds are exact:
%   P(W > x) = (1 - 2 theta^2/(mu(mu+theta)) + x(mu-theta)theta/(mu+theta)) e^{-theta x}
%   P(S > x) = (1 + theta x) e^{-theta x} .
% Away from it the bound stays sharp: against an exact CTMC reference for the
% Erlang(2)/M/1 -> ./M/1 tandem it is within 2% at P(S>x) = 1e-2 and within
% 0.6% at 5e-10, with the correct asymptotic slope theta^2/(mu(1-alpha mu)).
% Accuracy degrades with service variability, to about a factor of two at
% CV(Y) = 2.
%
% Returns a struct with fields:
%   S      - upper bound on P(S > x), one entry per threshold, capped at 1
%   W      - upper bound on P(W > x), NaN unless the service is exponential
%   theta  - the tail decay rate, positive root of E[e^{theta (Y-X)}] = 1
%   alpha  - E[X e^{-theta X}]
%   A,B,C,D- the coefficients of gamma fixed by Lemma 4
%   analyzer - 'qsys_tandem_ub_ciucu'
%
% Example, a D/M/1 -> ./M/1 tandem at utilization 3/4 with unit service rate:
%   r = qsys_tandem_ub_ciucu([5 10], @(s) exp(-s*4/3), 1, 1, @(s) (4/3)*exp(-s*4/3));
%
% Reference: F. Ciucu, S. Mehri, "On the Distribution of Sojourn Times in Tandem
% Queues", Proc. ACM Meas. Anal. Comput. Syst. 9(2), Article 27, 2025 (ACM
% SIGMETRICS 2025). Registered in .citations() as 'tandemub'.

if nargin < 5, dlst = []; end
x = x(:)';
p = p(:)'; mu = mu(:)';
if numel(p) ~= numel(mu)
    line_error(mfilename, 'p and mu must have the same number of phases.');
end
if any(x < 0)
    line_error(mfilename, 'The thresholds x must be nonnegative.');
end
if any(p < 0) || abs(sum(p) - 1) > 1e-10
    line_error(mfilename, 'The phase probabilities p must be nonnegative and sum to one.');
end
if any(mu <= 0)
    line_error(mfilename, 'The service rates mu must be positive.');
end

mgfY = @(t) sum(p.*mu./(mu - t));                  % E[e^{t Y}], t < min(mu)
mu1 = min(mu);
% Stability: E[X] > E[Y] is what makes E[e^{t(Y-X)}] - 1 cross zero on (0,mu1).
res = @(t) mgfY(t).*lst(t) - 1;
hi = mu1*(1 - 1e-12);
if res(hi) <= 0
    line_error(mfilename, 'No positive root of E[e^{theta(Y-X)}]=1 below min(mu): the tandem is unstable or the service is not the lighter tail.');
end
lo = mu1*1e-12;
while res(lo) >= 0 && lo > mu1*1e-16
    lo = lo/10;                                    % walk below the root at zero
end
if res(lo) >= 0
    % E[e^{t(Y-X)}]-1 is convex and vanishes at t=0, so it stays positive on the
    % whole of (0,mu1) exactly when its slope E[Y]-E[X] there is nonnegative.
    line_error(mfilename, 'The tandem is unstable, E[X] <= E[Y]: theta = 0 is the only root of E[e^{theta(Y-X)}]=1.');
end
theta = fzero(res, [lo, hi], optimset('TolX', 1e-14));

if isempty(dlst)
    alpha = local_dlst(lst, theta);
else
    alpha = dlst(theta);
end

EexpZ = mgfY(theta);                               % E[e^{theta Z}]
EZexp = sum(p.*mu./(mu - theta).^2);               % E[Z e^{theta Z}]
EY = sum(p./mu);
A = 1;
C = theta*sum(p./(mu - theta))/EZexp;
EUeV = EY - alpha*EexpZ;                           % E[U e^{theta V}]
EVeV = EZexp/EexpZ - alpha*EexpZ;                  % E[V e^{theta V}] > 0
D = -C*EUeV/EVeV;
if ~(D > 0), D = 0; end
if D > 0
    B = (C + D)/(mu1 - theta) - theta/mu1;
else
    B = C*(1/mu1 - alpha*EexpZ);
end

S = zeros(size(x));
for i = 1:numel(p)
    m = mu(i);
    S = S + p(i)*(exp(-m*x) ...
        + m/(m - theta)*(A + B)*(exp(-theta*x) - exp(-m*x)) ...
        + m/(m - theta)^2*(C + D)*(((m - theta)*x - 1).*exp(-theta*x) + exp(-m*x)));
end
S = min(S, 1);

if isscalar(p)
    beta = lst(mu1);                               % E[e^{-mu X}]
    if D == 0
        W = (1 - 2*theta^2/(mu1*(mu1 + theta)) + theta*(mu1 - theta)/(mu1 + theta)*x).*exp(-theta*x) ...
            + beta*(theta*mu1*alpha/(2*(mu1 - theta)) - theta/(2*mu1))*exp(-mu1*x);
    else
        W = (1 - 2*theta/mu1 + 2*theta^2*(2 - alpha*mu1)/((mu1 + theta)^2*(1 - alpha*mu1)) ...
            + theta^2*(mu1 - theta)/(mu1*(mu1 + theta)*(1 - alpha*mu1))*x).*exp(-theta*x);
    end
    W = min(W, 1);
else
    W = nan(size(x));                              % (23) for W is Exp-service only
end

result = struct('S', S, 'W', W, 'theta', theta, 'alpha', alpha, ...
    'A', A, 'B', B, 'C', C, 'D', D, 'analyzer', mfilename);
end

function a = local_dlst(lst, s)
% Richardson-extrapolated central difference of -LST at s, i.e. E[X e^{-s X}].
h = 1e-3*(1 + s);
if h > s, h = s/2; end
d1 = (lst(s - h) - lst(s + h))/(2*h);
d2 = (lst(s - h/2) - lst(s + h/2))/h;
a = (4*d2 - d1)/3;
end
