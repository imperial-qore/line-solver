%{ @file fj_ism_green.m
 %  @brief Green's independent server model of simultaneous server requests
 %
 %  @author LINE Development Team
%}

%{
 % @brief Green's independent server model of simultaneous server requests
 %
 % @details
 % In the independent server model a customer requires j servers simultaneously
 % to begin service, with probability c(j), and then releases them
 % asynchronously as each of its j tasks completes at rate mu. Servers can sit
 % idle while customers wait, which is what separates the model from M/G/s, and
 % customer service ends when the last of its tasks does, so the customer
 % service time is the maximum of j exponentials, H_j/mu.
 %
 % The cycle decomposition of Green splits time into a queueing period, during
 % which a queue exists, and a nonqueue period. Everything follows from two
 % quantities:
 %
 %   E[B] = sum_j c(j) * sum_{i=0..j-1} 1/((s-i)*mu),
 %
 % the interservice time, which is the j-th order statistic of s exponentials
 % because all s servers are busy whenever a customer enters service during a
 % queueing period, and
 %
 %   E[D] = sum_i sum_{k=1..i} [ sum_{m=0..k-1} 1/((i-m)*mu) ] * q(i)*c(s-i+k)/p_d,
 %
 % the initial delay of the customer that starts a queueing period, which finds
 % i servers busy and must wait for k of them to free. The busy-server
 % distribution q during a nonqueue period and its mean length come from the
 % embedded chain of arrivals and completions absorbed when a queue forms:
 %
 %   V = (I - T)^-1,   E[Qbar] = sum_j v(s,j)/(lambda + j*mu),
 %   q(i) = v(s,i)/((lambda + i*mu) * E[Qbar]),
 %   E[Q] = E[D]/(1 - lambda*E[B]),   p_q = E[Q]/(E[Q] + E[Qbar]).
 %
 % The waiting time transform of Eq. (61) factors into the equilibrium
 % transform of D and the Pollaczek-Khinchine transform of an M/G/1 queue with
 % service B, so its mean is
 %
 %   E[W] = (1 - pi0) * [ E[D^2]/(2*E[D]) + lambda*E[B^2]/(2*(1-lambda*E[B])) ],
 %   pi0  = (1 - lambda*E[B]) / (1 - lambda*(E[B] - E[D])).
 %
 % The second moments of B and D are exact, each stage of the order statistic
 % being an independent exponential.
 %
 % Note that Eq. (65) of the survey prints the inner sum of E[D] as starting at
 % 1/(s*mu) even though only i servers are busy; the sum is started at
 % 1/(i*mu) here, which is what the accompanying text prescribes and what makes
 % E[D] reduce to E[B] when i = s.
 %
 % @par Syntax:
 % @code
 % [W, R, out] = fj_ism_green(lambda, mu, s, c)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>lambda<td>Customer arrival rate
 % <tr><td>mu<td>Per-task service rate
 % <tr><td>s<td>Number of servers
 % <tr><td>c<td>Vector of s probabilities, c(j) = P(a customer needs j servers)
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>W<td>Mean waiting time before service starts
 % <tr><td>R<td>Mean response time, waiting plus the maximum of the j tasks
 % <tr><td>out<td>Struct with EB, EB2, ED, ED2, EQ, EQbar, pq, pd, pi0, rho, q, ES
 % </table>
 %
 % @par Reference:
 % A. Thomasian, "Analysis of Fork/Join and Related Queueing Systems",
 % ACM Computing Surveys, Vol. 47, No. 2, Article 17, July 2014, Section 6.4,
 % Eqs. (59)-(65).
 %
 % Original: L. Green, "A Queueing System in Which Customers Require a Random
 % Number of Servers", Operations Research 28(6), 1980.
%}
function [W, R, out] = fj_ism_green(lambda, mu, s, c)

c = c(:)';
if s < 1 || s ~= round(s)
    line_error(mfilename, 's must be a positive integer. Got s=%g.', s);
end
if numel(c) ~= s
    line_error(mfilename, 'c must have s=%d entries, one per server requirement. Got %d.', s, numel(c));
end
if any(c < 0)
    line_error(mfilename, 'The server-requirement probabilities must be non-negative.');
end
if abs(sum(c) - 1) > 1e-9
    line_error(mfilename, 'The server-requirement probabilities must sum to one. Got %g.', sum(c));
end
if lambda <= 0 || mu <= 0
    line_error(mfilename, 'lambda and mu must be positive. Got lambda=%g, mu=%g.', lambda, mu);
end

% Interservice time B: the j-th order statistic of s exponentials of rate mu,
% a sum of independent exponentials of rates (s-i)*mu for i = 0..j-1
EB = 0;
EB2 = 0;
for j = 1:s
    stages = 1 ./ (((s - (0:(j - 1)))) * mu);
    mj = sum(stages);
    vj = sum(stages.^2);
    EB = EB + c(j) * mj;
    EB2 = EB2 + c(j) * (vj + mj^2);
end

rho = lambda * EB;
if rho >= 1
    line_error(mfilename, 'System is unstable: rho = lambda*E[B] = %.4f >= 1.', rho);
end

% Embedded chain over the number of busy servers 0..s, absorbed when an
% arrival needs more servers than are free
T = zeros(s + 1, s + 1);
for i = 0:s
    tot = lambda + i * mu;
    if i > 0
        T(i + 1, i) = i * mu / tot;
    end
    for j = 1:(s - i)
        T(i + 1, i + j + 1) = T(i + 1, i + j + 1) + lambda * c(j) / tot;
    end
end

V = (eye(s + 1) - T) \ eye(s + 1);
% A nonqueue period starts with all s servers busy
v = V(s + 1, :);
holding = 1 ./ (lambda + (0:s) * mu);
EQbar = sum(v .* holding);
q = (v .* holding) / EQbar;

% A customer arriving during a nonqueue period is delayed when it needs more
% servers than the s-i free ones
pd = 0;
for i = 0:s
    free = s - i;
    if free < s
        pd = pd + q(i + 1) * sum(c((free + 1):s));
    end
end
if pd <= 0
    line_error(mfilename, 'No arrival can ever be delayed; the model degenerates to M/M/%d.', s);
end

% Initial delay D: i servers busy, the customer needs k = j-(s-i) more to free
ED = 0;
ED2 = 0;
for i = 1:s
    for k = 1:i
        j = s - i + k;
        if j < 1 || j > s
            continue
        end
        wgt = q(i + 1) * c(j) / pd;
        if wgt == 0
            continue
        end
        stages = 1 ./ (((i - (0:(k - 1)))) * mu);
        mk = sum(stages);
        vk = sum(stages.^2);
        ED = ED + wgt * mk;
        ED2 = ED2 + wgt * (vk + mk^2);
    end
end

EQ = ED / (1 - rho);
pq = EQ / (EQ + EQbar);
pi0 = (1 - rho) / (1 - lambda * (EB - ED));

% Eq. (61) factors into the equilibrium transform of D and the M/G/1 waiting
% time transform with service B, so the means add
Weq = ED2 / (2 * ED);
Wmg1 = lambda * EB2 / (2 * (1 - rho));
W = (1 - pi0) * (Weq + Wmg1);

% Customer service time: the maximum of the j tasks it holds
ES = 0;
for j = 1:s
    ES = ES + c(j) * fj_harmonic(j) / mu;
end
R = W + ES;

out = struct('EB', EB, 'EB2', EB2, 'ED', ED, 'ED2', ED2, 'EQ', EQ, ...
    'EQbar', EQbar, 'pq', pq, 'pd', pd, 'pi0', pi0, 'rho', rho, ...
    'q', q, 'ES', ES);

end
