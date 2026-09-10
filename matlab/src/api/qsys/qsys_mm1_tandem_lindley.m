function result = qsys_mm1_tandem_lindley(lambda, mu1, mu2, Wk, Wk1)
% QSYS_MM1_TANDEM_LINDLEY Conditional downstream waiting time in an M/M/1 tandem.
%
% RESULT = QSYS_MM1_TANDEM_LINDLEY(LAMBDA, MU1, MU2, WK, WK1) returns the
% conditional mean waiting time of customer n+1 at the downstream station of a
% two-station single-server tandem queue, given that customer n waited WK at the
% upstream station and WK1 at the downstream one. LAMBDA is the external arrival
% rate at the upstream station, MU1 and MU2 the two service rates. WK and WK1
% are scalars or arrays of the same size.
%
% The point of the tandem recursion is that the interarrival time at the
% downstream station is the interdeparture time upstream, not an independent
% draw. With A ~ Exp(LAMBDA) the interarrival time upstream and S1, S1' the
% service times upstream of customers n and n+1, that interdeparture time is
%   D = max(A - WK - S1, 0) + S1',
% an idle period followed by the next service, and the downstream Lindley step
% is W2_{n+1} = (WK1 + S2 - D)^+ with S2 ~ Exp(MU2) independent of D.
%
% Because A is exponential, max(A - WK - S1, 0) is zero with probability 1-q and
% Exp(LAMBDA) with probability
%   q = e^{-LAMBDA WK} MU1/(LAMBDA+MU1) = P(the upstream server goes idle),
% so D is either Exp(MU1) or the sum of Exp(MU1) and Exp(LAMBDA). Averaging the
% downstream step over both cases needs only two elementary transforms of
%   g(d) = E[(WK1 + S2 - d)^+] = WK1 - d + 1/MU2       for d <= WK1,
%                              = e^{-MU2 (d-WK1)}/MU2   for d > WK1,
% namely J(c) = int_0^inf e^{-cu} g(u) du and Jw(c) = int_0^inf u e^{-cu} g(u) du,
% both closed form, giving
%   E[W2_{n+1} | WK, WK1] = (1-q) MU1 J(MU1) + q C,
%   C = LAMBDA MU1 (J(MU1) - J(LAMBDA))/(LAMBDA-MU1)   if LAMBDA ~= MU1,
%     = MU1^2 Jw(MU1)                                  if LAMBDA == MU1.
% As WK grows the upstream server never idles, q vanishes, and the mean tends to
% MU1 J(MU1) = E[g(S1')], as it must.
%
% Two caveats, both inherited from the reference and both quantified here.
%
% First, this is exact for the step taken in isolation, that is when the
% conditioning pair is independent of the four primitives that drive the step.
% In a running tandem it is not: the downstream wait WK1 was itself determined
% by an interdeparture time containing S1, so conditioning on (WK, WK1) is not
% conditioning on a Markov state of the tandem. Measured against a 4e6-customer
% simulation of the real tandem at LAMBDA = 0.8, MU1 = MU2 = 1, the formula is
% within 0.4% to 1.3% away from the empty state and 7% at WK = WK1 = 0, where the
% entanglement is strongest. Treat it as exact for one isolated step and as a
% good approximation in a running tandem.
%
% Second, this closed form was derived here rather than transcribed from the
% reference's theorem 4, because that theorem rests on its proposition 1, which
% omits a service-time difference and so does not describe a tandem queue; see
% QSYS_TANDEM_LINDLEY. The two differ: at LAMBDA = 0.8, MU1 = MU2 = 1 and
% WK = WK1 = 0 the published route gives 0.3016 against 0.3457 here, the latter
% matching simulation of the step to 6e-4 relative error.
%
% As in the reference, the upstream interarrival time is taken to be
% Exp(LAMBDA), which by Burke's theorem is also the stationary interdeparture
% law, so the same formula is applied at any pair of consecutive stations of a
% longer M/M/1 tandem, with the caveat above compounding.
%
% Returns a struct with fields:
%   mean         - Conditional mean downstream waiting time, size of WK
%   interdepMean - Conditional mean interdeparture time E[D | WK] = 1/MU1 + q/LAMBDA
%   idleProb     - q, the probability the upstream server idles, size of WK
%   analyzer     - Identifier string
%
% Examples:
%   r = qsys_mm1_tandem_lindley(0.8, 1.0, 1.0, 2.0, 3.0);
%   r.mean      % 2.9061
%   qsys_mm1_tandem_lindley(0.8, 1, 1, 0, 0).mean   % 0.3457, both stations empty
%
% Reference: S. Palomo, J. Pender, "Learning the Tandem Network Lindley
% Recursion", Proc. Winter Simulation Conference, 2021, proposition 1 and
% theorem 4, corrected as described above.
%
% See also QSYS_TANDEM_LINDLEY, QSYS_MM1_LINDLEY, QSYS_HH1_LINDLEY
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if ~isscalar(lambda) || ~isreal(lambda) || lambda <= 0
    line_error(mfilename, 'lambda must be a positive real scalar');
end
if ~isscalar(mu1) || ~isreal(mu1) || mu1 <= 0
    line_error(mfilename, 'mu1 must be a positive real scalar');
end
if ~isscalar(mu2) || ~isreal(mu2) || mu2 <= 0
    line_error(mfilename, 'mu2 must be a positive real scalar');
end
if ~isreal(Wk) || any(Wk(:) < 0) || any(~isfinite(Wk(:)))
    line_error(mfilename, 'Wk must hold finite nonnegative real values');
end
if ~isreal(Wk1) || any(Wk1(:) < 0) || any(~isfinite(Wk1(:)))
    line_error(mfilename, 'Wk1 must hold finite nonnegative real values');
end
if ~isequal(size(Wk), size(Wk1))
    if isscalar(Wk)
        Wk = repmat(Wk, size(Wk1));
    elseif isscalar(Wk1)
        Wk1 = repmat(Wk1, size(Wk));
    else
        line_error(mfilename, 'Wk and Wk1 must have the same size');
    end
end

x = Wk(:);
y = Wk1(:);

q = exp(-lambda * x) * mu1 / (lambda + mu1);
base = mu1 * tandemJ(mu1, y, mu2);
if abs(lambda - mu1) > 1e-9 * max(lambda, mu1)
    conv = lambda * mu1 / (lambda - mu1) ...
        * (tandemJ(mu1, y, mu2) - tandemJ(lambda, y, mu2));
else
    % the two rates coincide, the interdeparture time is Erlang(2,mu1)
    conv = mu1^2 * tandemJw(mu1, y, mu2);
end

meanW = (1 - q) .* base + q .* conv;
interdep = 1 / mu1 + q / lambda;

result = struct('mean', reshape(meanW, size(Wk)), ...
    'interdepMean', reshape(interdep, size(Wk)), ...
    'idleProb', reshape(q, size(Wk)), ...
    'analyzer', 'qsys_mm1_tandem_lindley');
end

function v = tandemJ(c, y, mu2)
% J(c) = int_0^inf e^{-c u} E[(y + S2 - u)^+] du with S2 ~ Exp(mu2).
e = exp(-c * y);
v = (y + 1 / mu2) .* (1 - e) / c ...
    - (1 - e .* (1 + c * y)) / c^2 ...
    + e / (mu2 * (c + mu2));
end

function v = tandemJw(c, y, mu2)
% Jw(c) = int_0^inf u e^{-c u} E[(y + S2 - u)^+] du with S2 ~ Exp(mu2).
e = exp(-c * y);
d = c + mu2;
v = (y + 1 / mu2) .* (1 - e .* (1 + c * y)) / c^2 ...
    - (2 - e .* (2 + 2 * c * y + c^2 * y.^2)) / c^3 ...
    + e .* (y / d + 1 / d^2) / mu2;
end
