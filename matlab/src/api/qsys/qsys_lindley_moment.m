function value = qsys_lindley_moment(lambda, mu, Wn, m)
% QSYS_LINDLEY_MOMENT One conditional Lindley moment for exponential primitives.
%
% VALUE = QSYS_LINDLEY_MOMENT(LAMBDA, MU, WN, M) returns
% E[max(WN + S - A, 0)^M] with A ~ Exp(LAMBDA) and S ~ Exp(MU), for the vector
% of current waiting times WN and the integer order M >= 1.
%
% This is the kernel shared by QSYS_MM1_LINDLEY, which calls it once per moment
% order, and QSYS_HH1_LINDLEY, which mixes it over the arrival and service
% phases. See QSYS_MM1_LINDLEY for the derivation and for why the upper
% incomplete gamma function reduces to a finite sum here.
%
% See also QSYS_MM1_LINDLEY, QSYS_HH1_LINDLEY
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

w = Wn(:);

% S = sum_{k=0}^{m} C(m,k) w^k (m-k)! / mu^(m-k+1)
sTerm = zeros(size(w));
for k = 0:m
    sTerm = sTerm + nchoosek(m, k) * w.^k * factorial(m - k) / mu^(m - k + 1);
end

% T = (-1)^m m! ( sum_{k=0}^{m} (-lambda w)^k/k! - e^{-lambda w} ) / lambda^(m+1)
x = -lambda * w;
inner = zeros(size(w));
for k = 0:m
    inner = inner + x.^k / factorial(k);
end
tTerm = (-1)^m * factorial(m) * (inner - exp(x)) / lambda^(m + 1);

value = lambda * mu / (lambda + mu) * (sTerm + tTerm);
end
