function [W,rho]=qsys_mmk(lambda,mu,k)
% W=QSYS_MMK(LAMBDA,MU,K)
rho = lambda/mu/k;
Q = rho/(1-rho)*ErlangC(k, rho) + k * rho;
W = Q / lambda;
end

function C = ErlangC(k, rho)
% Erlang-C formula
% The probability that an arriving customer is forced to join the queue
% (all servers are occupied)
S = 0;
for j=0:(k-1)
    S = S + (k*rho)^j/factorial(j);
end
C = 1 / (1 + (1-rho)*factorial(k)/(k*rho)^k*S);
end
