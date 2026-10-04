function PB = me_gegecn_pb(p, K, N, c, Cs, Ca)
%ME_GEGECN_PB Blocking probability seen by one arrival stream of a censored
%GE/GE/c/K;N queue
%
% Evaluates equation (4.3) of Kouvatsos (1994) on the ME queue length
% distribution P returned by ME_GEGECN. Because a GE arrival process is a
% batch (bulk) process, an arrival can be blocked while the queue holds
% fewer than N jobs: the term (1-tau)^(N-n) is the probability that the
% batch overflows the residual room, and the extra factor of the first sum
% accounts for the servers that are still idle. With a Poisson stream
% (Ca = 1, tau = 1) every term but n = N vanishes and PB reduces to the
% PASTA value p(N).
%
% The stream scv CA is a per-stream quantity, so the same node solution P
% yields a different blocking probability for each of the flows merging
% into the queue: the external arrivals, the flow from each upstream
% station, and the flow released by each holding node. That is exactly how
% PBe_j, PB^i_j and PB^{h_ij}_j are obtained in the transfer-blocking
% algorithm of Tahilramani, Manjunath and Bose (1999).
%
% INPUTS:
%   p  - ME queue length distribution, p(idx) = Pr{n = K+idx-1}
%   K  - Minimum number of jobs in the queue
%   N  - Buffer capacity in jobs
%   c  - Number of servers
%   Cs - Squared coefficient of variation of the service times
%   Ca - Squared coefficient of variation of the interarrival times OF THE
%        STREAM whose blocking probability is requested
%
% OUTPUT:
%   PB - Probability that an arrival of this stream finds the queue full
%
% References:
%   D.D. Kouvatsos, "Entropy Maximisation and Queueing Network Models",
%   Annals of Operations Research, 48:63-126, 1994, equation (4.3).
%   H. Tahilramani, D. Manjunath, S.K. Bose, "Approximate analysis of open
%   network of GE/GE/m/N queues with transfer blocking", MASCOTS 1999,
%   164-171, equations (12) and (14).
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

tau = 2 / (Ca + 1);
sigma = 2 / (Cs + 1);
n = K:N;
nn = numel(n);
w = (1 - tau) .^ (N - n);

PB = 0;
% Jobs arriving while some servers are still idle: n = K,...,c-1
if K < c
    last = min(c - K, nn);
    idx = 1:last;
    fac = (sigma / (sigma * (1 - tau) + tau)) .^ (c - n(idx));
    PB = PB + sum(w(idx) .* fac .* p(idx));
end
% Jobs arriving with all servers busy: n = max(c,K),...,N
lo = max(c, K);
if lo <= N
    idx2 = (lo - K + 1):nn;
    PB = PB + sum(w(idx2) .* p(idx2));
end
end
