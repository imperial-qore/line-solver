function [pi_t, t] = solver_ctmc_chain_transient(chain, pi0, timespan)
% [PI_T, T] = SOLVER_CTMC_CHAIN_TRANSIENT(CHAIN, PI0, TIMESPAN)
%
% Transient distribution of a user-supplied Markov chain over TIMESPAN =
% [t0,t1]. For a MarkovProcess the Kolmogorov forward equations are integrated
% from PI0; for a MarkovChain the distribution is advanced one step per unit
% of time, so the returned T holds the integer steps within TIMESPAN.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if isa(chain,'MarkovChain')
    n = size(chain.getTransMat(),1);
elseif isa(chain,'MarkovProcess')
    n = size(chain.getGenerator(),1);
else
    line_error(mfilename,'solver_ctmc_chain_transient requires a MarkovProcess or a MarkovChain.');
end

if isempty(pi0)
    pi0 = ones(1,n)/n;
end
pi0 = reshape(pi0, 1, n);
if abs(sum(pi0) - 1) > GlobalConstants.FineTol
    line_error(mfilename,'The initial distribution must sum to one.');
end

t0 = timespan(1);
t1 = timespan(2);
if ~isfinite(t1)
    line_error(mfilename,'A finite timespan is required, e.g., SolverCTMC(chain,''timespan'',[0,T]).');
end
if ~isfinite(t0)
    t0 = 0;
end

if isa(chain,'MarkovChain')
    P = chain.getTransMat();
    k0 = ceil(t0);
    k1 = floor(t1);
    if k1 < k0
        line_error(mfilename,'The timespan [%g,%g] contains no integer step of the DTMC.', t0, t1);
    end
    t = (k0:k1)';
    pi_t = zeros(length(t), n);
    pik = pi0 * P^k0;
    for k = 1:length(t)
        pi_t(k,:) = pik;
        pik = pik * P;
    end
else
    [pi_t, t] = ctmc_transient(chain.getGenerator(), pi0, t0, t1);
    t = t(:);
end
end
