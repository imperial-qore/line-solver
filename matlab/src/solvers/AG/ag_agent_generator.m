function Qk = ag_agent_generator(k, x, Aa, Pb, L, ACT, PSV, A, N)
% AG_AGENT_GENERATOR Agent k's generator at the current reversed rates.
%
% QK = AG_AGENT_GENERATOR(K, X, AA, PB, L, ACT, PSV, A, N)
%
% SPLIT OUT FROM THE SOLVE BECAUSE THE TWO HALVES COST DIFFERENT ORDERS.
% Assembling the generator is O(N^2) and solving it is O(N^3), so the 'cluster'
% backend ships only the stationary vector back and rebuilds the generator on
% the coordinator: sending an N-by-N matrix per agent per sweep to save the
% cheaper half would spend more on the wire than it saves on the worker.

% Start with local/hidden rates
Qmat = L{k} - diag(L{k} * ones(N(k), 1));

% Add contributions from each action
for c = 1:A
    if PSV(c) == k
        % Process k is passive for action c: add x(c) * Pb{c}
        Qmat = Qmat + x(c) * Pb{c} - diag(Pb{c} * ones(N(k), 1));
    elseif ACT(c) == k
        % Process k is active for action c: add Aa{c}
        Qmat = Qmat + Aa{c} - diag(Aa{c} * ones(N(k), 1));
    end
end

% Convert to valid infinitesimal generator
Qk = ctmc_makeinfgen(Qmat);

end

%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
