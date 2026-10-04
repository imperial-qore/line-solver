function n = roundMarginalPreservingChains(n, sn)
% N = ROUNDMARGINALPRESERVINGCHAINS(N, SN)
%
% Round a fractional marginal queue-length matrix (station x class) to
% integers with the largest remainder method, so that every closed chain
% keeps exactly its own population. Plain element-wise rounding does not: it
% can move a job between two classes of the same chain or lose one
% altogether, which yields a state outside the state space, or inside a
% different chain population and hence a different steady state. Open chains
% are rounded element-wise, their population being unbounded.
%
% Twin of jline.lang.state.State.roundMarginalPreservingChains (JAR) and of
% line_solver.api.state.marginal.roundMarginalPreservingChains (python).

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

for c = 1:sn.nchains
    chainClasses = find(sn.chains(c,:) > 0);
    if isempty(chainClasses)
        continue
    end
    njobs_chain = sum(sn.njobs(chainClasses));
    if isinf(njobs_chain)
        % Open chain: element-wise rounding
        n(:, chainClasses) = round(n(:, chainClasses));
        continue
    end
    % Closed chain: largest remainder method over the chain's whole block.
    % vals is laid out STATION-MAJOR, class-minor, which is the order the JAR
    % and python twins build it in. sort is stable in all three, so the layout
    % is what breaks a tie between two cells with the same remainder: a
    % column-major vals would hand the spare job to a different station.
    blockT = n(:, chainClasses).';
    vals = blockT(:);
    floored = floor(vals);
    remainders = vals - floored;
    deficit = round(njobs_chain - sum(floored));
    if deficit > 0
        [~, sortIdx] = sort(remainders, 'descend');
        for d = 1:min(deficit, numel(sortIdx))
            floored(sortIdx(d)) = floored(sortIdx(d)) + 1;
        end
    end
    n(:, chainClasses) = reshape(floored, size(blockT)).';
end
end
