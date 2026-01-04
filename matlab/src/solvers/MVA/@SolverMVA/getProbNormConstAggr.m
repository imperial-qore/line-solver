function [lNormConst] = getProbNormConstAggr(self)
% [LNORMCONST] = GETPROBNORMCONST() Returns normalizing constant for state probabilities
%
% @brief Computes the logarithm of the normalizing constant for closed queueing networks
%
% This method returns the logarithm of the normalizing constant G, which is essential
% for computing state probabilities in closed queueing networks. The normalizing
% constant ensures that probability distributions sum to 1 across all possible states.
%
% For closed networks, state probabilities are computed as:
%   Pr[state] = (product of terms) / G
% where G is the normalizing constant. Computing log(G) avoids numerical underflow
% for large systems with many jobs.
%
% The MVA solver computes this using the exact MVA algorithm on the network lattice,
% which is numerically stable and efficient for closed networks.
%
% @param self SolverMVA instance
%
% @return lNormConst Logarithm of the normalizing constant (base e)
%         - Single real number (can be very large for large systems)
%         - Positive value representing log(G) where G > 1
%         - Used internally for probability computation
%         - Exponentiating gives actual normalizing constant: G = exp(lNormConst)
%
% @note Only defined for closed queueing networks with fixed populations.
%       For open or mixed networks, returns empty or NaN.
%       Used internally by getProb_aggr() for probability calculations.
%
% @warning Computing the normalizing constant for very large systems may be slow.
%          This is a fundamental computation in closed network analysis and cannot
%          be significantly accelerated without changing the algorithm.
%
% @see getProbAggr - Uses normalization constant to return state probabilities
% @see getProbSysAggr - Returns system-level probabilities
% @see getAvg - Get average metrics (alternative to probability computation)
%
% Example:
% @code
% model = Network('ClosedQN');
% % ... setup closed network ...
% solver = SolverMVA(model, 'method', 'exact');
%
% % Compute normalizing constant
% log_G = solver.getProbNormConstAggr();
% G = exp(log_G);
%
% fprintf('Normalizing constant G = %.4f\\n', G);
% fprintf('Log normalizing constant = %.4f\\n', log_G);
% @endcode

if ~isempty(self.result)
    lNormConst = self.result.Prob.logNormConstAggr;
else
    optnc = self.options;
    optnc.method = 'exact';
    [~,~,~,~,~,~,lNormConst] = solver_mva_analyzer(self.getStruct, optnc);
    self.result.Prob.logNormConstAggr = lNormConst;
end
end