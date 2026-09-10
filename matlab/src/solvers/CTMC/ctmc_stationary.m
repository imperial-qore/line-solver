function pi = ctmc_stationary(Q, StateSpace, sn, options)
% PI = CTMC_STATIONARY(Q, STATESPACE, SN, OPTIONS)
%
% Single entry point for the stationary distribution of a CTMC generated from a
% NetworkStruct.
%
% All the stationary mass of a reducible chain lives in its bottom strongly
% connected components, each weighted by the probability of being absorbed in it
% from the declared initial state; every other state is transient and carries
% zero. The block decomposition handles the irreducible case as the degenerate
% one BSCC / no transient states, so every CTMC solve goes through it and no
% dispatch can disagree with the algorithm about whether a chain is reducible.
%
% See _kb/11-conventions-and-gotchas.md.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 4
    options = [];
end
if nargin < 3
    sn = [];
end
if nargin < 2
    StateSpace = [];
end

pi0 = ctmc_initial_distribution(Q, StateSpace, sn, options);
if isempty(pi0)
    warn_if_unseeded_mixture(Q);
end
pi = ctmc_solve_reducible_blkdecomp(Q, pi0, options);
end

function warn_if_unseeded_mixture(Q)
% Warn when a reducible generator is solved with no initial distribution.
%
% Without a seed the block decomposition invents a start distribution -- here a
% uniform one over the SCCs with no incoming transition -- and no property of the
% model implies it: on a reducible chain the stationary distribution is fixed
% only by the initial state. Nor is it the product-form weighting, which weights
% the recurrent classes by their unnormalized Kelly mass. The invented start is
% order-independent and therefore looks more reproducible than the answer the
% declared initial state selects, but that reproducibility is bought by
% discarding the one input that makes the problem well posed.
%
% The fallback itself differs across codebases (python weights ALL the SCCs
% equally, this one the source SCCs), which is a second reason not to read the
% number as the model's answer. See G60 in line-gaps.md.
n = size(Q,1);
if n < 2
    return
end
A = abs(Q);
A(1:n+1:end) = 0;
[~, recurrent] = stronglyconncomp(A > GlobalConstants.ArcTol);
nbscc = sum(recurrent);
if nbscc > 1
    line_warning_always(mfilename, sprintf(['The generator has %d closed communicating classes and the ' ...
        'declared initial state could not be located in the enumerated state space, so the solve ' ...
        'starts from a distribution the model never stated (uniform over the SCCs with no incoming ' ...
        'transition). On a reducible chain the stationary distribution is determined only by the ' ...
        'initial state, so this answer is not the model''s. Call setState on the stations so the ' ...
        'class the model actually starts in is the one solved.'], nbscc));
end
end

function pi0 = ctmc_initial_distribution(Q, StateSpace, sn, options) %#ok<INUSD>
% Point mass at the initial state of SN, or empty when that state is absent from
% STATESPACE: stochastic complementation may have removed it (an SPN whose
% immediate ENABLE states were eliminated, for one). An empty seed makes the
% block decomposition start in the SCCs with no incoming transition.
pi0 = [];
if isempty(StateSpace) || isempty(sn) || ~isfield(sn,'state') || isempty(sn.state)
    return
end
if any(cellfun(@isempty, sn.state))
    return
end
initRow = [];
for isf = 1:sn.nstateful
    row_isf = sn.state{isf};
    row_isf = row_isf(1,:);
    w_isf = size(sn.space{isf},2);
    if numel(row_isf) < w_isf
        % the initial row carries only the buffer slots the population needs, while
        % the enumerated local space is sized for the full capacity; pad on the LEFT
        row_isf = [zeros(1, w_isf-numel(row_isf)), row_isf]; %#ok<AGROW>
    end
    initRow = [initRow, row_isf]; %#ok<AGROW>
end
if size(StateSpace,2) ~= numel(initRow)
    return
end
initState = matchrow(StateSpace, initRow);
if initState > 0
    pi0 = zeros(1, size(Q,1));
    pi0(initState) = 1.0;
    line_debug('Seeding absorption from initial state %d', initState);
else
    line_debug('Initial state absent from the state space, absorption seeded by the source SCCs');
end
end
