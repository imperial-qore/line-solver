function assignment = perm_maxweight_assignment(W)
% PERM_MAXWEIGHT_ASSIGNMENT  Maximum-weight perfect assignment of a square matrix.
%
% ASSIGNMENT = PERM_MAXWEIGHT_ASSIGNMENT(W) returns a permutation ASSIGNMENT of
% 1..N maximizing SUM_I W(I, ASSIGNMENT(I)).
%
% This replaces the row-by-row greedy that PERM_HUBERLAW called
% SUB_GREEDY_ASSIGNMENT. That routine was greedy despite the name its twins
% carry (_hungarian_assignment in python, hungarianAssignment in the JAR) and
% it returns a zero-weight assignment on inputs that admit a positive one: on
% [1 2; 0 3] row 1 takes the larger entry in column 2, leaving row 2 with the
% zero in column 1.
%
% That matters because the weight is ALPHA3 in PERM_HUBERLAW, and ALPHA3 sets
% the flooring level ALPHA1 = ALPHA3*DELTA/(3*N!) of the Huber-Law sampler --
% the one principled zero-handling in the whole family. A suboptimal assignment
% understates ALPHA3, hence ALPHA1, and weakens the method's own guarantee. A
% correct assignment makes ALPHA3 > 0 whenever the matrix has a positive
% permanent.
%
% W may carry -Inf for a forbidden pair. The three twins solve the same problem
% and agree on the optimal VALUE, which is all that ALPHA3 depends on; they need
% not return the same permutation when the optimum is degenerate.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

n = size(W, 1);
if size(W, 2) ~= n
    line_error(mfilename, 'The weight matrix must be square.');
end
assignment = zeros(1, n);
if n == 0
    return
end

% matchpairs minimizes, so negate. costUnmatched must exceed any cost a perfect
% matching could incur, so that leaving a row unmatched is never preferred.
finite = W(isfinite(W));
if isempty(finite)
    span = 1;
else
    span = max(1, max(finite(:)) - min(finite(:)));
end
costUnmatched = n * span + 1;

M = matchpairs(-W, costUnmatched);
if size(M, 1) < n
    line_error(mfilename, ['No perfect assignment exists: only %d of %d rows ' ...
        'could be matched to a finite-weight column.'], size(M, 1), n);
end
assignment(M(:,1)) = M(:,2).';
end
