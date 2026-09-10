%{
%{
 % @file pfqn_cntol.m
 % @brief Chandy-Neuse population-scaled termination cutoff for approximate MVA.
%}
%}

%{
%{
 % @brief Chandy-Neuse population-scaled termination cutoff for approximate MVA.
 % @fn pfqn_cntol(N)
 % @param N Population vector.
 % @return tol Termination cutoff 1/(4000+16*sum(N)).
%}
%}
function tol = pfqn_cntol(N)
% Published cutoff of the Linearizer termination test: K. M. Chandy, D. Neuse,
% "Linearizer: A Heuristic Algorithm for Queuing Network Models of Computing
% Systems", Commun. ACM 25(2):126-134, 1982, p.129 and appendix. The iteration
% continues while
%   max_{i,r} |Q^I(i,r) - Q^{I-1}(i,r)| / N_r > 1/(4000 + 16*|N|),
% |N| = sum(N). The paper motivates the scaling with |N|: at large populations
% removing one job changes the queue lengths very little, so a fixed cutoff
% would terminate the iteration prematurely. It also notes that the expression
% stays below 0.00025 even at very small populations.
%
% The same expression is what LQNS uses as its termination test, set in the
% SchweitzerCommon constructor of libmva/src/mva.cc; that code carries no
% citation, and the paper above is its source.
%
% Pass the result as the tol argument of pfqn_bs / pfqn_egflinearizer only if
% the plain cutoff is wanted with those functions' own convergence metric.
% Passing tol = 'cn' (or NaN) instead selects BOTH this cutoff and the
% normalized-maximum metric of the paper, which is the published test.
tol = 1 / (4000 + 16*sum(N(:)));
end
