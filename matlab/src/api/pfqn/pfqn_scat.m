%{
%{
 % @file pfqn_scat.m
 % @brief Neuse-Chandy SCAT (Self-Correcting Approximation Technique) AMVA.
%}
%}

%{
%{
 % @brief Neuse-Chandy SCAT approximate mean value analysis.
 %
 % @details
 % SCAT shares the Linearizer fixed point. It carries the mean queue lengths at
 % the target population N and at the R reduced populations N-e_s, and corrects
 % the Bard-Schweitzer proportionality assumption with the fraction difference
 %
 %   Delta(i,r,s) = Q(i,r|N-e_s)/(N-e_s)_r - Q(i,r|N)/N_r,
 %
 % held fixed while an inner MVA fixed point is iterated. It differs from
 % Linearizer in that this correction is refreshed ONCE: SCAT stops after the
 % first pass, where Linearizer performs the fixed three passes of Chandy and
 % Neuse (1982), Sec. 4. Cost is therefore about one third of Linearizer's, and
 % accuracy sits between Bard-Schweitzer (pfqn_bs, the Delta=0 special case) and
 % Linearizer.
 %
 % SCAT's second departure from Linearizer, fitting a probability mass function
 % centred on the mean queue length at queue-dependent centres instead of
 % propagating the MVA distribution recursion (Krzesinski and Greyling 1984,
 % Sec. 4), does not arise here: this entry point covers single-server and delay
 % stations only, exactly as pfqn_linearizer does. That mass function is
 % available separately as the 'scat' marginal rule of pfqn_ab_amva.
 %
 % Reference: D. Neuse, K. M. Chandy, "SCAT: A Heuristic Algorithm for Queueing
 % Network Models of Computing Systems", ACM SIGMETRICS Perform. Eval. Rev.
 % 10(3), 1981; K. M. Chandy, D. Neuse, "Linearizer: A Heuristic Algorithm for
 % Queuing Network Models of Computing Systems", Commun. ACM 25(2), 1982.
 %
 % @fn pfqn_scat(L, N, Z, type, tol, maxiter, QN0)
 % @param L Service demand matrix (M x R).
 % @param N Population vector (1 x R).
 % @param Z Think time vector (1 x R) or matrix summed over rows.
 % @param type Scheduling strategy per station; accepted for interface parity
 %        with pfqn_linearizer, but the residence-time recursion is
 %        discipline-independent.
 % @param tol Convergence tolerance (default: 1e-8); 'cn' or NaN selects the Chandy-Neuse (1982) population-scaled termination test, see pfqn_cntol.
 % @param maxiter Maximum inner iterations (default: 1000).
 % @param QN0 (M x R) queue lengths that warm-start the Bard-Schweitzer
 %        initialization; empty for the default cold start.
 % @return Q Mean queue lengths.
 % @return U Utilization.
 % @return W Residence times.
 % @return C Cycle times.
 % @return X Class throughputs.
 % @return totiter Total iterations performed.
%}
%}
function [Q,U,W,C,X,totiter] = pfqn_scat(L,N,Z,type,tol,maxiter,QN0)
if nargin<7
    QN0 = [];
end
if nargin<6 || isempty(maxiter)
    maxiter = 1000;
end
if nargin<5 || isempty(tol)
    tol = 1e-8;
end
if nargin<4
    type = [];
end
alpha = ones(size(N));
% npasses = 1 is what separates SCAT from Linearizer: one Delta refresh, not three
[Q,U,W,C,X,totiter] = pfqn_egflinearizer(L,N,Z,type,tol,maxiter,alpha,QN0,1);
end
