%{ @file ba_ignores_blocking.m
 %  @brief Whether a SolverBA bound method is blind to finite-buffer blocking
 %
 %  @author LINE Development Team
%}

%{
 % @brief True when METHOD bounds a model as if its buffers were unbounded
 %
 % @details
 % Every bound family in SolverBA except the QRF BLOCKING bounds is
 % parameterized by demands (visits x service time) and a population alone,
 % which is the BCMP parameterization: the buffers are unbounded and the
 % equilibrium distribution factorizes. A finite buffer that BINDS breaks
 % both premises -- the truncation couples the station occupancies -- so the
 % resulting numbers do not bracket the blocked model. On cqn_bas_blocking
 % (Queue2 capped at 1, N = 2) gb.upper reports QLen 1.28 at a station that
 % can never hold more than one job, which is not a loose bound but a wrong
 % one, and the lower sides are not one-sided at all.
 %
 % The exceptions are 'qrf.bas*' and 'qrf.rsrd', the Quadratic Reduction
 % Framework bounds that carry the blocking tables (MM, MM1, ZZ, ZM, BB, F)
 % explicitly and therefore model the finite buffer rather than ignore it,
 % and 'spnlp.*', whose polytope is indexed by the marking itself: a bounded
 % place enters it as a variable upper bound and as the P-invariant equality
 % that produced the bound, so the buffer is modelled and not assumed away.
 %
 % METHOD is the RESOLVED name (after the 'default'/'auto'/'qr'/'lr'
 % aliases), so a caller must resolve before asking.
 %
 % @par Syntax:
 % @code
 % bool = ba_ignores_blocking(method)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>method<td>Resolved SolverBA bound method name
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>bool<td>True if the method assumes unbounded buffers
 % </table>
%}
function bool = ba_ignores_blocking(method)
bool = ~(startsWith(method,'qrf.bas') || startsWith(method,'qrf.rsrd') || ...
    startsWith(method,'spnlp'));
end
