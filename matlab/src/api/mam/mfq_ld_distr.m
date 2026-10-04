%{ @file mfq_ld_distr.m
 %  @brief Stationary distribution of a level-dependent fluid queue
 %
 %  @author LINE Development Team
%}

%{
 % @brief Stationary fluid-level distribution from mfq_ld_solve output.
 %
 % @details
 % Thin LINE wrapper around LevelDependentFluidStationaryDistr. Evaluates the
 % stationary density or distribution of a first/second-order level-dependent
 % fluid queue at the requested points, using the matrix-exponential building
 % blocks returned by mfq_ld_solve.
 %
 % @par Syntax:
 % @code
 % res = mfq_ld_distr(masses,iniF,KF,cloF,iniB,KB,cloB,T,what,points)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>masses,iniF..cloB<td>Building blocks returned by mfq_ld_solve
 % <tr><td>T<td>Vector of regime thresholds (length K)
 % <tr><td>what<td>'pdf', 'pdfd' (density derivative), 'cdf' P(X<p), or 'cdfm' P(X<=p)
 % <tr><td>points<td>Fluid levels at which to evaluate
 % </table>
 %
 % @par Returns:
 % res: (numel(points), N) matrix with the per-state values at each point.
%}
function res = mfq_ld_distr(masses,iniF,KF,cloF,iniB,KB,cloB,T,what,points)
res = LevelDependentFluidStationaryDistr(masses,iniF,KF,cloF,iniB,KB,cloB,T,what,points);
end
