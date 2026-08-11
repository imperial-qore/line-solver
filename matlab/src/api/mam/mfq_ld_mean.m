%{ @file mfq_ld_mean.m
 %  @brief Stationary mean fluid level of a level-dependent fluid queue
 %
 %  @author LINE Development Team
%}

%{
 % @brief Stationary mean fluid level from mfq_ld_solve output.
 %
 % @details
 % Thin LINE wrapper around LevelDependentFluidStationaryMean. Returns the
 % scalar mean fluid level E[X] of a first/second-order level-dependent fluid
 % queue in closed form, using the matrix-exponential building blocks returned
 % by mfq_ld_solve.
 %
 % @par Syntax:
 % @code
 % res = mfq_ld_mean(masses,iniF,KF,cloF,iniB,KB,cloB,T)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>masses,iniF..cloB<td>Building blocks returned by mfq_ld_solve
 % <tr><td>T<td>Vector of regime thresholds (length K)
 % </table>
 %
 % @par Returns:
 % res: scalar mean fluid level E[X].
%}
function res = mfq_ld_mean(masses,iniF,KF,cloF,iniB,KB,cloB,T)
res = LevelDependentFluidStationaryMean(masses,iniF,KF,cloF,iniB,KB,cloB,T);
end
