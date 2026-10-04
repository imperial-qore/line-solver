%{ @file dtmc_isfeasible.m
 %  @brief Checks feasibility of a stochastic matrix
 %
 %  @author LINE Development Team
%}

%{
 % @brief Checks feasibility of a stochastic matrix
 %
 % @details
 % Verifies if row sums are close to 1 and elements are non-negative.
 %
 % @par Syntax:
 % @code
 % res = dtmc_isfeasible(P)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>P<td>Matrix to check
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>res<td>Precision level (1-15) if feasible, or 0 if not feasible
 % </table>
%}
function res=dtmc_isfeasible(P)
sP = sum(P,2);
res = 0;
for tol=1:15
    if min(sP) > 1-10^-tol && max(sP) < 1+10^-tol && min(P(:)) > -10^-tol
        res = tol;
    end
end

end