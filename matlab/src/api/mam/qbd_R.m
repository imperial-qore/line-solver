%{ @file qbd_R.m
 %  @brief Computes the rate matrix R using successive substitutions
 %
 %  @author LINE Development Team
%}

%{
 % @brief Computes QBD rate matrix R via successive substitutions
 %
 % @details
 % This function computes the rate matrix R for a Quasi-Birth-Death (QBD)
 % process using the successive substitutions method.
 %
 % @par Syntax:
 % @code
 % R = qbd_R(B, L, F)
 % R = qbd_R(B, L, F, iter_max)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>B<td>Backward transition block A_(-1)
 % <tr><td>L<td>Local transition block A_0
 % <tr><td>F<td>Forward transition block A_1
 % <tr><td>iter_max<td>(Optional) Maximum iterations (default: 100000)
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>R<td>Rate matrix R
 % </table>
%}
function R=qbd_R(B,L,F,iter_max)
% Successive substitutions method
if nargin<4
    iter_max = 100000;
end
Fil = F * inv(L);
BiL = B * inv(L);
R = -Fil;
Rprime = -Fil - R^2*BiL;
for iter=1:iter_max
    R = Rprime;
    Rprime = -Fil -R^2*BiL;
    if norm(R-Rprime,1)<=1e-12
        break
    end
end
R = Rprime;
end
