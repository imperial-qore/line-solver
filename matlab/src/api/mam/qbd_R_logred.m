%{ @file qbd_R_logred.m
 %  @brief Computes the rate matrix R using logarithmic reduction
 %
 %  @author LINE Development Team
%}

%{
 % @brief Computes QBD rate matrix R via logarithmic reduction
 %
 % @details
 % This function computes the rate matrix R for a Quasi-Birth-Death (QBD)
 % process using the logarithmic reduction method.
 %
 % @par Syntax:
 % @code
 % R = qbd_R_logred(B, L, F)
 % R = qbd_R_logred(B, L, F, iter_max)
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
function R=qbd_R_logred(B,L,F,iter_max)
% Logarithmic reduction method
if nargin<4
    iter_max = 100000;
end
r = size(L,1);
iLF = -inv(L)*F;
iLB = -inv(L)*B;
T = iLF;
S = iLB;
for iter=1:iter_max
    D = iLF*iLB + iLB*iLF;
    iLF = inv(eye(r)-D) *iLF*iLF;
    iLB = inv(eye(r)-D) *iLB*iLB;
    S = S + T*iLB;
    T = T*iLF;
    if norm(ones(r,1)-S*ones(r,1) ,1) <= 1e-12
        break
    end
end
U = L + F*S;
R = -F * inv(U);
end
