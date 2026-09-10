%{ @file qbd_depproc_jointmom.m
 %  @brief Computes joint moments of consecutive inter-departure times
 %
 %  @author LINE Development Team
%}

%{
 % @brief Computes joint moments E[X_0^i * X_1^j] of consecutive inter-departure times
 %
 % @details
 % Given arrival and service MAPs for a MAP/MAP/1 queue, this function
 % computes joint moments of consecutive inter-departure times using the
 % QBD structure. The initial vector is constructed from the stationary
 % distribution at departure epochs via ETAQA.
 %
 % The departure SCV and lag-1 ACF can be obtained as:
 %   E1  = qbd_depproc_jointmom(MAPa, MAPs, [1,0])  % E[X]
 %   E2  = qbd_depproc_jointmom(MAPa, MAPs, [2,0])  % E[X^2]
 %   E11 = qbd_depproc_jointmom(MAPa, MAPs, [1,1])  % E[X_0*X_1]
 %   SCV = (E2 - E1^2) / E1^2
 %   ACF = (E11 - E1^2) / (E2 - E1^2)
 %
 % @par Syntax:
 % @code
 % JM = qbd_depproc_jointmom(MAPa, MAPs, iset)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>MAPa<td>Arrival process in MAP format {D0, D1}
 % <tr><td>MAPs<td>Service process in MAP format {D0, D1}
 % <tr><td>iset<td>Matrix of moment orders [i1,j1; i2,j2; ...] where each
 %                  row specifies the exponents for E[X_0^i * X_1^j]
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>JM<td>Vector of joint moments, one per row of iset
 % </table>
%}
function [JM]=qbd_depproc_jointmom(MAPa,MAPs,iset)

na = length(MAPa{1});
ns = length(MAPs{1});
lvlsz = ns*na;

F = kron(MAPa{2},eye(ns));
L = krons(MAPa{1},MAPs{1});
B = kron(eye(na),MAPs{2});
L0 = kron(MAPa{1},eye(ns));
Z = 0*F;

[~,R,~] = QBD_CR(B,L,F);
pi = QBD_pi(B,L0,R);
v0  = pi(1,1:lvlsz); % QBD_pi returns every level in one row, only level 0 is the boundary vector

lambdaS = map_lambda(MAPs);
% departure epochs are the B transitions, so the embedded vector weighs the level
% probabilities by B and not by the arrival matrix F
v0D = 1/lambdaS*v0 *R * B;
v1D = 1/lambdaS*v0 *R^2 * B;
v2Dp = 1/lambdaS*v0 * R^3*inv(eye(size(R))-R)*B;
z = [v0D, v1D, v2Dp];
z = z / sum(z); % normalize to probability distribution

M0 = [L0,F,Z; Z,L,F; Z,Z,L+F];
M1 = [Z,Z,Z; B,Z,Z; Z,B,Z];

for k=1:size(iset,1)
    JM(k) = z * factorial(iset(k,1)) * (-M0)^(-iset(k,1)-1) * M1 * factorial(iset(k,2)) * (-M0)^(-iset(k,2)) * ones(length(M0),1);
end

end
