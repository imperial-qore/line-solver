%{ @file dtmc_stochcomp.m
 %  @brief Computes the stochastic complement of a DTMC partition
 %
 %  @author LINE Development Team
%}

%{
 % @brief Computes the stochastic complement of a DTMC partition
 %
 % @details
 % Performs stochastic complementation on a partitioned DTMC.
 %
 % @par Syntax:
 % @code
 % S = dtmc_stochcomp(P)
 % [S, P11, P12, P21, P22] = dtmc_stochcomp(P, I)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>P<td>Stochastic transition matrix
 % <tr><td>I<td>(Optional) Vector of indices representing the subset of states to retain
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>S<td>Stochastic complement matrix for the subset I
 % <tr><td>P11<td>Submatrix corresponding to transition within I
 % <tr><td>P12<td>Submatrix corresponding to transitions from I to complement of I
 % <tr><td>P21<td>Submatrix corresponding to transitions from complement of I to I
 % <tr><td>P22<td>Submatrix corresponding to transitions within complement of I
 % </table>
%}
function [S,P11,P12,P21,P22]=dtmc_stochcomp(P,I)
if nargin==1
    I=1:ceil(length(P)/2);
end
Iall = 1:length(P);
Ic = Iall(~ismember(Iall,I)); % slightly faster than setdiff
sparsity = nnz(P)/prod(size(P));
if sparsity <0.05
    P=sparse(P);
    P11=P(I,I);
    P12=P(I,Ic);
    P21=P(Ic,I);
    P22=P(Ic,Ic);
    S2 = sparse(eye(size(P22))-P22);
    S=P11+P12*local_solve(S2,P21);
    S=full(S);
else
    P11=P(I,I);
    P12=P(I,Ic);
    P21=P(Ic,I);
    P22=P(Ic,Ic);
    S2 = eye(size(P22))-P22;
    S=P11+P12*local_solve(S2,P21);
end
end

function T = local_solve(S2,P21)
% Iterative path. The backslash below factorizes the whole complement block,
% so its fill-in is what limits the model above the dispatch threshold rather
% than the arithmetic. One ILUT factorization serves every column of P21. The
% direct solve stays the default and remains the fallback.
GMRES_MIN_STATES = 6000;
T = [];
if size(S2,1) > GMRES_MIN_STATES
    [T,gflag] = ctmc_gmres_multi(sparse(S2), P21);
    if gflag ~= 0
        T = [];
    end
end
if isempty(T)
    T = S2 \ P21;
    if ~all(isfinite(T(:)))
        % Backslash returns NaN on a rank-deficient complement, which then
        % contaminates every downstream metric. Fall back to the minimum-norm
        % least-squares solution, as the Java and Python kernels do.
        T = lsqminnorm(full(S2), full(P21));
    end
end
end