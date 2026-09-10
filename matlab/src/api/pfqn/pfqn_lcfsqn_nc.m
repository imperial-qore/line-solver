%{
 % @file pfqn_lcfsqn_nc.m
 % @brief Normalizing constant for LCFS queueing networks
 %
 % @author LINE Development Team
%}

%{
 % @brief Computes the normalizing constant for LCFS queueing networks
 % @fn pfqn_lcfsqn_nc(alpha, beta, N)
 % @param alpha Service rates at LCFS station (1xR vector).
 % @param beta Service rates at LCFS-PR station (1xR vector).
 % @param N Population vector (default: ones(1,R)).
 % @return G Normalizing constant.
 % @return Ax Cell array of A matrices for each state.
%}
function [G,Ax] = pfqn_lcfsqn_nc(alpha,beta,N)
% [G,AX] = PFQN_LCFSQN_NC(ALPHA, BETA, N)
% Normalizing constant for multiclass LCFS queueing networks
%
% This function computes the normalizing constant for a 2-station closed
% queueing network with:
%   - Station 1: LCFS (Last-Come-First-Served, non-preemptive)
%   - Station 2: LCFS-PR (LCFS with Preemption-Resume)
%
% Parameters:
%   alpha - vector of inverse service rates at station 1 (LCFS)
%           alpha(r) = 1/mu(1,r) for class r
%   beta  - vector of inverse service rates at station 2 (LCFS-PR)
%           beta(r) = 1/mu(2,r) for class r
%   N     - population vector, N(r) = number of jobs of class r
%
% Returns:
%   G  - normalizing constant
%   Ax - cell array of A matrices for each state x=0:K
%
% Reference:
%   G. Casale, "A family of multiclass LCFS queueing networks with
%   order-dependent product-form solutions", QUESTA 2026.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

K = sum(N);
R = length(N);
G = 0;
Ax=cell(1,K+1);
for x=0:K
    Ax{1+x} = make_A(alpha, beta, x, K, R);
    G = G + perm(Ax{1+x}, N);
end
% The permanent counts the N(r)! orderings of the identical class-r jobs, so it
% overstates G by prod(N!). pfqn_joint.m divides by exactly this factor for the
% same reason. Without it this routine disagreed with pfqn_lcfsqn_ca, which the
% NC solver actually uses as the normalizing constant.
G = G / prod(factorial(N));
end


function A = make_A(alpha, beta, x, K, R)
% alpha : vector of length R
% beta  : vector of length R
% x     : integer
% K     : matrix size
% perm(A,N) requires A to be (sum(N) x R): rows are job slots, columns are
% class groups (see matlab/util/perm.m). This built the TRANSPOSE, class on the
% row, in a K x K buffer, so rows R+1..K stayed zero and every Ryser term
% carried a zero factor: G came back as 0 for every population other than the
% all-ones default, where the matrix happens to be square and symmetric under
% transposition.
if issym(alpha)
    A = sym(zeros(K, R));
else
    A = zeros(K, R);
end
for i = 1:R
    for j = 1:x
        A(j, i) = alpha(i)^j;
    end
    for j = 1:(K-x)
        A(x+j, i) = alpha(i)^(x+j-1) * beta(i);
    end
end
end
