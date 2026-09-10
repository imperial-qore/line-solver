%{
 % @file pfqn_lcfsqn_ca.m
 % @brief Convolution algorithm for 2-station LCFS queueing networks
 %
 % @author LINE Development Team
%}

%{
 % @brief Convolution algorithm for 2-station LCFS queueing networks
 % @fn pfqn_lcfsqn_ca(alpha, beta, N)
 % @param alpha Service rates at LCFS station (1xR vector).
 % @param beta Service rates at LCFS-PR station (1xR vector).
 % @param N Population vector (default: ones(1,R)).
 % @return G Normalizing constant.
 % @return V Auxiliary normalization term.
%}
function [G,V] = pfqn_lcfsqn_ca(alpha,beta,N)
% [G,V] = PFQN_LCFSQN_CA(ALPHA, BETA, N)
% Convolution algorithm for multiclass LCFS queueing networks
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
%           (default: ones(1,R) - one job per class)
%
% Returns:
%   G - normalizing constant
%   V - auxiliary normalization term
%
% Reference:
%   G. Casale, "A family of multiclass LCFS queueing networks with
%   order-dependent product-form solutions", QUESTA 2026.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

R = length(alpha);
if nargin<3
    N = ones(1,R);
end
K = sum(N);

if K ==0
    G = 1;
    V = 1;
    return
else
    prods = zeros(1,R);
    for r=1:R
        prods(r)=prod(N(1:r-1)+1);
    end
    G = zeros(prod(N+1),1);
    V = zeros(prod(N+1),1);
    n = pprod(N);
    while n>=0
        idx = hashpop(n,N,R,prods);
        if sum(n)==0
            G(idx) = 1;
            V(idx) = 1;
        else
            V(idx) = 0;
            G(idx) = 0;
            for r=1:R
                if n(r)>0
                    idx_r = hashpop(oner(n,r),N,R,prods);
                    % recursive definition of V
                    V(idx) = V(idx) + V(idx_r);
                    % recursive definition of G uses V and beta
                    G(idx) = G(idx) + alpha(r)^(sum(n)-1)*beta(r) * G(idx_r);
                end
            end
            V(idx) = prod(alpha.^n) * V(idx);
            G(idx) = G(idx) + V(idx);
        end
        n = pprod(n,N);
    end
    G=G(end);
    V=V(end);
end
end
