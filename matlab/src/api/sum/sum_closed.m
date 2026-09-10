function [XN,QN,UN,RN,it] = sum_closed(L,N,Z,mi,scv,tol,maxiter)
% [XN,QN,UN,RN,IT] = SUM_CLOSED(L,N,Z,MI,SCV,TOL,MAXITER)
%
% Summation method (SUM) for closed queueing networks, including the
% extended SUM (ESUM) node functions for non-product-form networks with
% generally distributed service times.
%
% The method expresses the mean queue length of each station as a function
% of its throughput, Ki = fi(lambdai), and solves the population constraint
% sum_i Ki = K. Single-class models are solved by bisection on the system
% throughput (Bolch et al., Sec. 9.2.1); multiclass models by fixed-point
% iteration on the class throughputs (Sec. 9.2.2, Eqs. 9.24-9.26).
%
% Node functions:
% - Product-form stations (scv=1, or insensitive disciplines PS/LCFS-PR,
%   for which the caller must pass scv=1): Eq. (9.15)/(9.19).
% - FCFS stations with general service (scv~=1): ESUM corrections,
%   Eq. (10.88) for -/G/1 and Eq. (10.89) for -/G/m, with
%   ai=(1+scv_i)/2 and Erlang-C waiting probability P_mi.
% - Infinite-server stations (mi=Inf) and think times Z: Ki = lambdai*Li.
%
% Input:
% L       - MxR service demand matrix, L(i,r) = e(i,r)/mu(i,r)
% N       - 1xR population vector
% Z       - 1xR think times (aggregated as a delay term)
% mi      - Mx1 number of servers (Inf for infinite-server stations)
% scv     - MxR squared coefficient of variation of service times
% tol     - convergence tolerance (default: 1e-6)
% maxiter - maximum number of iterations (default: 10000)
%
% Output:
% XN - 1xR class throughputs
% QN - MxR mean queue lengths
% UN - MxR utilizations (per-server for queueing stations, X.*L for IS)
% RN - MxR residence times, RN=QN./XN
% it - number of iterations
%
% Reference: G. Bolch, S. Greiner, H. de Meer, K.S. Trivedi, Queueing
% Networks and Markov Chains, 2nd ed., Wiley, 2006, Secs. 9.2 and 10.1.4.4.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

[M,R] = size(L);
if nargin<3 || isempty(Z), Z = zeros(1,R); end
if nargin<4 || isempty(mi), mi = ones(M,1); end
if nargin<5 || isempty(scv), scv = ones(M,R); end
if nargin<6 || isempty(tol), tol = 1e-6; end
if nargin<7 || isempty(maxiter), maxiter = 10000; end
mi = mi(:);
Z = Z(:)';
N = N(:)';
K = sum(N(isfinite(N)));

XN = zeros(1,R);
QN = zeros(M,R);
UN = zeros(M,R);
RN = zeros(M,R);
it = 0;

if K==0
    return
end

if R==1
    % single class: bisection on the throughput (Sec. 9.2.1)
    lambda_l = 0;
    lambda_u = Inf;
    for i=1:M
        if L(i)>0
            if isinf(mi(i))
                lambda_u = min(lambda_u, K/L(i));
            else
                lambda_u = min(lambda_u, mi(i)/L(i));
            end
        end
    end
    if Z>0
        lambda_u = min(lambda_u, K/Z);
    end
    if isinf(lambda_u)
        line_error(mfilename,'All service demands are zero.');
    end
    lambda = lambda_u;
    for it=1:maxiter
        lambda = (lambda_l+lambda_u)/2;
        g = lambda*Z + sum(sum_node_qlen(L,lambda,mi,scv,K));
        if abs(g-K)<=tol || (lambda_u-lambda_l)<=tol*lambda_u
            break
        end
        if g>K
            lambda_u = lambda;
        else
            lambda_l = lambda;
        end
    end
    XN = lambda;
else
    % multiclass: Gauss-Seidel sweeps of per-class bisections on the
    % population constraints sum_i fir(lambda_r)+lambda_r*Zr = Nr. This is
    % a robust alternative to the successive substitution of Sec. 9.2.2,
    % which can overshoot the saturation polytope for large populations.
    XN = zeros(1,R);
    for it=1:maxiter
        delta = 0;
        for r=1:R
            if N(r)==0
                continue
            end
            % upper bound for lambda_r given the other class throughputs
            ub = Inf;
            for i=1:M
                if L(i,r)>0
                    if isinf(mi(i))
                        ub = min(ub, K/L(i,r));
                    else
                        rem = mi(i) - (XN*L(i,:)' - XN(r)*L(i,r));
                        ub = min(ub, max(rem,0)/L(i,r));
                    end
                end
            end
            if Z(r)>0
                ub = min(ub, N(r)/Z(r));
            end
            if isinf(ub)
                line_error(mfilename,'All service demands are zero.');
            end
            lambda_old = XN(r);
            lambda_l = 0;
            lambda_u = ub;
            while (lambda_u-lambda_l) > tol*max(ub,1)/1e3
                lambda = (lambda_l+lambda_u)/2;
                XN(r) = lambda;
                Qir = sum_node_qlen(L,XN,mi,scv,K);
                g = lambda*Z(r) + sum(Qir(:,r));
                if g>N(r)
                    lambda_u = lambda;
                else
                    lambda_l = lambda;
                end
            end
            XN(r) = (lambda_l+lambda_u)/2;
            delta = max(delta, abs(XN(r)-lambda_old));
        end
        if delta<=tol
            break
        end
    end
end

QN = sum_node_qlen(L,XN,mi,scv,K);
for i=1:M
    for r=1:R
        if isinf(mi(i))
            UN(i,r) = XN(r)*L(i,r);
        else
            UN(i,r) = XN(r)*L(i,r)/mi(i);
        end
        if XN(r)>0
            RN(i,r) = QN(i,r)/XN(r);
        end
    end
end
end

function Qir = sum_node_qlen(L,XN,mi,scv,K)
% per-station per-class mean queue lengths Ki_r = fir(lambda_r)
[M,R] = size(L);
Qir = zeros(M,R);
for i=1:M
    if isinf(mi(i))
        Qir(i,:) = XN.*L(i,:); % Type 3, Eq. (9.15)
        continue
    end
    m = mi(i);
    Uir = XN.*L(i,:); % class offered loads lambda_r*e_ir/mu_ir
    rho = min(sum(Uir)/m, 1); % per-server utilization; the correction
    % factors keep the node functions finite at rho=1 (Ki(1)<=K)
    if sum(Uir)==0
        continue
    end
    % demand-weighted node service SCV
    ci2 = sum(Uir.*scv(i,:))/sum(Uir);
    ai = (1+ci2)/2;
    if K<=m
        % never more than m jobs at a m-server node: no queueing
        Qir(i,:) = Uir;
        continue
    end
    if m==1
        if ci2==1 || K<=1
            % Type 1,2,4 with mi=1, Eq. (9.15)/(9.19)
            Qir(i,:) = Uir/(1-(K-1)/K*rho);
        else
            % -/G/1 FCFS, Eq. (10.88)
            den = 1-(K-1-ai)/(K-1)*rho;
            Qir(i,:) = Uir.*(1+rho*ai/den);
        end
    else
        Pm = sum_erlangc(m,rho);
        if ci2==1
            % Type 1 with mi>1, Eq. (9.15)/(9.19)
            den = 1-(K-m-1)/(K-m)*rho;
            Qir(i,:) = Uir + (Uir/m)*Pm/den;
        else
            % -/G/m FCFS, Eq. (10.89)
            den = 1-(K-m-ai)/(K-m)*rho;
            Qir(i,:) = Uir + (Uir/m)*ai*Pm/den;
        end
    end
end
end

function Pm = sum_erlangc(m,rho)
% Erlang-C probability of waiting for an M/M/m queue (Eq. 6.28)
if rho>=1
    Pm = 1;
    return
end
a = m*rho;
s = 0;
for k=0:(m-1)
    s = s + a^k/factorial(k);
end
last = a^m/(factorial(m)*(1-rho));
Pm = last/(s+last);
end
