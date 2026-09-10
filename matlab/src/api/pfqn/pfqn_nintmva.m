%{
%{
 % @file pfqn_nintmva.m
 % @brief Mean value analysis at a nonintegral population (fractional-base aMVA).
%}
%}

function [X,Q,U,R] = pfqn_nintmva(L,N,Z)
%{
%{
 % @brief Exact MVA recursion started from the FRACTIONAL base
 %        n_0 = N - floor(N) instead of from the empty network, giving mean
 %        performance measures of a single-class closed product-form network at
 %        a real-valued population (Dowdy and Gordon 1984, "aMVA").
 %
 %        The recursion is the standard Reiser-Lavenberg one,
 %          R_i(n) = D_i (1 + Q_i(n-1)),  X(n) = n/(Z + sum_i R_i(n)),
 %          Q_i(n) = X(n) R_i(n),
 %        stepped in unit increments from n = n_0 (where the arrival theorem
 %        term Q_i(n_0 - 1) is taken as 0, the network below the base being
 %        empty) up to n = N. At integer N the base is 0 and the recursion is
 %        bit-identical to exact MVA; at fractional N it interpolates smoothly
 %        through the integral points, which is what a nonintegral degree of
 %        multiprogramming (a time-average over a measurement window) calls for.
 %
 %        Unlike pfqn_dnc this accepts a think time, since the delay enters the
 %        recursion and not a partial-fraction continuation. It is single-class:
 %        the multiclass recursion has no one-dimensional step. For fractional
 %        multiclass populations use pfqn_bs, which accepts them directly.
 %
 %        Reference: L. W. Dowdy, K. D. Gordon, "Algorithms for Nonintegral
 %        Degrees of Multiprogramming in Closed Queuing Networks", Performance
 %        Evaluation 4(1):19-28, 1984.
 % @fn pfqn_nintmva(L, N, Z)
 % @param L Service demand vector (M x 1) of the queueing stations.
 % @param N Population (real nonnegative scalar; may be fractional).
 % @param Z Think time (scalar, default 0).
 % @return X Throughput at population N.
 % @return Q Mean queue lengths (M x 1).
 % @return U Utilizations (M x 1).
 % @return R Residence times (M x 1).
%}
%}
if nargin < 3 || isempty(Z), Z = 0; end
if size(L,2) > 1 && size(L,1) > 1
    line_error(mfilename,'pfqn_nintmva is a single-class method, but the demand matrix has more than one class. Use pfqn_bs for fractional multiclass populations.');
end
if numel(N) > 1
    line_error(mfilename,'pfqn_nintmva is a single-class method, but the population vector has more than one entry.');
end
L = L(:);
Z = sum(Z(:));
if N < 0
    line_error(mfilename,'pfqn_nintmva requires a nonnegative population.');
end

M = numel(L);
Q = zeros(M,1);
X = 0;
R = zeros(M,1);
if N == 0
    U = zeros(M,1);
    return
end

n = N - floor(N);
if n == 0, n = 1; end
while n <= N + 1e-12
    R = L.*(1+Q);
    X = n/(Z + sum(R));
    Q = X*R;
    n = n + 1;
end
U = X*L;
end
