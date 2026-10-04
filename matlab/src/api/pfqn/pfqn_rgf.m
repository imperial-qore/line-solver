%{
%{
 % @file pfqn_rgf.m
 % @brief Recursion by Generating Functions (RGF) for the normalizing constant
 %        of single-class closed product-form networks with replicated stations.
%}
%}

function [G,lG,lg] = pfqn_rgf(L,N,Z)
%{
%{
 % @brief Exact normalizing constant of a single-class closed product-form
 %        network obtained by convolving the per-node generating-function
 %        sequences of Coury and Harrison (1997), Property 1, instead of the
 %        per-station Buzen recursion.
 %
 %        Each node contributes the coefficient sequence of its own generating
 %        function, and a GROUP of m stations sharing the same demand p is
 %        collapsed into the single negative-binomial sequence
 %          r(k) = C(k+m-1,k) p^k,
 %        i.e. the whole group costs one sequence rather than m convolution
 %        passes. The delay contributes the Poisson sequence r(k) = Z^k/k!.
 %        Convolving the G distinct sequences gives g(0..N) exactly, and
 %        lG = log g(N), X(n) = g(n-1)/g(n).
 %
 %        Cost O(G N^2) against Buzen's O(M N), so RGF is the cheaper route
 %        precisely when the model is heavily replicated and the population is
 %        moderate (G N < M); it is otherwise kept for its exact group closed
 %        form. The whole recursion runs in the log domain, so no intermediate
 %        overflows or underflows are possible.
 %
 %        Reference: J. Coury, P. G. Harrison, "Asymptotic properties of
 %        queuing networks", IEE Proc.-Comput. Digit. Tech. 144(5):247-254,
 %        1997 (Property 1 and the five-sequence decomposition of Sec. 4).
 % @fn pfqn_rgf(L, N, Z)
 % @param L Service demand vector (M x 1) of the queueing stations.
 % @param N Population (nonnegative integer scalar).
 % @param Z Think time (scalar, default 0). Aggregated delay demand.
 % @return G Normalizing constant.
 % @return lG Logarithm of the normalizing constant.
 % @return lg Logarithms of g(0), g(1), ..., g(N) (1 x (N+1) vector).
%}
%}
if nargin < 3 || isempty(Z), Z = 0; end
if size(L,2) > 1 && size(L,1) > 1
    line_error(mfilename,'pfqn_rgf is a single-class method, but the demand matrix has more than one class. Use pfqn_ca or pfqn_nc for multiclass models.');
end
if numel(N) > 1
    line_error(mfilename,'pfqn_rgf is a single-class method, but the population vector has more than one entry.');
end
L = L(:);
Z = sum(Z(:));
if N < 0 || abs(N-round(N)) > 0
    line_error(mfilename,'pfqn_rgf requires a nonnegative integer population.');
end
N = round(N);
if any(L < 0) || Z < 0
    line_error(mfilename,'pfqn_rgf requires nonnegative demands and think time.');
end

kk = 0:N;
lg = zeros(1,N+1);              % g(k) = 1 for the empty network
lg(2:end) = -Inf;

% Delay node: Poisson sequence Z^k/k!
if Z > 0
    lg = logconv(lg, kk*log(Z) - gammaln(kk+1));
end

% Queueing stations grouped by identical demand: one sequence per group
L = L(L > 0);
if ~isempty(L)
    [p, ~, grp] = unique(L);
    m = accumarray(grp, 1);
    for i = 1:numel(p)
        if m(i) == 1
            lr = kk*log(p(i));                                  % 1/(1-p u)
        else
            lr = gammaln(kk+m(i)) - gammaln(kk+1) - gammaln(m(i)) + kk*log(p(i));
        end
        lg = logconv(lg, lr);
    end
end

lG = lg(end);
G = exp(lG);
end

function c = logconv(u,v)
% Log-domain linear convolution truncated at the common length.
n = numel(u);
c = -inf(1,n);
for k = 1:n
    t = u(1:k) + v(k:-1:1);
    vm = max(t);
    if isinf(vm)
        c(k) = vm;
    else
        c(k) = vm + log(sum(exp(t-vm)));
    end
end
end
