%{
%{
 % @file pfqn_dnc.m
 % @brief Distinct-load Normalizing Constant (DNC) at a nonintegral population.
%}
%}

function [X,G,lG] = pfqn_dnc(L,N)
%{
%{
 % @brief Normalizing constant and throughput of a single-class closed
 %        product-form network at a REAL-VALUED population, by partial-fraction
 %        inversion of the network generating function (Dowdy and Gordon 1984).
 %
 %        With distinct loads x_1..x_G of multiplicities m_1..m_G the generating
 %        function prod_g (1-x_g u)^{-m_g} expands as
 %          G(n) = sum_g sum_{j=1..m_g} A_gj C(n+j-1,j-1) x_g^n,
 %        every term of which is an analytic function of n. Evaluating it at a
 %        real n therefore interpolates the integral normalizing constants
 %        exactly (it reproduces them at every integer) and gives a smooth
 %        throughput curve X(N) = G(N-1)/G(N) through the integral points,
 %        rather than the rounding or linear interpolation the paper compares
 %        against. For all-distinct loads the coefficients have the closed form
 %          A_g = prod_{l~=g} x_g/(x_g - x_l),
 %        used directly; with repeated loads they are recovered from the M
 %        integral constants G(0..M-1), which determine them uniquely.
 %
 %        Only the queueing part admits this continuation: the delay sequence
 %        Z^n/n! is entire and has no partial-fraction expansion, so a think
 %        time is rejected here. Use pfqn_nintmva for nonintegral populations
 %        with a delay.
 %
 %        Reference: L. W. Dowdy, K. D. Gordon, "Algorithms for Nonintegral
 %        Degrees of Multiprogramming in Closed Queuing Networks", Performance
 %        Evaluation 4(1):19-28, 1984.
 % @fn pfqn_dnc(L, N)
 % @param L Service demand vector (M x 1) of the queueing stations.
 % @param N Population (real nonnegative scalar; may be fractional).
 % @return X Throughput G(N-1)/G(N).
 % @return G Normalizing constant at population N.
 % @return lG Logarithm of the normalizing constant.
%}
%}
if size(L,2) > 1 && size(L,1) > 1
    line_error(mfilename,'pfqn_dnc is a single-class method, but the demand matrix has more than one class.');
end
if numel(N) > 1
    line_error(mfilename,'pfqn_dnc is a single-class method, but the population vector has more than one entry.');
end
L = L(:);
L = L(L > 0);
if isempty(L)
    line_error(mfilename,'pfqn_dnc requires at least one station with positive demand.');
end
if N < 0
    line_error(mfilename,'pfqn_dnc requires a nonnegative population.');
end

M = numel(L);
xmax = max(L);
y = L/xmax;                                     % scaled loads, max(y) = 1

% Distinct loads and multiplicities, merged under a relative tolerance so that
% numerically coincident loads are handled by the multiplicity branch rather
% than by a near-singular partial-fraction denominator.
ys = sort(y,'descend');
u = ys(1); mult = 1;
for i = 2:M
    if ys(i) > u(end)*(1-1e-9)
        mult(end) = mult(end)+1; %#ok<AGROW>
    else
        u(end+1) = ys(i); %#ok<AGROW>
        mult(end+1) = 1; %#ok<AGROW>
    end
end
u = u(:); mult = mult(:);
Gd = numel(u);

if all(mult == 1)
    % Closed form: A_g = prod_{l~=g} u_g/(u_g - u_l)
    A = ones(Gd,1);
    for g = 1:Gd
        idx = [1:g-1, g+1:Gd];
        A(g) = prod(u(g)./(u(g)-u(idx)));
    end
    j = ones(Gd,1);
    node = (1:Gd)';
else
    % Repeated loads: recover the coefficients from G(0..M-1), computed by
    % convolution on the scaled loads (bounded by construction since max u = 1).
    gint = 1;
    for i = 1:M
        gi = y(i).^(0:M-1);
        gint = conv(gint, gi);
        gint = gint(1:M);
    end
    node = zeros(M,1); j = zeros(M,1); c = 0;
    for g = 1:Gd
        for jj = 1:mult(g)
            c = c+1; node(c) = g; j(c) = jj;
        end
    end
    F = zeros(M,M);
    for n = 0:M-1
        F(n+1,:) = exp(gammaln(n+j') - gammaln(j') - gammaln(n+1) + n*log(u(node))');
    end
    A = F\gint(:);
end

GN  = dncEval(N,   A, u, node, j);
GN1 = dncEval(N-1, A, u, node, j);

lG = log(GN) + N*log(xmax);
G = exp(lG);
if N <= 0 || GN <= 0 || isnan(GN1)
    X = NaN;
else
    X = (GN1/GN)/xmax;
end
end

function g = dncEval(n, A, u, node, j)
% Partial-fraction series evaluated at a real population n. The continuation is
% analytic for n > -1; below that the binomial factor changes sign and the
% log-domain evaluation would lose it, so it is not extended there.
if n <= -1
    g = NaN;
    return
end
g = sum(A .* exp(gammaln(n+j) - gammaln(j) - gammaln(n+1) + n*log(u(node))));
end
