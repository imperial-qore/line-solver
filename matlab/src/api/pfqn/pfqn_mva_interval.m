%{
%{
 % @file pfqn_mva_interval.m
 % @brief Exact interval-valued MVA for single-class closed product-form networks.
%}
%}

function [X,Q,U,R,Rtot,Qtot] = pfqn_mva_interval(L,N,Z)
%{
%{
 % @brief Exact output intervals of single-class MVA when the service demands,
 %        the think time and the population are known only up to intervals.
 %
 %        Single-class MVA is monotone in every input: the throughput decreases
 %        in each demand and in the think time and increases in the population,
 %        the per-station queue length and residence time increase in the own
 %        demand and in the population and decrease in the other demands and in
 %        the think time, and the totals increase in every demand and in the
 %        population and decrease in the think time (Luthi and Haring 1998,
 %        Theorems 2-5, Table 1). By Theorem 1 the exact range of a function
 %        monotone in each argument is attained at the endpoints of the input
 %        box, so each bound below is one ordinary MVA call at the corner that
 %        the sign pattern selects. This is the algorithm of their Fig. 2 and it
 %        costs 2*(m+2) MVA calls, m being the number of thick demand intervals;
 %        evaluating the MVA recursion in interval arithmetic instead would be a
 %        valid but far wider enclosure, since every input recurs at each step
 %        (the dependency problem, 14x too wide on the paper's own example).
 %
 %        The returned interval is the exact hull of MVA over the input box, not
 %        a bound on the true network: it holds conditionally on the demands
 %        lying in the box, and says nothing about the accuracy of MVA itself.
 %        It must therefore not be composed with the brackets of SolverBA, which
 %        bracket the exact solution of a model whose demands are known.
 %
 %        Delay stations are folded into Z, exactly as in pfqn_mva: a delay
 %        demand interval enters as a term of the think-time interval, and the
 %        hull of the sum is the sum of the hulls when the delays vary
 %        independently. Load-independent single-server queueing stations only,
 %        one class only; the monotonicity theorems cover no other case.
 %
 %        Reference: J. Luthi, G. Haring, "Mean value analysis for queueing
 %        network models with intervals as input parameters", Performance
 %        Evaluation 32(3):185-215, 1998.
 % @fn pfqn_mva_interval(L, N, Z)
 % @param L Service demand intervals (M x 2), column 1 lower, column 2 upper.
 % @param N Population interval ([nlo nup], or a scalar for a thin population).
 % @param Z Think time interval ([zlo zup], or a scalar; default 0).
 % @return X Throughput interval (1 x 2).
 % @return Q Mean queue-length intervals per station (M x 2).
 % @return U Utilization enclosures per station (M x 2), exact where the demand
 %           interval is thin, capped at 1 elsewhere.
 % @return R Residence-time intervals per station (M x 2).
 % @return Rtot Total response-time interval (1 x 2).
 % @return Qtot Interval of the total number of jobs at the stations (1 x 2).
%}
%}
if nargin < 3 || isempty(Z), Z = 0; end
if isempty(L)
    line_error(mfilename,'pfqn_mva_interval requires at least one queueing station.');
end
if size(L,2) == 1
    L = [L(:), L(:)];
elseif size(L,2) ~= 2
    line_error(mfilename,'pfqn_mva_interval is a single-class method: L must be M x 2, one [lower upper] demand interval per station.');
end
N = N(:)'; if isscalar(N), N = [N N]; end
Z = Z(:)'; if isscalar(Z), Z = [Z Z]; end
if numel(N) ~= 2 || numel(Z) ~= 2
    line_error(mfilename,'population and think time must be given as scalars or as [lower upper] intervals.');
end
if any(L(:) < 0) || any(Z < 0)
    line_error(mfilename,'demands and think times must be nonnegative.');
end
if any(L(:,1) > L(:,2)) || N(1) > N(2) || Z(1) > Z(2)
    line_error(mfilename,'interval lower endpoints must not exceed the upper endpoints.');
end
if N(1) < 1
    line_error(mfilename,'pfqn_mva_interval requires a population interval with at least one job; the monotonicity theorems assume n >= 1.');
end
if any(abs(N - round(N)) > GlobalConstants.Zero)
    line_error(mfilename,'population endpoints must be integers; use pfqn_nintmva for a nonintegral population.');
end

M = size(L,1);
Llo = L(:,1); Lup = L(:,2);
nlo = round(N(1)); nup = round(N(2));
zlo = Z(1); zup = Z(2);
thick = Lup > Llo + GlobalConstants.Zero;

X = zeros(1,2); Q = zeros(M,2); U = zeros(M,2); R = zeros(M,2);
Rtot = zeros(1,2); Qtot = zeros(1,2);

% S1: throughput upper bound, and the upper bounds of the stations whose demand
% is thin (their own demand is fixed, so lowering the others maximizes them).
[x1,q1,~,c1] = pfqn_mva(Llo, nup, zlo);
X(2) = x1;
Q(~thick,2) = q1(~thick);
R(~thick,2) = c1(~thick);

% S2: the same quantities at the opposite corner, giving the lower bounds.
[x2,q2,~,c2] = pfqn_mva(Lup, nlo, zup);
X(1) = x2;
Q(~thick,1) = q2(~thick);
R(~thick,1) = c2(~thick);

% S3/S4: the totals increase in the demands and the population and decrease in
% the think time, so their corners differ from those of the throughput.
[~,q3,~,c3] = pfqn_mva(Llo, nlo, zup);
Rtot(1) = sum(c3); Qtot(1) = sum(q3);
[~,q4,~,c4] = pfqn_mva(Lup, nup, zlo);
Rtot(2) = sum(c4); Qtot(2) = sum(q4);

% S5/S6: one pair of calls per thick station, its own demand at the endpoint
% that maximizes (minimizes) it and the others at the opposite endpoint.
for k = find(thick)'
    d = Llo; d(k) = Lup(k);
    [~,q5,~,c5] = pfqn_mva(d, nup, zlo);
    Q(k,2) = q5(k); R(k,2) = c5(k);
    d = Lup; d(k) = Llo(k);
    [~,q6,~,c6] = pfqn_mva(d, nlo, zup);
    Q(k,1) = q6(k); R(k,1) = c6(k);
end

% U = X*D is not covered by the monotonicity table, so it is enclosed by the
% product of the two intervals, intersected with the range of a single-server
% utilization. Where the demand is thin the product is already exact.
U(:,1) = X(1)*Llo;
U(:,2) = min(1, X(2)*Lup);
end
