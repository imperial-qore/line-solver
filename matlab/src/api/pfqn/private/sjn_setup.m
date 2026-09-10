%{
%{
 % @file sjn_setup.m
 % @brief Grid and precomputed size-distribution integrals for an SJN station.
%}
%}

function G = sjn_setup(S, scv, ns, Lfactor)
%{
%{
 % @brief Build the job-size grid of one shortest-job-next station and the
 %        quantities of the size distribution that do not change with the
 %        population: the density, theta(x) = int_0^x t f(t) dt, and the two
 %        tail masses beyond the grid. The grid spans [0, Lfactor * max_r s_r]
 %        because the conditional waiting time has flattened out well before
 %        that point, its remainder being carried by the analytic tail of the
 %        recursion rather than by quadrature.
 % @fn sjn_setup(S, scv, ns, Lfactor)
 % @param S Mean service times at the station, one per class (1 x R).
 % @param scv Squared coefficients of variation, one per class (1 x R).
 % @param ns Number of grid subdivisions, even.
 % @param Lfactor Grid extent in units of the largest mean service time.
 % @return G Struct with fields x, dx, Lx, f, theta, tail0, tail1 and fit.
%}
%}
R = length(S);
smax = max(S);
if smax <= 0
    line_error(mfilename,'the station has zero service demand in every class');
end
G.Lx = Lfactor * smax;
G.x = linspace(0, G.Lx, ns+1);
G.dx = G.x(2) - G.x(1);
G.f = zeros(ns+1,R);
G.theta = zeros(ns+1,R);
G.tail0 = zeros(1,R);
G.tail1 = zeros(1,R);
G.fit = cell(1,R);
for r = 1:R
    G.fit{r} = sjn_fit(S(r), scv(r));
    if isempty(G.fit{r}.w)
        continue
    end
    G.f(:,r) = sjn_quad('pdf', G.fit{r}, G.x)';
    G.theta(:,r) = sjn_quad('theta', G.fit{r}, G.x)';
    G.tail0(r) = sjn_quad('ccdf', G.fit{r}, G.Lx);
    G.tail1(r) = S(r) - G.theta(end,r);
end
end
