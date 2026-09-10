%{
%{
 % @file sjn_station.m
 % @brief Conditional waiting time recursion at one shortest-job-next station.
%}
%}

function [Cmr, W, phi, phiinf, tailpar] = sjn_station(caller, m, r, G, S, scv, V, st, beta, useprio, prio)
%{
%{
 % @brief One evaluation of the SJN conditional waiting time equation for a
 %        tagged customer of class r at station m, given the state of the
 %        station at the reference population and the deflation factors that
 %        map that state onto the population the tagged customer sees.
 %
 %        The quantity lam_k W_k(x) f_k(x) is the density, in the job size x,
 %        of the queued class-k customers, so deflating it by beta_k is what
 %        turns the same equation into either the exact recursion (beta = 1,
 %        the state already being the one at n - e_r) or the Schweitzer
 %        closure of pfqn_amvasjn (beta_r = (N_r-1)/N_r, beta_k = 1 for
 %        k /= r, the state being the one at N). The tail beyond the grid is
 %        closed analytically by W(x) = a - b exp(-c (x - Lx)).
 % @fn sjn_station(caller, m, r, G, S, scv, V, st, beta, useprio, prio)
 % @param caller Name of the calling function, used in error messages.
 % @param m Station index, used in error messages.
 % @param r Index of the tagged class.
 % @param G Grid struct of the station, as returned by sjn_setup.
 % @param S Mean service times at the station, one per class (1 x R).
 % @param scv Squared coefficients of variation at the station (1 x R).
 % @param V Visit ratios at the station (1 x R).
 % @param st Struct with the station state at the reference population:
 %        lam, U, Q (1 x R each), W and phi (ngrid x R) and phiinf (1 x R).
 % @param beta Deflation factors applied to the queued and arriving work (1 x R).
 % @param useprio True for the priority reading, false for the pooled one.
 % @param prio Priority levels, one per class, lower is higher priority.
 % @return Cmr Residence time of class r at the station.
 % @return W Conditional waiting time profile on the grid (ngrid x 1).
 % @return phi The primitive int_0^x W(t) t f_r(t) dt on the grid (ngrid x 1).
 % @return phiinf The same primitive at infinity, tail included.
 % @return tailpar Tail parameters [a, b, c] of the profile.
%}
%}
R = length(S);
lamb = beta .* st.lam;
Ub = beta .* st.U;
Qb = beta .* st.Q;
RL = sum((1 + scv) .* S .* Ub) / 2;
if useprio
    hi = prio < prio(r);
    base = RL + sum(S(hi) .* (Qb(hi) - Ub(hi)));
    num = base + lamb(r) * st.phi(:,r)';
    den = 1 - sum(Ub(hi)) - lamb(r) * G.theta(:,r)';
    numinf = base + lamb(r) * st.phiinf(r);
    deninf = 1 - sum(Ub(hi)) - lamb(r) * S(r);
else
    num = RL + lamb * st.phi';
    den = 1 - lamb * G.theta';
    numinf = RL + lamb * st.phiinf';
    deninf = 1 - sum(lamb .* S);
end
if any(den <= 0) || deninf <= 0
    line_error(caller,sprintf(['the SJN recursion at station %d has no solution: the work brought by jobs no longer\n' ...
        'than the tagged one saturates the server, at which point long jobs starve and the arrival theorem\n' ...
        'no longer holds. Reduce the load at that station or model it with SolverCTMC or SolverLDES.'],m));
end
W = num ./ den;
Winf = numinf / deninf;
% eq. (14) of the reference, generalised by differentiating the recursion at the grid edge
if useprio
    slope = G.Lx * lamb(r) * G.f(end,r) * (st.W(end,r) + W(end)) / den(end);
else
    slope = 0;
    for k = 1:R
        slope = slope + lamb(k) * G.f(end,k) * (st.W(end,k) + W(end));
    end
    slope = G.Lx * slope / den(end);
end
a = Winf;
b = Winf - W(end);
if b <= 0
    b = 0;
    c = 0;
elseif slope < 0
    line_error(caller,sprintf(['the conditional waiting time at SJN station %d decreases in the job size, which the\n' ...
        'discipline forbids: the recursion has become numerically unstable.'],m));
else
    c = slope / b;
end
tailpar = [a, b, c];
phi = sjn_quad('cumsimpson', W .* (G.x .* G.f(:,r)'), G.dx);
phiinf = phi(end) + a * G.tail1(r) - b * sjn_quad('tailmom', G.fit{r}, G.Lx, c, 1);
Wbar = sjn_quad('simpson', W .* G.f(:,r)', G.dx) + a * G.tail0(r) - b * sjn_quad('tailmom', G.fit{r}, G.Lx, c, 0);
Cmr = V(r) * (S(r) + Wbar);
W = W(:);
phi = phi(:);
end
