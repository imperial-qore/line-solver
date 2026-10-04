%{ @file map2_fit_idc.m
 %  @brief Fits a MAP(2) to three moments and an index of dispersion
 %
 %  @author LINE Development Team
%}

%{
 % @brief Fits a second-order MAP matching the first three moments and the
 % asymptotic index of dispersion
 %
 % @details
 % A MAP(2) has a geometrically decaying autocorrelation, so its index of
 % dispersion obeys
 %
 %   I = SCV + (SCV-1)*g2/(1-g2)
 %
 % as reported in Section 5.2.2 of Casale, Mi, Cherkasova and Smirni, IEEE
 % Trans. Soft. Eng. 37(5), 2011. The relation is inverted in closed form
 % as g2 = (I-SCV)/(I-1) and the resulting decay rate is passed to
 % map2_fit, which is the explicit inverse characterization of Heindl,
 % Horvath and Gross. A third moment outside the feasible region is
 % replaced by its lower limit (3/2)*e2^2/e1, the largest heavy-tail decay
 % a MAP(2) admits.
 %
 % The paper returns an exponential whenever SCV <= 1 or I < SCV, on the
 % grounds that burstiness is then negligible. The rule does more than
 % avoid an infeasible fit and must not be relaxed: a flow-equivalent
 % server whose service is exponential and load dependent is exact for a
 % product-form subnetwork by Norton's theorem, whereas any MAP(2) fitted
 % to the marginal inter-departure statistics is not, because the departure
 % stream of the subnetwork is not independent of the rest of the model.
 % Fitting the sub-exponential SCV of a non-bursty aggregate was measured
 % to cost up to 2.2% of throughput on a three-station exponential network
 % that the exponential fallback reproduces exactly.
 %
 % @par Syntax:
 % @code
 % [MAP,status] = map2_fit_idc(e1,e2,e3,I)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>e1<td>Mean inter-arrival time
 % <tr><td>e2<td>Second moment of the inter-arrival times
 % <tr><td>e3<td>Third moment of the inter-arrival times
 % <tr><td>I<td>Asymptotic index of dispersion
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>MAP<td>Fitted process in the form {D0,D1}
 % <tr><td>status<td>0 all four descriptors matched, 1 exponential as
 %                   burstiness is not representable, 2 third moment
 %                   clamped, 3 third moment selected automatically,
 %                   4 fit failed and an exponential is returned
 % </table>
 %
 % @see map2_fit, map_idc, fes_map_moments
%}
function [MAP,status] = map2_fit_idc(e1,e2,e3,I)

scv = (e2-e1^2)/e1^2;

if scv <= 1 + GlobalConstants.FineTol || I < scv
    MAP = map_exponential(e1);
    status = 1;
    return
end

g2 = (I-scv)/(I-1);

[MAP,ERR] = map2_fit(e1,e2,e3,g2);
if ERR == 0
    status = 0;
    return
end

e3min = (3/2 + 1e-6)*e2^2/e1;
if e3 < e3min
    [MAP,ERR] = map2_fit(e1,e2,e3min,g2);
    if ERR == 0
        status = 2;
        return
    end
end

[MAP,ERR] = map2_fit(e1,e2,-1,g2);
if ERR == 0
    status = 3;
    return
end

MAP = map_exponential(e1);
status = 4;
end
