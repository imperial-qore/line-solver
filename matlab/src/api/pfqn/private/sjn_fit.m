%{
%{
 % @file sjn_fit.m
 % @brief Two-moment Erlang-mixture fit of a service time distribution.
%}
%}

function f = sjn_fit(s, cv2)
%{
%{
 % @brief Fit an Erlang mixture to a given mean and squared coefficient of
 %        variation, as the SJN response time equations require: a branching
 %        Erlang (Erlang(k-1) and Erlang(k) sharing a rate) when CV^2 < 1 and
 %        a balanced-means hyperexponential when CV^2 > 1. The mixture form
 %        makes theta(x) = int_0^x t f(t) dt and the tail integrals closed
 %        form, which is what lets the SJN quadrature run on a fixed grid.
 % @fn sjn_fit(s, cv2)
 % @param s Mean service time.
 % @param cv2 Squared coefficient of variation of the service time.
 % @return f Struct with the mixture weights w, phase counts k and rates mu.
%}
%}
tol = 1e-8;
if s <= 0
    f = struct('w',[],'k',[],'mu',[]);
    return
end
if cv2 < 0
    line_error(mfilename,'negative squared coefficient of variation');
end
if abs(cv2 - 1) < tol
    f = struct('w',1,'k',1,'mu',1/s);
elseif cv2 < 1
    k = ceil(1/cv2);
    p = (k*cv2 - sqrt(k*(1+cv2) - k^2*cv2)) / (1 + cv2);
    mu = (k - p) / s;
    f = struct('w',[p, 1-p],'k',[k-1, k],'mu',[mu, mu]);
else
    p = 0.5 * (1 + sqrt((cv2-1)/(cv2+1)));
    f = struct('w',[p, 1-p],'k',[1, 1],'mu',[2*p/s, 2*(1-p)/s]);
end
end
