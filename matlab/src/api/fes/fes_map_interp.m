%{ @file fes_map_interp.m
 %  @brief Shape-preserving interpolation of the flow-equivalent descriptors
 %
 %  @author LINE Development Team
%}

%{
 % @brief Monotone piecewise cubic Hermite interpolation
 %
 % @details
 % Interpolates the descriptors of a MAP flow-equivalent server between the
 % populations at which they were evaluated. Fritsch and Carlson slopes are
 % used, with the noncentered three-point endpoint rule of de Boor, so the
 % interpolant never overshoots and a monotone sequence of throughputs
 % stays monotone. The algorithm is written out rather than delegated to
 % the built-in pchip so that the MATLAB, Java, Python and C++ ports return
 % identical values.
 %
 % @par Syntax:
 % @code
 % yq = fes_map_interp(x, y, xq)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>x<td>Sample abscissae, strictly increasing
 % <tr><td>y<td>Sample values, one row per abscissa
 % <tr><td>xq<td>Query abscissae
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>yq<td>Interpolated values, one row per query point
 % </table>
 %
 % @see fes_map_aggregate
%}
function yq = fes_map_interp(x, y, xq)

x = x(:);
xq = xq(:);
n = numel(x);
if size(y,1) ~= n
    y = y.';
end
ncol = size(y,2);
yq = zeros(numel(xq), ncol);

if n == 1
    yq = repmat(y, numel(xq), 1);
    return
end

h = diff(x);
for c = 1:ncol
    v = y(:,c);
    delta = diff(v)./h;
    d = zeros(n,1);

    if n == 2
        d(:) = delta(1);
    else
        for i = 2:n-1
            if delta(i-1)*delta(i) > 0
                w1 = 2*h(i) + h(i-1);
                w2 = h(i) + 2*h(i-1);
                d(i) = (w1+w2)/(w1/delta(i-1) + w2/delta(i));
            end
        end
        d(1) = fes_map_interp_edge(h(1), h(2), delta(1), delta(2));
        d(n) = fes_map_interp_edge(h(n-1), h(n-2), delta(n-1), delta(n-2));
    end

    for q = 1:numel(xq)
        t = xq(q);
        if t <= x(1)
            i = 1;
        elseif t >= x(n)
            i = n-1;
        else
            i = find(x <= t, 1, 'last');
            i = min(i, n-1);
        end
        s = t - x(i);
        hi = h(i);
        c2 = (3*delta(i) - 2*d(i) - d(i+1))/hi;
        c3 = (d(i) - 2*delta(i) + d(i+1))/hi^2;
        yq(q,c) = v(i) + s*(d(i) + s*(c2 + s*c3));
    end
end
end
