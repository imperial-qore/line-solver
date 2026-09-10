%{ @file cache_rrm_meanfield_ode.m
 %  @brief ODE function for RRM mean-field cache dynamics
 %
 %  @author LINE Development Team
%}

%{
 % @brief Computes ODE dynamics for RRM mean-field cache model
 %
 % @details
 % This function defines the ODE system for the Random Replacement Model (RRM)
 % mean-field cache dynamics.
 %
 % @par Syntax:
 % @code
 % dxdt = cache_rrm_meanfield_ode(t, x, lambda, m, n, h)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>t<td>Time variable (unused)
 % <tr><td>x<td>State vector
 % <tr><td>lambda<td>Arrival rates per item
 % <tr><td>m<td>Cache capacity vector
 % <tr><td>n<td>Number of items
 % <tr><td>h<td>Number of cache levels
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>dxdt<td>Time derivative of state vector
 % </table>
%}
function dxdt = cache_rrm_meanfield_ode(~, x, lambda, m, n, h)
% Reshape state vector to matrix form: x_{k,s}
x = reshape(x, [n, 1+ h]);

% TODO check list 0 dxdt
dxdt = zeros(n, 1+ h);

for k = 1:n
    for s = 1:h

        % First term: promotion from list s-1
        sum1 = 0;
        for k1 = 1:n
            sum1 = sum1 + lambda(k1) / m(s) * x(k1,1+  s-1) * x(k, 1+ s);
        end

        % Second term: demotion from list s+1
        if s < h
            sum2 = 0;
            for k1 = 1:n
                sum2 = sum2 + lambda(k1) / m(s+1) * x(k1,1+  s) * x(k, 1+ s+1);
            end
            sum2 = sum2 - lambda(k) * x(k,1+  s);
        else
            sum2 = 0;
        end

        % Drift component
        dxdt(k, 1+s) = lambda(k) * x(k, 1+ s-1) - sum1 + sum2;
    end
    % case s=0
    dxdt(k, 1) = - sum(dxdt(k, 2:(1+h)));
end

dxdt = dxdt(:);  % Flatten to vector
end
