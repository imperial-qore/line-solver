%{ @file dist_fit_me.m
 %  @brief Fits a matrix exponential to a given mean and SCV
 %
 %  @author LINE Development Team
%}

%{
 % @brief Fits a matrix exponential to a mean and squared coefficient of variation
 %
 % @details
 % For SCV < 1 the fit is the convolution X = c*Y + Z of a scaled concentrated
 % matrix exponential Y (unit mean, minimal SCV sY for its order) with an
 % independent exponential Z. Writing c + d = mean and c^2*sY + d^2 = scv*mean^2,
 %
 %   c = mean*(1 - sqrt(1 - (1+sY)*(1-scv)))/(1 + sY),   d = mean - c,
 %
 % so every target in [sY/(1+sY), 1] is matched EXACTLY in 2n+2 phases. The
 % exponential tail is what makes the convolution reach up to SCV 1; the
 % concentrated part is what makes it reach far below the Erlang bound 1/order at
 % the same order.
 %
 % The order is the smallest tabulated one that reaches the target, capped by
 % maxPhases when given: with a phase budget an Erlang can only reach
 % 1/maxPhases, while this construction reaches O(1/maxPhases^2), and the
 % residual SCV is then the closest achievable from below.
 %
 % SCV >= 1 is outside the range of a concentrated ME (its SCV never exceeds
 % 0.34), and the caller keeps its own hyperexponential fit there.
 %
 % Mirrors native Python fit_me_mean_scv and jline.lang.processes.MEFit.
 %
 % @par Syntax:
 % @code
 % me = dist_fit_me(mean, scv)
 % me = dist_fit_me(mean, scv, maxPhases)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>mean<td>Target mean, positive
 % <tr><td>scv<td>Target squared coefficient of variation, in (0,1)
 % <tr><td>maxPhases<td>Optional cap on the number of phases, 0 for no cap
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>me<td>ME distribution with the requested mean and, budget permitting, SCV
 % </table>
%}
function me = dist_fit_me(mean, scv, maxPhases)

if nargin < 3 || isempty(maxPhases)
    maxPhases = 0;
end

if ~isscalar(mean) || ~isfinite(mean) || mean <= 0
    line_error(mfilename, 'dist_fit_me mean must be a positive finite number.');
end
if ~isscalar(scv) || ~isfinite(scv) || scv <= 0 || scv >= 1
    line_error(mfilename, 'dist_fit_me requires 0 < scv < 1; use a hyperexponential for scv >= 1 and a CME for scv = 0.');
end

% Smallest tabulated order whose convolution range covers the target, subject to
% the phase budget. reach(order) = sY/(1+sY) is the minimum SCV of the
% convolution, which sits just below the sY of the CME alone.
orders = CME.getSupportedOrders();
bestOrder = [];
for i = 1:numel(orders)
    order = orders(i);
    if maxPhases > 0 && order + 1 > maxPhases
        continue;
    end
    sY = CME.getMinSCV(order);
    bestOrder = order; % budget-limited: keep the most concentrated one that fits
    if sY/(1+sY) <= scv
        break;
    end
end
if isempty(bestOrder)
    line_error(mfilename, sprintf('No CME order fits a budget of %d phases; the smallest is 3 phases plus one exponential.', maxPhases));
end

[alphaY, AY, sY] = CME.representation(bestOrder);
reach = sY/(1+sY);
if scv < reach
    % Budget-limited: the target is below what this order can reach, so the most
    % concentrated member of the family is returned and the caller gets the
    % closest achievable SCV rather than a silent Erlang truncation.
    c = mean/(1+sY);
else
    c = mean*(1 - sqrt(1 - (1+sY)*(1-scv)))/(1+sY);
end
d = mean - c;

n = numel(alphaY);
if d <= mean*1e-12
    me = CME(mean, bestOrder);
    return;
end
if c <= mean*1e-12
    me = ME(1.0, -1.0/mean);
    return;
end

% Convolution of two matrix exponentials: the exit flow of the first block feeds
% the entry of the second, exactly as for a phase-type.
alpha = zeros(1, n+1);
alpha(1:n) = alphaY;
A = zeros(n+1, n+1);
A(1:n, 1:n) = AY/c;
A(1:n, n+1) = -(AY/c)*ones(n,1);
A(n+1, n+1) = -1/d;
me = ME(alpha, A, false);
end
