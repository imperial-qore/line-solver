function w = npfqn_rqna_weight(t)
% w = NPFQN_RQNA_WEIGHT(t) - Canonical RBM correlation weight function w*(t)
% used by the Robust Queueing Network Analyzer (RQNA).
%
% w*(t) = 1 - (1 - c*(t))/(2 t), where c*(t) is the correlation function of the
% stationary version of canonical reflected Brownian motion (drift -1, diffusion
% coefficient 1),
%
%   c*(t) = 2(1 - 2t - t^2) Phi^c(sqrt(t)) + 2 sqrt(t) phi(sqrt(t)) (1 + t),
%
% with Phi^c the standard-normal complementary cdf and phi its density. The
% weight is monotonically increasing with w*(0)=0 and w*(Inf)=1.
% Reference: W. Whitt and W. You (2018), "A Robust Queueing Network Analyzer
% Based on Indices of Dispersion", eqs. (24)-(25).
%
% Input:  t - array of nonnegative time arguments
% Output: w - array of weights w*(t), same shape as t

w = zeros(size(t));
for idx = 1:numel(t)
    ti = t(idx);
    if ti <= 0
        w(idx) = 0;
        continue;
    end
    if ~isfinite(ti)
        w(idx) = 1;
        continue;
    end
    st = sqrt(ti);
    Phic = 0.5*erfc(st/sqrt(2));            % complementary standard-normal cdf
    phi  = exp(-ti/2)/sqrt(2*pi);           % standard-normal density at sqrt(t)
    cstar = 2*(1 - 2*ti - ti^2)*Phic + 2*st*phi*(1 + ti);
    if ti < 1e-6
        % limit w*(t) -> 0 as t -> 0; use series to avoid catastrophic cancellation
        w(idx) = 0;
    else
        w(idx) = 1 - (1 - cstar)/(2*ti);
    end
    % numerical guard: w* is a weight in [0,1]
    if w(idx) < 0, w(idx) = 0; end
    if w(idx) > 1, w(idx) = 1; end
end
end
