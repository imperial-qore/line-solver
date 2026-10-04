function [meanAoI, lstAoI, peakAoI] = aoi_lcfsd_gim1(Y_lst, mu, E_Y, E_Y2)
%AOI_LCFSD_GIM1 Mean AoI for GI/M/1 non-preemptive LCFS-D queue
%
% [meanAoI, lstAoI, peakAoI] = aoi_lcfsd_gim1(Y_lst, mu, E_Y, E_Y2)
%
% Computes the Age of Information metrics for a GI/M/1 queue with
% non-preemptive Last-Come First-Served with Discarding (LCFS-D) discipline.
%
% In LCFS-D, when a new update arrives while the server is busy:
%   - If there's an update waiting in queue, it is discarded
%   - The new update takes its place in the queue
%   - When service completes, the waiting update (if any) is served
%
% Parameters:
%   Y_lst (function_handle): LST of interarrival time, @(s) -> complex
%   mu (double): Service rate (exponential service)
%   E_Y (double): Mean interarrival time (first moment)
%   E_Y2 (double): Second moment of interarrival time
%
% Returns:
%   meanAoI (double): Mean (average) Age of Information
%   lstAoI (function_handle): LST of AoI distribution
%   peakAoI (double): Mean Peak Age of Information
%
% Formulas (Inoue et al., IEEE Trans. IT, 2019, Section 3.3, NP-LCFS (C)),
% with G*(s) the interarrival LST and rho = 1/(E_Y*mu):
%   E[A]     = 1/mu + E[G^2]/(2E[G])
%              + rho*(-G*'(mu) + mu*G*(mu)*G*''(mu)/(1 + mu*G*'(mu)))   (eq. 69)
%   E[Apeak] = P0*(E[G] + (1+G*(mu))/mu) + Pw*(E[G]/(1-G*(mu)) + 1/mu),
%              P0 = q/(q+G*(mu)), Pw = 1-P0, q = 1 + mu*G*'(mu)/(1-G*(mu)) (eq. 91)
%   A*(s)    = eq. (65)
% Discarding keeps the system stable for every rho, so rho >= 1 is accepted.
%
% Reference:
%   Y. Inoue, H. Masuyama, T. Takine, T. Tanaka, "A General Formula for
%   the Stationary Distribution of the Age of Information and Its
%   Application to Single-Server Queues," IEEE Trans. Information Theory,
%   vol. 65, no. 12, pp. 8305-8324, 2019.
%
% See also: aoi_lcfss_gim1, aoi_fcfs_gim1, aoi_lcfspr_gim1

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

% Validate inputs
if mu <= 0
    line_error(mfilename, 'Service rate mu must be positive');
end
if E_Y <= 0
    line_error(mfilename, 'Mean interarrival time E_Y must be positive');
end

if E_Y2 < E_Y^2
    line_error(mfilename, 'Second moment E_Y2 must be >= E_Y^2');
end

% Discarding keeps at most one update waiting, so the age process is
% regenerative for every rho (Inoue et al. 2019, Section 3.3): no rho < 1 gate.
lambda = 1 / E_Y;
rho = lambda / mu;

% Y*(mu) = P(no service completion within an interarrival time), its first two
% derivatives at mu by central differences
gM = Y_lst(mu);
hstep = 1e-6 * max(1, mu);
mdY = -(Y_lst(mu + hstep) - Y_lst(mu - hstep)) / (2 * hstep);
% five-point stencil, O(h^4): truncation and roundoff both near 1e-10 at this step
h2 = 5e-3 * max(1, mu);
d2Y = (-Y_lst(mu + 2*h2) + 16 * Y_lst(mu + h2) - 30 * gM + 16 * Y_lst(mu - h2) ...
    - Y_lst(mu - 2*h2)) / (12 * h2^2);

% Mean AoI, Inoue et al. 2019, Corollary 36(ii), eq. (69)
meanAoI = 1 / mu + E_Y2 / (2 * E_Y) + rho * (mdY + mu * gM * d2Y / (1 - mu * mdY));

% Mean peak AoI from the two-state chain "no wait"/"wait" of informative updates
% (eq. 91): an update arriving to an empty system peaks at E[max(H,G)] + E[H];
% one that waited a residual service H<G peaks at E[G]/(1 - Y*(mu)) + E[H].
q = 1 - mu * mdY / (1 - gM);
P0 = q / (q + gM);
Pw = gM / (q + gM);
peakAoI = P0 * (E_Y + (1 + gM) / mu) + Pw * (E_Y / (1 - gM) + 1 / mu);

% LST of the AoI, Theorem 35(ii), eq. (65)
lstAoI = @(s) aoi_lcfsd_gim1_lst_eval(s, Y_lst, mu, E_Y, rho, gM, mdY);

end

function val = aoi_lcfsd_gim1_lst_eval(s, Y_lst, mu, E_Y, rho, gM, mdY)
%AOI_LCFSD_GIM1_LST_EVAL Evaluate eq. (65) at point(s) s
    val = zeros(size(s));
    for i = 1:numel(s)
        si = s(i);
        if abs(si) < 1e-12
            val(i) = 1; % removable singularity of the residual-interarrival LST at s = 0
            continue
        end
        ys = si + mu;
        hs = 1e-6 * max(1, abs(ys));
        mdYs = -(Y_lst(ys + hs) - Y_lst(ys - hs)) / (2 * hs);
        gres = (1 - Y_lst(si)) / (si * E_Y);            % residual interarrival LST
        val(i) = (gres + rho * mu / ys * (Y_lst(ys) ...
            - gM * (1 - mu * mdYs) / (1 - mu * mdY))) * mu / ys;
    end
end
