function [Z, W, Q, X] = qsys_gig1_rq(rho, mu, cs2, IaFun)
% [Z,W,Q,X] = QSYS_GIG1_RQ(rho,mu,cs2,IaFun) - Robust Queueing (RQ)
% approximation for a single G/GI/1 queue partially characterized by its
% arrival rate, index of dispersion for counts (IDC) and the first two moments
% of the service time. Implements the mean steady-state workload
%
%   Z* = sup_{x>=0} { -(1-rho) x + sqrt( 2 rho x (I_a(x) + c2_s) / mu ) }
%
% and the derived steady-state performance measures. Reference: W. Whitt and
% W. You (2018), "A Robust Queueing Network Analyzer Based on Indices of
% Dispersion", eqs. (13),(16)-(18).
%
% Inputs:
%   rho   : traffic intensity lambda/mu (0<rho<1)
%   mu    : service rate
%   cs2   : service SCV c2_s
%   IaFun : handle, IaFun(x) -> arrival IDC I_a(x) at time argument x>0
% Outputs:
%   Z : mean steady-state workload E[Z]
%   W : mean steady-state waiting time E[W]
%   Q : mean steady-state queue length E[Q] (number waiting + in service)
%   X : mean number in system (= Q here for single class), kept for interface

if rho <= 0
    Z = 0; W = 0; Q = 0; X = 0; return;
end
if rho >= 1
    Z = Inf; W = Inf; Q = Inf; X = Inf; return;
end
lambda = rho * mu;

    function val = negf(x)
        if x <= 0
            val = 0; return;
        end
        ia = IaFun(x);
        val = -( -(1-rho)*x + sqrt(max(0, 2*rho*x*(ia + cs2)/mu)) );
    end

% The objective is a unimodal-in-practice but occasionally multimodal function
% of x. Optimize with a coarse log-spaced scan for bracketing followed by a
% local golden-section refinement (fminbnd) on the best bracket.
xs = logspace(-6, 8, 200);
fv = arrayfun(@(x) -negf(x), xs);
[~, imax] = max(fv);
lo = xs(max(1,imax-1));
hi = xs(min(numel(xs),imax+1));
opt = optimset('TolX',1e-10);
[xopt, nfval] = fminbnd(@negf, lo, hi, opt);
Zscan = fv(imax);
Z = max(Zscan, -nfval);
Z = max(Z, 0);

% derived measures (eqs. 16-18)
W = max(0, Z/rho - (cs2 + 1)/(2*mu));
Q = lambda * W;      % E[Q] waiting (Little's law on the waiting time)
X = Q + rho;         % E[X] number in system including the one in service
end
