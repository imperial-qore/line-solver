function c = lqn_setup_charge(self, lqn, tidx)
% LQN_SETUP_CHARGE(SELF, LQN, TIDX) Mean cold start one request of task TIDX pays
%
% A SetupTask powers a thread down when it goes idle and pays a setup before it
% can serve again. The thread is released at a reply and starts a delay-off
% countdown D of mean d; it powers off only if D expires before the next request
% arrives, and a request arriving first cancels the countdown and pays nothing.
% With the idle interval I seen by one thread and exponential D,
%
%   p = P(D < I) = E[I] / (E[I] + d),   and the charge is  p * s.
%
% E[I] comes from the current iterate. Admission takes an ACTIVE idle thread
% before it wakes a sleeping one, so the pool that actually cycles is only as
% large as the load needs: with offered load b = X*S = RHO*MULT threads, about
% max(1,b) stay hot, each seeing arrivals at rate X/max(1,b) and busy S per
% arrival, so
%
%   E[I] = max(1,b)/X - S = (max(1,b) - b) / X.
%
% At MULT = 1 this is (1-RHO)/X and is EXACT given p: with one customer, caller
% mean a and callee mean b it returns p = a/(a+d), the closed form the LDES
% engine is checked against in SolverLDESLayeredSetupTest and
% test_ldes_ln_engine.cpp. Above one thread it is an approximation, the exact
% answer for c servers with setup being matrix-analytic (Gandhi, Harchol-Balter
% and Adan, Performance Evaluation 67(11), 2010).
%
% Shared by method 'srvn.cs' (updateMetricsDefault) and method 'srvn.ph'
% (phComposeEntryLaws/phSetupProb), so the two charge the same thing.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

c = 0;
if ~isfield(lqn,'hassetup') || isempty(lqn.hassetup) || tidx > numel(lqn.hassetup) ...
        || ~full(lqn.hassetup(tidx))
    return
end
s = setupmean(lqn, 'setuptime', tidx);
d = setupmean(lqn, 'delayofftime', tidx);
if ~(s > GlobalConstants.FineTol) || ~(d > GlobalConstants.FineTol)
    return
end
mult = full(lqn.mult(tidx));
if ~isfinite(mult) || mult <= 0
    return % an infinite-server task holds no thread to power down
end
if isempty(self.tput) || numel(self.tput) < tidx || isempty(self.util) || numel(self.util) < tidx
    c = s; % nothing has arrived yet, so the thread is down when the first does
    return
end
X = self.tput(tidx);
if ~isfinite(X) || X <= GlobalConstants.FineTol
    c = s;
    return
end
rho = self.util(tidx);
if ~isfinite(rho) || rho < 0
    rho = 0;
end
rho = min(rho, 1 - GlobalConstants.FineTol);
b = rho * mult;                 % offered load, in threads
EI = (max(1, b) - b) / X;       % idle interval of a thread in the hot pool
c = s * EI / (EI + d);
end

function m = setupmean(lqn, fieldname, tidx)
% Mean of task TIDX's setup or delay-off time, 0 when it declares none.
m = 0;
if ~isfield(lqn, fieldname) || isempty(lqn.(fieldname)) || tidx > numel(lqn.(fieldname))
    return
end
proc = lqn.(fieldname){tidx};
if isempty(proc) || ~isa(proc,'Distribution')
    return
end
m = proc.getMean();
if ~isfinite(m)
    m = 0;
end
end
