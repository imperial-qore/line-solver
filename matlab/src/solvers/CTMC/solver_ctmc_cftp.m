%{
%{
 % @file solver_ctmc_cftp.m
 % @brief Stationary analysis of a closed single-class product-form network by
 %        perfect sampling (Coupling From The Past) instead of state-space
 %        enumeration.
%}
%}

%{
%{
 % @brief Draws iid states from the exact stationary distribution with
 %        pfqn_cftp and reduces them to the standard mean performance metrics.
 % @fn solver_ctmc_cftp(sn, options)
 % @param sn Network structure.
 % @param options Solver options (method, samples, seed).
 % @return QN Mean queue length per station and class.
 % @return UN Utilization per station and class.
 % @return RN Response time per station and class.
 % @return TN Throughput per station and class.
 % @return CN System response time per class.
 % @return XN System throughput per class.
 % @return Xs Sampled states, one per row (samples x stations).
 % @return Ts Per-sample coalescence horizon or mixing steps.
 % @return pAggr Empirical probability of each distinct sampled state.
 % @return SSq Distinct sampled states, aligned with pAggr.
%}
%}
function [QN,UN,RN,TN,CN,XN,Xs,Ts,pAggr,SSq] = solver_ctmc_cftp(sn, options)
% [QN,UN,RN,TN,CN,XN,XS,TS,PAGGR,SSQ] = SOLVER_CTMC_CFTP(SN, OPTIONS)
%
% Perfect-sampling steady-state analysis of a closed single-class
% product-form network. States are drawn iid from the exact stationary
% distribution by monotone Coupling From The Past, so the estimator carries
% Monte Carlo error O(samples^(-1/2)) but never enumerates the state space.
%
% Reference: S. Kijima and T. Matsui, "Approximate/Perfect Samplers for
% Closed Jackson Networks", Proc. Winter Simulation Conference, 2005.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

M = sn.nstations;
K = sn.nclasses;

sampler = solver_ctmc_cftp_sampler(options.method);
solver_ctmc_cftp_assert(sn, options);

[Lchain,STchain,Vchain,~,Nchain,~,refstatchain] = sn_get_demands_chain(sn);
L = Lchain(:,1)';
N = Nchain(1);
S = sn.nservers(:)';
S(sn.sched == SchedStrategy.INF) = Inf;

nsamples = options.samples;
if isempty(nsamples) || ~isfinite(nsamples) || nsamples < 1
    line_error(mfilename,'The cftp method requires a finite positive options.samples.');
end
nsamples = round(nsamples);

[Qs,Xs,Ts] = pfqn_cftp(L, N, S, nsamples, sampler);

QN = zeros(M,K);
UN = zeros(M,K);
RN = zeros(M,K);
TN = zeros(M,K);
CN = zeros(1,K);
XN = zeros(1,K);

busy = zeros(1,M);
for i=1:M
    busy(i) = mean(min(Xs(:,i), S(i)));
end

% The utilization law X = mu_i*E[min(n_i,c_i)]/V_i holds at every station, but
% each station estimates it with its own Monte Carlo error. Taking the estimate
% at the reference station and propagating it through the visit ratios matches
% the CTMC convention (XN is the arrival rate at the reference station) and
% keeps flow balance, Little's law and C = N/X exact in the reported table.
iref = refstatchain(1);
if STchain(iref,1) > 0 && Vchain(iref,1) > 0
    XN(1) = busy(iref)/STchain(iref,1)/Vchain(iref,1);
end
if XN(1) > 0
    CN(1) = N/XN(1);
end

for i=1:M
    QN(i,1) = Qs(i);
    TN(i,1) = Vchain(i,1)*XN(1);
    % Utilization keeps its own estimator E[min(n_i,c_i)]/c_i: it is unbiased
    % and confined to [0,1] by construction, whereas deriving it from the
    % reference-station throughput lets Monte Carlo error push a saturated
    % station above 1.
    if sn.sched(i) == SchedStrategy.INF
        UN(i,1) = QN(i,1);
    else
        UN(i,1) = busy(i)/S(i);
    end
    if TN(i,1) > 0
        RN(i,1) = QN(i,1)/TN(i,1);
    end
end

[SSq,~,idx] = unique(Xs,'rows');
pAggr = accumarray(idx,1)/nsamples;
end

% ---- method string -> pfqn_cftp sampler ----------------------------------
function sampler = solver_ctmc_cftp_sampler(method)
switch lower(method)
    case {'cftp','cftp.exact'}
        sampler = 'cftp';
    case 'cftp.approx'
        sampler = 'approx';
    otherwise
        line_error(mfilename,'Unknown cftp variant ''%s''. Use ''cftp'' or ''cftp.approx''.',method);
end
end

% ---- model class gate ----------------------------------------------------
function solver_ctmc_cftp_assert(sn, options)
% The rules live in SOLVER_CTMC_CFTP_SUPPORTS, which SolverCTMC.supportsModelMethod
% also asks: the analyzer must refuse exactly what the report refuses, and one
% predicate with two callers is what keeps the two from drifting apart.
[bool, reason] = solver_ctmc_cftp_supports(sn, options);
if ~bool
    line_error(mfilename, reason);
end
end
