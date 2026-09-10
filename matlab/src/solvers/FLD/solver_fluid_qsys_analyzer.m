function [Q,U,R,T,C,X,Qt,Ut,Tt,runtime,actualmethod] = solver_fluid_qsys_analyzer(sn, options)
% [Q,U,R,T,C,X,QT,UT,TT,RUNTIME,ACTUALMETHOD] = SOLVER_FLUID_QSYS_ANALYZER(SN, OPTIONS)
%
% The single-station fluid limits: a Source -> Queue -> Sink model with one
% class, answered by a closed-form fluid or Gaussian limit rather than by
% integrating the network drift.
%
% WHY THESE ARE FLUID METHODS AND NOT MVA ONES. Each depends on the service or
% patience law BEYOND ITS MEAN -- the stationary point of the Liu-Whitt model
% is where the patience ccdf crosses 1/rho, the Mt/G/Inf mean is a convolution
% with the service ccdf -- and each is the limit of a sequence of systems, not
% an approximation to a fixed one. That is the fluid solver's contract.
%
% METHODS
%   'ggisgi.fluid' - stationary point of the G/GI/s+GI fluid model
%                    (Liu and Whitt, Operations Research 60(5), 2012)
%   'ggingi.tga'   - truncated Gaussian approximation, which adds the
%                    O(sqrt(n)) fluctuation around that point
%                    (Liu, Whitt and Yu, Naval Research Logistics 63(3), 2016)
%   'tvms'         - the Gt/Mt/st+GI many-server fluid queue at CONSTANT
%                    staffing (Liu and Whitt, INFORMS J. Computing 26(1), 2014)
%   'mtginf'       - the exact Mt/G/Inf mean (Eick, Massey and Whitt,
%                    Management Science 39(2), 1993)
%   'mol'          - the modified-offered-load approximation for the finite
%                    server count (Massey and Whitt, Ann. Appl. Prob. 4(4), 1994)
%
% The transient tables Qt, Ut, Tt are returned for the three time-varying
% methods and are empty for the two stationary ones.
%
% See also SOLVER_FLUID_ANALYZER, SN_ARRIVAL_RATE_FUN, SN_PATIENCE_HANDLES.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

T0 = tic;
M = sn.nstations; K = sn.nclasses;
Q = zeros(M,K); U = zeros(M,K); R = zeros(M,K); T = zeros(M,K);
C = zeros(1,K); X = zeros(1,K);
Qt = {}; Ut = {}; Tt = {};
method = options.method;
actualmethod = method;

% THE SHAPE THESE LIMITS ARE STATED FOR, refused by name rather than answered
% on a model they do not describe: one open class through one queueing station,
% a patience law for the abandonment limits, a finite server count for three of
% them. The MVA qsys analyzer is reached by a structural dispatch that
% guarantees it; these methods are selected by NAME, so the check lives in
% FLUID_QSYS_ADMITS, which SolverFLD.supportsModelMethod asks as well.
[shapeOk, shapeWhy] = fluid_qsys_admits(sn, method);
if ~shapeOk
    line_error(mfilename, shapeWhy);
end
source_ist = sn.nodeToStation(sn.nodetype == NodeType.Source);
queue_ist = sn.nodeToStation(sn.nodetype == NodeType.Queue);
if isempty(queue_ist)
    queue_ist = sn.nodeToStation(sn.nodetype == NodeType.Delay);
end
Vq = sn.visits{1}(sn.stationToStateful(queue_ist));
lambda = sn.rates(source_ist)*Vq;
mu = sn.rates(queue_ist);
s = sn.nservers(queue_ist);
ca = sqrt(sn.scv(source_ist));
cs = sqrt(sn.scv(queue_ist));
hpat = sn_patience_handles(sn, queue_ist, 1);

% The service ccdf, needed by the two Mt/G methods: they are exact in the
% service DISTRIBUTION, not in its mean, which is the whole point of the
% Eick-Massey-Whitt lag.
svcMAP = [];
if iscell(sn.proc{queue_ist}) && numel(sn.proc{queue_ist}) >= 1
    pr = sn.proc{queue_ist}{1};
    if iscell(pr) && numel(pr) >= 2 && all(size(pr{1}) == size(pr{2}))
        svcMAP = {pr{1}, pr{2}};
    end
end
if isempty(svcMAP)
    serviceCcdf = @(x) exp(-mu*x);
else
    serviceCcdf = @(x) 1 - map_cdf(svcMAP, x);
end
ES = 1/mu;
ES2 = (1 + sn.scv(queue_ist)) * ES^2;

switch method
    case {'ggisgi.fluid','fluid.ggisgi','ggisgi'}
        res = qsys_ggisgi_fluid(lambda, mu, s, hpat.ccdf);
        [Q,U,R,T,C,X] = local_stationary(Q,U,R,T,C,X, queue_ist, source_ist, Vq, ...
            lambda, res.meanNumber, res.throughput, res.utilization);
        actualmethod = 'ggisgi.fluid';

    case {'ggingi.tga','fluid.tga','tga'}
        res = qsys_ggingi_tga(lambda, mu, s, ca, cs, hpat.ccdf, ...
            'patiencePdf', hpat.pdf, 'serviceCcdf', serviceCcdf);
        Tq = lambda*(1 - res.probAbandon);
        [Q,U,R,T,C,X] = local_stationary(Q,U,R,T,C,X, queue_ist, source_ist, Vq, ...
            lambda, res.meanNumber, Tq, min(res.meanNumberInService/s, 1));
        actualmethod = 'ggingi.tga';

    case {'tvms','fluid.tvms'}
        [lambdaFun, ~] = sn_arrival_rate_fun(sn, source_ist, 1);
        % CONSTANT STAFFING. Nothing in a Network declares a time-varying
        % server count, so s(t) is the station's own s; the time variation the
        % method is for enters through lambda(t) alone. A staffing schedule
        % would need a model feature that does not exist, and inventing one
        % here would make the solver answer a model the user did not build.
        [t0, t1] = local_horizon(options);
        res = qsys_gtmtst_fluid(lambdaFun, @(t) s*ones(size(t)), @(t) mu*ones(size(t)), ...
            hpat.ccdf, t1 - t0, 'patiencePdf', hpat.pdf);
        tt = t0 + res.times;
        Xtot = res.X;
        Btot = res.B;
        Tser = mu * Btot;
        [Q,U,R,T,C,X,Qt,Ut,Tt] = local_transient(Q,U,R,T,C,X, queue_ist, source_ist, Vq, M, K, ...
            tt, Xtot, res.utilization, Tser, res.arrivalRate);
        actualmethod = 'tvms';

    case {'mtginf','fluid.mtginf'}
        [lambdaFun, ~] = sn_arrival_rate_fun(sn, source_ist, 1);
        [t0, t1] = local_horizon(options);
        tt = linspace(t0, t1, 200);
        res = qsys_mtginf(lambdaFun, serviceCcdf, ES, tt, 'ES2', ES2);
        % An infinite-server station serves everything that arrives, so the
        % throughput is the arrival rate and the utilization is not a fraction
        % of anything; report the busy servers, as the fluid analyzer does for
        % a Delay.
        [Q,U,R,T,C,X,Qt,Ut,Tt] = local_transient(Q,U,R,T,C,X, queue_ist, source_ist, Vq, M, K, ...
            res.times, res.meanNumber, res.meanNumber, res.arrivalRate, res.arrivalRate);
        actualmethod = 'mtginf';

    case {'mol','fluid.mol'}
        [lambdaFun, ~] = sn_arrival_rate_fun(sn, source_ist, 1);
        [t0, t1] = local_horizon(options);
        tt = linspace(t0, t1, 200);
        % A finite buffer beyond the servers is not part of the loss model the
        % approximation is for; only s servers and no waiting room is.
        useDelay = isfinite(sn.cap(queue_ist)) && sn.cap(queue_ist) > s;
        res = qsys_mtgs0_mol(lambdaFun, serviceCcdf, ES, s, tt, 'delay', useDelay, 'ES2', ES2);
        busy = res.meanBusyMOL;
        Tser = mu * busy;
        [Q,U,R,T,C,X,Qt,Ut,Tt] = local_transient(Q,U,R,T,C,X, queue_ist, source_ist, Vq, M, K, ...
            res.times, busy, busy/s, Tser, res.arrivalRate);
        actualmethod = 'mol';

    otherwise
        line_error(mfilename,sprintf('The ''%s'' method is not a single-station fluid limit.',method));
end

runtime = toc(T0);
end

function [t0, t1] = local_horizon(options)
% The integration window. The rule lives in FLUID_QSYS_HORIZON, which
% SolverFLD.supportsModelMethod also asks: a report that offers 'mol' or
% 'mtginf' on a solver with no finite timespan offers a run that stops here,
% and one predicate with two callers is what keeps the two answers the same.
[bool, reason, t0, t1] = fluid_qsys_horizon(options);
if ~bool
    line_error(mfilename, reason);
end
end

function [Q,U,R,T,C,X] = local_stationary(Q,U,R,T,C,X, qi, si, Vq, lambda, Lsys, Tq, Uq)
% Little's law on the CARRIED rate, as every LINE solver reports a station that
% loses work.
if Tq > 0
    R(qi,1) = Lsys/Tq;
else
    R(qi,1) = 0;
end
Q(qi,1) = Lsys;
U(qi,1) = Uq;
T(qi,1) = Tq;
T(si,1) = lambda/Vq;
X(1) = Tq;
C(1) = R(qi,1)*Vq;
end

function [Q,U,R,T,C,X,Qt,Ut,Tt] = local_transient(Q,U,R,T,C,X, qi, si, Vq, M, K, ...
        tt, Lt, Ut_, Tt_, arrival)
% The steady-state row of a time-varying model is the TIME AVERAGE over the
% horizon, which is what a stationary reader of a periodic system measures; the
% trajectory itself is returned beside it.
tt = tt(:); Lt = Lt(:); Ut_ = Ut_(:); Tt_ = Tt_(:); arrival = arrival(:);
span = tt(end) - tt(1);
if span <= 0
    Lbar = Lt(1); Ubar = Ut_(1); Tbar = Tt_(1); Abar = arrival(1);
else
    Lbar = trapz(tt, Lt)/span;
    Ubar = trapz(tt, Ut_)/span;
    Tbar = trapz(tt, Tt_)/span;
    Abar = trapz(tt, arrival)/span;
end
Q(qi,1) = Lbar;
U(qi,1) = Ubar;
T(qi,1) = Tbar;
T(si,1) = Abar/Vq;
if Tbar > 0
    R(qi,1) = Lbar/Tbar;
end
X(1) = Tbar;
C(1) = R(qi,1)*Vq;
Qt = cell(M,K); Ut = cell(M,K); Tt = cell(M,K);
for i = 1:M
    for r = 1:K
        Qt{i,r} = [zeros(numel(tt),1), tt];
        Ut{i,r} = [zeros(numel(tt),1), tt];
        Tt{i,r} = [zeros(numel(tt),1), tt];
    end
end
Qt{qi,1} = [Lt, tt];
Ut{qi,1} = [Ut_, tt];
Tt{qi,1} = [Tt_, tt];
Tt{si,1} = [arrival/Vq, tt];
end
