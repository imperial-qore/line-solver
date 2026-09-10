function [SensTable, sens] = getSensitivityTable(self, varargin)
% GETSENSITIVITYTABLE Performance sensitivities with respect to service rates.
%
% [SENSTABLE, SENS] = GETSENSITIVITYTABLE(SELF) returns a table with one row
% per (Station, JobClass) giving the derivative of that row's mean performance
% measures with respect to that station-class service RATE:
%   dTput_dRate, dRespT_dRate, dQLen_dRate, dUtil_dRate.
%
% Two branches produce the derivatives, selected automatically:
%
%   'exact'  Analytic differentiation of a product-form recursion, exact to
%            machine precision and cheaper than a single extra solve. Closed
%            networks use pfqn_sens (differentiated MVA), open networks the
%            closed-form BCMP sensitivities (their stations decouple). Rate
%            derivatives follow the chain rule d(.)/d(rate) = -(L/rate) d(.)/dL,
%            since L(i,r) = visits(i,r)/rate(i,r). Available only on the solvers
%            that evaluate that recursion, SolverMVA and SolverNC, and only for
%            single-server queues plus an optional delay, with mixed
%            (open+closed) models excluded.
%
%   'fd'     Forward or central finite differences on the CALLING solver's own
%            predictions: the service process at (station,class) is rate-scaled
%            by (1+h), the same solver with the same options is re-run, and the
%            difference quotient is formed. This costs 1+M*R solves (forward) or
%            2*M*R (central), and it is the only branch that applies to
%            non-product-form models, so it is what every solver other than
%            SolverMVA and SolverNC uses.
%
% [...] = GETSENSITIVITYTABLE(SELF, 'method', M, 'step', H, 'scheme', S) sets
%   M in {'auto','exact','fd'}   default 'auto': exact where available and in
%                                scope, finite differences otherwise
%   H  relative step of the rate perturbation, default 1e-4 for deterministic
%      solvers and 1e-2 for the simulators, whose Monte Carlo error would
%      otherwise dominate the difference quotient
%   S in {'forward','central'}   default 'forward'
%
% SENS is the pfqn_sens struct on the closed exact branch, and empty on the
% open exact branch and on the finite-difference branch. The branch actually
% taken is reported in SENSTABLE.Properties.UserData.method.
%
% Simulation solvers must be run with common random numbers for the difference
% quotient to be meaningful: the same options, and hence the same seed, are
% reused for the base and the perturbed runs. An unset seed is pinned before
% the sweep so that the runs remain paired.

options = struct('method', 'auto', 'step', [], 'scheme', 'forward');
if mod(numel(varargin), 2) ~= 0
    line_error(mfilename, 'Options must be given as name-value pairs.');
end
for a = 1:2:numel(varargin)
    name = lower(varargin{a});
    if ~isfield(options, name)
        line_error(mfilename, sprintf('Unknown option ''%s''.', varargin{a}));
    end
    options.(name) = varargin{a+1};
end
if ischar(options.method), options.method = lower(options.method); end
if ischar(options.scheme), options.scheme = lower(options.scheme); end
if ~ismember(options.method, {'auto', 'exact', 'fd'})
    line_error(mfilename, 'The method must be one of ''auto'', ''exact'', ''fd''.');
end
if ~ismember(options.scheme, {'forward', 'central'})
    line_error(mfilename, 'The scheme must be ''forward'' or ''central''.');
end

sn = self.model.getStruct();
R = sn.nclasses;
queueIndices = find(sn.nodetype == NodeType.Queue);
Mq = numel(queueIndices);

exactAvailable = self.supportsExactSensitivity();
[inScope, scopeMsg] = i_exactScope(sn);
switch options.method
    case 'exact'
        if ~exactAvailable
            line_error(mfilename, sprintf(['Exact analytic sensitivities differentiate a ', ...
                'product-form recursion and are available on SolverMVA and SolverNC only; ', ...
                '%s must use ''fd''.'], class(self)));
        end
        if ~inScope
            line_error(mfilename, scopeMsg);
        end
        useExact = true;
    case 'fd'
        useExact = false;
    otherwise
        useExact = exactAvailable && inScope;
end

sens = [];
if useExact
    [dTput_self, dRespT_self, dQLen_self, dUtil_self, mask, sens] = ...
        i_sensitivityExact(sn, queueIndices, R);
    methodUsed = 'exact';
else
    [dTput_self, dRespT_self, dQLen_self, dUtil_self, mask] = ...
        i_sensitivityFD(self, sn, queueIndices, R, options);
    methodUsed = 'fd';
end

Station = {}; JobClass = {};
dTput = []; dRespT = []; dQLen = []; dUtil = [];
for ist = 1:Mq
    node = queueIndices(ist);
    for r = 1:R
        if ~mask(ist, r), continue; end
        Station{end+1, 1} = sn.nodenames{node}; %#ok<AGROW>
        JobClass{end+1, 1} = sn.classnames{r};   %#ok<AGROW>
        dTput(end+1, 1)  = dTput_self(ist, r);   %#ok<AGROW>
        dRespT(end+1, 1) = dRespT_self(ist, r);  %#ok<AGROW>
        dQLen(end+1, 1)  = dQLen_self(ist, r);   %#ok<AGROW>
        dUtil(end+1, 1)  = dUtil_self(ist, r);   %#ok<AGROW>
    end
end

SensTable = table(Station, JobClass, dTput, dRespT, dQLen, dUtil, ...
    'VariableNames', {'Station', 'JobClass', 'dTput_dRate', ...
    'dRespT_dRate', 'dQLen_dRate', 'dUtil_dRate'});
SensTable.Properties.UserData = struct('method', methodUsed);
end

function [inScope, msg] = i_exactScope(sn)
% Scope of the analytic branch: single-server stations (a delay, with infinite
% servers, is allowed) and not a mixed open-and-closed model. Class switching is
% supported: the branch aggregates classes into chains before differentiating.
inScope = false;
msg = '';
try
    [~, ~, ~, ~, ~, S, ~] = sn_get_product_form_params(sn);
catch
    msg = 'getSensitivityTable could not extract the product-form parameters of this model.';
    return;
end
if any(S(isfinite(S)) > 1)
    msg = 'getSensitivityTable supports single-server stations only.';
    return;
end
N = sn.njobs;
if any(isinf(N)) && any(isfinite(N))
    msg = 'getSensitivityTable does not yet support mixed (open+closed) networks.';
    return;
end
inScope = true;
end

function [dTput_self, dRespT_self, dQLen_self, dUtil_self, mask, sens] = ...
    i_sensitivityExact(sn, queueIndices, R)
% Analytic branch: differentiated MVA (closed) or closed-form BCMP (open), both
% evaluated at CHAIN level and then disaggregated back to the classes.
%
% A product-form model is solved per chain, not per class: a chain carries the
% population and the arrival rate, and a class is a share of it. With class
% switching the two differ, the whole chain population sitting on the reference
% class, so the recursion must see the chain demands
%   Dc(i,c) = sum_{r in c} D(i,r),  Nc(c) = sum_{r in c} N(r),  Zc likewise.
% The parameter of the table is still the per-class rate mu(i,r), which enters
% exactly one chain demand, giving the chain rule
%   d(.)/dmu(i,r) = -(D(i,r)/mu(i,r)) d(.)/dDc(i,c).
% The per-class measures are then composed from the chain ones, which is also
% what fixes the visit ratios: a class throughput at a station is X_c*v(i,r) and
% a per-visit response time is Q(i,r)/T(i,r), whereas the chain quantities are
% per chain and per visit-chain respectively.
[lambda, D, Np, Z, ~, ~, ~] = sn_get_product_form_params(sn);
N = sn.njobs;
isOpen = any(isinf(N));
Mq = numel(queueIndices);
C = sn.nchains;

sens = [];
dRespT_self = zeros(Mq, R);
dQLen_self = zeros(Mq, R);
dUtil_self = zeros(Mq, R);
dTput_self = zeros(Mq, R);
mask = false(Mq, R);
rates = zeros(Mq, R);
for ist = 1:Mq
    sIdx = sn.nodeToStation(queueIndices(ist));
    for r = 1:R
        rates(ist, r) = sn.rates(sIdx, r);
        mask(ist, r) = isfinite(rates(ist, r)) && rates(ist, r) > 0 && D(ist, r) > 0;
    end
end

% chainOf(r) is the chain class r belongs to; Dc, Zc, Nc and lambdac aggregate
% the per-class quantities over the classes of each chain.
chainOf = zeros(1, R);
Dc = zeros(Mq, C);
Zc = zeros(1, C);
Nc = zeros(1, C);
lambdac = zeros(1, C);
Zrow = sum(Z, 1);
for c = 1:C
    cls = find(sn.chains(c, :));
    chainOf(cls) = c;
    Dc(:, c) = sum(D(:, cls), 2);
    Zc(c) = sum(Zrow(cls));
    Nc(c) = sum(Np(cls(isfinite(Np(cls)))));
    lambdac(c) = sum(lambda(cls));
end

if ~isOpen
    % ---- closed branch: differentiated MVA at chain level -------------
    sens = pfqn_sens(Dc, Nc, Zc);
    for ist = 1:Mq
        for r = 1:R
            if ~mask(ist, r), continue; end
            c = chainOf(r);
            if Dc(ist, c) <= 0, continue; end
            rate = rates(ist, r);
            Dir = D(ist, r);
            visits = Dir * rate;              % chain-normalized visit ratio
            p = (ist-1)*C + c;
            chain = -Dir / rate;              % dDc(i,c)/dmu(i,r)
            Xc = sens.X(c);
            Qc = sens.Q(ist, c);
            dXc = sens.dX(c, p) * chain;
            dQc = sens.dQ(ist, c, p) * chain;
            % Class share of the chain queue at this station, and its own
            % dependence on the rate.
            alpha = Dir / Dc(ist, c);
            dalpha = chain * (Dc(ist, c) - Dir) / Dc(ist, c)^2;
            Qir = alpha * Qc;
            dQir = dalpha * Qc + alpha * dQc;
            Tir = Xc * visits;
            dTir = dXc * visits;
            dTput_self(ist, r) = dTir;
            dQLen_self(ist, r) = dQir;
            dUtil_self(ist, r) = dXc * Dir + Xc * chain;
            if Tir > 0
                % Per-visit response time by Little's law, R = Q/T.
                dRespT_self(ist, r) = (dQir * Tir - Qir * dTir) / Tir^2;
            end
        end
    end
else
    % ---- open branch: exact closed-form BCMP --------------------------
    % rho(i,r) = lambdac(c)*D(i,r) with c the chain of r; U(i) = sum_r rho(i,r).
    % The stations decouple, so only the own service rate mu(i,r) moves the
    % measures at (i,r). The throughput T(i,r) = lambdac(c)*v(i,r) is fixed by
    % the arrival rate, hence its rate derivative is exactly zero.
    rho = zeros(Mq, R);
    for ist = 1:Mq
        for r = 1:R
            if D(ist, r) > 0
                rho(ist, r) = lambdac(chainOf(r)) * D(ist, r);
            end
        end
    end
    Ui = sum(rho, 2);
    for ist = 1:Mq
        denom = 1 - Ui(ist);
        for r = 1:R
            if ~mask(ist, r), continue; end
            rate = rates(ist, r);
            svct = 1 / rate;                  % per-visit service time
            drho = -rho(ist, r) / rate;
            dU = drho;                        % own class only
            dsvct = -svct / rate;
            dRespT_self(ist, r) = (dsvct * denom + svct * dU) / denom^2;
            dQLen_self(ist, r)  = (drho * denom + rho(ist, r) * dU) / denom^2;
            dUtil_self(ist, r)  = drho;
            dTput_self(ist, r)  = 0;   % open throughput = lambda*visits (fixed)
        end
    end
end
end

function [dTput_self, dRespT_self, dQLen_self, dUtil_self, mask] = ...
    i_sensitivityFD(self, sn, queueIndices, R, options)
% Finite-difference branch: re-run the calling solver on rate-perturbed copies
% of the model. The perturbation is a pure time scaling of the service process,
% so the shape of the distribution, and in particular its SCV, is preserved and
% only the rate moves.
Mq = numel(queueIndices);
central = strcmp(options.scheme, 'central');

h = options.step;
if isempty(h)
    if i_isSimulation(self)
        h = 1e-2;
    else
        h = 1e-4;
    end
end
if ~isnumeric(h) || ~isscalar(h) || ~isfinite(h) || h <= 0 || h >= 1
    line_error(mfilename, 'The finite-difference step must be a scalar in (0,1).');
end
if i_isSimulation(self)
    % Common random numbers: the base and perturbed runs must share a seed,
    % otherwise the difference quotient measures Monte Carlo noise.
    if ~isfield(self.options, 'seed') || isempty(self.options.seed) || ...
            ~isfinite(self.options.seed)
        self.options.seed = 23000;
    end
end

mask = false(Mq, R);
visited = i_visitMask(sn);
for ist = 1:Mq
    sIdx = sn.nodeToStation(queueIndices(ist));
    for r = 1:R
        rate = sn.rates(sIdx, r);
        mask(ist, r) = isfinite(rate) && rate > 0 && visited(sIdx, r);
    end
end

dRespT_self = zeros(Mq, R);
dQLen_self = zeros(Mq, R);
dUtil_self = zeros(Mq, R);
dTput_self = zeros(Mq, R);

[Q0, U0, R0, T0] = i_solveOnce(self);

for ist = 1:Mq
    sIdx = sn.nodeToStation(queueIndices(ist));
    station = self.model.stations{sIdx};
    for r = 1:R
        if ~mask(ist, r), continue; end
        jobclass = self.model.classes{r};
        rate = sn.rates(sIdx, r);
        base = station.getService(jobclass);
        restore = onCleanup(@() i_restoreService(self.model, station, jobclass, base));

        i_setService(self.model, station, jobclass, dist_scale_rate(base, 1 + h));
        [Qp, Up, Rp, Tp] = i_solveOnce(self);

        if central
            i_setService(self.model, station, jobclass, dist_scale_rate(base, 1 - h));
            [Qm, Um, Rm, Tm] = i_solveOnce(self);
            denom = 2 * rate * h;
        else
            Qm = Q0; Um = U0; Rm = R0; Tm = T0;
            denom = rate * h;
        end

        dTput_self(ist, r)  = (Tp(sIdx, r) - Tm(sIdx, r)) / denom;
        dRespT_self(ist, r) = (Rp(sIdx, r) - Rm(sIdx, r)) / denom;
        dQLen_self(ist, r)  = (Qp(sIdx, r) - Qm(sIdx, r)) / denom;
        dUtil_self(ist, r)  = (Up(sIdx, r) - Um(sIdx, r)) / denom;

        clear restore;   % restores the unperturbed service process
    end
end
end

function [QN, UN, RN, TN] = i_solveOnce(self)
% Solve with a fresh instance of the calling solver class, carrying its options
% over so that seeds, tolerances and method selection are those of the caller.
% A fresh instance is used because a solver caches its results, and reset()
% alone would not discard a warm start.
solver = feval(class(self), self.model, self.options);
[QN, UN, RN, TN] = solver.getAvg();
end

function i_setService(model, station, jobclass, distrib)
station.setService(jobclass, distrib);
model.refreshProcesses();
end

function i_restoreService(model, station, jobclass, distrib)
station.setService(jobclass, distrib);
model.refreshProcesses();
end

function tf = i_isSimulation(self)
tf = ismember(class(self), {'SolverSSA', 'SolverJMT', 'SolverLDES'});
end

function visited = i_visitMask(sn)
% A (nstations x nclasses) mask of the station-class pairs that carry visits,
% used in place of the demand matrix, which a non-product-form model may not
% admit.
visited = false(sn.nstations, sn.nclasses);
for c = 1:sn.nchains
    V = sn.visits{c};
    if isempty(V), continue; end
    visited = visited | (V > GlobalConstants.Zero);
end
end
