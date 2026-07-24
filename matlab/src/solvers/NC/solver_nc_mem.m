function [QN,UN,RN,TN,CN,XN,totiter,actualmethod] = solver_nc_mem(sn, options)
% [QN,UN,RN,TN,CN,XN,TOTITER,ACTUALMETHOD] = SOLVER_NC_MEM(SN, OPTIONS)
%
% Maximum Entropy Method (MEM) for Open and Closed Queueing Networks
%
% Implements the ME algorithms from Kouvatsos (1994) for analyzing
% queueing networks with general arrival and service processes under
% non-priority scheduling disciplines. Open models (Section 3.2) are
% decomposed into GE/GE/1, GE/GE/c and GE/GE/inf building blocks; closed
% models (Section 3.3) are solved by the two-stage pseudo-open network
% plus convolution algorithm on G/G/1 and G/G/inf building blocks; mixed
% models compose the two algorithms by product-form-style conditioning
% (open classes reduce the server capacity seen by the closed classes,
% closed occupancy inflates the open queue lengths).
%
% Parameters:
%   sn - Network structure from getStruct()
%   options - Solver options with MEM-specific fields:
%             - config.mem_tol: convergence tolerance (default 1e-6)
%             - config.mem_maxiter: maximum iterations (default 1000)
%             - config.mem_verbose: print iteration info (default false)
%
% Returns:
%   QN - Mean queue lengths [M x R matrix]
%   UN - Utilizations [M x R matrix]
%   RN - Mean response times [M x R matrix]
%   TN - Throughputs [M x R matrix]
%   CN - System response times per class [1 x R], by Little's law
%   XN - System throughputs per class [1 x R]
%   totiter - Number of iterations until convergence
%   actualmethod - 'mem', or 'mem.blocking' when the model carries a finite
%             station buffer and was solved by the censored GE/GE/c/0;N
%             building blocks of Section 4.1
%
% Reference:
%   D.D. Kouvatsos, "Entropy Maximisation and Queueing Network Models",
%   Annals of Operations Research, 48:63-126, 1994.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

actualmethod = 'mem';

% Validate the model against the MEM feature set
[memok, memreason, memblocking] = solver_nc_mem_supports(sn);
if ~memok
    line_error(mfilename, memreason);
end

M = sn.nstations;
R = sn.nclasses;

% Get MEM options from config
if isfield(options, 'config') && isfield(options.config, 'mem_tol')
    tol = options.config.mem_tol;
else
    tol = 1e-6;
end

if isfield(options, 'config') && isfield(options.config, 'mem_maxiter')
    maxiter = options.config.mem_maxiter;
else
    maxiter = 1000;
end

if isfield(options, 'config') && isfield(options.config, 'mem_verbose')
    verbose = options.config.mem_verbose;
else
    verbose = false;
end

% Create MEM options structure
mem_options = struct();
mem_options.tol = tol;
mem_options.maxiter = maxiter;
mem_options.verbose = verbose;

% Closed models: two-stage pseudo-open + convolution algorithm (Section 3.3)
if sn_is_closed_model(sn)
    N = zeros(1, R);
    refstat = zeros(1, R);
    for r = 1:R
        N(r) = sn.njobs(r);
        refstat(r) = sn.refstat(r);
    end
    mu = zeros(M, R);
    Cs = ones(M, R);
    nservers = ones(M, 1);
    for ist = 1:M
        nservers(ist) = sn.nservers(ist);
        for r = 1:R
            if isfinite(sn.rates(ist, r)) && sn.rates(ist, r) > 0
                mu(ist, r) = sn.rates(ist, r);
                if isfinite(sn.scv(ist, r)) && sn.scv(ist, r) > 0
                    Cs(ist, r) = sn.scv(ist, r);
                end
            end
        end
    end
    % Routing probabilities between stations
    P = zeros(M, M, R);
    for r = 1:R
        for j = 1:M
            jNode = sn.stationToNode(j);
            for k = 1:M
                iNode = sn.stationToNode(k);
                P(j, k, r) = sn.rtnodes((jNode-1)*R + r, (iNode-1)*R + r);
            end
        end
    end
    insens = false(M, 1);
    for ist = 1:M
        insens(ist) = any(sn.sched(ist) == [SchedStrategy.PS, SchedStrategy.LCFSPR]);
    end
    [L, W, ~, ~, lam, rho, X, totiter] = me_cqn(M, R, N, mu, Cs, P, nservers, refstat, insens, mem_options);
    QN = L;
    UN = rho;
    RN = W;
    TN = lam;
    XN = X;
    CN = zeros(1, R);
    for r = 1:R
        if X(r) > 0
            CN(r) = N(r) / X(r); % class cycle time by Little's law
        end
    end
    return
end

% Mixed models: composition of the open (Section 3.2) and closed
% (Section 3.3) algorithms with product-form-style conditioning
if ~sn_is_open_model(sn)
    openCls = false(1, R);
    N = zeros(1, R);
    for r = 1:R
        N(r) = sn.njobs(r);
        openCls(r) = isinf(sn.njobs(r));
    end
    sourceIdx = 0;
    for ind = 1:sn.nnodes
        if sn.nodetype(ind) == NodeType.Source
            sourceIdx = sn.nodeToStation(ind);
            break;
        end
    end
    qs = setdiff(1:M, sourceIdx); % queueing/delay stations
    Mq = length(qs);
    mu = zeros(Mq, R);
    Cs = ones(Mq, R);
    nservers = ones(Mq, 1);
    refstat = zeros(1, R);
    for k = 1:Mq
        ist = qs(k);
        nservers(k) = sn.nservers(ist);
        for r = 1:R
            if isfinite(sn.rates(ist, r)) && sn.rates(ist, r) > 0
                mu(k, r) = sn.rates(ist, r);
                if isfinite(sn.scv(ist, r)) && sn.scv(ist, r) > 0
                    Cs(k, r) = sn.scv(ist, r);
                end
            end
        end
    end
    for r = 1:R
        if ~openCls(r)
            refstat(r) = find(qs == sn.refstat(r), 1);
        end
    end
    % External arrivals of the open classes along the source routing
    lambda0 = zeros(Mq, R);
    Ca0 = zeros(Mq, R);
    sourceNode = sn.stationToNode(sourceIdx);
    for r = 1:R
        if openCls(r) && isfinite(sn.rates(sourceIdx, r)) && sn.rates(sourceIdx, r) > 0
            extRate = sn.rates(sourceIdx, r);
            Ca_ext = 1.0;
            if isfinite(sn.scv(sourceIdx, r)) && sn.scv(sourceIdx, r) > 0
                Ca_ext = sn.scv(sourceIdx, r);
            end
            for k = 1:Mq
                destNode = sn.stationToNode(qs(k));
                routeProb = sn.rtnodes((sourceNode-1)*R + r, (destNode-1)*R + r);
                if routeProb > 0
                    lambda0(k, r) = extRate * routeProb;
                    Ca0(k, r) = Ca_ext;
                end
            end
        end
    end
    % Routing probabilities between queueing stations
    P = zeros(Mq, Mq, R);
    for r = 1:R
        for j = 1:Mq
            jNode = sn.stationToNode(qs(j));
            for k = 1:Mq
                iNode = sn.stationToNode(qs(k));
                P(j, k, r) = sn.rtnodes((jNode-1)*R + r, (iNode-1)*R + r);
            end
        end
    end
    insens = false(Mq, 1);
    for k = 1:Mq
        insens(k) = any(sn.sched(qs(k)) == [SchedStrategy.PS, SchedStrategy.LCFSPR]);
    end
    [L, W, ~, ~, lam, rho, X, totiter] = me_mqn(Mq, R, openCls, lambda0, Ca0, N, mu, Cs, P, nservers, refstat, insens, mem_options);
    QN = zeros(M, R);
    UN = zeros(M, R);
    RN = zeros(M, R);
    TN = zeros(M, R);
    QN(qs, :) = L;
    UN(qs, :) = rho;
    RN(qs, :) = W;
    TN(qs, :) = lam;
    % Cap utilization of unstable stations at 1 (LINE convention)
    for k = 1:Mq
        ist = qs(k);
        if isfinite(sn.nservers(ist))
            Utot = sum(UN(ist, :));
            if Utot > 1
                UN(ist, :) = UN(ist, :) / Utot;
            end
        end
    end
    XN = zeros(1, R);
    CN = zeros(1, R);
    for r = 1:R
        if openCls(r)
            TN(sourceIdx, r) = X(r);
            XN(r) = X(r);
            if X(r) > 0
                CN(r) = sum(QN(qs, r)) / X(r);
            end
        else
            XN(r) = X(r);
            if X(r) > 0
                CN(r) = N(r) / X(r); % class cycle time
            end
        end
    end
    return
end

% Locate the source station (external arrivals)
sourceIdx = 0;
for ind = 1:sn.nnodes
    if sn.nodetype(ind) == NodeType.Source
        sourceIdx = sn.nodeToStation(ind);
        break;
    end
end
if sourceIdx <= 0
    line_error(mfilename, 'MEM requires a Source node.');
end

qs = setdiff(1:M, sourceIdx); % queueing/delay stations
Mq = length(qs);

% Service rates, scvs and server counts (rates/scv are NaN for classes
% not served at a station)
mu = zeros(Mq, R);
Cs = ones(Mq, R);
nservers = ones(Mq, 1);
for k = 1:Mq
    ist = qs(k);
    nservers(k) = sn.nservers(ist);
    for r = 1:R
        if isfinite(sn.rates(ist, r)) && sn.rates(ist, r) > 0
            mu(k, r) = sn.rates(ist, r);
            if isfinite(sn.scv(ist, r)) && sn.scv(ist, r) > 0
                Cs(k, r) = sn.scv(ist, r);
            end
        end
    end
end

% External arrivals: distribute the source output along its routing
lambda0 = zeros(Mq, R);
Ca0 = zeros(Mq, R);
sourceNode = sn.stationToNode(sourceIdx);
for r = 1:R
    if isfinite(sn.rates(sourceIdx, r)) && sn.rates(sourceIdx, r) > 0
        extRate = sn.rates(sourceIdx, r);
        Ca_ext = 1.0;
        if isfinite(sn.scv(sourceIdx, r)) && sn.scv(sourceIdx, r) > 0
            Ca_ext = sn.scv(sourceIdx, r);
        end
        for k = 1:Mq
            destNode = sn.stationToNode(qs(k));
            routeProb = sn.rtnodes((sourceNode-1)*R + r, (destNode-1)*R + r);
            if routeProb > 0
                lambda0(k, r) = extRate * routeProb;
                Ca0(k, r) = Ca_ext;
            end
        end
    end
end

% Routing probabilities between queueing stations
P = zeros(Mq, Mq, R);
for r = 1:R
    for j = 1:Mq
        jNode = sn.stationToNode(qs(j));
        for k = 1:Mq
            iNode = sn.stationToNode(qs(k));
            P(j, k, r) = sn.rtnodes((jNode-1)*R + r, (iNode-1)*R + r);
        end
    end
end

% Finite buffers: censored GE/GE/c/0;N building blocks, with the
% holding-node expansion where the drop rule is BAS (Section 4.1 of the
% source plus Tahilramani, Manjunath and Bose 1999)
if memblocking
    actualmethod = 'mem.blocking';
    Nbuf = zeros(Mq, 1);
    blockrule = zeros(Mq, 1);
    for k = 1:Mq
        ist = qs(k);
        Nbuf(k) = sn_get_buffer_size(sn, ist);
        dr = DropStrategy.DROP;
        if ~isempty(sn.droprule) && size(sn.droprule, 1) >= ist
            dr = sn.droprule(ist, 1);
        end
        if dr == DropStrategy.BAS
            blockrule(k) = 1;
        end
    end
    [L, W, Tk, rho, ~, ~, PBa, ~, totiter] = me_oqn_blk(Mq, lambda0(:, 1), Ca0(:, 1), ...
        mu(:, 1), Cs(:, 1), P(:, :, 1), nservers, Nbuf, blockrule, mem_options);
    QN = zeros(M, R);
    UN = zeros(M, R);
    RN = zeros(M, R);
    TN = zeros(M, R);
    QN(qs, 1) = L;
    UN(qs, 1) = rho;
    RN(qs, 1) = W;
    TN(qs, 1) = Tk;
    % The source emits at its nominal rate; the jobs lost at a full buffer
    % never reach a station, so the carried flow reported per station is
    % below it by the loss probability at the entry stations.
    XN = zeros(1, R);
    CN = zeros(1, R);
    if isfinite(sn.rates(sourceIdx, 1)) && sn.rates(sourceIdx, 1) > 0
        TN(sourceIdx, 1) = sn.rates(sourceIdx, 1);
        XN(1) = sn.rates(sourceIdx, 1);
        accepted = XN(1);
        for k = 1:Mq
            if lambda0(k, 1) > 0
                accepted = accepted - lambda0(k, 1) * PBa(k) * (blockrule(k) == 0);
            end
        end
        if accepted > 0
            CN(1) = sum(QN(qs, 1)) / accepted;
        end
    end
    return
end

% Run the Maximum Entropy fixed-point algorithm
insens = false(Mq, 1);
for k = 1:Mq
    insens(k) = any(sn.sched(qs(k)) == [SchedStrategy.PS, SchedStrategy.LCFSPR]);
end
[L, W, ~, ~, lam, rho, totiter] = me_oqn(Mq, R, lambda0, Ca0, mu, Cs, P, nservers, insens, mem_options);

% Map results back to station-indexed LINE outputs
QN = zeros(M, R);
UN = zeros(M, R);
RN = zeros(M, R);
TN = zeros(M, R);
QN(qs, :) = L;
UN(qs, :) = rho;
RN(qs, :) = W;
TN(qs, :) = lam;

% Cap utilization of unstable stations at 1 (LINE convention)
for k = 1:Mq
    ist = qs(k);
    if isfinite(sn.nservers(ist))
        Utot = sum(UN(ist, :));
        if Utot > 1
            UN(ist, :) = UN(ist, :) / Utot;
        end
    end
end

% Source station: report the external arrival rates as throughputs and
% derive system metrics by Little's law
XN = zeros(1, R);
CN = zeros(1, R);
for r = 1:R
    if isfinite(sn.rates(sourceIdx, r)) && sn.rates(sourceIdx, r) > 0
        TN(sourceIdx, r) = sn.rates(sourceIdx, r);
        XN(r) = sn.rates(sourceIdx, r);
        CN(r) = sum(QN(qs, r)) / XN(r);
    end
end

end
