function [QN,UN,RN,TN,CN,XN,totiter] = solver_nc_mem(sn, options)
% [QN,UN,RN,TN,CN,XN,TOTITER] = SOLVER_NC_MEM(SN, OPTIONS)
%
% Maximum Entropy Method (MEM) for Open Queueing Networks
%
% Implements the ME algorithm from Kouvatsos (1994) for analyzing open
% queueing networks with general arrival and service processes.
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
%   TN - Throughputs (arrival rates) [M x R matrix]
%   CN - Completion times (not computed by MEM, returns zeros)
%   XN - System throughputs (not computed by MEM, returns zeros)
%   totiter - Number of iterations until convergence
%
% Reference:
%   D.D. Kouvatsos, "Entropy Maximisation and Queueing Network Models",
%   Annals of Operations Research, 48:63-126, 1994.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

% Validate that model is open
if ~sn_is_open_model(sn)
    line_error(mfilename, 'MEM only supports open queueing networks (all classes must be open).');
end

if sn_is_closed_model(sn)
    line_error(mfilename, 'MEM does not support models with closed classes.');
end

% Extract dimensions
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

% Initialize parameter matrices for me_oqn
lambda0 = zeros(M, R);  % External arrival rates
Ca0 = zeros(M, R);      % External arrival SCV
mu = zeros(M, R);       % Service rates
Cs = zeros(M, R);       % Service SCV
P = zeros(M, M, R);     % Routing probabilities

% Extract service rates and service SCVs
for i = 1:M
    for r = 1:R
        if sn.rates(i, r) > 0
            mu(i, r) = sn.rates(i, r);

            % Extract service SCV from service process
            if ~isempty(sn.proc) && i <= length(sn.proc) && ~isempty(sn.proc{i}) && ...
               r <= length(sn.proc{i}) && ~isempty(sn.proc{i}{r})
                MAP = sn.proc{i}{r};
                if iscell(MAP) && length(MAP) >= 2
                    Cs(i, r) = map_scv(MAP);
                else
                    Cs(i, r) = 1.0;  % Default to exponential
                end
            else
                Cs(i, r) = 1.0;  % Default to exponential
            end
        end
    end
end

% Find source node for external arrivals
sourceIdx = [];
for i = 1:sn.nnodes
    if sn.nodetype(i) == NodeType.Source
        sourceIdx = sn.nodeToStation(i);
        break;
    end
end

% Extract external arrival rates and SCVs
if ~isempty(sourceIdx) && sourceIdx > 0
    for r = 1:R
        if sn.rates(sourceIdx, r) > 0
            externalRate = sn.rates(sourceIdx, r);
            Ca_external = 1.0;

            % Extract arrival SCV from arrival process
            if ~isempty(sn.proc) && sourceIdx <= length(sn.proc) && ~isempty(sn.proc{sourceIdx}) && ...
               r <= length(sn.proc{sourceIdx}) && ~isempty(sn.proc{sourceIdx}{r})
                MAP = sn.proc{sourceIdx}{r};
                if iscell(MAP) && length(MAP) >= 2
                    Ca_external = map_scv(MAP);
                end
            end

            % Distribute external arrivals based on routing from source
            for i = 1:M
                if i ~= sourceIdx
                    % Check if queue i receives external arrivals
                    if ~isempty(sn.rt)
                        sourceNode = sn.stationToNode(sourceIdx);
                        queueNode = sn.stationToNode(i);

                        if sourceNode > 0 && queueNode > 0
                            routeIdx = (sourceNode-1)*R + r;
                            destIdx = (queueNode-1)*R + r;

                            if routeIdx <= size(sn.rt, 1) && destIdx <= size(sn.rt, 2)
                                routeProb = sn.rt(routeIdx, destIdx);
                                if routeProb > 0
                                    lambda0(i, r) = externalRate * routeProb;
                                    Ca0(i, r) = Ca_external;
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

% Extract routing probabilities between queues
if ~isempty(sn.rt)
    for r = 1:R
        for j = 1:M
            for i = 1:M
                jNode = sn.stationToNode(j);
                iNode = sn.stationToNode(i);

                if jNode > 0 && iNode > 0
                    % Skip Source and Sink nodes
                    if sn.nodetype(jNode) ~= NodeType.Source && ...
                       sn.nodetype(jNode) ~= NodeType.Sink && ...
                       sn.nodetype(iNode) ~= NodeType.Source && ...
                       sn.nodetype(iNode) ~= NodeType.Sink

                        routeIdx = (jNode-1)*R + r;
                        destIdx = (iNode-1)*R + r;

                        if routeIdx <= size(sn.rt, 1) && destIdx <= size(sn.rt, 2)
                            P(j, i, r) = sn.rt(routeIdx, destIdx);
                        end
                    end
                end
            end
        end
    end
end

% Call Maximum Entropy algorithm
[L, W, Ca, Cd, lambda, rho] = me_oqn(M, R, lambda0, Ca0, mu, Cs, P, mem_options);

% Map to LINE standard output format
QN = L;          % Queue lengths
UN = rho;        % Utilizations
RN = W;          % Response times (waiting times)
TN = lambda;     % Throughputs (arrival rates)
CN = zeros(1, R);  % Completion times (not computed by MEM)
XN = zeros(1, R);  % System throughputs (not computed by MEM)

% Extract iteration count if available
% (me_oqn doesn't currently return this, so we return maxiter as placeholder)
totiter = mem_options.maxiter;

end
