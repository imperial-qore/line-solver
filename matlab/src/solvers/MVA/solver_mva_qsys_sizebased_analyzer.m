function [Q,U,R,T,C,X,lG,runtime,totiter] = solver_mva_qsys_sizebased_analyzer(sn, options, schedType)
% [Q,U,R,T,C,X,LG,RUNTIME,ITER] = SOLVER_MVA_QSYS_SIZEBASED_ANALYZER(SN, OPTIONS, SCHEDTYPE)
%
% Analyzer for M/G/1 queueing systems with size-based scheduling.
% Supports: SRPT, PSJF, FB/LAS, LRPT, SETF
%
% This function handles multiclass open queueing systems with size-based
% scheduling policies using the analytical formulas from:
%   A. Wierman and M. Harchol-Balter, "Classifying scheduling policies with
%   respect to unfairness in an M/GI/1", SIGMETRICS 2003.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

T0 = tic;
Q = []; U = [];
R = []; T = [];
C = []; X = [];
totiter = 1;

% Extract model parameters
source_ist = sn.nodeToStation(sn.nodetype == NodeType.Source);
queue_ist = sn.nodeToStation(sn.nodetype == NodeType.Queue);

K = sn.nclasses;
lambda = zeros(K, 1);
mu = zeros(K, 1);
cs = zeros(K, 1);

for k = 1:K
    lambda(k) = sn.rates(source_ist, k) * sn.visits{source_ist}(sn.stationToStateful(queue_ist), k);
    mu(k) = sn.rates(queue_ist, k);
    cs(k) = sqrt(sn.scv(queue_ist, k));
end

% Check for valid parameters
if any(lambda <= 0) || any(mu <= 0)
    line_error(mfilename, 'Invalid arrival or service rates (must be positive).');
end

% Compute total utilization
rho = sum(lambda ./ mu);
if rho >= 1
    line_warning(mfilename, 'System is unstable (rho = %.4f >= 1).', rho);
end

% Call appropriate qsys function based on scheduling type
switch schedType
    case SchedStrategy.SRPT
        line_debug('Using M/G/1/SRPT exact solution');
        [W, ~] = qsys_mg1_srpt(lambda, mu, cs);
    case SchedStrategy.PSJF
        line_debug('Using M/G/1/PSJF exact solution');
        [W, ~] = qsys_mg1_psjf(lambda, mu, cs);
    case SchedStrategy.FB
        line_debug('Using M/G/1/FB (LAS) exact solution');
        [W, ~] = qsys_mg1_fb(lambda, mu, cs);
    case SchedStrategy.LRPT
        line_debug('Using M/G/1/LRPT exact solution');
        [W, ~] = qsys_mg1_lrpt(lambda, mu, cs);
    case SchedStrategy.SETF
        line_debug('Using M/G/1/SETF (non-preemptive FB) exact solution');
        [W, ~] = qsys_mg1_setf(lambda, mu, cs);
    otherwise
        line_error(mfilename, 'Unsupported scheduling type for size-based analyzer.');
end

% Initialize result matrices
M = sn.nstations;
R = zeros(M, K);
Q = zeros(M, K);
U = zeros(M, K);
T = zeros(M, K);
C = zeros(M, K);
X = zeros(M, K);

% Populate results
for k = 1:K
    visits = sn.visits{source_ist}(sn.stationToStateful(queue_ist), k);
    R(queue_ist, k) = W(k) * visits;
    T(source_ist, k) = lambda(k);
    T(queue_ist, k) = lambda(k);
    X(queue_ist, k) = lambda(k);
    U(queue_ist, k) = lambda(k) / mu(k);
    Q(queue_ist, k) = lambda(k) * W(k);
    C(queue_ist, k) = R(queue_ist, k);
end

lG = 0;
runtime = toc(T0);
end
