function [QN,UN,RN,TN,xvec_it,QNt,UNt,TNt,xvec_t,t,iters,runtime] = solver_mfq_prio(sn, options)
%SOLVER_MFQ_PRIO Single-queue open priority model via the fluid priority queue.
%
% [QN,UN,RN,TN,xvec_it,QNt,UNt,TNt,xvec_t,t,iters,runtime] = SOLVER_MFQ_PRIO(sn, options)
%
% Analyzes a single-queue open system with class priorities using the fluid
% priority queue mfq_prio_queue (G. Horvath, "Efficient analysis of the
% MMAP[K]/PH[K]/1 priority queue", EJOR 246(1):128-139, 2015). The per-class
% Markovian arrival processes are superposed into a joint background CTMC that
% modulates the per-class fluid input rates; the single server drains fluid at
% the constant rate d = mu, serving higher-priority fluid first (preemptive).
%
% The fluid level is interpreted as the (per-class) queue length and the fluid
% sojourn time as the response time, consistent with solver_mfq. The method
% applies when arrivals are Markov-modulated (MAP) and the service rate is
% class-independent; simple-exponential arrivals or class-dependent service
% degenerate the fluid model and fall back to the matrix fluid method.
%
% See also: solver_mfq, mfq_prio_queue, fluid_is_single_queue

runtime = tic;
iters = 1;

M = sn.nstations;
K = sn.nclasses;
QN = zeros(M,K); UN = zeros(M,K); RN = zeros(M,K); TN = zeros(M,K);

t = [0; options.timespan(2)];
QNt = cell(M,K); UNt = cell(M,K); TNt = cell(M,K);
for ist = 1:M
    for k = 1:K
        QNt{ist,k} = zeros(2,1); UNt{ist,k} = zeros(2,1); TNt{ist,k} = zeros(2,1);
    end
end
xvec_t = [];
xvec_it = {zeros(size(sn.state{1}, 2), 1)};

    function [QN,UN,RN,TN,xvec_it,QNt,UNt,TNt,xvec_t,t] = fallback_matrix(reason)
        % The documented fallback. FLUID_MFQ_ADMITS states the same conditions
        % ahead of the run, so SolverFLD.resolveMethod labels the pair 'matrix'
        % wherever this branch would take it.
        line_warning(mfilename, 'MFQ-prio not applicable (%s); falling back to matrix method.', reason);
        opts = options; opts.method = 'matrix';
        [QN,UN,RN,TN,xvec_it,QNt,UNt,TNt,xvec_t,t] = solver_fluid_matrix(sn, opts);
    end

% ---- topology ----
[isSingleQueue, fluidInfo] = fluid_is_single_queue(sn);
if ~isSingleQueue
    line_error(mfilename, 'MFQ-prio requires single-queue topology: %s', fluidInfo.errorMsg);
end
sourceIdx = fluidInfo.sourceStation;
queueIdx  = fluidInfo.queueStation;

% ---- classes ordered for mfq_prio_queue (LINE classprio: lower value = higher
% priority; FluidPrioQueue uses a higher row index for higher priority, so the
% highest-priority class must occupy the last row) ----
openClasses = find(isinf(sn.njobs));
[~, ord] = sort(sn.classprio(openClasses), 'descend');
ordC = openClasses(ord);
Kc = numel(ordC);

% ---- constant service rate (class-independent required) ----
mu = sn.rates(queueIdx, ordC);
if Kc == 0 || any(~isfinite(mu)) || any(mu <= 0)
    [QN,UN,RN,TN,xvec_it,QNt,UNt,TNt,xvec_t,t] = fallback_matrix('invalid service rates'); return;
end
if any(abs(mu - mu(1)) > GlobalConstants.FineTol * max(1, mu(1)))
    [QN,UN,RN,TN,xvec_it,QNt,UNt,TNt,xvec_t,t] = fallback_matrix('class-dependent service'); return;
end
d = mu(1);

% ---- per-class arrival fluid representations (Q_k = D0+D1, R_k = D1*1) ----
Qk = cell(1,Kc); Rk = cell(1,Kc); Nk = zeros(1,Kc); lambda = zeros(1,Kc);
for i = 1:Kc
    k = ordC(i);
    proc = sn.proc{sourceIdx}{k};
    if ~iscell(proc) || numel(proc) < 2
        [QN,UN,RN,TN,xvec_it,QNt,UNt,TNt,xvec_t,t] = fallback_matrix('non-MAP arrival'); return;
    end
    D0 = proc{1}; D1 = proc{2};
    Qk{i} = D0 + D1;
    Rk{i} = sum(D1, 2);
    Nk(i) = size(D0, 1);
    lambda(i) = sn.rates(sourceIdx, k);
end

if prod(Nk) < 2
    [QN,UN,RN,TN,xvec_it,QNt,UNt,TNt,xvec_t,t] = fallback_matrix('non-modulated (exponential) arrivals'); return;
end

% ---- joint background CTMC (Kronecker sum) and per-class fluid rates ----
Njoint = prod(Nk);
Qjoint = zeros(Njoint);
for i = 1:Kc
    term = 1;
    for j = 1:Kc
        if j == i, term = kron(term, Qk{i}); else, term = kron(term, eye(Nk(j))); end
    end
    Qjoint = Qjoint + term;
end
Rjoint = zeros(Kc, Njoint);
for i = 1:Kc
    v = 1;
    for j = 1:Kc
        if j == i, v = kron(v, Rk{i}); else, v = kron(v, ones(Nk(j),1)); end
    end
    Rjoint(i,:) = v(:)';
end

prec = 1e-14;
if isfield(options, 'tol') && ~isempty(options.tol), prec = options.tol; end

% ---- fluid priority queue (last row = highest priority) ----
fl = cell(1,Kc); st = cell(1,Kc);
try
    [fl{:}] = mfq_prio_queue(Qjoint, Rjoint, d, 'classes', 1:Kc, 'flMoms', 1, 'prec', prec);
    [st{:}] = mfq_prio_queue(Qjoint, Rjoint, d, 'classes', 1:Kc, 'stMoms', 1, 'prec', prec);
catch ME
    [QN,UN,RN,TN,xvec_it,QNt,UNt,TNt,xvec_t,t] = fallback_matrix(ME.message); return;
end

% ---- map fluid measures to LINE metrics ----
for i = 1:Kc
    k = ordC(i);
    QN(queueIdx, k) = fl{i}(1);        % mean queue length = mean fluid level
    RN(queueIdx, k) = st{i}(1);        % mean response time = mean fluid sojourn
    TN(queueIdx, k) = lambda(i);       % throughput = arrival rate
    TN(sourceIdx, k) = lambda(i);
    UN(queueIdx, k) = min(1, lambda(i) / d);
    QNt{queueIdx, k} = [0; QN(queueIdx, k)];
    UNt{queueIdx, k} = [0; UN(queueIdx, k)];
    TNt{queueIdx, k} = [0; TN(queueIdx, k)];
    QNt{sourceIdx, k} = [0; 0];
    UNt{sourceIdx, k} = [0; 0];
    TNt{sourceIdx, k} = [0; TN(sourceIdx, k)];
end

runtime = toc(runtime);
end
