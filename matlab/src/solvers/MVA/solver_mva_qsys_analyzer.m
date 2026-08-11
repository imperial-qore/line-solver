function [Q,U,R,T,C,X,lG,runtime,totiter,actualmethod] = solver_mva_qsys_analyzer(sn, options)
% [Q,U,R,T,C,X,LG,RUNTIME,ITER,ACTUALMETHOD] = SOLVER_MVA_QSYS_ANALYZER(QN, OPTIONS)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

T0=tic;
M = sn.nstations;
K = sn.nclasses;
Q = zeros(M,K); U = zeros(M,K);
R = zeros(M,K); T = zeros(M,K);
C = zeros(1,K); X = zeros(M,K);
totiter = 1;

method = options.method;
actualmethod = method; % resolved algorithm name, reported in the solver banner
source_ist = sn.nodeToStation(sn.nodetype == NodeType.Source);
queue_ist = sn.nodeToStation(sn.nodetype == NodeType.Queue);
lambda = sn.rates(source_ist)*sn.visits{source_ist}(sn.stationToStateful(queue_ist));
k = sn.nservers(queue_ist);
mu = sn.rates(queue_ist);
ca = sqrt(sn.scv(source_ist));
cs = sqrt(sn.scv(queue_ist));

line_debug('MVA qsys analyzer starting: method=%s, lambda=%g, mu=%g, k=%d', method, lambda, mu, k);

% Finite-capacity loss branch (M/M/1/K with tail drop). Uses the moment-based
% (MacGregor Smith) qsys_mg1k_loss_mgs, exact only at scv=1; queue-length
% metrics come from the truncated M/M/1/K distribution. Being an approximation
% in general, it is not offered under method='exact'.
if sn_is_mm1k_loss(sn)
    if strcmp(method, 'exact')
        line_error(mfilename, 'M/M/1/K tail-drop is solved by the approximate ''mg1k.mgs'' method (MacGregor Smith); it is not available under method=''exact''. Use the default method, or SolverCTMC/SolverNC for an exact result.');
    end
    Kcap = sn.cap(queue_ist);
    rho = lambda/mu;
    Ploss = qsys_mg1k_loss_mgs(lambda, mu, cs^2, Kcap);
    Tq = lambda*(1-Ploss);              % carried throughput
    Uq = Tq/mu;                         % single-server utilization
    if abs(rho-1) < 1e-10
        Lsys = Kcap/2;                  % L'Hopital limit at rho=1
    else
        Lsys = rho/(1-rho) - (Kcap+1)*rho^(Kcap+1)/(1-rho^(Kcap+1));
    end
    Vq = sn.visits{1}(sn.stationToStateful(queue_ist));
    R(queue_ist,1) = Lsys/Tq;           % per-visit response time (Little)
    Q(queue_ist,1) = Lsys;
    U(queue_ist,1) = Uq;
    T(queue_ist,1) = Tq;                % carried (effective) rate
    T(source_ist,1) = lambda;          % offered arrival rate
    X(queue_ist,1) = Tq;               % system throughput = carried rate
    C(1,1) = R(queue_ist,1)*Vq;
    actualmethod = 'mg1k.mgs';
    lG = 0; totiter = 1; runtime = toc(T0);
    return
end

% Check for BMAP arrivals (batch Markovian)
if sn.procid(source_ist) == ProcessType.BMAP
    line_debug('BMAP arrival process detected');

    % Check if service is exponential
    if sn.procid(queue_ist) == ProcessType.EXP
        line_debug('Service is exponential, using MX/M/1 queue model');
        % sn.proc holds {D0, D1, D_batch1, ..., D_batchK} (JAR MatrixCell
        % layout); batch rates are weighted by the stationary vector of
        % the underlying CTMC, as in BMAP.getBatchRates
        proc_src = sn.proc{source_ist};
        if iscell(proc_src) && ~isempty(proc_src) && iscell(proc_src{1})
            proc_src = proc_src{1};
        end
        pie_b = ctmc_solve(proc_src{1} + proc_src{2});
        nbatch = length(proc_src) - 2;
        batch_rates = zeros(1, nbatch);
        for bsize = 1:nbatch
            batch_rates(bsize) = pie_b * proc_src{2+bsize} * ones(size(proc_src{1},1),1);
        end
        total_rate = sum(batch_rates);
        lambda_batch = total_rate; % batch event rate
        if total_rate > 0
            E_X = sum((1:nbatch) .* batch_rates) / total_rate;
            E_X2 = sum(((1:nbatch).^2) .* batch_rates) / total_rate;
        else
            E_X = 1;
            E_X2 = 1;
        end

        [W, Wq, U_mxm1, Q_mxm1] = qsys_mxm1(lambda_batch, mu, E_X, E_X2);
        R(queue_ist,1) = W * sn.visits{1}(sn.stationToStateful(queue_ist));
        C(1,1) = R(queue_ist,1);
        X(queue_ist,1) = lambda_batch * E_X;  % Job arrival rate
        U(queue_ist,1) = U_mxm1;
        T(source_ist,1) = lambda_batch * E_X;
        T(queue_ist,1) = lambda_batch * E_X;
        Q(queue_ist,1) = Q_mxm1;
        lG = 0;
        runtime=toc(T0);
        actualmethod = 'mxm1';
        return;
    end
end

if strcmpi(method,'exact')
    if ca == 1 && cs == 1 && k==1
        method = 'mm1';
        line_debug('Exact method selected: M/M/1 (ca=1, cs=1, k=1)');
    elseif ca == 1 && cs == 1 && k>1
        method = 'mmk';
        line_debug('Exact method selected: M/M/k (ca=1, cs=1, k=%d)', k);
    elseif ca == 1 && k==1
        method = 'mg1';
        line_debug('Exact method selected: M/G/1 (ca=1, k=1)');
    elseif cs == 1 && k==1
        method = 'gm1';
        line_debug('Exact method selected: G/M/1 (cs=1, k=1)');
    else
        line_error(mfilename,'MVA exact method unavailable for this model.');
    end
end

switch method
    case 'default'
        if ca == 1 && cs == 1 && k == 1
            method = 'mm1';
            line_debug('Default method: using M/M/1 exact solution\n');
        elseif ca == 1 && cs == 1 && k > 1
            method = 'mmk';
            line_debug('Default method: using M/M/k exact solution (k=%d)\n', k);
        elseif ca == 1 && k == 1
            method = 'mg1';
            line_debug('Default method: using M/G/1 exact solution\n');
        elseif cs == 1 && k == 1
            method = 'gm1';
            line_debug('Default method: using G/M/1 exact solution\n');
        elseif k > 1
            method = 'gigk';
            line_debug('Default method: using G/G/k approximation (k=%d)\n', k);
        else
            method = 'gig1.klb';
            line_debug('Default method: using G/G/1 KLB approximation\n');
        end
end

switch method
    case 'mm1'
        line_debug('Using M/M/1 exact solution');
        R = qsys_mm1(lambda,mu);
    case 'mmk'
        line_debug('Using M/M/k exact solution (k=%d)', k);
        R = qsys_mmk(lambda,mu,k);
    case {'rqna'}
        line_debug('Using RQNA (robust queueing) single-queue solution');
        arvMAP = sn.proc{source_ist}{1};
        rho1 = lambda/mu;
        IaFun1 = @(x) map_count_idc(arvMAP, x);
        [~, W1] = qsys_gig1_rq(rho1, mu, cs^2, IaFun1);
        R = W1 + 1/mu;
    case {'mg1', 'mgi1'}  % verified
        line_debug('Using M/G/1 exact solution');
        R = qsys_mg1(lambda,mu,cs);
    case {'gigk'}
        line_debug('Using G/G/k approximation (k=%d)', k);
        R = qsys_gigk_approx(lambda,mu,ca,cs,k);
    case {'gigk.kingman_approx'}
        line_debug('Using G/G/k Kingman approximation (k=%d)', k);
        R = qsys_gigk_approx_kingman(lambda,mu,ca,cs,k);
    case 'gig1.kingman'  % verified
        line_debug('Using G/G/1 Kingman upper bound');
        R = qsys_gig1_ubnd_kingman(lambda,mu,ca,cs);
    case 'gig1.heyman'
        line_debug('Using G/G/1 Heyman approximation');
        R = qsys_gig1_approx_heyman(lambda,mu,ca,cs);
    case {'gig1', 'gig1.allen'}
        line_debug('Using G/G/1 Allen-Cunneen approximation');
        R = qsys_gig1_approx_allencunneen(lambda,mu,ca,cs);
    case 'gig1.kobayashi'
        line_debug('Using G/G/1 Kobayashi approximation');
        R = qsys_gig1_approx_kobayashi(lambda,mu,ca,cs);
    case 'gig1.klb'
        line_debug('Using G/G/1 KLB approximation');
        R = qsys_gig1_approx_klb(lambda,mu,ca,cs);
    case 'gig1.marchal' % verified
        line_debug('Using G/G/1 Marchal approximation');
        R = qsys_gig1_approx_marchal(lambda,mu,ca,cs);
    case {'gm1', 'gim1'}
        line_debug('Using G/M/1 exact solution');
        mu = sn.rates(queue_ist);
        % Prefer the exact PH/M/1 sigma-root, then LST-based fzero, then the
        % two-moment qsys_gg1 fit. The PH path is exact only for a Markovian
        % arrival law; see _kb/06-solver-catalog.md (MVA section, gm1 note).
        R = [];
        if ProcessType.isMarkovian(sn.procid(source_ist))
            try
                phPair = sn.proc{source_ist}{1};
                if iscell(phPair) && numel(phPair) >= 2
                    D0_ph = phPair{1};
                    D1_ph = phPair{2};
                    if size(D0_ph, 1) == size(D0_ph, 2) && all(size(D0_ph) == size(D1_ph))
                        pie_src = map_pie({D0_ph, D1_ph});
                        res_phm1 = qsys_phm1(pie_src, D0_ph, mu);
                        R = res_phm1.meanSojournTime;
                    end
                end
            catch
                R = [];
            end
        end
        if isempty(R)
            try
                LA = @(s) sn.lst{source_ist}{1}(s);
                sigma = fzero(@(x) LA(mu-mu*x)-x, 0.5);
                R = qsys_gm1(sigma, mu);
            catch
                R = qsys_gg1(lambda, mu, ca^2, 1);
            end
        end
    otherwise
        line_error(mfilename,'Unsupported method for a model with 1 station and 1 class.');
end
actualmethod = method;

Rscalar = R;  % save scalar from qsys function
% per-visit vs per-job: RespT/QLen are per-visit, ResidT and system response C are
% per-job (=per-visit * Vq); see _kb/06-solver-catalog.md (MVA section, gm1 note).
Vq = sn.visits{1}(sn.stationToStateful(queue_ist));
srcRate = sn.rates(source_ist);
R = zeros(M,K);
R(queue_ist,1) = Rscalar;              % per-visit response time
C(1,1) = Rscalar * Vq;                 % system (per-job) response time
X(queue_ist,1) = srcRate;              % system throughput = external arrival rate
U(queue_ist,1) = lambda/mu/k;
T(source_ist,1) = srcRate;             % source throughput = external arrival rate
T(queue_ist,1) = lambda;               % queue throughput = effective arrival rate
Q(queue_ist,1) = lambda * Rscalar;     % Little's law at the queue (per-visit)
lG = 0;
runtime=toc(T0);
end
