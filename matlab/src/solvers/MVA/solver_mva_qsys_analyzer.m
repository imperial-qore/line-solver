function [Q,U,R,T,C,X,lG,runtime,totiter] = solver_mva_qsys_analyzer(sn, options)
% [Q,U,R,T,C,X,LG,RUNTIME,ITER] = SOLVER_MVA_QSYS_ANALYZER(QN, OPTIONS)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

T0=tic;
Q = []; U = [];
R = []; T = [];
C = []; X = [];
totiter = 1;

method = options.method;
source_ist = sn.nodeToStation(sn.nodetype == NodeType.Source);
queue_ist = sn.nodeToStation(sn.nodetype == NodeType.Queue);
lambda = sn.rates(source_ist)*sn.visits{source_ist}(sn.stationToStateful(queue_ist));
k = sn.nservers(queue_ist);
mu = sn.rates(queue_ist);
ca = sqrt(sn.scv(source_ist));
cs = sqrt(sn.scv(queue_ist));

line_debug('MVA qsys analyzer starting: method=%s, lambda=%g, mu=%g, k=%d', method, lambda, mu, k);

% Check for BMAP arrivals (batch Markovian)
arvproc = sn.proc{source_ist};
if iscell(arvproc) && length(arvproc) >= 1
    arvproc = arvproc{1};
end

if isa(arvproc, 'BMAP')
    line_debug('BMAP arrival process detected');
    srvproc = sn.proc{queue_ist};
    if iscell(srvproc) && length(srvproc) >= 1
        srvproc = srvproc{1};
    end

    % Check if service is exponential
    if isa(srvproc, 'Exp')
        line_debug('Service is exponential, using MX/M/1 queue model');
        lambda_batch = arvproc.getRate();
        E_X = arvproc.getMeanBatchSize();

        % Compute second moment of batch size
        % For BMAP, we need to extract this from the D matrices
        batch_rates = arvproc.getBatchRates();
        total_rate = sum(batch_rates);
        if total_rate > 0
            E_X2 = 0;
            for bsize = 1:length(batch_rates)
                E_X2 = E_X2 + (bsize^2) * batch_rates(bsize) / total_rate;
            end
        else
            E_X2 = E_X^2;  % Fallback to deterministic batch size
        end

        [W, Wq, U_mxm1, Q_mxm1] = qsys_mxm1(lambda_batch, mu, E_X, E_X2);
        R(queue_ist,1) = W * sn.visits{1}(sn.stationToStateful(queue_ist));
        C(queue_ist,1) = R(1,1);
        X(queue_ist,1) = lambda_batch * E_X;  % Job arrival rate
        U(queue_ist,1) = U_mxm1;
        T(source_ist,1) = lambda_batch * E_X;
        T(queue_ist,1) = lambda_batch * E_X;
        Q(queue_ist,1) = Q_mxm1;
        lG = 0;
        runtime=toc(T0);
        method = 'mxm1';
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
    case {'mg1', 'mgi1'}  % verified
        line_debug('Using M/G/1 exact solution');
        R = qsys_mg1(lambda,mu,cs);
    case {'gigk'}
        line_debug('Using G/G/k approximation (k=%d)', k);
        R = qsys_gigk_approx(lambda,mu,ca,cs,k);
    case {'gigk.kingman_approx'}
        line_debug('Using G/G/k Kingman approximation (k=%d)', k);
        R = qsys_gigk_approx_kingman(lambda,mu,ca,cs,k);
    case {'gig1', 'gig1.kingman'}  % verified
        line_debug('Using G/G/1 Kingman upper bound');
        R = qsys_gig1_ubnd_kingman(lambda,mu,ca,cs);
    case 'gig1.heyman'
        line_debug('Using G/G/1 Heyman approximation');
        R = qsys_gig1_approx_heyman(lambda,mu,ca,cs);
    case 'gig1.allen'
        line_debug('Using G/G/1 Allen-Cunneen approximation');
        R = qsys_gig1_approx_allencunneen(lambda,mu,ca,cs);
    case 'gig1.kobayashi'
        line_debug('Using G/G/1 Kobayashi approximation');
        R = qsys_gig1_approx_kobayashi(lambda,mu,ca,cs);
    case 'gig1.klb'
        line_debug('Using G/G/1 KLB approximation');
        R = qsys_gig1_approx_klb(lambda,mu,ca,cs);
        if strcmpi(options.method,'default')
            method = sprintf('default [%s]','gig1.klb');
        end
    case 'gig1.marchal' % verified
        line_debug('Using G/G/1 Marchal approximation');
        R = qsys_gig1_approx_marchal(lambda,mu,ca,cs);
    case {'gm1', 'gim1'}
        line_debug('Using G/M/1 exact solution');
        % sigma = Load at arrival instants (Laplace transform of the inter-arrival times)
        LA = @(s) sn.lst{source_ist,1}(s);
        mu = sn.rates(queue_ist);
        sigma = fzero(@(x) LA(mu-mu*x)-x,0.5);
        R = qsys_gm1(sigma,mu);
    otherwise
        line_error(mfilename,'Unsupported method for a model with 1 station and 1 class.');
end

R(queue_ist,1) = R *sn.visits{1}(sn.stationToStateful(queue_ist));
C(queue_ist,1) = R(1,1);
X(queue_ist,1) = lambda;
U(queue_ist,1) = lambda/mu/k;
T(source_ist,1) = lambda;
T(queue_ist,1) = lambda;
Q(queue_ist,1) = X(queue_ist,1) * R(queue_ist,1);
lG = 0;
runtime=toc(T0);
end
