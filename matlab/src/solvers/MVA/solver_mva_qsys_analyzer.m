function [Q,U,R,T,C,X,lG,runtime,totiter] = solver_mva_qsys_analyzer(sn, options)
% [Q,U,R,T,C,X,LG,RUNTIME,ITER] = SOLVER_MVA_QSYS_ANALYZER(QN, OPTIONS)

% Copyright (c) 2012-2024, Imperial College London
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
if strcmpi(method,'exact')
    if ca == 1 && cs == 1 && k==1
        method = 'mm1';
    elseif ca == 1 && cs == 1 && k>1
        method = 'mmk';
    elseif ca == 1 && k==1
        method = 'mg1';
    elseif cs == 1 && k==1
        method = 'gm1';
    else
        line_error(mfilename,'MVA exact method unavailable for this model.');
    end
end

switch method
    case 'default'
        if k>1
            method = 'gigk';
        else
            method = 'gig1.klb';
        end
end

switch method
    case 'mm1'
        R = qsys_mm1(lambda,mu);
    case 'mmk'
        R = qsys_mmk(lambda,mu,k);
    case {'mg1', 'mgi1'}  % verified
        R = qsys_mg1(lambda,mu,cs);
    case {'gigk'}
        R = qsys_gigk_approx(lambda,mu,ca,cs,k);
    case {'gigk.kingman_approx'}
        R = qsys_gigk_approx_kingman(lambda,mu,ca,cs,k);
    case {'gig1', 'gig1.kingman'}  % verified
        R = qsys_gig1_ubnd_kingman(lambda,mu,ca,cs);
    case 'gig1.heyman'
        R = qsys_gig1_approx_heyman(lambda,mu,ca,cs);
    case 'gig1.allen'
        R = qsys_gig1_approx_allencunneen(lambda,mu,ca,cs);
    case 'gig1.kobayashi'
        R = qsys_gig1_approx_kobayashi(lambda,mu,ca,cs);
    case 'gig1.klb'
        R = qsys_gig1_approx_klb(lambda,mu,ca,cs);
        if strcmpi(options.method,'default')
            method = sprintf('default [%s]','gig1.klb');
        end
    case 'gig1.marchal' % verified
        R = qsys_gig1_approx_marchal(lambda,mu,ca,cs);
    case {'gm1', 'gim1'}
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
