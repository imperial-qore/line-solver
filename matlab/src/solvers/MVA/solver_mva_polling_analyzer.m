function [Q,U,R,T,C,X,lG,runtime,totiter] = solver_mva_polling_analyzer(sn, options)
% [Q,U,R,T,C,X,LG,RUNTIME,ITER] = SOLVER_MVA_POLLING_ANALYZER(QN, OPTIONS)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

T0=tic;
Q = []; U = [];
R = []; T = [];
C = []; X = [];
totiter = 1;

method = options.method;

line_debug('MVA polling analyzer starting: method=%s, nclasses=%d', method, sn.nclasses);

source_ist = sn.nodeToStation(sn.nodetype == NodeType.Source);
queue_ist = sn.nodeToStation(sn.nodetype == NodeType.Queue);
lambda = sn.rates(source_ist,:)*sn.visits{source_ist}(sn.stationToStateful(queue_ist));
k = sn.nservers(queue_ist);
mu = sn.rates(queue_ist,:);
ca = sqrt(sn.scv(source_ist,:));
cs = sqrt(sn.scv(queue_ist,:));
switchover = cell(1,sn.nclasses);
for r=1:sn.nclasses
    if isnan(mu(r))
        polling_type(r) = NaN;
        switchover{r} = Immediate();
    else
        polling_type(r) = PollingType.toId(sn.nodeparam{queue_ist}{r}.pollingType); % we assume this is identical across all buffers
        switchover{r} = sn.nodeparam{queue_ist}{r}.switchoverTime;
    end
end

if strcmpi(method,'exact')
    if all(ca == 1) && k==1 && ((polling_type == PollingType.EXHAUSTIVE) || (polling_type == PollingType.GATED))
        method = 'stationtime';
    else
        line_error(mfilename,'MVA exact method unavailable for this model.');
    end
end

switch method
    case {'default','stationtime'}
        if strcmpi(method, 'default')
            line_debug('Default method: using stationtime analysis for polling\n');
        end
        line_debug('Using stationtime method for polling analysis');
        switch max(polling_type,[], 'omitnan') % we assume polling types to be identical
            case PollingType.EXHAUSTIVE
                W = polling_qsys_exhaustive(sn.proc{source_ist},sn.proc{queue_ist},switchover);
            case PollingType.GATED
                W = polling_qsys_gated(sn.proc{source_ist},sn.proc{queue_ist},switchover);
            case PollingType.KLIMITED
                K = sn.nodeparam{sn.stationToNode(queue_ist)}{1}.pollingPar;
                if K==1
                    W = polling_qsys_1limited(sn.proc{source_ist},sn.proc{queue_ist},switchover);
                else
                    line_error(mfilename,'MVA method unavailable for K-limited polling with K>1.');
                end
            otherwise
                line_error(mfilename,'Unsupported polling type.');
        end
    otherwise
        line_error(mfilename,'Unsupported polling solution method.');
end

R = zeros(1,sn.nclasses);
for r=1:sn.nclasses
    R(r) = W(r) + 1/mu(r);
end

R(queue_ist,1:sn.nclasses) = R .*sn.visits{1}(sn.stationToStateful(queue_ist));
C(queue_ist,1:sn.nclasses) = R(1,:);
X(queue_ist,1:sn.nclasses) = lambda;
U(queue_ist,1:sn.nclasses) = lambda./mu/k;
T(source_ist,1:sn.nclasses) = lambda;
T(queue_ist,1:sn.nclasses) = lambda;
Q(queue_ist,1:sn.nclasses) = X(queue_ist,1:sn.nclasses) .* R(queue_ist,1:sn.nclasses);
lG = 0;
runtime=toc(T0);
end
