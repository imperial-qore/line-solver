function [Q,U,R,T,C,X,lG,totiter,method] = solver_amva(sn,options)
% [Q,U,R,T,C,X,lG,ITER] = SOLVER_AMVA(SN, OPTIONS)
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
if nargin < 2
    options = SolverMVA.defaultOptions;
end
%% aggregate chains
[Lchain,STchain,Vchain,alpha,Nchain,SCVchain,refstatchain] = sn_get_demands_chain(sn);

%% check options
if ~isfield(options.config,'np_priority')
    options.config.np_priority = 'default';
end
if ~isfield(options.config,'multiserver')
    options.config.multiserver = 'default';
end
if ~isfield(options.config,'highvar')
    options.config.highvar = 'default';
end

switch options.method
    case 'amva.qli'
        options.method = 'qli';
    case {'amva.qd', 'amva.qdamva', 'qdamva'}
        options.method = 'qd';
    case 'amva.aql'
        options.method = 'aql';
    case 'amva.qdaql'
        options.method = 'qdaql';
    case 'amva.lin'
        options.method = 'lin';
    case 'amva.qdlin'
        options.method = 'qdlin';
    case 'amva.fli'
        options.method = 'fli';
    case 'amva.bs'
        options.method = 'bs';
    case 'amva.ab'
        options.method = 'ab';
    case 'amva.schmidt'
        options.method = 'schmidt';
    case 'amva.schmidt-ext'
        options.method = 'schmidt-ext';
    case {'default','amva'}
        if (sum(Nchain)<=2 || any(Nchain<1))
            options.method = 'qd'; % changing to bs degrades accuracy
        else
            if max(sn.nservers(isfinite(sn.nservers)))==1
                options.method = 'egflin'; % if single server model
            else
                options.method = 'lin'; % lin seems way worse than aql in test_LQN_8.xml
            end
        end
end
method = options.method;

%% trivial models
if sn_has_homogeneous_scheduling(sn,SchedStrategy.INF)
    options.config.multiserver = 'default';
    [Q,U,R,T,C,X,lG,totiter] = solver_amvald(sn,Lchain,STchain,Vchain,alpha,Nchain,SCVchain,refstatchain,options);
    return
end

sourceIdx = sn.nodetype == NodeType.Source;
queueIdx = sn.nodetype == NodeType.Queue;
delayIdx = sn.nodetype == NodeType.Delay;

%% run amva method
M = sn.nstations;
%K = sn.nclasses;
C = sn.nchains;
V = zeros(M,C);
for i=sn.nodeToStation(sourceIdx)
    for c=1:sn.nchains
        inchain = find(sn.chains(c,:));
        if sum(sn.rates(sn.nodeToStation(sourceIdx),inchain), 'omitnan') > 0
            V(i,c) = 1;
        end
    end
end
Q = zeros(M,C);
U = zeros(M,C);

if sn_has_product_form_not_het_fcfs(sn) && ~sn_has_load_dependence(sn) && (~sn_has_open_classes(sn) || (sn_has_product_form(sn) && sn_has_open_classes(sn) && strcmpi(options.method,'lin')))
    % we can use linearizer only if the open model is not heterfcfs as that approximation is not supported, so strict product-form is required
    [lambda,L0,N,Z0,~,nservers,V(sn.nodeToStation(queueIdx|delayIdx),:)] = sn_get_product_form_chain_params(sn);
    L = L0;
    Z = Z0;
    switch options.config.multiserver
        case {'default','seidmann'}
            % apply seidmann
            L = L ./ repmat(nservers(:),1,C);
            for j=1:size(L,1) % move component of queue j to the first delay
                Z(1,:) = Z(1,:) + L0(j,:) .* (repmat(nservers(j),1,C) - 1)./ repmat(nservers(j),1,C);
            end
        case 'softmin'
            [Q,U,R,T,C,X,lG,totiter] = solver_amvald(sn,Lchain,STchain,Vchain,alpha,Nchain,SCVchain,refstatchain,options);
            return
        otherwise
            %no-op
    end

    switch options.method
        case 'sqni' % square root non-iterative approximation
            if sn.nstations==2 && sum(sn.sched==SchedStrategy.INF)==1
                [Q(sn.nodeToStation(queueIdx),:),U(sn.nodeToStation(queueIdx),:),X] = pfqn_sqni(N,L,Z);
            else
                [Q,U,R,T,C,X,lG,totiter,method] = deal([],[],[],[],[],[],[],0,options.method);
                return
            end
            totiter=1;
        case 'bs'
            [X,Q(sn.nodeToStation(queueIdx),:),U(sn.nodeToStation(queueIdx),:),~,totiter] = pfqn_bs(L,N,Z,options.tol,options.iter_max,[],sn.sched(queueIdx));
        case 'aql'
            if sn_has_multi_server(sn)
                line_error(mfilename,'AQL cannot handle multi-server stations. Try with the ''default'' or ''lin'' methods.');
            end
            [X,Q(sn.nodeToStation(queueIdx),:),U(sn.nodeToStation(queueIdx),:),~,totiter] = pfqn_aql(L,N,Z,options.tol,options.iter_max);
        case 'ab'
            % Akyildiz-Bolch AMVA for multi-server networks
            % Use L0/Z0 (original demands) since ab handles multi-server directly
            % without needing Seidmann transformation
            nDelays = sum(delayIdx);  % Count actual delay stations, not size(Z)
            if nDelays > 0
                nservers_full = [Inf*ones(nDelays,1); nservers];
                D_full = [Z0; L0];
                V_full = [ones(nDelays,C); V(sn.nodeToStation(queueIdx),:)];
                sched_full = [SchedStrategy.INF*ones(nDelays,1); sn.sched(sn.nodeToStation(queueIdx))];
            else
                nservers_full = nservers;
                D_full = L0;
                V_full = V(sn.nodeToStation(queueIdx),:);
                sched_full = sn.sched(sn.nodeToStation(queueIdx));
            end
            [Q_tmp,U_tmp,~,~,X,totiter] = pfqn_ab_amva(D_full,N,V_full,nservers_full,sched_full,false,'ab');
            if nDelays > 0
                Q(sn.nodeToStation(delayIdx),:) = Q_tmp(1:nDelays,:);
                U(sn.nodeToStation(delayIdx),:) = Q_tmp(1:nDelays,:); % delay utilization = queue length
            end
            Q(sn.nodeToStation(queueIdx),:) = Q_tmp(nDelays+1:end,:);
            U(sn.nodeToStation(queueIdx),:) = U_tmp(nDelays+1:end,:);
        case 'schmidt'
            % Schmidt's exact MVA for class-dependent FCFS
            % Use L0/Z0 (original demands) since schmidt handles multi-server directly
            % without needing Seidmann transformation
            nDelays = sum(delayIdx);  % Count actual delay stations, not size(Z)
            if nDelays > 0
                D_full = [Z0; L0];
                S_full = [Inf*ones(nDelays,1); nservers];
                sched_full = [SchedStrategy.INF*ones(nDelays,1); sn.sched(sn.nodeToStation(queueIdx))];
                V_full = [ones(nDelays,C); V(sn.nodeToStation(queueIdx),:)];
            else
                D_full = L0;
                S_full = nservers;
                sched_full = sn.sched(sn.nodeToStation(queueIdx));
                V_full = V(sn.nodeToStation(queueIdx),:);
            end
            [X_tmp,Q_tmp,~,~,~] = pfqn_schmidt(D_full,N,S_full,sched_full,V_full);
            if nDelays > 0
                Q(sn.nodeToStation(delayIdx),:) = Q_tmp(1:nDelays,:);
                U(sn.nodeToStation(delayIdx),:) = Q_tmp(1:nDelays,:); % delay utilization = queue length
            end
            Q(sn.nodeToStation(queueIdx),:) = Q_tmp(nDelays+1:end,:);
            X = X_tmp(1,:);
            % Compute utilization from throughput: U = X * D / nservers
            U(sn.nodeToStation(queueIdx),:) = repmat(X,size(L0,1),1) .* L0 ./ repmat(nservers,1,C);
            totiter = 1;
        case 'schmidt-ext'
            % Extended Schmidt MVA with alpha corrections
            % Use L0/Z0 (original demands) since schmidt-ext handles multi-server directly
            % without needing Seidmann transformation
            nDelays = sum(delayIdx);  % Count actual delay stations, not size(Z)
            if nDelays > 0
                D_full = [Z0; L0];
                S_full = [Inf*ones(nDelays,1); nservers];
                sched_full = [SchedStrategy.INF*ones(nDelays,1); sn.sched(sn.nodeToStation(queueIdx))];
            else
                D_full = L0;
                S_full = nservers;
                sched_full = sn.sched(sn.nodeToStation(queueIdx));
            end
            [X_tmp,Q_tmp,~,~,~] = pfqn_schmidt_ext(D_full,N,S_full,sched_full);
            if nDelays > 0
                Q(sn.nodeToStation(delayIdx),:) = Q_tmp(1:nDelays,:);
                U(sn.nodeToStation(delayIdx),:) = Q_tmp(1:nDelays,:); % delay utilization = queue length
            end
            Q(sn.nodeToStation(queueIdx),:) = Q_tmp(nDelays+1:end,:);
            X = X_tmp(1,:);
            % Compute utilization from throughput: U = X * D / nservers
            U(sn.nodeToStation(queueIdx),:) = repmat(X,size(L0,1),1) .* L0 ./ repmat(nservers,1,C);
            totiter = 1;
        case {'lin','gflin','egflin'}
            % For models with joint dependence (LJD/LJCD), use solver_amvald which handles it
            if sn_has_joint_dependence(sn)
                [Q,U,R,T,C,X,lG,totiter] = solver_amvald(sn,Lchain,STchain,Vchain,alpha,Nchain,SCVchain,refstatchain,options);
                return
            elseif max(nservers)==1
                % remove sources from L
                [Q(sn.nodeToStation(queueIdx),:),U(sn.nodeToStation(queueIdx),:),~,~,X,totiter] = pfqn_linearizermx(lambda,L,N,Z,nservers,sn.sched(sn.nodeToStation(queueIdx | delayIdx)),options.tol,options.iter_max, options.method);
            else
                switch options.config.multiserver
                    case 'conway'
                        [Q(sn.nodeToStation(queueIdx),:),U(sn.nodeToStation(queueIdx),:),~,~,X,totiter] = pfqn_conwayms(L,N,Z,nservers,sn.sched(queueIdx),options.tol,options.iter_max);
                    case 'krzesinski'
                        [Q(sn.nodeToStation(queueIdx),:),U(sn.nodeToStation(queueIdx),:),~,~,X,totiter] = pfqn_linearizermx(lambda,L,N,Z,nservers,sn.sched(sn.nodeToStation(queueIdx | delayIdx)),options.tol,options.iter_max, options.method);
                    case {'default', 'softmin', 'seidmann'}
                        [Q,U,R,T,C,X,lG,totiter] = solver_amvald(sn,Lchain,STchain,Vchain,alpha,Nchain,SCVchain,refstatchain,options);
                        return
                end
            end
        otherwise
            switch options.config.multiserver
                case {'conway','erlang','krzesinski'}
                    options.config.multiserver = 'default';
            end
            [Q,U,R,T,C,X,lG,totiter] = solver_amvald(sn,Lchain,STchain,Vchain,alpha,Nchain,SCVchain,refstatchain,options);
            return
    end

    % compute performance at delay, then unapply seidmann if needed
    for i=1:size(Z0,1)
        % For ab/schmidt methods, Q for delays was already set correctly by the algorithm
        % using original demands Z0. For other methods, use Seidmann-modified Z.
        if strcmp(options.method,'ab') || startsWith(options.method,'schmidt')
            Q(sn.nodeToStation(delayIdx),:) = repmat(X,sum(delayIdx),1) .* Z0;
        else
            Q(sn.nodeToStation(delayIdx),:) = repmat(X,sum(delayIdx),1) .* Z;
        end
        U(sn.nodeToStation(delayIdx),:) = Q(sn.nodeToStation(delayIdx),:);
        switch options.config.multiserver
            case {'default','seidmann'}
                % Skip Seidmann un-apply for ab and schmidt methods as it removes queue length
                if ~strcmp(options.method,'ab') && ~startsWith(options.method,'schmidt')
                    for j=1:size(L,1)
                        if i == 1 && nservers(j)>1
                            % un-apply seidmann from first delay    and move it to
                            % the origin queue
                            jq = find(queueIdx,j);
                            Q(jq,:) = Q(jq,:) + (L0(j,:) .* (repmat(nservers(j),1,C) - 1)./ repmat(nservers(j),1,C)) .* X;
                        end
                    end
                end
        end
    end
    T = V .* repmat(X,M,1);
    R = Q ./ T;
    % For ab/schmidt methods, use original delay demands Z0
    if strcmp(options.method,'ab') || startsWith(options.method,'schmidt')
        C = N ./ X - Z0;
    else
        C = N ./ X - Z;
    end
    lG = NaN;
    if sn_has_class_switching(sn)
        [Q,U,R,T,C,X] = sn_deaggregate_chain_results(sn, Lchain, [], STchain, Vchain, alpha, [], [], R, T, [], X);
    end
else
    switch options.config.multiserver
        case {'conway','erlang','krzesinski'}
            options.config.multiserver = 'default';
    end
    [Q,U,R,T,C,X,lG,totiter] = solver_amvald(sn,Lchain,STchain,Vchain,alpha,Nchain,SCVchain,refstatchain,options);
end
end

