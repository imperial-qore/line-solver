function [Q,U,R,T,C,X,lG,totiter] = solver_amva(sn,options)
% [Q,U,R,T,C,X,lG,ITER] = SOLVER_AMVA(SN, OPTIONS)
%
% Copyright (c) 2012-2024, Imperial College London
% All rights reserved.
if nargin < 2
    options = SolverMVA.defaultOptions;
end
%% aggregate chains
[Lchain,STchain,Vchain,alpha,Nchain,SCVchain,refstatchain] = snGetDemandsChain(sn);

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
    case 'default'
        if (sum(Nchain)<=2 || any(Nchain<1))
            options.method = 'qd'; % changing to bs degrades accuracy
        else
            if max(sn.nservers(isfinite(sn.nservers)))==1
                options.method = 'egflin'; % lin seems way worse than aql in test_LQN_8.xml
            else
                options.method = 'lin'; % lin seems way worse than aql in test_LQN_8.xml
            end
        end
end

%% trivial models
if snHasHomogeneousScheduling(sn,SchedStrategy.INF)
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
        if nansum(sn.rates(sn.nodeToStation(sourceIdx),inchain)) > 0
            V(i,c) = 1;
        end
    end
end
Q = zeros(M,C);
U = zeros(M,C);

if snHasProductFormExceptMultiClassHeterExpFCFS(sn) && ~snHasLoadDependence(sn) && (~snHasOpenClasses(sn) || (snHasProductForm(sn) && snHasOpenClasses(sn) && strcmpi(options.method,'lin')))
    % we can use linearizer only if the open model is not heterfcfs as that approximation is not supported, so strict product-form is required
    [lambda,L0,N,Z0,~,nservers,V(sn.nodeToStation(queueIdx|delayIdx),:)] = snGetProductFormChainParams(sn);
    L = L0;
    Z = Z0;
    switch options.config.multiserver
        case {'default','seidmann'}
            % apply seidmann
            L = L ./ repmat(nservers(:),1,C);
            for j=1:size(L,1) % move componnet of queue j to the first delay
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
            if sn.nstations==2
                Nvec = N;
                Nt = sum(sn.njobs);
                X = zeros(1,C);
                for r=1:C
                    Nr = Nvec(r);
                    Lr = L(r);
                    Zr = Z(r);
                    Nvec_1r = Nvec; Nvec_1r(r)=Nvec_1r(r)-1;
                    Br = Nvec./(Z+L+L.*(sum(Nvec)-1-sum(Z.*Nvec_1r./(Z+L+L*(sum(Nvec)-2))))).*Z;
                    Br = Lr*sum(Br(setdiff(1:C,r)));
                    X(r) = (Zr - (Br^2 - 2*Br*Lr*Nt - 2*Br*Zr + Lr^2*Nt^2 + 2*Lr*Nt*Zr - 4*Nr*Lr*Zr + Zr^2)^(1/2) - Br + Lr*Nt)/(2*Lr*Zr);
                    U(queueIdx,r) = X(r)*L(r);
                    Q(queueIdx,r) = N(r)-X(r)*Z(r);
                    totiter=1;
                end
            end
        case 'bs'
            [X,Q(sn.nodeToStation(queueIdx),:),U(sn.nodeToStation(queueIdx),:),~,totiter] = pfqn_bs(L,N,Z,options.tol,options.iter_max,[],sn.schedid(queueIdx));
        case 'aql'
            if snHasMultiServer(sn)
                line_error(mfilename,'AQL cannot handle multi-server stations. Try with the ''default'' or ''lin'' methods.');
            end
            [X,Q(sn.nodeToStation(queueIdx),:),U(sn.nodeToStation(queueIdx),:),~,totiter] = pfqn_aql(L,N,Z,options.tol,options.iter_max);
        case {'lin','gflin','egflin'}
            if max(nservers)==1
                % remove sources from L
                [Q(sn.nodeToStation(queueIdx),:),U(sn.nodeToStation(queueIdx),:),~,~,X,totiter] = pfqn_linearizermx(lambda,L,N,Z,nservers,sn.schedid(sn.nodeToStation(queueIdx | delayIdx)),options.tol,options.iter_max, options.method);
            else
                switch options.config.multiserver
                    case 'conway'
                        [Q(sn.nodeToStation(queueIdx),:),U(sn.nodeToStation(queueIdx),:),~,~,X,totiter] = pfqn_conwayms(L,N,Z,nservers,sn.schedid(queueIdx),options.tol,options.iter_max);
                    case 'erlang'
                        [Q(sn.nodeToStation(queueIdx),:),U(sn.nodeToStation(queueIdx),:),~,~,X,totiter] = pfqn_conwayms_heur(L,N,Z,nservers,sn.schedid(queueIdx),options.tol,options.iter_max);
                    case 'krzesinski'
                        [Q(sn.nodeToStation(queueIdx),:),U(sn.nodeToStation(queueIdx),:),~,~,X,totiter] = pfqn_linearizermx(lambda,L,N,Z,nservers,sn.schedid(sn.nodeToStation(queueIdx | delayIdx)),options.tol,options.iter_max, options.method);
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
        Q(sn.nodeToStation(delayIdx),:) = repmat(X,sum(delayIdx),1) .* Z;
        U(sn.nodeToStation(delayIdx),:) = Q(sn.nodeToStation(delayIdx),:);                     
        switch options.config.multiserver
            case {'default','seidmann'}
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
    T = V .* repmat(X,M,1);
    R = Q ./ T;
    C = N ./ X - Z;
    lG = NaN;
    if snHasClassSwitching(sn)
        [Q,U,R,T,C,X] = snDeaggregateChainResults(sn, Lchain, [], STchain, Vchain, alpha, [], [], R, T, [], X);
    end
else
    switch options.config.multiserver
        case {'conway','erlang','krzesinski'}
            options.config.multiserver = 'default';
    end
    [Q,U,R,T,C,X,lG,totiter] = solver_amvald(sn,Lchain,STchain,Vchain,alpha,Nchain,SCVchain,refstatchain,options);
end
end
    
