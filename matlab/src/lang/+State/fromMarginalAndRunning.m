function space = fromMarginalAndRunning(sn, ind, n, s, options)
% SPACE = FROMMARGINALANDRUNNING(QN, IND, N, S, OPTIONS)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin<5 %~exist('options','var')
    options.force = false;
end
if isa(sn,'Network')
    sn=sn.getStruct();
end
% ind: node index
ist = sn.nodeToStation(ind);
isf = sn.nodeToStateful(ind);

% generate one initial state such that the marginal queue-lengths are as in vector n
% n(r): number of jobs at the station in class r
% s(r): jobs of class r that are running
R = sn.nclasses;
S = sn.nservers;
K = zeros(1,R);
for r=1:R
    if isempty(sn.proc{ist}{r})
        K(r) = 0;
    else
        K(r) = length(sn.proc{ist}{r}{1});
    end
end
state = [];
space = [];
if any(n>sn.classcap(ist,:))
    exceeded = n>sn.classcap(ist,:);
    for r=find(exceeded)
        if ~isempty(sn.proc) && ~isempty(sn.proc{ist}{r}) && any(any(isnan(sn.proc{ist}{r}{1})))
            line_warning(mfilename,sprintf('State vector at station %d (n=%s) exceeds the class capacity (classcap=%s). Some service classes are disabled.\n',ist,mat2str(n),mat2str(sn.classcap(ist,:))));
        else
            line_warning(mfilename,sprintf('State vector at station %d (n=%s) exceeds the class capacity (classcap=%s).\n',ist,mat2str(n),mat2str(sn.classcap(ist,:))));
        end
    end
    return
end
if (sn.nservers(ist)>0 && sum(s) > sn.nservers(ist))
    return
end
% generate local-state space
switch sn.nodetype(ind)
    case {NodeType.Queue, NodeType.Delay, NodeType.Source}
        switch sn.sched(ist)
            case SchedStrategy.EXT
                for r=1:R
                    init = State.spaceClosedSingle(K(r),0);
                    if isinf(sn.njobs(r))
                        if isnan(sn.rates(ist,r))
                            init(1) = 0; % class is not processed at this source
                        else
                            init(1) = 1;
                        end
                    end
                    state = State.cartesian(state,init);
                end
                space = State.cartesian(space,state);
                space = [Inf*ones(size(space,1),1),space];
            case {SchedStrategy.INF, SchedStrategy.PS, SchedStrategy.DPS, SchedStrategy.GPS, SchedStrategy.PSPRIO, SchedStrategy.DPSPRIO, SchedStrategy.GPSPRIO, SchedStrategy.LPS}
                % in these policies we only track the jobs in the servers
                for r=1:R
                    init = State.spaceClosedSingle(K(r),n(r));
                    state = State.cartesian(state,init);
                end
                space = State.cartesian(space,state);
            case {SchedStrategy.SIRO, SchedStrategy.LEPT, SchedStrategy.SEPT, SchedStrategy.SRPT, SchedStrategy.SRPTPRIO}
                % in these policies we track an un-ordered buffer and
                % the jobs in the servers
                % build list of job classes in the node, with repetition
                if sum(n) <= S(ist)
                    for r=1:R
                        init = State.spaceClosedSingle(K(r),n(r));
                        state = State.cartesian(state,init);
                    end
                    space = State.cartesian(space,[zeros(size(state,1),R),state]);
                else
                    %            si = multichoosecon(n,S(i)); % jobs of class r that are running
                    si = s;
                    mi_buf = repmat(n,size(si,1),1) - si; % jobs of class r in buffer
                    for k=1:size(si,1)
                        % determine number of classes r jobs running in phase j
                        kstate=[];
                        for r=1:R
                            init = State.spaceClosedSingle(K(r),si(k,r));
                            kstate = State.cartesian(kstate,init);
                        end
                        state = [repmat(mi_buf(k,:),size(kstate,1),1), kstate];
                        space = [space; state];
                    end
                end
            case {SchedStrategy.FCFS, SchedStrategy.HOL, SchedStrategy.LCFS}
                sizeEstimator = multinomialln(n-s);
                sizeEstimator = round(sizeEstimator/log(10));
                if sizeEstimator > 2
                    if ~isfield(options,'force') || options.force == false
                        line_warning(mfilename,sprintf('State space size is very large: 1e%d states. Stopping execution. Set options.force=true to bypass this control.\n',round(sizeEstimator/log(10))));
                    end
                end
                if sum(n) == 0
                    space = zeros(1,1+sum(K));
                    return
                end
                % in these policies we track an ordered buffer and
                % the jobs in the servers
                
                % build list of job classes in the buffer, with repetition
                inbuf = [];
                for r=1:R
                    if n(r)>s(r)
                        inbuf=[inbuf, r*ones(1,n(r)-s(r))];
                    end
                end
                
                % gen permutation of their positions in the fcfs buffer
                mi = uniqueperms(inbuf); %unique(perms(vi),'rows) is notoriously slow
                if isempty(mi)
                    mi_buf = zeros(1,max(0,sum(n)-S(ist)));
                    state = zeros(1,R);
                    state = State.cartesian(state,[mi_buf,state]);
                else
                    % mi_buf: class of job in buffer position i (0=empty)
                    if sum(n)>sum(s)
                        mi_buf = mi(:,1:(sum(n)-sum(s)));
                    else % set an empty buffer
                        mi_buf = 0;
                    end
                    % si: number of class r jobs that are running
                    si = s;
                    %si = unique(si,'rows');
                    for b=1:size(mi_buf,1)
                        for k=1:size(si,1)
                            % determine number of classs r jobs running in phase
                            % j in server state mi_srv(kjs,:) and build
                            % state
                            kstate=[];
                            for r=1:R
                                % kstate = State.cartesian(kstate,State.spaceClosedSingle(K(r),si(k,r)));
                                init = State.spaceClosedSingle(K(r),si(k,r));
                                kstate = State.cartesian(kstate,init);
                            end
                            state = [state; repmat(mi_buf(b,:),size(kstate,1),1), kstate];
                        end
                    end
                end
                space = state;
            case SchedStrategy.LCFSPR
                %% TODO
                sizeEstimator = multinomialln(n-s);
                sizeEstimator = round(sizeEstimator/log(10));
                if sizeEstimator > 2
                    if ~isfield(options,'force') || options.force == false
                        line_warning(mfilename,sprintf('State space size is very large: 1e%d states. Stopping execution. Set options.force=true to bypass this control.\n',round(sizeEstimator/log(10))));
                    end
                end
                if sum(n) == 0
                    space = zeros(1,1+sum(K));
                    return
                end
                % in these policies we track an ordered buffer and
                % the jobs in the servers
                
                % build list of job classes in the buffer, with repetition
                inbuf = [];
                for r=1:R
                    if n(r)>s(r)
                        inbuf=[inbuf, r*ones(1,n(r)-s(r))];
                    end
                end
                
                % gen permutation of their positions in the fcfs buffer
                mi = uniqueperms(inbuf); %unique(perms(vi),'rows) is notoriously slow
                if isempty(mi)
                    mi_buf = zeros(1,max(0,sum(n)-S(ist)));
                    state = zeros(1,R);
                    state = State.cartesian(state,[mi_buf,state]);
                else
                    % mi_buf: class of job in buffer position i (0=empty)
                    if sum(n)>sum(s)
                        mi_buf = mi(:,1:(sum(n)-sum(s)));
                    else % set an empty buffer
                        mi_buf = 0;
                    end
                    % si: number of class r jobs that are running
                    si = s;
                    %si = unique(si,'rows');
                    for b=1:size(mi_buf,1)
                        for k=1:size(si,1)
                            % determine number of classs r jobs running in phase
                            % j in server state mi_srv(kjs,:) and build
                            % state
                            kstate=[];
                            for r=1:R
                                % kstate = State.cartesian(kstate,State.spaceClosedSingle(K(r),si(k,r)));
                                init = State.spaceClosedSingle(K(r),si(k,r));
                                kstate = State.cartesian(kstate,init);
                            end
                            bkstate = [];
                            for j=mi_buf(b,:) % for each job in the buffer
                                if j>0
                                    bkstate = State.cartesian(bkstate,[1:K(j)]');
                                else
                                    bkstate = 0;
                                end
                            end
                            bufstate_tmp = State.cartesian(mi_buf(b,:), bkstate);
                            % here interleave positions of class and phases in
                            % buf
                            bufstate = zeros(size(bufstate_tmp));
                            bufstate(:,1:2:end)=bufstate_tmp(:,1:size(mi_buf,2));
                            bufstate(:,2:2:end)=bufstate_tmp(:,(size(mi_buf,2)+1):end);
                            state = [state; State.cartesian(bufstate, kstate)];
                        end
                    end
                end
                space = state;
            case {SchedStrategy.SJF, SchedStrategy.LJF}
                % in these policies the state space includes continuous
                % random variables for the service times
                line_warning(mfilename,'The scheduling policy does not admit a discrete state space.\n');
        end
    case NodeType.Cache
        switch sn.sched(ist)
            case SchedStrategy.INF
                % in this policies we only track the jobs in the servers
                for r=1:R
                    init = State.spaceClosedSingle(K(r),n(r));
                    state = State.cartesian(state,init);
                end
                space = State.cartesian(space,state);
        end
end
space = unique(space,'rows'); % do not comment, required to sort empty state as first
space = space(end:-1:1,:); % so that states with jobs in phase 1 comes earlier
end
