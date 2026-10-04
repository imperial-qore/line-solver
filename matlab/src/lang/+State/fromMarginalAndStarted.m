function space = fromMarginalAndStarted(sn, ind, n, s, options)
% Wrapper: the discipline branches below return early from several places, so
% the synchronous-call (REPLY) counter columns are appended here, once, for
% every exit path. Without them the initial state is narrower than the
% enumerated local space, matchrow fails, and solver_ctmc silently skips its
% unreachable-state pruning -- leaving the enumerated-but-unreachable
% "counter set while every job is here" states as a second absorbing class.
if nargin < 5
    space = sub_fromMarginalAndStarted(sn, ind, n, s);
else
    space = sub_fromMarginalAndStarted(sn, ind, n, s, options);
end
if isfield(sn,'replyblock') && ~isempty(sn.replyblock) && size(sn.replyblock,1) >= ind ...
        && any(sn.replyblock(ind,:) > 0) && ~isempty(space)
    space = [space, zeros(size(space,1), sum(sn.replyblock(ind,:) > 0))];
end
end

function space = sub_fromMarginalAndStarted(sn, ind, n, s, options)
% SPACE = FROMMARGINALANDSTARTED(QN, IND, N, S, OPTIONS)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin<5 %~exist('options','var')
    options.force = true;
end
if isa(sn,'Network')
    sn = sn.getStruct();
end
% generate one initial state such that the marginal queue-lengths are as in vector n
% n(r): number of jobs at the station in class r
% s(r): jobs of class r that are running
R = sn.nclasses;
S = sn.nservers;

% ind: node index
ist = sn.nodeToStation(ind);
%isf = sn.nodeToStateful(ind);
if ~isnan(ist)
    K = zeros(1,R);
    for r=1:R
        if isempty(sn.proc{ist}{r})
            K(r) = 0;
        else
            K(r) = length(sn.proc{ist}{r}{1});
        end
    end
else % node or stateful
    if sn.nodetype(ind) == NodeType.Transition
        K = zeros(1,sn.nodeparam{ind}.nmodes);
        for m=1:sn.nodeparam{ind}.nmodes
            if isempty(sn.nodeparam{ind}.firingproc{m})
                K(m) = 0;
            else
                K(m) = length(sn.nodeparam{ind}.firingproc{m}{1});
            end
        end
    end
end

if sn.isstation(ind) && any(sn.procid(ist,:)==ProcessType.MAP | sn.procid(ist,:)==ProcessType.MMPP2)
    if sn.sched(ist) ~= SchedStrategy.FCFS && sn.nodetype(ind) ~= NodeType.Source
        % Non-FCFS MAP/MMPP2 stations still get a valid marginal state (jobs in first phase) for sim backends; analytical solvers reject via feature set
    end
end

state = [];
space = [];
if sn.isstation(ind) && any(n>sn.classcap(ist,:))
    exceeded = n>sn.classcap(ist,:);
    for r=find(exceeded)
        if ~isempty(sn.proc) && ~isempty(sn.proc{ist}{r}) && any(any(isnan(sn.proc{ist}{r}{1})))
            line_warning(mfilename,'State vector at station %d (n=%s) exceeds the class capacity (classcap=%s). Some service classes are disabled.\n',ist,mat2str(n(ist,:)),mat2str(sn.classcap(ist,:)));
        else
            line_warning(mfilename,'State vector at station %d (n=%s) exceeds the class capacity (classcap=%s).\n',ist,mat2str(n(ist,:)),mat2str(sn.classcap(ist,:)));
        end
    end
    return
end
if sn.isstation(ind) && (sn.nservers(ist)>0 && sum(s) > sn.nservers(ist))
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
                        if ~isempty(sn.proc) && ~isempty(sn.proc{ist}{r}) && any(any(isnan(sn.proc{ist}{r}{1})))
                            init(1) = 0; % class is not processed at this source
                        else
                            % init the job generation
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
                    init = State.spaceClosedSingle(K(r),0);
                    init(1) = n(r);
                    state = State.cartesian(state,init);
                end
                space = State.cartesian(space,state);
            case {SchedStrategy.SIRO, SchedStrategy.LEPT, SchedStrategy.SEPT, SchedStrategy.SRPT, SchedStrategy.SRPTPRIO, SchedStrategy.PSJF, SchedStrategy.FB, SchedStrategy.LRPT, SchedStrategy.POLLING, SchedStrategy.SETF, SchedStrategy.FSP}
                % Unordered buffer + in-server jobs -- see _kb/04-networkstruct.md
                % build list of job classes in the node, with repetition
                if sum(n) <= S(ist)
                    for r=1:R
                        init = State.spaceClosedSingle(K(r),0);
                        init(1) = n(r);
                        state = State.cartesian(state,init);
                    end
                    space = State.cartesian(space,[zeros(size(state,1),R),state]);
                else
                    % si = multichoosecon(n,S(i)); % jobs of class r that are running
                    si = s;
                    mi_buf = repmat(n,size(si,1),1) - si; % jobs of class r in buffer
                    for k=1:size(si,1)
                        % determine number of classes r jobs running in phase j
                        kstate=[];
                        for r=1:R
                            init = State.spaceClosedSingle(K(r),0);
                            init(1) = si(k,r);
                            kstate = State.cartesian(kstate,init);
                        end
                        state = [repmat(mi_buf(k,:),size(kstate,1),1), kstate];
                        space = [space; state];
                    end
                end
            case {SchedStrategy.FCFS, SchedStrategy.HOL,SchedStrategy.FCFSPRIO, SchedStrategy.LCFS,SchedStrategy.LCFSPRIO, SchedStrategy.EDD}
                if sum(n) == 0
                    space = zeros(1,1+max(R,sum(K)));
                    return
                end
                % Ordered buffer + in-server jobs -- see _kb/04-networkstruct.md

                % build list of job classes in the buffer, with repetition
                inbuf = [];
                for r=1:R
                    if n(r)>0
                        inbuf=[inbuf, r*ones(1,n(r)-s(r))];
                    end
                end

                % One buffer ordering is enough: this function returns an INITIAL
                % state, and the trailing unique/flip keeps only the lexicographic
                % maximum, i.e. the descending-sorted buffer. Enumerating every
                % permutation first was factorial in the buffer content.
                if sum(n)>sum(s)
                    mi_buf = sort(inbuf,'descend'); % class of job in buffer position i
                else % set an empty buffer
                    mi_buf = 0;
                end
                % kstate: number of class r jobs running in each service phase
                kstate=[];
                for r=1:R
                    init = State.spaceClosedSingle(K(r),0);
                    init(1) = s(r);
                    kstate = State.cartesian(kstate,init);
                end
                space = [repmat(mi_buf,size(kstate,1),1), kstate];
            case {SchedStrategy.FCFSPR, SchedStrategy.FCFSPI, SchedStrategy.FCFSPRPRIO, SchedStrategy.FCFSPIPRIO, SchedStrategy.LCFSPR, SchedStrategy.LCFSPI, SchedStrategy.LCFSPRPRIO, SchedStrategy.LCFSPIPRIO, SchedStrategy.EDF}
                %% TODO
                if sum(n) == 0
                    % Preempt-resume buffer is (class,phase) pairs, must be even-width -- see _kb/04-networkstruct.md
                    space = zeros(1,2+sum(K));
                    return
                end
                % Ordered buffer + in-server jobs -- see _kb/04-networkstruct.md

                % build list of job classes in the buffer, with repetition
                inbuf = [];
                for r=1:R
                    if n(r)>0
                        inbuf=[inbuf, r*ones(1,n(r)-s(r))];
                    end
                end

                % As in the FCFS branch, only the lexicographic maximum survives the
                % trailing unique/flip, so it is built directly: buffer classes in
                % descending order, each buffered job in its last service phase.
                if sum(n)>sum(s)
                    mi_buf = sort(inbuf,'descend'); % class of job in buffer position i
                else % set an empty buffer
                    mi_buf = 0;
                end
                % kstate: number of class r jobs running in each service phase
                kstate=[];
                for r=1:R
                    init = State.spaceClosedSingle(K(r),0);
                    init(1) = s(r);
                    kstate = State.cartesian(kstate,init);
                end
                % buffer slots interleave (class, phase); an empty slot carries phase 0
                bkstate = zeros(1,numel(mi_buf));
                for jj=1:numel(mi_buf)
                    if mi_buf(jj)>0
                        bkstate(jj) = K(mi_buf(jj));
                    end
                end
                bufstate = zeros(1,2*numel(mi_buf));
                bufstate(1:2:end) = mi_buf;
                bufstate(2:2:end) = bkstate;
                space = State.cartesian(bufstate, kstate);
            case SchedStrategy.PAS
                % PAS/OI: local state is the ordered class-index list; "started" counts s are immaterial -- see _kb/04-networkstruct.md
                W = sn.cap(ist);
                if isinf(W)
                    line_error(mfilename,'PAS stations require finite capacity for state-space generation.');
                end
                if sum(n) == 0
                    space = zeros(1, W);
                elseif sum(n) > W
                    space = zeros(0, W); % infeasible: exceeds total capacity
                else
                    vi = [];
                    for r=1:R
                        if n(r)>0
                            vi=[vi, r*ones(1,n(r))];
                        end
                    end
                    mi = multiset_perms(vi);
                    space = [mi, zeros(size(mi,1), W - size(mi,2))];
                end
            case {SchedStrategy.SJF, SchedStrategy.LJF}
                % in these policies the state space includes continuous
                % random variables for the service times
                % in these policies we only track the jobs in the servers

                for r=1:R
                    init = State.spaceClosedSingle(K(r),0);
                    init(1) = n(r);
                    state = State.cartesian(state,init);
                end
                space = State.cartesian(space,state);
                % this is not casted as an error since this function is
                % IS

                % called to initial models with SJF and LJF
                line_warning(mfilename,'The scheduling policy does not admit a discrete state space.\n');
        end
        % True BAS blocked marker: pin blocked=0 (an initial state has no held job); must still be appended -- see _kb/04-networkstruct.md
        % Col 2*R+1 is shared with the polling controller, so gate on the dedicated
        % sn.isbasblocking field (set for the blocking station under BOTH declaration
        % forms) rather than the station's own drop rule, which fails for a
        % destination-declared BAS. Mirrors State.fromMarginal. See BUG-83.
        if ~isempty(sn.isbasblocking) && numel(sn.isbasblocking) >= ind ...
                && sn.isbasblocking(ind) == 1 && ~isempty(space)
            space = State.cartesian(space, 0);
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
    case NodeType.Join
        if isfield(sn,'isfjaugmented') && sn.isfjaugmented
            % FJ-augmented struct: join state is the per-class buffered-jobs count, deterministic given the marginals
            space = n(:)';
        else
            space = 0;
        end
    case NodeType.Transition
        line_error(mfilename, 'fromMarginalAndStarted cannot be used on Petri net elements');
    case NodeType.Place
        line_error(mfilename, 'fromMarginalAndStarted cannot be used on Petri net elements');
end
space = unique(space,'rows'); % do not comment, required to sort empty state as first
space = space(end:-1:1,:); % this ensures that states where jobs start in phase 1 are first, which is used eg in SSA
end
