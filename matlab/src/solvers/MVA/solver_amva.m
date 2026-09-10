function [Q,U,R,T,C,X,lG,totiter,method,converged] = solver_amva(sn,options)
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

% [] when the handler reports no flag (single-loop handlers, where the count is
% sound); see _kb/06-solver-catalog.md (MVA section) for AMVA convergence
converged = [];

% Apply hard cap on iter_max for stability (Python parity)
max_iter_cap = 10000;
options.iter_max = min(options.iter_max, max_iter_cap);

switch options.method
    case 'amva.qli'
        options.method = 'qli';
    case {'amva.qd', 'amva.qdamva', 'qdamva'}
        options.method = 'qd';
    case 'amva.aql'
        options.method = 'aql';
    case 'amva.qsa'
        options.method = 'qsa';
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
    % 2026-07-29: listValidMethods (SolverMVA.m:63) advertises 'amva.tay' but
    % this switch had no arm for it, so the name survived unchanged and fell
    % through to the qd-equivalent path: on a Delay(Z=1)+FCFS(D=0.5) N=3 model
    % 'tay' gave Qd=1.585456789 in 8 iterations while 'amva.tay' gave
    % Qd=1.499999844 in 31, exactly the 'qd' figures. An advertised method must
    % not silently compute a different one, and the other ten amva.X spellings
    % all alias, so this follows them.
    case 'amva.tay'
        options.method = 'tay';
    case 'amva.scat'
        options.method = 'scat';
    case 'amva.lcp'
        options.method = 'lcp';
    case 'amva.chow'
        options.method = 'chow';
    case 'amva.pamb'
        options.method = 'pamb';
    case 'amva.pami'
        options.method = 'pami';
    case 'amva.pamt'
        options.method = 'pamt';
    case 'amva.clust'
        options.method = 'clust';
    case 'amva.dmlin'
        options.method = 'dmlin';
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

% The closed-population AMVA family (Bard-Schweitzer, SQNI, Tay, SCAT, AQL,
% QSA, Bard LCP, Chow SA, Hsieh-Lam PAM, clustering, Improved Linearizer,
% Akyildiz-Bolch, Schmidt) lives ONLY in the product-form branch below. The
% same predicate the report gates on decides here, so a name the report offers
% is a name that runs and a name it withholds errors rather than falling
% through to solver_amvald and returning the qd-family answer -- or a table of
% zeros -- under a method the caller did not ask for.
[amvaOk, amvaReason] = SolverMVA.supportsClosedPopulation(sn, options.method);
if ~amvaOk
    line_error(mfilename, amvaReason);
end
% The preemptive-resume priority arm of solver_amvald_forward is a single-server
% one; the same predicate the report gates on refuses a multiserver FCFSPRPRIO
% station here, before the forward step where the refusal used to surface.
[prsOk, prsReason] = SolverMVA.supportsPreemptivePriority(sn, options.method);
if ~prsOk
    line_error(mfilename, prsReason);
end

%% trivial models
if sn_has_homogeneous_scheduling(sn,SchedStrategy.INF)
    options.config.multiserver = 'default';
    [Q,U,R,T,C,X,lG,totiter,converged] = solver_amvald(sn,Lchain,STchain,Vchain,alpha,Nchain,SCVchain,refstatchain,options);
    return
end

%% interlocked flow (Franks 1999, Eq. 4.7)
% options.config.interlock arrives CLASS-indexed, and it is translated to the
% chain basis the handlers work in without overwriting it: this function
% re-enters itself on the Conway fallback below, where the class-level matrix
% has to survive. Only SOLVER_AMVALD carries the correction, so an interlocked
% model goes there rather than to the closed-form or linearizermx branches,
% which have no interlock term and would drop it silently.
options.config.interlock_chain = [];
if isfield(options.config,'interlock') && ~isempty(options.config.interlock)
    options.config.interlock_chain = sn_interlock_chain(sn, options.config.interlock);
    if ~isempty(options.config.interlock_chain)
        [Q,U,R,T,C,X,lG,totiter,converged] = solver_amvald(sn,Lchain,STchain,Vchain,alpha,Nchain,SCVchain,refstatchain,options);
        return
    end
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

% ab / schmidt / schmidt-ext ARE the class-dependent FCFS algorithms and live
% only in the product-form branch below, so the het-FCFS exclusion must not
% divert them: doing so would return the qd-family answer under their name.
cond1 = sn_has_product_form_not_het_fcfs(sn) || ...
    (any(strcmpi(options.method,{'ab','schmidt','schmidt-ext'})) && ...
     sn_has_product_form_not_het_fcfs(sn, false));
cond2 = ~sn_has_load_dependence(sn);
cond3 = ~sn_has_open_classes(sn);
if cond1 && cond2 && (cond3 || (sn_has_product_form(sn) && sn_has_open_classes(sn) && strcmpi(options.method,'lin')))
    % we can use linearizer only if the open model is not heterfcfs as that approximation is not supported, so strict product-form is required
    [lambda,L0,N,Z0,~,nservers,V(sn.nodeToStation(queueIdx|delayIdx),:)] = sn_get_product_form_chain_params(sn);
    if size(L0,1) == 0
        % Nothing to correct: with no queueing station the arrival-instant queue
        % length is identically zero, so every AMVA approximation coincides with
        % the exact delay solution and the name a caller passed selects nothing.
        % SOLVER_AMVALD computes it; the kernels below were instead handed a
        % zero-row demand matrix and returned a table of zeros under whatever
        % name was asked for. The C++ port already guards this
        % (solver_mva.h, pf.queue_stations.empty()); the JAR reaches it through
        % SnHasHomogeneousScheduling, which does not reproduce the MATLAB
        % findstring quirk that reduces the trivial-models exit above to
        % nstations == 1.
        [Q,U,R,T,C,X,lG,totiter,converged] = solver_amvald(sn,Lchain,STchain,Vchain,alpha,Nchain,SCVchain,refstatchain,options);
        return
    end
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
            [Q,U,R,T,C,X,lG,totiter,converged] = solver_amvald(sn,Lchain,STchain,Vchain,alpha,Nchain,SCVchain,refstatchain,options);
            return
        otherwise
            %no-op
    end

    % Warm-start queue lengths from a supplied initial solution (options.init_sol),
    % aligned to the algorithm's queueing-station rows. The AMVA fixed point starts
    % from Q0 instead of the default N/M guess, which sharply cuts iterations when
    % the seed is near the fixed point (e.g., the previous SolverLN outer iterate).
    Q0 = [];
    if isfield(options,'init_sol') && ~isempty(options.init_sol)
        Qinit = options.init_sol;
        stationsQ = sort(sn.nodeToStation(queueIdx));
        stationsQ = stationsQ(stationsQ>0);
        if size(Qinit,2) == C && size(Qinit,1) == sn.nstations && numel(stationsQ) == size(L,1)
            Q0 = Qinit(stationsQ, :);
        elseif isequal(size(Qinit), size(L))
            Q0 = Qinit;
        end
        if ~isempty(Q0) && (any(~isfinite(Q0(:))) || any(Q0(:) < 0))
            Q0 = [];   % reject malformed seed
        end
    end

    switch options.method
        case 'sqni' % square root non-iterative approximation
            % pfqn_sqni is a closed form for one queueing station with a delay;
            % SUPPORTSCLOSEDPOPULATION above refuses any other shape (returning
            % empty here silently reported zeros for every metric).
            [Q(sn.nodeToStation(queueIdx),:),U(sn.nodeToStation(queueIdx),:),X] = pfqn_sqni(N,L,Z);
            totiter=1;
        case 'bs'
            [X,Q(sn.nodeToStation(queueIdx),:),U(sn.nodeToStation(queueIdx),:),~,totiter] = pfqn_bs(L,N,Z,options.tol,options.iter_max,Q0,sn.sched(queueIdx));
        case 'lcp'
            % Bard LCP: the Schweitzer proportional term set to zero
            [X,Q(sn.nodeToStation(queueIdx),:),U(sn.nodeToStation(queueIdx),:),~,totiter] = pfqn_lcp(L,N,Z,options.tol,options.iter_max,Q0,sn.sched(queueIdx));
        case 'chow'
            % Chow Second Approximation: theta-terms taken off the LCP solution
            [X,Q(sn.nodeToStation(queueIdx),:),U(sn.nodeToStation(queueIdx),:),~,totiter] = pfqn_chow(L,N,Z,options.tol,options.iter_max,Q0,sn.sched(queueIdx));
        case {'pamb','pami','pamt'}
            % Hsieh-Lam proportional approximations, noniterative
            [X,Q(sn.nodeToStation(queueIdx),:),U(sn.nodeToStation(queueIdx),:)] = pfqn_pam(L,N,Z,options.method);
            totiter=1;
        case 'clust'
            % de Souza e Silva-Lavenberg-Muntz clustering approximation
            [X,Q(sn.nodeToStation(queueIdx),:),U(sn.nodeToStation(queueIdx),:),~,totiter] = pfqn_clust(L,N,Z,[],[],'lin',options.tol,options.iter_max);
        case 'dmlin'
            % de Souza e Silva-Muntz Improved Linearizer
            [Q(sn.nodeToStation(queueIdx),:),U(sn.nodeToStation(queueIdx),:),~,~,X,totiter] = pfqn_dmlin(L,N,Z,sn.sched(sn.nodeToStation(queueIdx)),options.tol,options.iter_max,Q0);
        case 'tay'
            % single-server recursion: SUPPORTSCLOSEDPOPULATION above refuses a
            % multiserver model, as it does for 'aql' and 'qsa'
            [X,Q(sn.nodeToStation(queueIdx),:),U(sn.nodeToStation(queueIdx),:),~,totiter] = pfqn_tay(L,N,Z,options.tol,options.iter_max,Q0);
        case 'scat'
            % Neuse-Chandy SCAT: the Linearizer fixed point with a single Delta
            % refresh. Multiserver stations arrive here already Seidmann-scaled,
            % as they do for 'bs', so no separate guard is needed.
            [Q(sn.nodeToStation(queueIdx),:),U(sn.nodeToStation(queueIdx),:),~,~,X,totiter] = pfqn_scat(L,N,Z,sn.sched(sn.nodeToStation(queueIdx)),options.tol,options.iter_max,Q0);
        case 'aql'
            [X,Q(sn.nodeToStation(queueIdx),:),U(sn.nodeToStation(queueIdx),:),~,totiter] = pfqn_aql(L,N,Z,options.tol,options.iter_max,Q0);
        case 'qsa'
            [Q(sn.nodeToStation(queueIdx),:),U(sn.nodeToStation(queueIdx),:),~,~,X,totiter] = pfqn_qsa(L,N,Z,sn.sched(queueIdx),options.tol,options.iter_max,3,Q0);
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
            % One predicate for the gate and the run, asked about the numbers
            % THIS arm passes: pfqn_schmidt_ext forms its alpha correction from
            % the network with one class-r customer tagged, and a chain holding
            % no customer has none to tag.
            [sxOk, sxReason] = SolverMVA.supportsSchmidtExt(N, ...
                sn.sched(sn.nodeToStation(queueIdx)) == SchedStrategy.FCFS, 'schmidt-ext');
            if ~sxOk
                line_error(mfilename, sxReason);
            end
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
            % For class-dependent models, use solver_amvald which handles them
            if ~isempty(sn.cdscaling) || ~isempty(sn.jdscaling)
                [Q,U,R,T,C,X,lG,totiter,converged] = solver_amvald(sn,Lchain,STchain,Vchain,alpha,Nchain,SCVchain,refstatchain,options);
                return
            elseif max(nservers)==1
                % remove sources from L
                [Q(sn.nodeToStation(queueIdx),:),U(sn.nodeToStation(queueIdx),:),~,~,X,totiter] = pfqn_linearizermx(lambda,L,N,Z,nservers,sn.sched(sn.nodeToStation(queueIdx)),options.tol,options.iter_max, options.method, Q0);
            else
                % 'erlang' names no algorithm of its own: SOLVER_AMVALD_FORWARD
                % implements default/softmin/seidmann/suri and rejects the rest,
                % so both other sites in this file (the non-lin arm and the
                % non-product-form tail) alias it to 'default' before calling in.
                if strcmp(options.config.multiserver,'erlang')
                    options.config.multiserver = 'default';
                end
                switch options.config.multiserver
                    case 'conway'
                        [Q(sn.nodeToStation(queueIdx),:),U(sn.nodeToStation(queueIdx),:),~,~,X,totiter] = pfqn_conwayms(L,N,Z,nservers,sn.sched(queueIdx),options.tol,options.iter_max,Q0);
                    case 'krzesinski'
                        [Q(sn.nodeToStation(queueIdx),:),U(sn.nodeToStation(queueIdx),:),~,~,X,totiter] = pfqn_linearizermx(lambda,L,N,Z,nservers,sn.sched(sn.nodeToStation(queueIdx)),options.tol,options.iter_max, options.method, Q0);
                    case {'default', 'softmin', 'seidmann', 'suri'}
                        [Q,U,R,T,C,X,lG,totiter,converged] = solver_amvald(sn,Lchain,STchain,Vchain,alpha,Nchain,SCVchain,refstatchain,options);
                        if ~isempty(converged) && ~converged && ~amvaConservesPopulation(sn,Q,Nchain)
                            % The amvald fixed point can oscillate without converging on
                            % het-FCFS multiserver layers (e.g. lqn_bpmn layer 10), and its
                            % unconverged iterate violates the closed populations. Conway's
                            % multiserver linearizer is stable there, so re-solve with it.
                            % The population test, not the convergence flag alone, gates the
                            % re-solve: exhausting the iteration budget is routine on a
                            % many-job multichain layer and leaves a feasible iterate, and
                            % pfqn_conwayms is the worse answer there (on lqn_basic layer
                            % P:P1 it returned U=1.0106 against amvald's 0.9959).
                            line_warning(mfilename,'AMVA (%s multiserver) did not converge and violates the closed populations; re-solving with the Conway multiserver linearizer.\n',options.config.multiserver);
                            options.config.multiserver = 'conway';
                            [Q,U,R,T,C,X,lG,totiter2,~,converged] = solver_amva(sn,options);
                            totiter = totiter + totiter2;
                        end
                        return
                    otherwise
                        % Falling through left Q, U and X unassigned and the
                        % caller died on 'Unrecognized function or variable X',
                        % naming nothing. An unhandled rule must name itself.
                        line_error(mfilename,sprintf(['Unrecognized multiserver approximation ''%s''. ' ...
                            'Supported: default, softmin, seidmann, suri, conway, erlang, krzesinski.'], ...
                            options.config.multiserver));
                end
            end
        otherwise
            switch options.config.multiserver
                case {'conway','erlang','krzesinski'}
                    options.config.multiserver = 'default';
            end
            [Q,U,R,T,C,X,lG,totiter,converged] = solver_amvald(sn,Lchain,STchain,Vchain,alpha,Nchain,SCVchain,refstatchain,options);
            return
    end

    % compute performance at delay, then unapply seidmann if needed
    % The delay is charged the ORIGINAL think time Z0 for every method. Seidmann
    % folds L(m-1)/m of each multiserver station into Z, but that population is
    % in service at the station and is given back to it below; charging Z here
    % as well counted it twice, so sum(Q) exceeded N.
    Q(sn.nodeToStation(delayIdx),:) = repmat(X,sum(delayIdx),1) .* Z0;
    U(sn.nodeToStation(delayIdx),:) = Q(sn.nodeToStation(delayIdx),:);
    switch options.config.multiserver
        case {'default','seidmann'}
            % Skip Seidmann un-apply for ab and schmidt methods: they were handed
            % the original demands, so the transform was never applied to them
            if ~strcmp(options.method,'ab') && ~startsWith(options.method,'schmidt')
                % station of algorithm row j, in the order the solves above wrote
                % Q(sn.nodeToStation(queueIdx),:); find(queueIdx,j) returned the
                % first j queues instead, spraying station j's term over all of them
                stationsQ = sn.nodeToStation(queueIdx);
                for j=1:size(L,1)
                    if nservers(j)>1
                        % move the folded delay component back to its own station
                        jq = stationsQ(j);
                        Q(jq,:) = Q(jq,:) + (L0(j,:) .* (repmat(nservers(j),1,C) - 1)./ repmat(nservers(j),1,C)) .* X;
                    end
                end
            end
    end
    T = V .* repmat(X,M,1);
    R = Q ./ T;
    % Cycle time excludes the think time actually spent at the delays, which is
    % Z0 summed over them: with the Seidmann Z it disagreed with sum(R.*V)
    C = N ./ X - sum(Z0,1);
    lG = NaN;
    if sn_has_class_switching(sn)
        [Q,U,R,T,C,X] = sn_deaggregate_chain_results(sn, Lchain, [], STchain, Vchain, alpha, [], [], R, T, [], X);
    end
else
    % Nothing of the closed-population family can reach here: it is refused
    % above by SUPPORTSCLOSEDPOPULATION, the one predicate the report also
    % gates on. Keeping a second copy of that rule here is what let the two
    % drift, so that 'bs', 'sqni', 'ab' and the two Schmidt arms fell through
    % to solver_amvald and returned the qd-family answer under their name while
    % 'aql', 'qsa', 'tay' and the Chapter-2 survey algorithms errored.
    switch options.config.multiserver
        case {'conway','erlang','krzesinski'}
            options.config.multiserver = 'default';
    end
    [Q,U,R,T,C,X,lG,totiter,converged] = solver_amvald(sn,Lchain,STchain,Vchain,alpha,Nchain,SCVchain,refstatchain,options);
end
end

function tf = amvaConservesPopulation(sn,Q,Nchain)
% True when the class-level queue lengths Q still add up, chain by chain, to
% the closed populations Nchain. Population is conserved per CHAIN, not per
% class, since class switching moves jobs between the classes of a chain.
tf = true;
if isempty(Q) || isempty(sn.chains)
    return
end
for c=1:size(sn.chains,1)
    if c > numel(Nchain) || ~isfinite(Nchain(c)) || Nchain(c) <= 0
        continue
    end
    inchain = sn.chains(c,:) > 0;
    qc = sum(sum(Q(:,inchain),'omitnan'),'omitnan');
    if ~isfinite(qc) || abs(qc - Nchain(c)) > 1e-3 * Nchain(c)
        tf = false;
        return
    end
end
end

