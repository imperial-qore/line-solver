function [outglspace, outrate, outprob, outcomp] =  afterGlobalEvent(sn, ind, glspace, glevent, isSimulation)
% [OUTGLSPACE, OUTRATE, OUTPROB, OUTCOMP] =  AFTERGLOBALEVENT(QN, IND, GLSPACE, GLEVENT, ISSIMULATION)
%
% OUTCOMP(io) is true when outcome io is a firing completion (a D1 firing of
% the active mode, which applies the PRE/POST place updates), and false when
% it is not (e.g. every ENABLE outcome). Callers must not re-derive this from
% the place markings: a transition whose firing outcome returns exactly what
% its enabling condition consumed leaves every marking invariant, yet still
% completes at a nonzero rate.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

outspace = []; % state space of ind node
outrate = []; % rate of synchronization
outprob = []; % probability of synchronization
outcomp = []; % true if the outcome is a firing completion (see above)
outglspace = glspace; % new global state after synchronization

R = sn.nclasses;
%phasessz = sn.phasessz;
%phaseshift = sn.phaseshift;
if sn.nodetype(ind) == NodeType.Transition % same isa(glevent,'ModeEvent')
    % ModeEvent
    isf = sn.nodeToStateful(ind);
    inspace = glspace{isf};
    V = sum(sn.nvars(ind,:));
    % in this case service is always immediate so sum(K)=1
    space_var = inspace(:,(end-V+1):end); % local state variables
    if sn.nodetype(ind) == NodeType.Transition
        fK = sn.nodeparam{ind}.firingphases;
        nmodes = sn.nodeparam{ind}.nmodes;
        % Handle NaN firingphases (non-phase-type distributions like Pareto)
        if any(isnan(fK))
            fK = zeros(1, nmodes);
            for m = 1:nmodes
                if iscell(sn.nodeparam{ind}.firingproc) && ~isempty(sn.nodeparam{ind}.firingproc{m})
                    fK(m) = size(sn.nodeparam{ind}.firingproc{m}{1}, 1);
                else
                    fK(m) = 1;
                end
            end
        end
        fKs = [0,cumsum(fK,2)];
        mode = glevent.active{1}.mode;
        % Transition state format: [idle_counts(nmodes), phase_counts(sum(fK)), fired_counts(nmodes)]
        space_buf = inspace(:,1:nmodes); % idle servers count put in buf
        space_srv = inspace(:,(nmodes+1):(nmodes+sum(fK))); % enabled servers' phases
        % Handle both state formats: with and without fired component
        expected_len_with_fired = 2*nmodes + sum(fK);
        expected_len_without_fired = nmodes + sum(fK);
        if size(inspace, 2) >= expected_len_with_fired
            space_fired = inspace(:,(nmodes+sum(fK)+1):(2*nmodes+sum(fK))); % servers that just fired
        elseif size(inspace, 2) == expected_len_without_fired
            % Legacy format without fired component - initialize to zeros
            space_fired = zeros(size(inspace,1), nmodes);
        else
            line_error(mfilename, 'Unexpected state vector length for Transition node');
        end

        switch glevent.active{1}.event
            case EventType.ENABLE
                enabling_m = sn.nodeparam{ind}.enabling{mode}; % enabling requirement for mode m
                inhibiting_m = sn.nodeparam{ind}.inhibiting{mode}; % inhibitor thresholds (Inf = no inhibition)
                ep_space = zeros(sn.nnodes,R);
                for j=1:length(glevent.passive)
                    ep_linidx = glevent.passive{j}.node; % linear index into (nnodes x nclasses) matrix
                    % Decode linear index to (node, class) - the passive node is a linear index from find() on enabling matrix
                    [ep_ind, ~] = ind2sub([sn.nnodes, R], ep_linidx);
                    if ep_ind > length(sn.nodeToStateful) || isnan(sn.nodeToStateful(ep_ind)) || sn.nodeToStateful(ep_ind) <= 0
                        continue; % skip non-stateful nodes
                    end
                    ep_isf = sn.nodeToStateful(ep_ind);
                    K = ones(1,R);
                    Ks = [0,cumsum(K,2)];
                    ep_space_buf = glspace{ep_isf};
                    ep_space_srv = zeros(1,R);
                    ep_space_var = [];
                    [~,ep_space(ep_ind,1:R)] = State.toMarginalAggr(sn,ep_ind, glspace{ep_isf},K,Ks,ep_space_buf,ep_space_srv,ep_space_var);
                end
                m  = mode;
                % Inhibitor arcs: mode is disabled while any inhibited input
                % place has reached its threshold (Inf default => never true).
                inhibited = any(ep_space(:) >= inhibiting_m(:));
                if any(ep_space(:) < enabling_m(:)) || inhibited
                    % check if any servers need to be disabled
                    new_space_buf = space_buf;
                    new_space_srv = space_srv;
                    % Cap nmodeservers to MaxInt for state (SSA needs finite states)
                    nmodeservers_m = sn.nodeparam{ind}.nmodeservers(m);
                    if isinf(nmodeservers_m)
                        nmodeservers_m = GlobalConstants.MaxInt();
                    end
                    new_space_buf(m) = nmodeservers_m;
                    new_space_srv((fKs(m)+1):(fKs(m)+fK(m))) = 0;
                    new_state = [new_space_buf, new_space_srv, space_fired, space_var];
                    old_state = [space_buf, space_srv, space_fired, space_var];
                    if ~isequal(new_state, old_state)
                        % State changes - return the disable action
                        outspace = new_state;
                        outrate = GlobalConstants.Immediate;
                        outprob = 1.0;
                    else
                        % Already disabled, no state change
                        outspace = [];
                        outrate = [];
                        outprob = 1.0;
                    end
                else
                    % find enabling degree, i.e., running servers
                    en_degree_m = 1;
                    while all(ep_space >= en_degree_m * enabling_m)
                        en_degree_m = en_degree_m + 1;
                    end
                    % Cap nmodeservers to MaxInt for enabling degree calculation
                    nmodeservers_m = sn.nodeparam{ind}.nmodeservers(m);
                    if isinf(nmodeservers_m)
                        nmodeservers_m = GlobalConstants.MaxInt();
                    end
                    en_degree_m = min(en_degree_m - 1, nmodeservers_m);
                    running_m = sum(space_srv((fKs(m)+1):(fKs(m)+fK(m))));
                    if running_m == en_degree_m
                        % all running as expected, do nothing
                        outspace = [];
                        outrate = [];
                        outprob = 1.0;
                    elseif running_m < en_degree_m
                        % fewer expected, start en_degree_m - running_m
                        space_buf(m) = space_buf(m) - (en_degree_m - running_m);
                        pentry = sn.nodeparam{ind}.firingpie{mode};
                        nadd = en_degree_m - running_m;
                        combs = sortrows(multichoose(fK(m), nadd),'descend');

                        for i = 1:size(combs,1)
                            comb = combs(i,:);
                            space_srv_k = space_srv;
                            for k = 1:fK(m)
                                space_srv_k(fKs(m)+k) = space_srv_k(fKs(m)+k) + comb(k);
                            end
                            outspace = [outspace; space_buf, space_srv_k, space_fired, space_var]; %#ok<AGROW>
                            outrate = [outrate; GlobalConstants.Immediate]; %#ok<AGROW>

                            % compute multinomial coefficient
                            % multinomial probability: n! / (k1! * k2! * ...) * p1^k1 * p2^k2 * ...
                            logprob = factln(nadd);
                            for k = 1:fK(m)
                                if pentry(k)>0
                                    logprob = logprob + comb(k)*log(pentry(k)) - factln(comb(k));                                
                                elseif (pentry(k)==0 && comb(k)==0)
                                    continue
                                else % (pentry(k)==0 && comb(k)>0)
                                    logprob = -Inf; % set zero outprob
                                end
                            end
                            outprob = [outprob; exp(logprob)]; %#ok<AGROW>
                        end
                        %outprob = outprob / sum(outprob);
                    else %running_m > en_degree_m
                        % more than expected, stop running_m - en_degree_m
                        % servers chosen uniformly at random across phases.
                        % We enumerate how many servers to stop per phase (comb),
                        % where 0 <= comb(k) <= srv_vec(k) and sum(comb) = ndiff.
                        % Probability of each comb follows the hypergeometric
                        % mixture: prod_k C(srv_vec(k), comb(k)) / C(running_m, ndiff).
                        ndiff = running_m - en_degree_m; % number of servers to stop
                        srv_vec = space_srv((fKs(m)+1):(fKs(m)+fK(m)));
                        n = numel(srv_vec);
                        % Enumerate all (comb) of length n with sum == ndiff and
                        % comb(k) <= srv_vec(k). multichoose gives sum==ndiff with
                        % entries >= 0 but without the per-phase cap, so filter.
                        all_combs = multichoose(n, ndiff);
                        valid = all(bsxfun(@le, all_combs, srv_vec(:)'), 2);
                        all_combs = all_combs(valid, :);
                        % Compute weights via multivariate hypergeometric
                        logW = zeros(size(all_combs,1),1);
                        for i = 1:size(all_combs,1)
                            for k = 1:n
                                logW(i) = logW(i) + factln(srv_vec(k)) ...
                                    - factln(all_combs(i,k)) ...
                                    - factln(srv_vec(k) - all_combs(i,k));
                            end
                        end
                        W = exp(logW - max(logW));
                        W = W / sum(W);
                        for i = 1:size(all_combs,1)
                            comb = all_combs(i,:);
                            space_srv_reduced = space_srv;
                            space_srv_reduced((fKs(m)+1):(fKs(m)+fK(m))) = srv_vec - comb;
                            space_buf_reduced = space_buf;
                            space_buf_reduced(m) = space_buf_reduced(m) + ndiff; % return stopped servers to idle pool
                            outspace = [outspace; space_buf_reduced, space_srv_reduced, space_fired, space_var]; %#ok<AGROW>
                            outrate = [outrate; GlobalConstants.Immediate]; %#ok<AGROW>
                            outprob = [outprob; W(i)]; %#ok<AGROW>
                        end
                    end
                end
                % update new state space
                outglspace{isf} = outspace;
            case EventType.FIRE
                %% Update transition servers
                fK = sn.nodeparam{ind}.firingphases;
                % Handle NaN firingphases (non-phase-type distributions like Pareto)
                if any(isnan(fK))
                    fK = zeros(1, nmodes);
                    for m = 1:nmodes
                        if iscell(sn.nodeparam{ind}.firingproc) && ~isempty(sn.nodeparam{ind}.firingproc{m})
                            fK(m) = size(sn.nodeparam{ind}.firingproc{m}{1}, 1);
                        else
                            fK(m) = 1;
                        end
                    end
                end
                fKs = [0,cumsum(fK,2)];
                mode = glevent.active{1}.mode;
                % find transition server counts
                [~,nim,~,kim] = State.toMarginal(sn,ind,inspace,fK,fKs,space_buf,space_srv,space_var);
                % state of transition mode
                enabling_m = sn.nodeparam{ind}.enabling{mode}; % enabling requirement for mode m
                inhibiting_m = sn.nodeparam{ind}.inhibiting{mode}; % inhibitor thresholds (Inf = no inhibition)
                firing_m = sn.nodeparam{ind}.firing{mode}; % firing requirement for mode m
                % find enabling degree, i.e., running servers
                ep_space = zeros(sn.nnodes,R);
                for j=1:length(glevent.passive)
                    ep_linidx = glevent.passive{j}.node; % linear index into (nnodes x nclasses) matrix
                    % Decode linear index to (node, class) - the passive node is a linear index from find() on enabling/firing matrix
                    [ep_ind, ~] = ind2sub([sn.nnodes, R], ep_linidx);
                    if ep_ind > length(sn.nodeToStateful) || isnan(sn.nodeToStateful(ep_ind)) || sn.nodeToStateful(ep_ind) <= 0
                        continue; % skip non-stateful nodes
                    end
                    ep_isf = sn.nodeToStateful(ep_ind);
                    K = ones(1,R);
                    Ks = [0,cumsum(K,2)];
                    ep_space_buf = outglspace{ep_isf}; % this equals glspace at this point
                    ep_space_srv = zeros(1,R);
                    ep_space_var = [];
                    [~,ep_space(ep_ind,1:R)] = State.toMarginalAggr(sn,ep_ind, glspace{ep_isf},K,Ks,ep_space_buf,ep_space_srv,ep_space_var);
                end
                if any(ep_space(:) >= inhibiting_m(:))
                    % Inhibitor arc active: mode cannot fire.
                    en_degree_m = 0;
                else
                    en_degree_m = 1;
                    while all(ep_space >= en_degree_m * enabling_m)
                        en_degree_m = en_degree_m + 1;
                    end
                    mark_degree_m = en_degree_m - 1; % concurrent firings the marking alone supports
                    %en_degree_m = min(en_degree_m - 1, sn.nodeparam{ind}.nmodeservers(m));
                    en_degree_m = min(mark_degree_m, nim(mode)); % the actual enabling degree depends on servers actually in execution, though this should coincide
                end

                % Get D0 and D1 matrices for phase transitions and completions
                D0 = sn.nodeparam{ind}.firingproc{mode}{1}; % Internal phase transitions
                D1 = sn.nodeparam{ind}.firingproc{mode}{2}; % Completions/firings

                % Track whether this is a completion (with place updates) or just phase transition
                completion_start_idx = []; % indices in outspace that are completions

                % see _kb/04-networkstruct.md (SPN/GSPN state generation) for rationale

                % Add completions (D1) - these update places
                % Track which outcomes are completions (vs phase transitions)
                is_completion = false(size(outspace,1), 1); % mark existing as phase transitions
                % see _kb/04-networkstruct.md (SPN/GSPN state generation) for rationale
                isImmediateMode = false;
                if isfield(sn.nodeparam{ind}, 'timing') && ~isempty(sn.nodeparam{ind}.timing) ...
                        && mode <= numel(sn.nodeparam{ind}.timing)
                    isImmediateMode = sn.nodeparam{ind}.timing(mode) == TimingStrategy.IMMEDIATE;
                end
                immediateWeight = 1.0;
                if isImmediateMode && isfield(sn.nodeparam{ind}, 'fireweight') ...
                        && ~isempty(sn.nodeparam{ind}.fireweight) ...
                        && mode <= numel(sn.nodeparam{ind}.fireweight)
                    immediateWeight = sn.nodeparam{ind}.fireweight(mode);
                end

                % see _kb/04-networkstruct.md (SPN/GSPN state generation) for rationale
                if isImmediateMode
                    imm_nservers = sn.nodeparam{ind}.nmodeservers(mode);
                    if isinf(imm_nservers)
                        imm_nservers = GlobalConstants.MaxInt();
                    end
                    imm_servers_m = min(mark_degree_m, imm_nservers);
                else
                    imm_servers_m = 0;
                end
                for k = 1:fK(mode) % source phase
                    if isImmediateMode
                        fires_k = (k == 1) && imm_servers_m >= 1;
                    else
                        fires_k = kim(:,mode,k) > 0 && en_degree_m >= 1; % servers in phase k, and the place can still pay the PRE
                    end
                    if fires_k
                        if isImmediateMode
                            rate_kd = GlobalConstants.Immediate * immediateWeight * imm_servers_m;
                        else
                            rate_kd = sum(D1(k,:)) * kim(:,mode,k);
                            % Marking-dependent firing-rate multiplier g_mode(marking):
                            % ep_space is the node-indexed input-place marking, the
                            % same object the JSON tabulation was built over. Exact
                            % because CTMC evaluates it per enumerated state.
                            if isfield(sn.nodeparam{ind}, 'firingdep') ...
                                    && numel(sn.nodeparam{ind}.firingdep) >= mode ...
                                    && ~isempty(sn.nodeparam{ind}.firingdep{mode})
                                rate_kd = rate_kd * double(sn.nodeparam{ind}.firingdep{mode}(ep_space));
                            end
                        end
                        if rate_kd > 0
                            space_buf_kd = space_buf;
                            space_srv_kd = space_srv;
                            if kim(:,mode,k) > 0
                                % Decrease by one the firing server
                                space_srv_kd(fKs(mode)+k) = space_srv_kd(fKs(mode)+k) - 1;
                                % Move the server back to the disabled pool
                                space_buf_kd(mode) = space_buf_kd(mode) + 1;
                            end
                            space_fired_kd = space_fired;
                            if isSimulation
                                space_fired_kd(mode) = space_fired_kd(mode) + 1; % increment fired count only for simulation
                            end
                            outspace = [outspace; space_buf_kd, space_srv_kd, space_fired_kd, space_var]; %#ok<AGROW>
                            outrate = [outrate; rate_kd]; %#ok<AGROW>
                            outprob = [outprob; 1.0]; %#ok<AGROW>
                            is_completion = [is_completion; true]; %#ok<AGROW>
                        end
                    end
                end
                outglspace{isf} = outspace;
                outcomp = is_completion;

                % For simulation: select outcome first, then apply PRE/POST only if completion
                if isSimulation && size(outspace,1) > 1 && ~isempty(outprob)
                    % Use effective rate (rate * prob) for selection
                    eff_rate = outrate .* outprob;
                    tot_rate = sum(eff_rate);
                    if tot_rate > 0
                        cum_rate = cumsum(eff_rate) / tot_rate;
                        firing_ctr = 1 + max([0,find( rand > cum_rate' )]);
                        selected_is_completion = is_completion(firing_ctr);
                        outcomp = selected_is_completion;
                        outspace = outspace(firing_ctr,:);
                        outrate = sum(eff_rate);  % return effective rate
                        outprob = 1.0;  % probability already incorporated
                        outglspace{isf} = outspace;

                        % Only process PRE/POST if this is a completion
                        if selected_is_completion
                            %% Process ID_PRE events, i.e., consume from all input places
                            for j=1:length(glevent.passive)
                                if glevent.passive{j}.event == EventType.PRE
                                    ep_linidx = glevent.passive{j}.node;
                                    [ep_ind, ep_class] = ind2sub([sn.nnodes, R], ep_linidx);
                                    if ep_ind > length(sn.nodeToStateful) || isnan(sn.nodeToStateful(ep_ind)) || sn.nodeToStateful(ep_ind) <= 0
                                        continue;
                                    end
                                    ep_isf = sn.nodeToStateful(ep_ind);
                                    ep_space_buf = outglspace{ep_isf};
                                    ep_ist = sn.nodeToStation(ep_ind);
                                    % see _kb/04-networkstruct.md (SPN/GSPN state generation) for rationale
                                    consume_count = glevent.passive{j}.weight;
                                    % Check scheduling strategy and state format to determine handling
                                    % Places with INF scheduling use queue-based state format
                                    is_queue_based = sn.nodetype(ep_ind) == NodeType.Place && ...
                                        length(ep_space_buf) ~= R && length(ep_space_buf) ~= 2*R;
                                    if is_queue_based || sn.sched(ep_ist) == SchedStrategy.FCFS
                                        % FCFS queue-based handling
                                        idx = find(ep_space_buf == ep_class);
                                        if length(idx) >= consume_count
                                            ep_space_buf(idx(end - consume_count + 1:end)) = [];
                                        end
                                    elseif sn.sched(ep_ist) == SchedStrategy.LCFS
                                        idx = find(ep_space_buf == ep_class);
                                        if length(idx) >= consume_count
                                            ep_space_buf(idx(1:consume_count)) = [];
                                        end
                                    elseif sn.sched(ep_ist) == SchedStrategy.SIRO
                                        ep_space_buf(ep_class) = ep_space_buf(ep_class) - consume_count;
                                    elseif sn.nodetype(ep_ind) == NodeType.Place
                                        % Place count-based state format: [buffer(R), server_phases(sum(K))]
                                        state_len = length(ep_space_buf);
                                        if state_len > R
                                            buf_jobs = ep_space_buf(ep_class);
                                            srv_jobs = ep_space_buf(R + ep_class);
                                            total_jobs = buf_jobs + srv_jobs;
                                            remaining = total_jobs - consume_count;
                                            ep_space_buf(ep_class) = remaining;
                                            ep_space_buf(R + ep_class) = 0;
                                        else
                                            ep_space_buf(ep_class) = ep_space_buf(ep_class) - consume_count;
                                        end
                                    else
                                        line_error(mfilename,sprintf('Scheduling strategy %s is unsupported at places.', SchedStrategy.toText(sn.sched(ep_ist))));
                                    end
                                    outglspace{ep_isf} = ep_space_buf;
                                end
                            end

                            %% Process ID_POST events, i.e., produce to all output places
                            for j=1:length(glevent.passive)
                                if glevent.passive{j}.event == EventType.POST
                                    fp_linidx = glevent.passive{j}.node;
                                    [fp_ind, fp_class] = ind2sub([sn.nnodes, R], fp_linidx);
                                    if fp_ind > length(sn.nodeToStateful) || isnan(sn.nodeToStateful(fp_ind)) || sn.nodeToStateful(fp_ind) <= 0
                                        continue;
                                    end
                                    fp_isf = sn.nodeToStateful(fp_ind);
                                    fp_ist = sn.nodeToStation(fp_ind);
                                    fp_space_buf = outglspace{fp_isf};
                                    produce_count = glevent.passive{j}.weight;
                                    % Check scheduling strategy and state format to determine handling
                                    % Use original glspace length for detection (before PRE modified it)
                                    orig_len = length(glspace{fp_isf});
                                    is_queue_based = sn.nodetype(fp_ind) == NodeType.Place && ...
                                        orig_len ~= R && orig_len ~= 2*R;
                                    if is_queue_based || sn.sched(fp_ist) == SchedStrategy.FCFS || sn.sched(fp_ist) == SchedStrategy.LCFS
                                        % Queue-based handling - prepend to buffer
                                        fp_space_buf = [repmat(fp_class, 1, produce_count), fp_space_buf];
                                    elseif sn.sched(fp_ist) == SchedStrategy.SIRO
                                        fp_space_buf(fp_class) = fp_space_buf(fp_class) + produce_count;
                                    elseif sn.nodetype(fp_ind) == NodeType.Place
                                        fp_space_buf(fp_class) = fp_space_buf(fp_class) + produce_count;
                                    else
                                        line_error(mfilename,sprintf('Scheduling strategy %s is unsupported at places.', SchedStrategy.toText(sn.sched(fp_ist))));
                                    end
                                    outglspace{fp_isf} = fp_space_buf;
                                end
                            end
                        end
                    else
                        % All effective rates are zero, pick first option with non-zero prob if any
                        valid_idx = find(outprob > 0, 1);
                        if isempty(valid_idx)
                            valid_idx = 1;
                        end
                        outspace = outspace(valid_idx,:);
                        outrate = 0;
                        outprob = outprob(valid_idx,:);
                        outglspace{isf} = outspace;
                    end
                elseif ~isSimulation
                    %% Non-simulation mode: process all completions for PRE/POST
                    % For state space analysis, we apply PRE/POST for completion outcomes
                    for j=1:length(glevent.passive)
                        if glevent.passive{j}.event == EventType.PRE
                            ep_linidx = glevent.passive{j}.node;
                            [ep_ind, ep_class] = ind2sub([sn.nnodes, R], ep_linidx);
                            if ep_ind > length(sn.nodeToStateful) || isnan(sn.nodeToStateful(ep_ind)) || sn.nodeToStateful(ep_ind) <= 0
                                continue;
                            end
                            ep_isf = sn.nodeToStateful(ep_ind);
                            ep_space_buf = outglspace{ep_isf};
                            ep_ist = sn.nodeToStation(ep_ind);
                            consume_count = glevent.passive{j}.weight; % one completion = one firing (see PRE note above)
                            % Check scheduling strategy and state format to determine handling
                            is_queue_based = sn.nodetype(ep_ind) == NodeType.Place && ...
                                length(ep_space_buf) ~= R && length(ep_space_buf) ~= 2*R;
                            if is_queue_based || sn.sched(ep_ist) == SchedStrategy.FCFS
                                % FCFS queue-based handling
                                idx = find(ep_space_buf == ep_class);
                                if length(idx) >= consume_count
                                    ep_space_buf(idx(end - consume_count + 1:end)) = [];
                                end
                            elseif sn.sched(ep_ist) == SchedStrategy.LCFS
                                idx = find(ep_space_buf == ep_class);
                                if length(idx) >= consume_count
                                    ep_space_buf(idx(1:consume_count)) = [];
                                end
                            elseif sn.sched(ep_ist) == SchedStrategy.SIRO
                                ep_space_buf(ep_class) = ep_space_buf(ep_class) - consume_count;
                            elseif sn.nodetype(ep_ind) == NodeType.Place
                                % Place count-based state format
                                state_len = length(ep_space_buf);
                                if state_len > R
                                    buf_jobs = ep_space_buf(ep_class);
                                    srv_jobs = ep_space_buf(R + ep_class);
                                    total_jobs = buf_jobs + srv_jobs;
                                    remaining = total_jobs - consume_count;
                                    ep_space_buf(ep_class) = remaining;
                                    ep_space_buf(R + ep_class) = 0;
                                else
                                    ep_space_buf(ep_class) = ep_space_buf(ep_class) - consume_count;
                                end
                            else
                                line_error(mfilename,sprintf('Scheduling strategy %s is unsupported at places.', SchedStrategy.toText(sn.sched(ep_ist))));
                            end
                            outglspace{ep_isf} = ep_space_buf;
                        end
                    end

                    for j=1:length(glevent.passive)
                        if glevent.passive{j}.event == EventType.POST
                            fp_linidx = glevent.passive{j}.node;
                            [fp_ind, fp_class] = ind2sub([sn.nnodes, R], fp_linidx);
                            if fp_ind > length(sn.nodeToStateful) || isnan(sn.nodeToStateful(fp_ind)) || sn.nodeToStateful(fp_ind) <= 0
                                continue;
                            end
                            fp_isf = sn.nodeToStateful(fp_ind);
                            fp_ist = sn.nodeToStation(fp_ind);
                            fp_space_buf = outglspace{fp_isf};
                            produce_count = glevent.passive{j}.weight;
                            % Check scheduling strategy and state format to determine handling
                            % Use original glspace length for detection (before PRE modified it)
                            orig_len = length(glspace{fp_isf});
                            is_queue_based = sn.nodetype(fp_ind) == NodeType.Place && ...
                                orig_len ~= R && orig_len ~= 2*R;
                            if is_queue_based || sn.sched(fp_ist) == SchedStrategy.FCFS || sn.sched(fp_ist) == SchedStrategy.LCFS
                                fp_space_buf = [repmat(fp_class, 1, produce_count), fp_space_buf];
                            elseif sn.sched(fp_ist) == SchedStrategy.SIRO
                                fp_space_buf(fp_class) = fp_space_buf(fp_class) + produce_count;
                            elseif sn.nodetype(fp_ind) == NodeType.Place
                                fp_space_buf(fp_class) = fp_space_buf(fp_class) + produce_count;
                            else
                                line_error(mfilename,sprintf('Scheduling strategy %s is unsupported at places.', SchedStrategy.toText(sn.sched(fp_ist))));
                            end
                            outglspace{fp_isf} = fp_space_buf;
                        end
                    end
                else
                    % Single outcome case in simulation mode
                    if ~isempty(is_completion) && is_completion(1)
                        %% Process ID_PRE events
                        for j=1:length(glevent.passive)
                            if glevent.passive{j}.event == EventType.PRE
                                ep_linidx = glevent.passive{j}.node;
                                [ep_ind, ep_class] = ind2sub([sn.nnodes, R], ep_linidx);
                                if ep_ind > length(sn.nodeToStateful) || isnan(sn.nodeToStateful(ep_ind)) || sn.nodeToStateful(ep_ind) <= 0
                                    continue;
                                end
                                ep_isf = sn.nodeToStateful(ep_ind);
                                ep_space_buf = outglspace{ep_isf};
                                ep_ist = sn.nodeToStation(ep_ind);
                                consume_count = glevent.passive{j}.weight; % one completion = one firing (see PRE note above)
                                % Check scheduling strategy and state format to determine handling
                                is_queue_based = sn.nodetype(ep_ind) == NodeType.Place && ...
                                    length(ep_space_buf) ~= R && length(ep_space_buf) ~= 2*R;
                                if is_queue_based || sn.sched(ep_ist) == SchedStrategy.FCFS
                                    % FCFS queue-based handling
                                    idx = find(ep_space_buf == ep_class);
                                    if length(idx) >= consume_count
                                        ep_space_buf(idx(end - consume_count + 1:end)) = [];
                                    end
                                elseif sn.sched(ep_ist) == SchedStrategy.LCFS
                                    idx = find(ep_space_buf == ep_class);
                                    if length(idx) >= consume_count
                                        ep_space_buf(idx(1:consume_count)) = [];
                                    end
                                elseif sn.sched(ep_ist) == SchedStrategy.SIRO
                                    ep_space_buf(ep_class) = ep_space_buf(ep_class) - consume_count;
                                elseif sn.nodetype(ep_ind) == NodeType.Place
                                    % Place count-based state format
                                    state_len = length(ep_space_buf);
                                    if state_len > R
                                        buf_jobs = ep_space_buf(ep_class);
                                        srv_jobs = ep_space_buf(R + ep_class);
                                        total_jobs = buf_jobs + srv_jobs;
                                        remaining = total_jobs - consume_count;
                                        ep_space_buf(ep_class) = remaining;
                                        ep_space_buf(R + ep_class) = 0;
                                    else
                                        ep_space_buf(ep_class) = ep_space_buf(ep_class) - consume_count;
                                    end
                                else
                                    line_error(mfilename,sprintf('Scheduling strategy %s is unsupported at places.', SchedStrategy.toText(sn.sched(ep_ist))));
                                end
                                outglspace{ep_isf} = ep_space_buf;
                            end
                        end

                        %% Process ID_POST events
                        for j=1:length(glevent.passive)
                            if glevent.passive{j}.event == EventType.POST
                                fp_linidx = glevent.passive{j}.node;
                                [fp_ind, fp_class] = ind2sub([sn.nnodes, R], fp_linidx);
                                if fp_ind > length(sn.nodeToStateful) || isnan(sn.nodeToStateful(fp_ind)) || sn.nodeToStateful(fp_ind) <= 0
                                    continue;
                                end
                                fp_isf = sn.nodeToStateful(fp_ind);
                                fp_ist = sn.nodeToStation(fp_ind);
                                fp_space_buf = outglspace{fp_isf};
                                produce_count = glevent.passive{j}.weight;
                                % Check scheduling strategy and state format to determine handling
                                % Use original glspace length for detection (before PRE modified it)
                                orig_len = length(glspace{fp_isf});
                                is_queue_based = sn.nodetype(fp_ind) == NodeType.Place && ...
                                    orig_len ~= R && orig_len ~= 2*R;
                                if is_queue_based || sn.sched(fp_ist) == SchedStrategy.FCFS || sn.sched(fp_ist) == SchedStrategy.LCFS
                                    fp_space_buf = [repmat(fp_class, 1, produce_count), fp_space_buf];
                                elseif sn.sched(fp_ist) == SchedStrategy.SIRO
                                    fp_space_buf(fp_class) = fp_space_buf(fp_class) + produce_count;
                                elseif sn.nodetype(fp_ind) == NodeType.Place
                                    fp_space_buf(fp_class) = fp_space_buf(fp_class) + produce_count;
                                else
                                    line_error(mfilename,sprintf('Scheduling strategy %s is unsupported at places.', SchedStrategy.toText(sn.sched(fp_ist))));
                                end
                                outglspace{fp_isf} = fp_space_buf;
                            end
                        end
                    end
                end
        end
    end
end

if isempty(outcomp) && ~isempty(outrate)
    outcomp = false(numel(outrate),1); % ENABLE outcomes are never completions
end

% TODO: if there are FIRE events with Immediate rate it
% is unclear if the ENABLE event should get priority.
% The order of choice of what fires could change the system behavior.
% On the other hand, this could cause perpetual loops if the enabling is
% applied with priority on the same state with a FIRE?
if isSimulation
    if size(outspace,1) > 1 && ~isempty(outprob)
        % Use effective rate (rate * prob) for selection to avoid selecting zero-prob outcomes
        eff_rate = outrate .* outprob;
        tot_rate = sum(eff_rate);
        if tot_rate > 0
            cum_rate = cumsum(eff_rate) / tot_rate;
            firing_ctr = 1 + max([0,find( rand > cum_rate' )]); % select action
            outspace = outspace(firing_ctr,:);
            outrate = sum(eff_rate);  % return effective rate
            outprob = 1.0;  % probability is already incorporated into rate
            if numel(outcomp) >= firing_ctr
                outcomp = outcomp(firing_ctr);
            end
        else
            % All effective rates are zero, pick first option with non-zero prob if any
            valid_idx = find(outprob > 0, 1);
            if isempty(valid_idx)
                valid_idx = 1;
            end
            outspace = outspace(valid_idx,:);
            outrate = 0;
            outprob = outprob(valid_idx,:);
            if numel(outcomp) >= valid_idx
                outcomp = outcomp(valid_idx);
            end
        end
        % Update the global state with the selected outcome for the active node
        if exist('isf', 'var') && isf > 0 && isf <= length(outglspace)
            outglspace{isf} = outspace;
        end
    end
end
end
