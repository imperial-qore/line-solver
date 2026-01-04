function [outspace, outrate, outprob, eventCache] = afterEventStation(sn, ind, inspace, event, class, isSimulation, eventCache, ...
    M, R, S, phasessz, phaseshift, pie, isf, ismkvmod, ismkvmodclass, lldscaling, lldlimit, cdscaling, ...
    hasOnlyExp, ist, K, Ks, mu, phi, proc, capacity, classcap, V, space_buf, space_srv, space_var, key)
outspace = [];
outrate = [];
outprob = 1;
switch event
    case EventType.ARV %% passive
        % return if there is no space to accept the arrival
        [ni,nir] = State.toMarginalAggr(sn,ind,inspace,K,Ks,space_buf,space_srv,space_var);
        % otherwise check scheduling strategy
        pentry = pie{ist}{class};
        % For Place nodes (INF scheduling with NaN service), use uniform entry probability
        if all(isnan(pentry))
            pentry = ones(size(pentry)) / length(pentry);
        end
        outprob = [];
        outprob_k = [];
        for kentry = 1:K(class)
            space_var_k = space_var;
            space_srv_k = space_srv;
            space_buf_k = space_buf;
            switch sn.sched(ist)
                case SchedStrategy.EXT % source, can receive any "virtual" arrival from the sink as long as it is from an open class
                    if isinf(sn.njobs(class))
                        outspace = inspace;
                        outrate = -1*zeros(size(outspace,1)); % passive action, rate is unspecified
                        outprob = ones(size(outspace,1));
                        break
                    end
                case {SchedStrategy.PS, SchedStrategy.INF, SchedStrategy.DPS, SchedStrategy.GPS, SchedStrategy.PSPRIO, SchedStrategy.DPSPRIO, SchedStrategy.GPSPRIO, SchedStrategy.LPS}
                    % job enters service immediately
                    if space_srv_k(:,Ks(class)+kentry) < classcap(ist,class)
                        space_srv_k(:,Ks(class)+kentry) = space_srv_k(:,Ks(class)+kentry) + 1;
                        outprob_k = pentry(kentry)*ones(size(space_srv_k,1));
                    else
                        outprob_k = pentry(kentry)*zeros(size(space_srv_k,1));
                    end
                case {SchedStrategy.SIRO, SchedStrategy.SEPT, SchedStrategy.LEPT}
                    if ni<S(ist)
                        space_srv_k(:,Ks(class)+kentry) = space_srv_k(:,Ks(class)+kentry) + 1;
                        outprob_k = pentry(kentry)*ones(size(space_srv_k,1));
                    else
                        space_buf_k(:,class) = space_buf_k(:,class) + 1;
                        outprob_k = pentry(kentry)*ones(size(space_srv_k,1));
                    end
                case {SchedStrategy.FCFS, SchedStrategy.HOL, SchedStrategy.LCFS, SchedStrategy.LCFSPRIO}
                    % find states with all servers busy - this
                    % needs not to be moved

                    % if MAP service, when empty restart from the phase
                    % stored in space_var for this class
                    if ~ismkvmodclass(class) || (ismkvmodclass(class) && kentry == space_var(sum(sn.nvars(ind,1:class))))
                        if ismkvmodclass(class)
                            pentry = zeros(size(pentry));
                            pentry(kentry) = 1.0;
                        end
                        all_busy_srv = sum(space_srv_k,2) >= S(ist);

                        % find and modify states with an idle server
                        idle_srv = sum(space_srv_k,2) < S(ist);
                        space_srv_k(idle_srv, end-sum(K)+Ks(class)+kentry) = space_srv_k(idle_srv,end-sum(K)+Ks(class)+kentry) + 1; % job enters service

                        % this section dynamically grows the number of
                        % elements in the buffer

                        if any(ni < capacity(ist))
                            if any(nir(:,class) < classcap(ist,class)) % if there is room
                                if ~any(space_buf_k(:)==0) % but the buffer has no empty slots
                                    % append job slot
                                    space_buf_k = [zeros(size(space_buf_k,1),1),space_buf_k];
                                end
                            end
                        end
                        %end
                        %get position of first empty slot
                        empty_slots = -1*ones(size(all_busy_srv,1),1);
                        if size(space_buf_k,2) == 0
                            empty_slots(all_busy_srv) = false;
                        elseif size(space_buf_k,2) == 1
                            empty_slots(all_busy_srv) = space_buf_k(all_busy_srv,:)==0;
                        else
                            empty_slots(all_busy_srv) = max(bsxfun(@times, space_buf_k(all_busy_srv,:)==0, [1:size(space_buf_k,2)]),[],2);
                        end

                        % ignore states where the buffer has no empty slots
                        wbuf_empty = empty_slots>0;
                        if any(wbuf_empty)
                            space_srv_k = space_srv_k(wbuf_empty,:);
                            space_buf_k = space_buf_k(wbuf_empty,:);
                            space_var_k = space_var_k(wbuf_empty,:);
                            empty_slots = empty_slots(wbuf_empty);
                            space_buf_k(sub2ind(size(space_buf_k),1:size(space_buf_k,1),empty_slots')) = class;
                            %outspace(all_busy_srv(wbuf_empty),:) = [space_buf, space_srv, space_var];
                        end
                        outprob_k = pentry(kentry)*ones(size(space_srv_k,1),1);
                    else
                        outprob_k = 0*ones(size(space_srv_k,1),1); % zero probability event
                    end
                case {SchedStrategy.FCFSPR,SchedStrategy.FCFSPI,SchedStrategy.FCFSPRPRIO,SchedStrategy.FCFSPIPRIO,SchedStrategy.LCFSPR,SchedStrategy.LCFSPI,SchedStrategy.LCFSPRPRIO,SchedStrategy.LCFSPIPRIO}
                    % find states with all servers busy - this
                    % must not be moved
                    all_busy_srv = sum(space_srv_k,2) >= S(ist);
                    % find states with an idle server
                    idle_srv = sum(space_srv_k,2) < S(ist);

                    % reorder states so that idle ones come first
                    space_buf_k_reord = space_buf_k(idle_srv,:);
                    space_srv_k_reord = space_srv_k(idle_srv,:);
                    space_var_k_reord = space_var_k(idle_srv,:);

                    % if idle, the job enters service in phase kentry
                    if any(idle_srv)
                        space_srv_k_reord(:, end-sum(K)+Ks(class)+kentry) = space_srv_k_reord(:,end-sum(K)+Ks(class)+kentry) + 1;
                        outprob_k = pentry(kentry);
                    else
                        % if all busy, expand output states for all possible choices of job class to preempt
                        psentry = ones(size(space_buf_k_reord,1),1); % probability scaling due to preemption
                        for classpreempt = 1:R
                            % For priority variants, only higher-priority (lower class number) jobs can preempt
                            if (sn.sched(ist) == SchedStrategy.FCFSPRPRIO || sn.sched(ist) == SchedStrategy.FCFSPIPRIO || sn.sched(ist) == SchedStrategy.LCFSPRPRIO || sn.sched(ist) == SchedStrategy.LCFSPIPRIO)
                                if class >= classpreempt % arriving job has same or lower priority, cannot preempt
                                    continue;
                                end
                            end
                            for phasepreempt = 1:K(classpreempt) % phase of job to preempt
                                si_preempt = space_srv_k(:, (end-sum(K)+Ks(classpreempt)+phasepreempt));
                                busy_preempt = si_preempt > 0; % states where there is at least on class-r job in execution
                                if any(busy_preempt)
                                    psentry = [psentry; si_preempt(busy_preempt) ./ sum(space_srv_k,2)];
                                    space_srv_k_preempt = space_srv_k(busy_preempt,:);
                                    space_buf_k_preempt = space_buf_k(busy_preempt,:);
                                    space_var_k_preempt = space_var_k(busy_preempt,:);
                                    space_srv_k_preempt(:, end-sum(K)+Ks(classpreempt)+phasepreempt) = space_srv_k_preempt(:,end-sum(K)+Ks(classpreempt)+phasepreempt) - 1; % remove preempted job
                                    space_srv_k_preempt(:, end-sum(K)+Ks(class)+kentry) = space_srv_k_preempt(:,end-sum(K)+Ks(class)+kentry) + 1;

                                    % dynamically grow buffer lenght in
                                    % simulation
                                    if isSimulation
                                        if ni < capacity(ist) && nir(class) < classcap(ist,class) % if there is room
                                            if ~any(space_buf_k_preempt(:)==0) % but the buffer has no empty slots
                                                % append job slot
                                                space_buf_k_preempt = [zeros(size(space_buf_k_preempt,1),2),space_buf_k]; % append two columns for (class, preempt-phase)
                                            end
                                        end
                                    end

                                    %get position of first empty slot
                                    empty_slots = -1*ones(sum(busy_preempt),1);
                                    if size(space_buf_k_preempt,2) == 0
                                        empty_slots(busy_preempt) = false;
                                    elseif size(space_buf_k_preempt,2) == 2 % 2 due to (class, preempt-phase) pairs
                                        empty_slots(busy_preempt) = space_buf_k_preempt(busy_preempt,1:2:end)==0;
                                    else
                                        empty_slots(busy_preempt) = max(bsxfun(@times, space_buf_k_preempt(busy_preempt,:)==0, [1:size(space_buf_k_preempt,2)]),[],2)-1; %-1 due to (class, preempt-phase) pairs
                                    end

                                    % ignore states where the buffer has no empty slots
                                    wbuf_empty = empty_slots>0;
                                    if any(wbuf_empty)
                                        space_srv_k_preempt = space_srv_k_preempt(wbuf_empty,:);
                                        space_buf_k_preempt = space_buf_k_preempt(wbuf_empty,:);
                                        space_var_k_preempt = space_var_k_preempt(wbuf_empty,:);
                                        empty_slots = empty_slots(wbuf_empty);
                                        if sn.sched(ist) == SchedStrategy.LCFSPR || sn.sched(ist) == SchedStrategy.LCFSPRPRIO || sn.sched(ist) == SchedStrategy.FCFSPR || sn.sched(ist) == SchedStrategy.FCFSPRPRIO % preempt-resume
                                            space_buf_k_preempt(sub2ind(size(space_buf_k_preempt),1:size(space_buf_k_preempt,1),empty_slots')+1) = phasepreempt;
                                        elseif sn.sched(ist) == SchedStrategy.LCFSPI || sn.sched(ist) == SchedStrategy.LCFSPIPRIO || sn.sched(ist) == SchedStrategy.FCFSPI || sn.sched(ist) == SchedStrategy.FCFSPIPRIO % preempt-independent
                                            space_buf_k_preempt(sub2ind(size(space_buf_k_preempt),1:size(space_buf_k_preempt,1),empty_slots')+1) = 1;
                                        end
                                        space_buf_k_preempt(sub2ind(size(space_buf_k_preempt),1:size(space_buf_k_preempt,1),empty_slots')) = classpreempt;
                                        %outspace(all_busy_srv(wbuf_empty),:) = [space_buf, space_srv, space_var];
                                    end
                                    space_srv_k_reord = [space_srv_k_reord; space_srv_k_preempt];
                                    space_buf_k_reord = [space_buf_k_reord; space_buf_k_preempt];
                                    space_var_k_reord = [space_var_k_reord; space_var_k_preempt];
                                end
                            end
                        end
                        outprob_k = pentry(kentry) * psentry .* ones(size(space_srv_k_reord,1),1);
                    end
                    space_buf_k = space_buf_k_reord; % save reordered output states
                    space_srv_k = space_srv_k_reord; % save reordered output states
                    space_var_k = space_var_k_reord; % save reordered output states
            end
            % form the new state
            outspace_k = [space_buf_k, space_srv_k, space_var_k];
            % remove states where new arrival violates capacity or cutoff constraints
            [oi,oir] = State.toMarginalAggr(sn,ind,outspace_k,K,Ks,space_buf_k,space_srv_k,space_var_k);
            en_o = classcap(ist,class)>= oir(:,class) | capacity(ist)*ones(size(oi,1),1) >= oi;

            if size(outspace,2)>size(outspace_k(en_o,:),2)
                outspace = [outspace; zeros(1,size(outspace,2)-size(outspace_k(en_o,:),2)),outspace_k(en_o,:)];
            elseif size(outspace,2)<size(outspace_k(en_o,:),2)
                outspace = [zeros(size(outspace,1),size(outspace_k(en_o,:),2)-size(outspace,2)), outspace; outspace_k(en_o,:)];
            else
                outspace = [outspace; outspace_k(en_o,:)];
            end
            outrate = [outrate; -1*ones(size(outspace_k(en_o,:),1),1)]; % passive action, rate is unspecified
            outprob = [outprob; outprob_k(en_o,:)];
        end
        if isSimulation
            if size(outprob,1) > 1
                cum_prob = cumsum(outprob) / sum(outprob);
                firing_ctr = 1 + max([0,find( rand > cum_prob' )]); % select action
                outspace = outspace(firing_ctr,:);
                outrate = -1;
                outprob = 1;
            end
        end
    case EventType.DEP
        if any(any(space_srv(:,(Ks(class)+1):(Ks(class)+K(class))))) % something is busy
            if hasOnlyExp && (sn.sched(ist) == SchedStrategy.PS || sn.sched(ist) == SchedStrategy.DPS || sn.sched(ist) == SchedStrategy.GPS || sn.sched(ist) == SchedStrategy.INF || sn.sched(ist) == SchedStrategy.PSPRIO || sn.sched(ist) == SchedStrategy.DPSPRIO || sn.sched(ist) == SchedStrategy.GPSPRIO || sn.sched(ist) == SchedStrategy.LPS)
                nir = space_srv;
                ni = sum(nir,2);
                sir = nir;
                kir = sir;
            else
                [ni,nir,sir,kir] = State.toMarginal(sn,ind,inspace,K,Ks,space_buf,space_srv,space_var);
            end
            switch sn.routing(ind,class)
                case RoutingStrategy.RROBIN
                    idx = find(space_var(sum(sn.nvars(ind,1:(R+class)))) == sn.nodeparam{ind}{class}.outlinks);
                    if idx < length(sn.nodeparam{ind}{class}.outlinks)
                        space_var(sum(sn.nvars(ind,1:(R+class)))) = sn.nodeparam{ind}{class}.outlinks(idx+1);
                    else
                        space_var(sum(sn.nvars(ind,1:(R+class)))) = sn.nodeparam{ind}{class}.outlinks(1);
                    end
            end
            if sir(class)>0 % is a job of class is in service
                outprob = [];
                for k=1:K(class)
                    space_srv = inspace(:,(end-sum(K)-V+1):(end-V)); % server state
                    space_buf = inspace(:,1:(end-sum(K)-V)); % buffer state
                    rate = zeros(size(space_srv,1),1);
                    en =  space_srv(:,Ks(class)+k) > 0;
                    if any(en)
                        switch sn.sched(ist)
                            case SchedStrategy.EXT % source, can produce an arrival from phase-k as long as it is from an open class
                                if isinf(sn.njobs(class))
                                    pentry = pie{ist}{class};
                                    for kentry = 1:K(class)
                                        space_srv = inspace(:,(end-sum(K)-V+1):(end-V)); % server state
                                        space_srv(en,Ks(class)+k) = space_srv(en,Ks(class)+k) - 1; % record departure
                                        space_srv(en,Ks(class)+kentry) = space_srv(en,Ks(class)+kentry) + 1; % new job
                                        outspace = [outspace; space_buf(en,:), space_srv(en,:), space_var(en,:)];
                                        if isinf(ni) % hit limited load-dependence
                                            outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,end).*pentry(kentry)*mu{ist}{class}(k)*phi{ist}{class}(k)*ones(size(inspace(en,:),1),1)];
                                        else
                                            outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,min(ni(en),lldlimit)).*pentry(kentry)*mu{ist}{class}(k)*phi{ist}{class}(k)*ones(size(inspace(en,:),1),1)];
                                        end
                                        outprob = [outprob; ones(size(space_buf(en,:),1),1)];
                                    end
                                end
                            case SchedStrategy.INF % move first job in service
                                space_srv(en,Ks(class)+k) = space_srv(en,Ks(class)+k) - 1; % record departure
                                rate(en) = mu{ist}{class}(k)*(phi{ist}{class}(k)).*kir(en,class,k); % assume active
                                % if state is unchanged, still add with rate 0
                                outspace = [outspace; space_buf(en,:), space_srv(en,:), space_var(en,:)];
                                if isinf(ni) % hit limited load-dependence
                                    outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,end).*rate(en,:)];
                                else
                                    outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,min(ni(en),lldlimit)).*rate(en,:)];
                                end
                                outprob = [outprob; ones(size(rate(en,:),1),1)];
                            case SchedStrategy.PS % move first job in service
                                space_srv(en,Ks(class)+k) = space_srv(en,Ks(class)+k) - 1; % record departure
                                rate(en) = mu{ist}{class}(k)*(phi{ist}{class}(k)).*(kir(en,class,k)./ni(en)).*min(ni(en),S(ist)); % assume active
                                % if state is unchanged, still add with rate 0
                                outspace = [outspace; space_buf(en,:), space_srv(en,:), space_var(en,:)];
                                if isinf(ni) % hit limited load-dependence
                                    outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,end).*rate(en,:)];
                                else
                                    outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,min(ni(en),lldlimit)).*rate(en,:)];
                                end
                                outprob = [outprob; ones(size(rate(en,:),1),1)];
                                %                                    end
                            case SchedStrategy.PSPRIO
                                % unclear if LD scaling should be with
                                % ni or with niprio, for now left as ni
                                % for consistency with HOL multiserver
                                if all(ni(en) <= S(ist))
                                    % n <= c: all jobs get service, priority doesn't matter
                                    space_srv(en,Ks(class)+k) = space_srv(en,Ks(class)+k) - 1; % record departure
                                    rate(en) = mu{ist}{class}(k)*(phi{ist}{class}(k)).*(kir(en,class,k)./ni(en)).*min(ni(en),S(ist)); % assume active
                                    outspace = [outspace; space_buf(en,:), space_srv(en,:), space_var(en,:)];
                                    if isinf(ni) % hit limited load-dependence
                                        outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,end).*rate(en,:)];
                                    else
                                        outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,min(ni(en),lldlimit)).*rate(en,:)];
                                    end
                                    outprob = [outprob; ones(size(rate(en,:),1),1)];
                                elseif sn.classprio(class) == min(sn.classprio(nir>0)) % if this class is in the most urgent priority group (lower value = higher priority in LINE)
                                    space_srv(en,Ks(class)+k) = space_srv(en,Ks(class)+k) - 1; % record departure
                                    niprio(en) = sum(nir(sn.classprio==sn.classprio(class))); % jobs at the same priority class
                                    rate(en) = mu{ist}{class}(k)*(phi{ist}{class}(k)).*(kir(en,class,k)./niprio(en)).*min(niprio(en),S(ist)); % assume active
                                    % if state is unchanged, still add with rate 0
                                    outspace = [outspace; space_buf(en,:), space_srv(en,:), space_var(en,:)];
                                    if isinf(ni) % hit limited load-dependence
                                        outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,end).*rate(en,:)];
                                    else
                                        outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,min(niprio(en),lldlimit)).*rate(en,:)];
                                    end
                                    outprob = [outprob; ones(size(rate(en,:),1),1)];
                                else % n > c and not highest priority: set rate to zero
                                    outspace = [outspace; space_buf(en,:), space_srv(en,:), space_var(en,:)];
                                    outrate = [outrate; zeros(size(rate(en,:),1),1)];
                                    outprob = [outprob; ones(size(rate(en,:),1),1)];
                                end
                            case SchedStrategy.DPS
                                space_srv(en,Ks(class)+k) = space_srv(en,Ks(class)+k) - 1; % record departure
                                if S(ist) > 1
                                    line_error(mfilename,'Multi-server DPS stations are not supported yet.');
                                end
                                % in GPS, the scheduling parameter are the weights
                                w_i = sn.schedparam(ist,:);
                                w_i = w_i / sum(w_i);
                                rate(en) = mu{ist}{class}(k)*(phi{ist}{class}(k))*(kir(en,class,k)/nir(class))*w_i(class)*nir(class)./(sum(repmat(w_i,sum(en),1)*nir',2));
                                % if state is unchanged, still add with rate 0
                                outspace = [outspace; space_buf(en,:), space_srv(en,:), space_var(en,:)];
                                if isinf(ni) % hit limited load-dependence
                                    outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,end).*rate(en,:)];
                                else
                                    outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,min(ni(en),lldlimit)).*rate(en,:)];
                                end
                                outprob = [outprob; ones(size(rate(en,:),1),1)];
                            case SchedStrategy.DPSPRIO
                                if all(ni(en) <= S(ist))
                                    % n <= c: all jobs get service, behave like regular DPS
                                    space_srv(en,Ks(class)+k) = space_srv(en,Ks(class)+k) - 1; % record departure
                                    if S(ist) > 1
                                        line_error(mfilename,'Multi-server DPS stations are not supported yet.');
                                    end
                                    w_i = sn.schedparam(ist,:);
                                    w_i = w_i / sum(w_i);
                                    rate(en) = mu{ist}{class}(k)*(phi{ist}{class}(k))*(kir(en,class,k)/nir(class))*w_i(class)*nir(class)./(sum(repmat(w_i,sum(en),1)*nir',2));
                                    outspace = [outspace; space_buf(en,:), space_srv(en,:), space_var(en,:)];
                                    if isinf(ni) % hit limited load-dependence
                                        outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,end).*rate(en,:)];
                                    else
                                        outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,min(ni(en),lldlimit)).*rate(en,:)];
                                    end
                                    outprob = [outprob; ones(size(rate(en,:),1),1)];
                                elseif sn.classprio(class) == min(sn.classprio(nir>0)) % if this class is in the most urgent priority group (lower value = higher priority in LINE)
                                    nirprio = nir;
                                    nirprio(sn.classprio~=sn.classprio(class)) = 0; % ignore jobs of lower priority
                                    niprio = sum(nirprio);
                                    space_srv(en,Ks(class)+k) = space_srv(en,Ks(class)+k) - 1; % record departure
                                    if S(ist) > 1
                                        line_error(mfilename,'Multi-server DPS stations are not supported yet.');
                                    end
                                    % in GPS, the scheduling parameter are the weights
                                    w_i = sn.schedparam(ist,:);
                                    w_i = w_i / sum(w_i);
                                    rate(en) = mu{ist}{class}(k)*(phi{ist}{class}(k))*(kir(en,class,k)/nirprio(class))*w_i(class)*nirprio(class)./(sum(repmat(w_i,sum(en),1)*nirprio',2));
                                    % if state is unchanged, still add with rate 0
                                    outspace = [outspace; space_buf(en,:), space_srv(en,:), space_var(en,:)];
                                    if isinf(ni) % hit limited load-dependence
                                        outrate = [outrate; cdscaling{ist}(nirprio).*lldscaling(ist,end).*rate(en,:)];
                                    else
                                        outrate = [outrate; cdscaling{ist}(nirprio).*lldscaling(ist,min(niprio(en),lldlimit)).*rate(en,:)];
                                    end
                                    outprob = [outprob; ones(size(rate(en,:),1),1)];
                                else % n > c and not most urgent priority: set rate to zero
                                    outspace = [outspace; space_buf(en,:), space_srv(en,:), space_var(en,:)];
                                    outrate = [outrate; zeros(size(rate(en,:),1),1)];
                                    outprob = [outprob; ones(size(rate(en,:),1),1)];
                                end
                            case SchedStrategy.GPS
                                space_srv(en,Ks(class)+k) = space_srv(en,Ks(class)+k) - 1; % record departure
                                if S(ist) > 1
                                    line_error(mfilename,'Multi-server GPS stations are not supported yet.');
                                end
                                % in GPS, the scheduling parameter are the weights
                                w_i = sn.schedparam(ist,:);
                                w_i = w_i / sum(w_i);
                                cir = min(nir,ones(size(nir)));
                                rate = mu{ist}{class}(k)*(phi{ist}{class}(k))*(kir(en,class,k)/nir(class))*w_i(class)/(w_i*cir(:)); % assume active
                                % if state is unchanged, still add with rate 0
                                outspace = [outspace; space_buf(en,:), space_srv(en,:), space_var(en,:)];
                                if isinf(ni) % hit limited load-dependence
                                    outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,end).*rate(en,:)];
                                else
                                    outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,min(ni(en),lldlimit)).*rate(en,:)];
                                end
                                outprob = [outprob; ones(size(rate(en,:),1),1)];
                            case SchedStrategy.GPSPRIO
                                if all(ni(en) <= S(ist))
                                    % n <= c: all jobs get service, behave like regular GPS
                                    space_srv(en,Ks(class)+k) = space_srv(en,Ks(class)+k) - 1; % record departure
                                    if S(ist) > 1
                                        line_error(mfilename,'Multi-server GPS stations are not supported yet.');
                                    end
                                    w_i = sn.schedparam(ist,:);
                                    w_i = w_i / sum(w_i);
                                    cir = min(nir,ones(size(nir)));
                                    rate = mu{ist}{class}(k)*(phi{ist}{class}(k))*(kir(en,class,k)/nir(class))*w_i(class)/(w_i*cir(:)); % assume active
                                    outspace = [outspace; space_buf(en,:), space_srv(en,:), space_var(en,:)];
                                    if isinf(ni) % hit limited load-dependence
                                        outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,end).*rate(en,:)];
                                    else
                                        outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,min(ni(en),lldlimit)).*rate(en,:)];
                                    end
                                    outprob = [outprob; ones(size(rate(en,:),1),1)];
                                elseif sn.classprio(class) == min(sn.classprio(nir>0)) % if this class is in the most urgent priority group (lower value = higher priority in LINE)
                                    nirprio = nir;
                                    nirprio(sn.classprio~=sn.classprio(class)) = 0; % ignore jobs of lower priority
                                    niprio = sum(nirprio);
                                    space_srv(en,Ks(class)+k) = space_srv(en,Ks(class)+k) - 1; % record departure
                                    if S(ist) > 1
                                        line_error(mfilename,'Multi-server DPS stations are not supported yet.');
                                    end
                                    % in GPS, the scheduling parameter are the weights
                                    w_i = sn.schedparam(ist,:);
                                    w_i = w_i / sum(w_i);
                                    cir = min(nirprio,ones(size(nirprio)));
                                    rate = mu{ist}{class}(k)*(phi{ist}{class}(k))*(kir(en,class,k)/nirprio(class))*w_i(class)/(w_i*cir(:)); % assume active
                                    % if state is unchanged, still add with rate 0
                                    outspace = [outspace; space_buf(en,:), space_srv(en,:), space_var(en,:)];
                                    if isinf(ni) % hit limited load-dependence
                                        outrate = [outrate; cdscaling{ist}(nirprio).*lldscaling(ist,end).*rate(en,:)];
                                    else
                                        outrate = [outrate; cdscaling{ist}(nirprio).*lldscaling(ist,min(niprio(en),lldlimit)).*rate(en,:)];
                                    end
                                    outprob = [outprob; ones(size(rate(en,:),1),1)];
                                else % n > c and not most urgent priority: set rate to zero
                                    outspace = [outspace; space_buf(en,:), space_srv(en,:), space_var(en,:)];
                                    outrate = [outrate; zeros(size(rate(en,:),1),1)];
                                    outprob = [outprob; ones(size(rate(en,:),1),1)];
                                end
                            case SchedStrategy.FCFS % move first job in service
                                space_srv(en,Ks(class)+k) = space_srv(en,Ks(class)+k) - 1; % record departure
                                en_wbuf = en & ni>S(ist); %states with jobs in buffer
                                for kdest=1:K(class) % new phase
                                    space_buf_kd = space_buf;
                                    space_var_kd = space_var;
                                    if ismkvmodclass(class)
                                        space_var_kd(en,sum(sn.nvars(ind,1:class))) = kdest;
                                    end
                                    rate_kd = rate;
                                    rate_kd(en) = proc{ist}{class}{2}(k,kdest).*kir(en,class,k); % assume active
                                    % first process states without jobs in buffer
                                    en_wobuf = ~en_wbuf;
                                    if any(en_wobuf) %any state without jobs in buffer
                                        outspace = [outspace; space_buf_kd(en_wobuf,:), space_srv(en_wobuf,:), space_var_kd(en_wobuf,:)];
                                        if isinf(ni) % hit limited load-dependence
                                            outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,end).*rate_kd(en_wobuf,:)];
                                        else
                                            outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,min(ni(en),lldlimit)).*rate_kd(en_wobuf,:)];
                                        end
                                    end
                                    % now process states with jobs in buffer
                                    outprob = [outprob; ones(size(rate_kd(en_wobuf,:),1),1)];
                                    if any(en_wbuf) %any state with jobs in buffer
                                        % get class of job at head
                                        start_svc_class = space_buf_kd(en_wbuf,end);
                                        if start_svc_class > 0 % redunant if?
                                            % update input buffer
                                            space_buf_kd(en_wbuf,:) = [zeros(sum(en_wbuf),1),space_buf_kd(en_wbuf,1:end-1)];
                                            % probability vector for the next job of starting in phase kentry
                                            if ismkvmodclass(start_svc_class) % if markov-modulated
                                                if start_svc_class==class % if successive service from the same class
                                                    kentry_range = kdest; % new job enters in phase left by departing job
                                                else % resume phase from local variables
                                                    kentry_range = space_var_kd(en,sum(sn.nvars(ind,1:start_svc_class)));
                                                end
                                                pentry_svc_class = 0*pie{ist}{start_svc_class};
                                                pentry_svc_class(kentry_range) = 1.0;
                                            else % if i.i.d.
                                                pentry_svc_class = pie{ist}{start_svc_class};
                                                kentry_range = 1:K(start_svc_class);
                                            end
                                            for kentry = kentry_range
                                                space_srv(en_wbuf,Ks(start_svc_class)+kentry) = space_srv(en_wbuf,Ks(start_svc_class)+kentry) + 1;
                                                outspace = [outspace; space_buf_kd(en,:), space_srv(en,:), space_var_kd(en,:)];
                                                rate_k = rate_kd;
                                                rate_k(en_wbuf,:) = rate_kd(en_wbuf,:)*pentry_svc_class(kentry);
                                                if isinf(ni) % use limited load-dependence at the latest user-provided level
                                                    outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,end).*rate_k(en,:)];
                                                else
                                                    outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,min(ni(en),lldlimit)).*rate_k(en,:)];
                                                end
                                                outprob_cur = ones(size(rate_kd(en,:),1),1);
                                                outprob_cur(outrate==0.0) = 0;
                                                outprob = [outprob; outprob_cur'];
                                                space_srv(en_wbuf,Ks(start_svc_class)+kentry) = space_srv(en_wbuf,Ks(start_svc_class)+kentry) - 1;
                                            end
                                        end
                                    end
                                end
                                % if state is unchanged, still add with rate 0
                            case SchedStrategy.HOL % FCFS priority
                                rate(en) = mu{ist}{class}(k)*(phi{ist}{class}(k)).*kir(:,class,k); % assume active
                                en_wbuf = en & ni>S(ist); %states with jobs in buffer
                                en_wobuf = ~en_wbuf;
                                space_srv(en,Ks(class)+k) = space_srv(en,Ks(class)+k) - 1; % record departure
                                priogroup = [Inf,sn.classprio]; % Inf for empty positions (lower value = higher priority)
                                space_buf_groupg = arrayfun(@(x) priogroup(1+x), space_buf);
                                start_classprio = min(space_buf_groupg(en_wbuf,:),[],2); % min finds highest priority
                                isrowmax = space_buf_groupg == repmat(start_classprio, 1, size(space_buf_groupg,2));
                                [~,rightmostMaxPosFlipped]=max(fliplr(isrowmax),[],2);
                                rightmostMaxPos = size(isrowmax,2) - rightmostMaxPosFlipped + 1;
                                start_svc_class = space_buf(en_wbuf, rightmostMaxPos);
                                outspace = [outspace; space_buf(en_wobuf,:), space_srv(en_wobuf,:), space_var(en_wobuf,:)];
                                if isinf(ni) % hit limited load-dependence
                                    outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,end).*rate(en_wobuf,:)];
                                else
                                    outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,min(ni(en),lldlimit)).*rate(en_wobuf,:)];
                                end
                                outprob = [outprob; ones(size(rate(en_wobuf,:),1),1)];
                                if start_svc_class > 0
                                    pentry_svc_class = pie{ist}{start_svc_class};
                                    for kentry = 1:K(start_svc_class)
                                        space_srv_k = space_srv;
                                        space_buf_k = space_buf;
                                        space_srv_k(en_wbuf,Ks(start_svc_class)+kentry) = space_srv_k(en_wbuf,Ks(start_svc_class)+kentry) + 1;
                                        for j=find(en_wbuf)'
                                            space_buf_k(j,:) = [0, space_buf_k(j,1:rightmostMaxPos(j)-1), space_buf_k(j,(rightmostMaxPos(j)+1):end)];
                                        end
                                        % if state is unchanged, still add with rate 0
                                        outspace = [outspace; space_buf_k(en_wbuf,:), space_srv_k(en_wbuf,:), space_var(en_wbuf,:)];
                                        rate_k = rate;
                                        rate_k(en_wbuf,:) = rate(en_wbuf,:) * pentry_svc_class(kentry);
                                        if isinf(ni) % hit limited load-dependence
                                            outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,end).*rate_k(en_wbuf,:)];
                                        else
                                            outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,min(ni(en),lldlimit)).*rate_k(en_wbuf,:)];
                                        end
                                        outprob = [outprob; ones(size(rate_k(en_wbuf,:),1),1)];
                                    end
                                end
                            case SchedStrategy.LCFSPRIO % LCFS priority - like HOL but LCFS order within priority groups
                                rate(en) = mu{ist}{class}(k)*(phi{ist}{class}(k)).*kir(:,class,k); % assume active
                                en_wbuf = en & ni>S(ist); %states with jobs in buffer
                                en_wobuf = ~en_wbuf;
                                space_srv(en,Ks(class)+k) = space_srv(en,Ks(class)+k) - 1; % record departure
                                priogroup = [Inf,sn.classprio]; % Inf for empty positions (lower value = higher priority)
                                space_buf_groupg = arrayfun(@(x) priogroup(1+x), space_buf);
                                start_classprio = min(space_buf_groupg(en_wbuf,:),[],2); % min finds highest priority
                                isrowmax = space_buf_groupg == repmat(start_classprio, 1, size(space_buf_groupg,2));
                                % LCFS: Find leftmost (first) position instead of rightmost for LCFS order
                                [~,leftmostMaxPos]=max(isrowmax,[],2);
                                start_svc_class = space_buf(en_wbuf, leftmostMaxPos);
                                outspace = [outspace; space_buf(en_wobuf,:), space_srv(en_wobuf,:), space_var(en_wobuf,:)];
                                if isinf(ni) % hit limited load-dependence
                                    outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,end).*rate(en_wobuf,:)];
                                else
                                    outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,min(ni(en),lldlimit)).*rate(en_wobuf,:)];
                                end
                                outprob = [outprob; ones(size(rate(en_wobuf,:),1),1)];
                                if start_svc_class > 0
                                    pentry_svc_class = pie{ist}{start_svc_class};
                                    for kentry = 1:K(start_svc_class)
                                        space_srv_k = space_srv;
                                        space_buf_k = space_buf;
                                        space_srv_k(en_wbuf,Ks(start_svc_class)+kentry) = space_srv_k(en_wbuf,Ks(start_svc_class)+kentry) + 1;
                                        for j=find(en_wbuf)'
                                            % LCFS: Remove from leftmost position instead of rightmost
                                            space_buf_k(j,:) = [0, space_buf_k(j,1:leftmostMaxPos(j)-1), space_buf_k(j,(leftmostMaxPos(j)+1):end)];
                                        end
                                        % if state is unchanged, still add with rate 0
                                        outspace = [outspace; space_buf_k(en_wbuf,:), space_srv_k(en_wbuf,:), space_var(en_wbuf,:)];
                                        rate_k = rate;
                                        rate_k(en_wbuf,:) = rate_k(en_wbuf,:)*pentry_svc_class(kentry);
                                        if isinf(ni) % hit limited load-dependence
                                            outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,end).*rate_k(en_wbuf,:)];
                                        else
                                            outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,min(ni(en),lldlimit)).*rate_k(en_wbuf,:)];
                                        end
                                        outprob = [outprob; ones(size(rate_k(en_wbuf,:),1),1)];
                                        space_srv_k(en_wbuf,Ks(start_svc_class)+kentry) = space_srv_k(en_wbuf,Ks(start_svc_class)+kentry) - 1;
                                    end
                                end
                            case SchedStrategy.LCFS % move last job in service
                                space_srv(en,Ks(class)+k) = space_srv(en,Ks(class)+k) - 1; % record departure
                                rate(en) = mu{ist}{class}(k)*(phi{ist}{class}(k)).*kir(:,class,k); % assume active
                                en_wbuf = en & ni>S(ist); %states with jobs in buffer
                                [~, colfirstnnz] = max( space_buf(en_wbuf,:) ~=0, [], 2 ); % find first nnz column
                                start_svc_class = space_buf(en_wbuf,colfirstnnz); % job entering service
                                space_buf(en_wbuf,colfirstnnz)=0;
                                if isempty(start_svc_class)
                                    outspace = [outspace; space_buf(en,:), space_srv(en,:), space_var(en,:)];
                                    if isinf(ni) % hit limited load-dependence
                                        outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,end).*rate(en,:)];
                                    else
                                        outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,min(ni(en),lldlimit)).*rate(en,:)];
                                    end
                                    outprob = [outprob; ones(size(rate(en,:),1),1)];
                                    if isSimulation && nargin>=7 && isobject(eventCache)
                                        eventCache(key) = {outprob, outspace,outrate};
                                    end
                                    return
                                end
                                for kentry = 1:K(start_svc_class)
                                    pentry_svc_class = pie{ist}{start_svc_class};
                                    space_srv(en_wbuf,Ks(start_svc_class)+kentry) = space_srv(en_wbuf,Ks(start_svc_class)+kentry) + 1;
                                    % if state is unchanged, still add with rate 0
                                    outspace = [outspace; space_buf(en,:), space_srv(en,:), space_var(en,:)];
                                    rate_k = rate;
                                    rate_k(en_wbuf,:) = rate(en_wbuf,:)*pentry_svc_class(kentry);
                                    if isinf(ni) % hit limited load-dependence
                                        outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,end).*rate_k(en,:)];
                                    else
                                        outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,min(ni(en),lldlimit)).*rate_k(en,:)];
                                    end
                                    outprob = [outprob; ones(size(rate(en,:),1),1)];
                                    space_srv(en_wbuf,Ks(start_svc_class)+kentry) = space_srv(en_wbuf,Ks(start_svc_class)+kentry) - 1;
                                end
                            case SchedStrategy.LCFSPR % move last job in service (preempt-resume)
                                space_srv(en,Ks(class)+k) = space_srv(en,Ks(class)+k) - 1; % record departure
                                rate(en) = mu{ist}{class}(k)*(phi{ist}{class}(k)).*kir(:,class,k); % assume active
                                en_wbuf = en & ni>S(ist); %states with jobs in buffer
                                [~, colfirstnnz] = max( space_buf(en_wbuf,:) ~=0, [], 2 ); % find first nnz column
                                start_svc_class = space_buf(en_wbuf,colfirstnnz); % job entering service
                                kentry = space_buf(en_wbuf,colfirstnnz+1); % entry phase of job resuming service
                                space_buf(en_wbuf,colfirstnnz)=0;% zero popped job
                                space_buf(en_wbuf,colfirstnnz+1)=0; % zero popped phase
                                if isempty(start_svc_class)
                                    outspace = [outspace; space_buf(en,:), space_srv(en,:), space_var(en,:)];
                                    if isinf(ni) % hit limited load-dependence
                                        outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,end).*rate(en,:)];
                                    else
                                        outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,min(ni(en),lldlimit)).*rate(en,:)];
                                    end
                                    outprob = [outprob; ones(size(rate(en,:),1),1)];
                                    if isSimulation && nargin>=7 && isobject(eventCache)
                                        eventCache(key) = {outprob, outspace,outrate};
                                    end
                                    return
                                end
                                space_srv(en_wbuf,Ks(start_svc_class)+kentry) = space_srv(en_wbuf,Ks(start_svc_class)+kentry) + 1;
                                % if state is unchanged, still add with rate 0
                                outspace = [outspace; space_buf(en,:), space_srv(en,:), space_var(en,:)];
                                if isinf(ni) % hit limited load-dependence
                                    outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,end).*rate(en,:)];
                                else
                                    outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,min(ni(en),lldlimit)).*rate(en,:)];
                                end
                                outprob = [outprob; ones(size(rate(en,:),1),1)];
                            case SchedStrategy.LCFSPI % move last job in service (preempt-independent)
                                space_srv(en,Ks(class)+k) = space_srv(en,Ks(class)+k) - 1; % record departure
                                rate(en) = mu{ist}{class}(k)*(phi{ist}{class}(k)).*kir(:,class,k); % assume active
                                en_wbuf = en & ni>S(ist); %states with jobs in buffer
                                [~, colfirstnnz] = max( space_buf(en_wbuf,:) ~=0, [], 2 ); % find first nnz column
                                start_svc_class = space_buf(en_wbuf,colfirstnnz); % job entering service
                                space_buf(en_wbuf,colfirstnnz)=0;% zero popped job
                                space_buf(en_wbuf,colfirstnnz+1)=0; % zero popped phase (ignored for LCFSPI)
                                if isempty(start_svc_class)
                                    outspace = [outspace; space_buf(en,:), space_srv(en,:), space_var(en,:)];
                                    if isinf(ni) % hit limited load-dependence
                                        outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,end).*rate(en,:)];
                                    else
                                        outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,min(ni(en),lldlimit)).*rate(en,:)];
                                    end
                                    outprob = [outprob; ones(size(rate(en,:),1),1)];
                                    if isSimulation && nargin>=7 && isobject(eventCache)
                                        eventCache(key) = {outprob, outspace,outrate};
                                    end
                                    return
                                end
                                % For LCFSPI, jobs restart from pie distribution instead of stored phase
                                pentry_svc_class = pie{ist}{start_svc_class};
                                for kentry = 1:K(start_svc_class)
                                    space_srv_k = space_srv;
                                    space_srv_k(en_wbuf,Ks(start_svc_class)+kentry) = space_srv_k(en_wbuf,Ks(start_svc_class)+kentry) + 1;
                                    % if state is unchanged, still add with rate 0
                                    outspace = [outspace; space_buf(en,:), space_srv_k(en,:), space_var(en,:)];
                                    rate_k = rate;
                                    rate_k(en_wbuf,:) = rate_k(en_wbuf,:) * pentry_svc_class(kentry);
                                    if isinf(ni) % hit limited load-dependence
                                        outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,end).*rate_k(en,:)];
                                    else
                                        outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,min(ni(en),lldlimit)).*rate_k(en,:)];
                                    end
                                    outprob = [outprob; ones(size(rate_k(en,:),1),1)];
                                end
                            case SchedStrategy.FCFSPR % FCFS preempt-resume (no priority)
                                space_srv(en,Ks(class)+k) = space_srv(en,Ks(class)+k) - 1; % record departure
                                rate(en) = mu{ist}{class}(k)*(phi{ist}{class}(k)).*kir(:,class,k); % assume active
                                en_wbuf = en & ni>S(ist); %states with jobs in buffer
                                % FCFS: Find rightmost (last) non-zero position (longest waiting job)
                                [~, colLastNnz] = max(fliplr(space_buf(en_wbuf,:) ~= 0), [], 2);
                                colLastNnz = size(space_buf,2) - colLastNnz; % convert to actual position (odd index for class)
                                start_svc_class = space_buf(en_wbuf, colLastNnz); % job entering service
                                kentry = space_buf(en_wbuf, colLastNnz+1); % entry phase of job resuming service
                                space_buf(en_wbuf, colLastNnz) = 0; % zero popped job
                                space_buf(en_wbuf, colLastNnz+1) = 0; % zero popped phase
                                if isempty(start_svc_class)
                                    outspace = [outspace; space_buf(en,:), space_srv(en,:), space_var(en,:)];
                                    if isinf(ni)
                                        outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,end).*rate(en,:)];
                                    else
                                        outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,min(ni(en),lldlimit)).*rate(en,:)];
                                    end
                                    outprob = [outprob; ones(size(rate(en,:),1),1)];
                                    if isSimulation && nargin>=7 && isobject(eventCache)
                                        eventCache(key) = {outprob, outspace, outrate};
                                    end
                                    return
                                end
                                space_srv(en_wbuf, Ks(start_svc_class)+kentry) = space_srv(en_wbuf, Ks(start_svc_class)+kentry) + 1;
                                outspace = [outspace; space_buf(en,:), space_srv(en,:), space_var(en,:)];
                                if isinf(ni)
                                    outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,end).*rate(en,:)];
                                else
                                    outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,min(ni(en),lldlimit)).*rate(en,:)];
                                end
                                outprob = [outprob; ones(size(rate(en,:),1),1)];
                            case SchedStrategy.FCFSPI % FCFS preempt-independent (no priority)
                                space_srv(en,Ks(class)+k) = space_srv(en,Ks(class)+k) - 1; % record departure
                                rate(en) = mu{ist}{class}(k)*(phi{ist}{class}(k)).*kir(:,class,k); % assume active
                                en_wbuf = en & ni>S(ist); %states with jobs in buffer
                                % FCFS: Find rightmost (last) non-zero position (longest waiting job)
                                [~, colLastNnz] = max(fliplr(space_buf(en_wbuf,:) ~= 0), [], 2);
                                colLastNnz = size(space_buf,2) - colLastNnz; % convert to actual position (odd index for class)
                                start_svc_class = space_buf(en_wbuf, colLastNnz); % job entering service
                                space_buf(en_wbuf, colLastNnz) = 0; % zero popped job
                                space_buf(en_wbuf, colLastNnz+1) = 0; % zero popped phase (ignored for FCFSPI)
                                if isempty(start_svc_class)
                                    outspace = [outspace; space_buf(en,:), space_srv(en,:), space_var(en,:)];
                                    if isinf(ni)
                                        outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,end).*rate(en,:)];
                                    else
                                        outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,min(ni(en),lldlimit)).*rate(en,:)];
                                    end
                                    outprob = [outprob; ones(size(rate(en,:),1),1)];
                                    if isSimulation && nargin>=7 && isobject(eventCache)
                                        eventCache(key) = {outprob, outspace, outrate};
                                    end
                                    return
                                end
                                % For FCFSPI, jobs restart from pie distribution instead of stored phase
                                pentry_svc_class = pie{ist}{start_svc_class};
                                for kentry = 1:K(start_svc_class)
                                    space_srv_k = space_srv;
                                    space_srv_k(en_wbuf, Ks(start_svc_class)+kentry) = space_srv_k(en_wbuf, Ks(start_svc_class)+kentry) + 1;
                                    outspace = [outspace; space_buf(en,:), space_srv_k(en,:), space_var(en,:)];
                                    rate_k = rate;
                                    rate_k(en_wbuf,:) = rate_k(en_wbuf,:) * pentry_svc_class(kentry);
                                    if isinf(ni)
                                        outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,end).*rate_k(en,:)];
                                    else
                                        outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,min(ni(en),lldlimit)).*rate_k(en,:)];
                                    end
                                    outprob = [outprob; ones(size(rate_k(en,:),1),1)];
                                end
                            case SchedStrategy.FCFSPRPRIO % FCFS preempt-resume with priority groups
                                space_srv(en,Ks(class)+k) = space_srv(en,Ks(class)+k) - 1; % record departure
                                rate(en) = mu{ist}{class}(k)*(phi{ist}{class}(k)).*kir(:,class,k); % assume active
                                en_wbuf = en & ni>S(ist); %states with jobs in buffer
                                en_wobuf = ~en_wbuf;
                                priogroup = [Inf,sn.classprio]; % Inf for empty positions (lower value = higher priority)
                                space_buf_groupg = arrayfun(@(x) priogroup(1+x), space_buf);
                                start_classprio = min(space_buf_groupg(en_wbuf,:),[],2); % min finds highest priority
                                isrowmax = space_buf_groupg == repmat(start_classprio, 1, size(space_buf_groupg,2));
                                % FCFS: Find rightmost (last) position for highest priority class (longest waiting)
                                [~,rightmostMaxPos]=max(fliplr(isrowmax),[],2);
                                rightmostMaxPos = size(space_buf_groupg,2) - rightmostMaxPos + 1;
                                start_svc_class = space_buf(en_wbuf, rightmostMaxPos); % job entering service
                                kentry = space_buf(en_wbuf, rightmostMaxPos+1); % entry phase of job resuming service (preempt-resume)

                                % Handle states without buffer jobs
                                outspace = [outspace; space_buf(en_wobuf,:), space_srv(en_wobuf,:), space_var(en_wobuf,:)];
                                if isinf(ni)
                                    outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,end).*rate(en_wobuf,:)];
                                else
                                    outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,min(ni(en),lldlimit)).*rate(en_wobuf,:)];
                                end
                                outprob = [outprob; ones(size(rate(en_wobuf,:),1),1)];

                                % Handle states with buffer jobs
                                if any(en_wbuf) && start_svc_class > 0
                                    space_srv_k = space_srv;
                                    space_buf_k = space_buf;
                                    space_srv_k(en_wbuf,Ks(start_svc_class)+kentry) = space_srv_k(en_wbuf,Ks(start_svc_class)+kentry) + 1;
                                    for j=find(en_wbuf)'
                                        % Remove both class and phase from rightmost position (preempt-resume)
                                        space_buf_k(j,:) = [0, space_buf_k(j,1:rightmostMaxPos(j)-1), space_buf_k(j,(rightmostMaxPos(j)+2):end), 0];
                                    end
                                    outspace = [outspace; space_buf_k(en_wbuf,:), space_srv_k(en_wbuf,:), space_var(en_wbuf,:)];
                                    if isinf(ni)
                                        outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,end).*rate(en_wbuf,:)];
                                    else
                                        outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,min(ni(en),lldlimit)).*rate(en_wbuf,:)];
                                    end
                                    outprob = [outprob; ones(size(rate(en_wbuf,:),1),1)];
                                end
                            case SchedStrategy.LCFSPRPRIO % LCFS preempt-resume with priority groups
                                space_srv(en,Ks(class)+k) = space_srv(en,Ks(class)+k) - 1; % record departure
                                rate(en) = mu{ist}{class}(k)*(phi{ist}{class}(k)).*kir(:,class,k); % assume active
                                en_wbuf = en & ni>S(ist); %states with jobs in buffer
                                en_wobuf = ~en_wbuf;
                                priogroup = [Inf,sn.classprio]; % Inf for empty positions (lower value = higher priority)
                                space_buf_groupg = arrayfun(@(x) priogroup(1+x), space_buf);
                                start_classprio = min(space_buf_groupg(en_wbuf,:),[],2); % min finds highest priority
                                isrowmax = space_buf_groupg == repmat(start_classprio, 1, size(space_buf_groupg,2));
                                % LCFS: Find leftmost (first) position for highest priority class (most recently preempted)
                                [~,leftmostMaxPos]=max(isrowmax,[],2);
                                start_svc_class = space_buf(en_wbuf, leftmostMaxPos); % job entering service
                                kentry = space_buf(en_wbuf, leftmostMaxPos+1); % entry phase of job resuming service (preempt-resume)
                                
                                % Handle states without buffer jobs
                                outspace = [outspace; space_buf(en_wobuf,:), space_srv(en_wobuf,:), space_var(en_wobuf,:)];
                                if isinf(ni)
                                    outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,end).*rate(en_wobuf,:)];
                                else
                                    outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,min(ni(en),lldlimit)).*rate(en_wobuf,:)];
                                end
                                outprob = [outprob; ones(size(rate(en_wobuf,:),1),1)];
                                
                                % Handle states with buffer jobs
                                if any(en_wbuf) && start_svc_class > 0
                                    space_srv_k = space_srv;
                                    space_buf_k = space_buf;
                                    space_srv_k(en_wbuf,Ks(start_svc_class)+kentry) = space_srv_k(en_wbuf,Ks(start_svc_class)+kentry) + 1;
                                    for j=find(en_wbuf)'
                                        % Remove both class and phase from leftmost position (preempt-resume)
                                        space_buf_k(j,:) = [0, space_buf_k(j,1:leftmostMaxPos(j)-1), space_buf_k(j,(leftmostMaxPos(j)+2):end), 0];
                                    end
                                    outspace = [outspace; space_buf_k(en_wbuf,:), space_srv_k(en_wbuf,:), space_var(en_wbuf,:)];
                                    if isinf(ni)
                                        outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,end).*rate(en_wbuf,:)];
                                    else
                                        outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,min(ni(en),lldlimit)).*rate(en_wbuf,:)];
                                    end
                                    outprob = [outprob; ones(size(rate(en_wbuf,:),1),1)];
                                end
                            case SchedStrategy.FCFSPIPRIO % FCFS preempt-independent with priority groups
                                space_srv(en,Ks(class)+k) = space_srv(en,Ks(class)+k) - 1; % record departure
                                rate(en) = mu{ist}{class}(k)*(phi{ist}{class}(k)).*kir(:,class,k); % assume active
                                en_wbuf = en & ni>S(ist); %states with jobs in buffer
                                en_wobuf = ~en_wbuf;
                                priogroup = [Inf,sn.classprio]; % Inf for empty positions (lower value = higher priority)
                                space_buf_groupg = arrayfun(@(x) priogroup(1+x), space_buf);
                                start_classprio = min(space_buf_groupg(en_wbuf,:),[],2); % min finds highest priority
                                isrowmax = space_buf_groupg == repmat(start_classprio, 1, size(space_buf_groupg,2));
                                % FCFS: Find rightmost (last) position for highest priority class
                                [~,rightmostMaxPos]=max(fliplr(isrowmax),[],2);
                                rightmostMaxPos = size(space_buf_groupg,2) - rightmostMaxPos + 1;
                                start_svc_class = space_buf(en_wbuf, rightmostMaxPos); % job entering service

                                % Handle states without buffer jobs
                                outspace = [outspace; space_buf(en_wobuf,:), space_srv(en_wobuf,:), space_var(en_wobuf,:)];
                                if isinf(ni)
                                    outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,end).*rate(en_wobuf,:)];
                                else
                                    outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,min(ni(en),lldlimit)).*rate(en_wobuf,:)];
                                end
                                outprob = [outprob; ones(size(rate(en_wobuf,:),1),1)];

                                % Handle states with buffer jobs
                                if any(en_wbuf) && start_svc_class > 0
                                    % For FCFSPIPRIO, jobs restart from pie distribution instead of stored phase
                                    pentry_svc_class = pie{ist}{start_svc_class};
                                    for kentry = 1:K(start_svc_class)
                                        space_srv_k = space_srv;
                                        space_buf_k = space_buf;
                                        space_srv_k(en_wbuf,Ks(start_svc_class)+kentry) = space_srv_k(en_wbuf,Ks(start_svc_class)+kentry) + 1;
                                        for j=find(en_wbuf)'
                                            % Remove both class and phase from rightmost position (preempt-independent ignores stored phase)
                                            space_buf_k(j,:) = [0, space_buf_k(j,1:rightmostMaxPos(j)-1), space_buf_k(j,(rightmostMaxPos(j)+2):end), 0];
                                        end
                                        outspace = [outspace; space_buf_k(en_wbuf,:), space_srv_k(en_wbuf,:), space_var(en_wbuf,:)];
                                        rate_k = rate;
                                        rate_k(en_wbuf,:) = rate_k(en_wbuf,:) * pentry_svc_class(kentry);
                                        if isinf(ni)
                                            outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,end).*rate_k(en_wbuf,:)];
                                        else
                                            outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,min(ni(en),lldlimit)).*rate_k(en_wbuf,:)];
                                        end
                                        outprob = [outprob; ones(size(rate_k(en_wbuf,:),1),1)];
                                    end
                                end
                            case SchedStrategy.LCFSPIPRIO % LCFS preempt-independent with priority groups
                                space_srv(en,Ks(class)+k) = space_srv(en,Ks(class)+k) - 1; % record departure
                                rate(en) = mu{ist}{class}(k)*(phi{ist}{class}(k)).*kir(:,class,k); % assume active
                                en_wbuf = en & ni>S(ist); %states with jobs in buffer
                                en_wobuf = ~en_wbuf;
                                priogroup = [Inf,sn.classprio]; % Inf for empty positions (lower value = higher priority)
                                space_buf_groupg = arrayfun(@(x) priogroup(1+x), space_buf);
                                start_classprio = min(space_buf_groupg(en_wbuf,:),[],2); % min finds highest priority
                                isrowmax = space_buf_groupg == repmat(start_classprio, 1, size(space_buf_groupg,2));
                                % LCFS: Find leftmost (first) position for highest priority class
                                [~,leftmostMaxPos]=max(isrowmax,[],2);
                                start_svc_class = space_buf(en_wbuf, leftmostMaxPos); % job entering service
                                
                                % Handle states without buffer jobs
                                outspace = [outspace; space_buf(en_wobuf,:), space_srv(en_wobuf,:), space_var(en_wobuf,:)];
                                if isinf(ni)
                                    outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,end).*rate(en_wobuf,:)];
                                else
                                    outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,min(ni(en),lldlimit)).*rate(en_wobuf,:)];
                                end
                                outprob = [outprob; ones(size(rate(en_wobuf,:),1),1)];
                                
                                % Handle states with buffer jobs
                                if any(en_wbuf) && start_svc_class > 0
                                    % For LCFSPIPRIO, jobs restart from pie distribution instead of stored phase
                                    pentry_svc_class = pie{ist}{start_svc_class};
                                    for kentry = 1:K(start_svc_class)
                                        space_srv_k = space_srv;
                                        space_buf_k = space_buf;
                                        space_srv_k(en_wbuf,Ks(start_svc_class)+kentry) = space_srv_k(en_wbuf,Ks(start_svc_class)+kentry) + 1;
                                        for j=find(en_wbuf)'
                                            % Remove both class and phase from leftmost position (preempt-independent ignores stored phase)
                                            space_buf_k(j,:) = [0, space_buf_k(j,1:leftmostMaxPos(j)-1), space_buf_k(j,(leftmostMaxPos(j)+2):end), 0];
                                        end
                                        outspace = [outspace; space_buf_k(en_wbuf,:), space_srv_k(en_wbuf,:), space_var(en_wbuf,:)];
                                        rate_k = rate;
                                        rate_k(en_wbuf,:) = rate_k(en_wbuf,:) * pentry_svc_class(kentry);
                                        if isinf(ni)
                                            outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,end).*rate_k(en_wbuf,:)];
                                        else
                                            outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,min(ni(en),lldlimit)).*rate_k(en_wbuf,:)];
                                        end
                                        outprob = [outprob; ones(size(rate_k(en_wbuf,:),1),1)];
                                    end
                                end
                            case SchedStrategy.SIRO
                                rate = zeros(size(space_srv,1),1);
                                rate(en) = mu{ist}{class}(k)*(phi{ist}{class}(k)).*kir(:,class,k); % this is for states not in en_buf
                                space_srv = inspace(:,end-sum(K)+1:end); % server state
                                space_srv(en,Ks(class)+k) = space_srv(en,Ks(class)+k) - 1; % record departure
                                % first record departure in states where the buffer is empty
                                en_wobuf = en & sum(space_buf(en,:),2) == 0;
                                outspace = [outspace; space_buf(en_wobuf,:), space_srv(en_wobuf,:), space_var(en_wobuf,:)];
                                if isinf(ni) % hit limited load-dependence
                                    outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,end).*rate(en_wobuf,:)];
                                else
                                    outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,min(ni(en),lldlimit)).*rate(en_wobuf,:)];
                                end
                                outprob = [outprob; ones(size(rate(en_wobuf,:),1),1)];
                                % let's go now to states where the buffer is non-empty
                                for r=1:R % pick a job of a random class
                                    rate_r = rate;
                                    space_buf = inspace(:,1:(end-sum(K))); % buffer state
                                    en_wbuf = en & space_buf(en,r) > 0; % states where the buffer is non-empty
                                    space_buf(en_wbuf,r) = space_buf(en_wbuf,r) - 1; % remove from buffer
                                    space_srv_r = space_srv;
                                    pentry_svc_class = pie{ist}{r};
                                    pick_prob = (nir(r)-sir(r)) / (ni-sum(sir));
                                    if pick_prob >= 0
                                        rate_r(en_wbuf,:) = rate_r(en_wbuf,:) * pick_prob;
                                    end
                                    for kentry=1:K(r)
                                        space_srv_r(en_wbuf,Ks(r)+kentry) = space_srv_r(en_wbuf,Ks(r)+kentry) + 1; % bring job in service
                                        outspace = [outspace; space_buf(en_wbuf,:), space_srv_r(en_wbuf,:), space_var(en_wbuf,:)];
                                        rate_k = rate_r;
                                        rate_k(en_wbuf,:) = rate_k(en_wbuf,:) * pentry_svc_class(kentry);
                                        if isinf(ni) % hit limited load-dependence
                                            outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,end).*rate_k(en_wbuf,:)];
                                        else
                                            outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,min(ni(en),lldlimit)).*rate_k(en_wbuf,:)];
                                        end
                                        outprob = [outprob; ones(size(rate(en_wbuf,:),1),1)];
                                        space_srv_r(en_wbuf,Ks(r)+kentry) = space_srv_r(en_wbuf,Ks(r)+kentry) - 1; % bring job in service
                                    end
                                end
                            case {SchedStrategy.SEPT,SchedStrategy.LEPT} % move last job in service
                                rate = zeros(size(space_srv,1),1);
                                rate(en) = mu{ist}{class}(k)*(phi{ist}{class}(k)).*kir(:,class,k); % this is for states not in en_buf
                                space_srv = inspace(:,end-sum(K)+1:end); % server state
                                space_srv(en,Ks(class)+k) = space_srv(en,Ks(class)+k) - 1; % record departure
                                space_buf = inspace(:,1:(end-sum(K))); % buffer state
                                % in SEPT, the scheduling parameter is the priority order of the class means
                                % en_wbuf: states where the buffer is non-empty
                                % sept_class: class to pick in service
                                [en_wbuf, first_class_inrow] = max(space_buf(:,sn.schedparam(ist,:))~=0, [], 2);
                                sept_class = sn.schedparam(ist,first_class_inrow); % this is different for sept and lept

                                space_buf(en_wbuf,sept_class) = space_buf(en_wbuf,sept_class) - 1; % remove from buffer
                                pentry = pie{ist}{sept_class};
                                for kentry=1:K(sept_class)
                                    space_srv(en_wbuf,Ks(sept_class)+kentry) = space_srv(en_wbuf,Ks(sept_class)+kentry) + 1; % bring job in service
                                    if isSimulation
                                        % break the tie
                                        outspace = [outspace; space_buf(en,:), space_srv(en,:), space_var(en,:)];
                                        rate_k = rate;
                                        rate_k(en,:) = rate_k(en,:) * pentry(kentry);
                                        if isinf(ni) % hit limited load-dependence
                                            outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,end).*rate_k(en,:)];
                                        else
                                            outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,min(ni(en),lldlimit)).*rate_k(en,:)];
                                        end
                                        outprob = [outprob; ones(size(rate(en,:),1),1)];
                                    else
                                        outspace = [outspace; space_buf(en,:), space_srv(en,:), space_var(en,:)];
                                        rate_k = rate;
                                        rate_k(en,:) = rate_k(en,:) * pentry(kentry);
                                        if isinf(ni) % hit limited load-dependence
                                            outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,end).*rate_k(en,:)];
                                        else
                                            outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,min(ni(en),lldlimit)).*rate_k(en,:)];
                                        end
                                        outprob = [outprob; ones(size(rate(en,:),1),1)];
                                    end
                                    space_srv(en_wbuf,Ks(sept_class)+kentry) = space_srv(en_wbuf,Ks(sept_class)+kentry) - 1; % bring job in service
                                end
                            otherwise
                                line_error(mfilename,sprintf('Scheduling strategy %s is not supported.', SchedStrategy.toText(sn.sched(ist))));
                        end
                    end
                end
                if isSimulation
                    if nargin>=7 && isobject(eventCache)
                        eventCache(key) = {outprob, outspace,outrate};
                    end

                    if size(outspace,1) > 1
                        tot_rate = sum(outrate);
                        cum_rate = cumsum(outrate) / tot_rate;
                        firing_ctr = 1 + max([0,find( rand > cum_rate' )]); % select action
                        outspace = outspace(firing_ctr,:);
                        outrate = sum(outrate);
                        outprob = outprob(firing_ctr,:);
                    end
                end
            end
        end
    case EventType.PHASE
        outspace = [];
        outrate = [];
        outprob = [];
        [ni,nir,~,kir] = State.toMarginal(sn,ind,inspace,K,Ks,space_buf,space_srv,space_var);
        if nir(class)>0
            for k=1:K(class)
                en = space_srv(:,Ks(class)+k) > 0;
                if any(en)
                    for kdest=setdiff(1:K(class),k) % new phase
                        rate = 0;
                        space_srv_k = space_srv(en,:);
                        space_buf_k = space_buf(en,:);
                        space_var_k = space_var(en,:);
                        if ismkvmodclass(class)
                            space_var_k(sum(sn.nvars(ind,1:class))) = kdest;
                        end
                        space_srv_k(:,Ks(class)+k) = space_srv_k(:,Ks(class)+k) - 1;
                        space_srv_k(:,Ks(class)+kdest) = space_srv_k(:,Ks(class)+kdest) + 1;
                        switch sn.sched(ist)
                            case SchedStrategy.EXT
                                rate = proc{ist}{class}{1}(k,kdest); % move next job forward
                            case SchedStrategy.INF
                                rate = proc{ist}{class}{1}(k,kdest)*kir(:,class,k); % assume active
                            case {SchedStrategy.PS, SchedStrategy.LPS}
                                rate = proc{ist}{class}{1}(k,kdest)*kir(:,class,k)./ni(:).*min(ni(:),S(ist)); % assume active
                            case SchedStrategy.PSPRIO
                                if all(ni <= S(ist)) || sn.classprio(class) == min(sn.classprio(nir>0))
                                    rate = proc{ist}{class}{1}(k,kdest)*kir(:,class,k)./ni(:).*min(ni(:),S(ist)); % assume active
                                else
                                    rate = 0; % not in most urgent priority group
                                end
                            case SchedStrategy.DPSPRIO
                                if all(ni <= S(ist)) || sn.classprio(class) == min(sn.classprio(nir>0))
                                    w_i = sn.schedparam(ist,:);
                                    w_i = w_i / sum(w_i);
                                    nirprio = nir;
                                    if ~all(ni <= S(ist))
                                        nirprio(sn.classprio~=sn.classprio(class)) = 0;
                                    end
                                    rate = proc{ist}{class}{1}(k,kdest)*kir(:,class,k)*w_i(class)./(sum(repmat(w_i,size(nirprio,1),1)*nirprio',2)); % assume active
                                else
                                    rate = 0; % not in most urgent priority group
                                end
                            case SchedStrategy.GPSPRIO
                                if all(ni <= S(ist)) || sn.classprio(class) == min(sn.classprio(nir>0))
                                    w_i = sn.schedparam(ist,:);
                                    w_i = w_i / sum(w_i);
                                    nirprio = nir;
                                    if ~all(ni <= S(ist))
                                        nirprio(sn.classprio~=sn.classprio(class)) = 0;
                                    end
                                    cir = min(nirprio,ones(size(nirprio)));
                                    rate = proc{ist}{class}{1}(k,kdest)*kir(:,class,k)/nirprio(class)*w_i(class)/(w_i*cir(:)); % assume active
                                else
                                    rate = 0; % not in most urgent priority group
                                end
                            case SchedStrategy.DPS
                                if S(ist) > 1
                                    line_error(mfilename,'Multi-server DPS not supported yet');
                                end
                                w_i = sn.schedparam(ist,:);
                                w_i = w_i / sum(w_i);
                                rate = proc{ist}{class}{1}(k,kdest)*kir(:,class,k)*w_i(class)./(sum(repmat(w_i,size(nir,1),1)*nir',2)); % assume active
                            case SchedStrategy.GPS
                                if S(ist) > 1
                                    line_error(mfilename,'Multi-server GPS not supported yet');
                                end
                                cir = min(nir,ones(size(nir)));
                                w_i = sn.schedparam(ist,:); w_i = w_i / sum(w_i);
                                rate = proc{ist}{class}{1}(k,kdest)*kir(:,class,k)/nir(class)*w_i(class)/(w_i*cir(:)); % assume active

                            case {SchedStrategy.FCFS, SchedStrategy.HOL, SchedStrategy.LCFS, SchedStrategy.LCFSPR, SchedStrategy.LCFSPI, SchedStrategy.LCFSPRIO, SchedStrategy.SIRO, SchedStrategy.SEPT, SchedStrategy.LEPT}
                                rate = proc{ist}{class}{1}(k,kdest)*kir(:,class,k); % assume active
                        end
                        % if the class cannot be served locally,
                        % then rate = NaN since mu{i,class}=NaN
                        if isinf(ni) % hit limited load-dependence
                            outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,end).*rate];
                        else
                            outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,min(ni(en),lldlimit)).*rate];
                        end
                        outspace = [outspace; space_buf_k, space_srv_k, space_var_k];
                        outprob = [outprob; ones(size(rate,1),1)];
                    end
                end
            end
            if isSimulation
                if nargin>=7 && isobject(eventCache)
                    eventCache(key) = {outprob, outspace,outrate};
                end

                if size(outspace,1) > 1
                    tot_rate = sum(outrate);
                    cum_rate = cumsum(outrate) / tot_rate;
                    firing_ctr = 1 + max([0,find( rand > cum_rate' )]); % select action
                    outspace = outspace(firing_ctr,:);
                    outrate = sum(outrate);
                    outprob = outprob(firing_ctr,:);
                end
            end
        end
end

end