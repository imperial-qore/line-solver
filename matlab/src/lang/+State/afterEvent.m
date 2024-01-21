function [outspace, outrate, outprob] =  afterEvent(sn, ind, inspace, event, class, isSimulation)
% [OUTSPACE, OUTRATE, OUTPROB] =  AFTEREVENT(QN, IND, INSPACE, EVENT, CLASS, ISSIMULATION)

% Copyright (c) 2012-2024, Imperial College London
% All rights reserved.

% outprob: event probability if the action is passive

%if ~exist('isSimulation','var')
%    isSimulation = false;
%end
%global GlobalConstants.Immediate

M = sn.nstations;
R = sn.nclasses;
S = sn.nservers;
phasessz = sn.phasessz;
phaseshift = sn.phaseshift;
pie = sn.pie;
outspace = [];
outrate = [];
outprob = 1;
% ind: node index
isf = sn.nodeToStateful(ind);

if sn.isstation(ind)
    ismkvmod = any(sn.procid(sn.nodeToStation(ind),:)==ProcessType.ID_MAP | sn.procid(sn.nodeToStation(ind),:)==ProcessType.ID_MMPP2);
    ismkvmodclass = zeros(R,1);
    for r=1:R
        ismkvmodclass(r) = any(sn.procid(sn.nodeToStation(ind),r)==ProcessType.ID_MAP | sn.procid(sn.nodeToStation(ind),r)==ProcessType.ID_MMPP2);
    end
end

lldscaling = sn.lldscaling;
if isempty(lldscaling)
    lldlimit = max(sum(sn.nclosedjobs),1);
    lldscaling = ones(M,lldlimit);
else
    lldlimit = size(lldscaling,2);
end

cdscaling = sn.cdscaling;
if isempty(cdscaling)
    cdscaling = cell(M,1);
    cdscaling{1} = @(ni) 1;
    for i=2:M % faster this way
        cdscaling{i} = cdscaling{1};
    end
end

hasOnlyExp = false; % true if all service processes are exponential
if sn.isstation(ind)
    ist = sn.nodeToStation(ind);
    K = phasessz(ist,:);
    Ks = phaseshift(ist,:);
    if max(K)==1
        hasOnlyExp = true;
    end
    mu = sn.mu;
    phi = sn.phi;
    proc = sn.proc;
    capacity = sn.cap;
    classcap = sn.classcap;
    if K(class) == 0 % if this class is not accepted at the resource
        return
    end
    V = sum(sn.nvars(ind,:));
    space_var = inspace(:,(end-V+1):end); % local state variables
    space_srv = inspace(:,(end-sum(K)-V+1):(end-V)); % server state
    space_buf = inspace(:,1:(end-sum(K)-V)); % buffer state
elseif sn.isstateful(ind)
    V = sum(sn.nvars(ind,:));
    % in this case service is always immediate so sum(K)=1
    space_var = inspace(:,(end-V+1):end); % local state variables
    space_srv = inspace(:,(end-R-V+1):(end-V)); % server state
    space_buf = []; % buffer state
else % stateless node
    space_var = [];
    space_srv = [];
    space_buf = [];
end

%switch sn.nodetype(ind)
%    case {NodeType.Queue, NodeType.Delay, NodeType.Source}
if sn.isstation(ind)
    switch event
        case EventType.ID_ARV %% passive
            % return if there is no space to accept the arrival
            [ni,nir] = State.toMarginalAggr(sn,ind,inspace,K,Ks,space_buf,space_srv,space_var);
            % otherwise check scheduling strategy
            pentry = pie{ist}{class};
            outprob = [];
            outprob_k = [];
            for kentry = 1:K(class)
                space_var_k = space_var;
                space_srv_k = space_srv;
                space_buf_k = space_buf;
                switch sn.schedid(ist)
                    case SchedStrategy.ID_EXT % source, can receive any "virtual" arrival from the sink as long as it is from an open class
                        if isinf(sn.njobs(class))
                            outspace = inspace;
                            outrate = -1*zeros(size(outspace,1)); % passive action, rate is unspecified
                            outprob = ones(size(outspace,1));
                            break
                        end
                    case {SchedStrategy.ID_PS, SchedStrategy.ID_INF, SchedStrategy.ID_DPS, SchedStrategy.ID_GPS}
                        % SPN code:
                        %                         if sn.nodetype(ind) == NodeType.Place
                        %                             indp = sn.varsparam{ind}.nodeToPlace;
                        %                             % Find all connecting transition into this place node
                        %                             ct = find(sn.rtnodes(:, ind));
                        %                             % Iterate through all transitions
                        %                             for itr=1:length(ct)
                        %                                 tr = ct(itr);
                        %                                 if sn.isstation(tr)
                        %                                     % space_srv_k(:,Ks(class)+kentry) = space_srv_k(:,Ks(class)+kentry) + 1;
                        %                                     outprob_k = pentry(kentry)*ones(size(space_srv_k,1));
                        %                                 else
                        %                                     nmodes = sn.varsparam{tr}.nmodes;
                        %                                     % Iterate through all modes of transition
                        %                                     for im=1:nmodes
                        %                                         space_srv_k(:,Ks(class)+kentry) = space_srv_k(:,Ks(class)+kentry) + sn.varsparam{tr}.forw(indp,nmodes);
                        %                                         en = classcap(ist,class)> nir(:,class) | capacity(ist)*ones(size(ni,1),1) > ni;
                        %                                         outspace = [outspace; space_srv_k(en,:)];
                        %                                         outprob_k = ones(size(space_srv_k,1));
                        %                                         outprob = [outprob; outprob_k];
                        %                                         space_srv_k(:,Ks(class)+kentry) = space_srv_k(:,Ks(class)+kentry) - sn.varsparam{tr}.forw(indp,nmodes);
                        %                                         trans = [trans, tr];
                        %                                         modes = [modes, im];
                        %                                     end
                        %                                 end
                        %                             end
                        %                             space_srv_k = [];
                        %                             outprob_k = [];
                        %                         else
                        % job enters service immediately
                        if space_srv_k(:,Ks(class)+kentry) < classcap(ist,class)
                            space_srv_k(:,Ks(class)+kentry) = space_srv_k(:,Ks(class)+kentry) + 1;
                            outprob_k = pentry(kentry)*ones(size(space_srv_k,1));
                        else
                            outprob_k = pentry(kentry)*zeros(size(space_srv_k,1));
                        end
                        %                        end
                    case {SchedStrategy.ID_SIRO, SchedStrategy.ID_SEPT, SchedStrategy.ID_LEPT}
                        if ni<S(ist)
                            space_srv_k(:,Ks(class)+kentry) = space_srv_k(:,Ks(class)+kentry) + 1;
                            outprob_k = pentry(kentry)*ones(size(space_srv_k,1));
                        else
                            space_buf_k(:,class) = space_buf_k(:,class) + 1;
                            outprob_k = pentry(kentry)*ones(size(space_srv_k,1));
                        end
                    case {SchedStrategy.ID_FCFS, SchedStrategy.ID_HOL, SchedStrategy.ID_LCFS}
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
                            %if isSimulation
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
                    case SchedStrategy.ID_LCFSPR
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
                        end

                        % if all busy, expand output states for all possible choices of job class to preempt
                        psentry = ones(size(space_buf_k_reord,1),1); % probability scaling due to preemption
                        for classpreempt = 1:R
                            for phasepreempt = 1:K(classpreempt) % phase of job to preempt
                                si_preempt = sum(space_srv_k(:, (end-sum(K)+Ks(classpreempt)+phasepreempt)),2);
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
                                        space_buf_k_preempt(sub2ind(size(space_buf_k_preempt),1:size(space_buf_k_preempt,1),empty_slots')+1) = phasepreempt;
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
        case EventType.ID_DEP
            if any(any(space_srv(:,(Ks(class)+1):(Ks(class)+K(class))))) % something is busy
                if hasOnlyExp && (sn.schedid(ist) == SchedStrategy.ID_PS || sn.schedid(ist) == SchedStrategy.ID_DPS || sn.schedid(ist) == SchedStrategy.ID_GPS || sn.schedid(ist) == SchedStrategy.ID_INF)
                    nir = space_srv;
                    ni = sum(nir,2);
                    sir = nir;
                    kir = sir;
                else
                    [ni,nir,sir,kir] = State.toMarginal(sn,ind,inspace,K,Ks,space_buf,space_srv,space_var);
                end
                switch sn.routing(ind,class)
                    case RoutingStrategy.ID_RROBIN
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
                            switch sn.schedid(ist)
                                case SchedStrategy.ID_EXT % source, can produce an arrival from phase-k as long as it is from an open class
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
                                                outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,min(ni,lldlimit)).*pentry(kentry)*mu{ist}{class}(k)*phi{ist}{class}(k)*ones(size(inspace(en,:),1),1)];
                                            end
                                            outprob = [outprob; ones(size(space_buf(en,:),1),1)];
                                        end
                                    end
                                case SchedStrategy.ID_INF % move first job in service
                                    space_srv(en,Ks(class)+k) = space_srv(en,Ks(class)+k) - 1; % record departure
                                    rate(en) = mu{ist}{class}(k)*(phi{ist}{class}(k)).*kir(en,class,k); % assume active
                                    % if state is unchanged, still add with rate 0
                                    outspace = [outspace; space_buf(en,:), space_srv(en,:), space_var(end,:)];
                                    if isinf(ni) % hit limited load-dependence
                                        outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,end).*rate(en,:)];
                                    else
                                        outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,min(ni,lldlimit)).*rate(en,:)];
                                    end
                                    outprob = [outprob; ones(size(rate(en,:),1),1)];
                                case SchedStrategy.ID_PS % move first job in service
                                    %                                     SPN code:
                                    %                                     if sn.nodetype(ind) == NodeType.Place
                                    %                                         indp = sn.varsparam{ind}.nodeToPlace;
                                    %                                         % Find all connecting transition from this place node
                                    %                                         ct = find(sn.rtnodes(ind,:));
                                    %                                         % Iterate through all transitions
                                    %                                         for itr=1:length(ct)
                                    %                                             tr = ct(itr);
                                    %                                             nmodes = sn.varsparam{tr}.nmodes;
                                    %                                             % Iterate through all modes of transition
                                    %                                             for im=1:nmodes
                                    %                                                 afterstate = space_srv(en,Ks(class)+k) - sn.varsparam{tr}.back(indp,nmodes);
                                    %                                                 if afterstate >= 0
                                    %                                                     % #input arcs x #output arcs
                                    %                                                     mult = length(find((sn.varsparam{tr}.back(:,nmodes)>0)))*length(find((sn.varsparam{tr}.forw(:,nmodes)>0)));
                                    %                                                     space_srv(en,Ks(class)+k) = afterstate;
                                    %                                                     rate(en) = 0;
                                    %                                                     % If transition is timed
                                    %                                                     if sn.varsparam{tr}.timingstrategies(im) == TimingStrategy.Timed
                                    %                                                         % mu and phi of mode
                                    %                                                         mu_l = sn.varsparam{tr}.mu{nmodes};
                                    %                                                         phi_l = sn.varsparam{tr}.phi{nmodes};
                                    %                                                         rate(en) = mu_l*(phi_l).*(kir(en,class,k)./ni(en)).*min(ni(en),S(ist))./mult;
                                    %                                                     end
                                    %                                                     outspace = [outspace; space_srv(en,:)];
                                    %                                                     outrate = [outrate; rate(en,:)];
                                    %                                                     outprob = [outprob; ones(size(rate(en,:),1),1)];
                                    %                                                     trans = [trans, tr];
                                    %                                                     modes = [modes, im];
                                    %                                                     space_srv(en,Ks(class)+k) = space_srv(en,Ks(class)+k) + sn.varsparam{tr}.back(indp,nmodes);
                                    %                                                 end
                                    %                                             end
                                    %                                         end
                                    %                                     else
                                    space_srv(en,Ks(class)+k) = space_srv(en,Ks(class)+k) - 1; % record departure
                                    rate(en) = mu{ist}{class}(k)*(phi{ist}{class}(k)).*(kir(en,class,k)./ni(en)).*min(ni(en),S(ist)); % assume active
                                    % if state is unchanged, still add with rate 0
                                    outspace = [outspace; space_buf(en,:), space_srv(en,:), space_var(en,:)];
                                    if isinf(ni) % hit limited load-dependence
                                        outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,end).*rate(en,:)];
                                    else
                                        outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,min(ni,lldlimit)).*rate(en,:)];
                                    end
                                    outprob = [outprob; ones(size(rate(en,:),1),1)];
                                    %                                    end
                                case SchedStrategy.ID_DPS
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
                                        outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,min(ni,lldlimit)).*rate(en,:)];
                                    end
                                    outprob = [outprob; ones(size(rate(en,:),1),1)];
                                case SchedStrategy.ID_GPS
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
                                        outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,min(ni,lldlimit)).*rate(en,:)];
                                    end
                                    outprob = [outprob; ones(size(rate(en,:),1),1)];
                                case SchedStrategy.ID_FCFS % move first job in service
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
                                                outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,min(ni,lldlimit)).*rate_kd(en_wobuf,:)];
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
                                                        outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,min(ni,lldlimit)).*rate_k(en,:)];
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
                                case SchedStrategy.ID_HOL % FCFS priority
                                    rate(en) = mu{ist}{class}(k)*(phi{ist}{class}(k)).*kir(:,class,k); % assume active
                                    en_wbuf = en & ni>S(ist); %states with jobs in buffer
                                    en_wobuf = ~en_wbuf;
                                    space_srv(en,Ks(class)+k) = space_srv(en,Ks(class)+k) - 1; % record departure
                                    priogroup = [0,sn.classprio];
                                    space_buf_groupg = arrayfun(@(x) priogroup(1+x), space_buf);
                                    start_classprio = max(space_buf_groupg(en_wbuf,:),[],2);
                                    isrowmax = space_buf_groupg == repmat(start_classprio, 1, size(space_buf_groupg,2));
                                    [~,rightmostMaxPosFlipped]=max(fliplr(isrowmax),[],2);
                                    rightmostMaxPos = size(isrowmax,2) - rightmostMaxPosFlipped + 1;
                                    start_svc_class = space_buf(en_wbuf, rightmostMaxPos);
                                    outspace = [outspace; space_buf(en_wobuf,:), space_srv(en_wobuf,:), space_var(en_wobuf,:)];
                                    if isinf(ni) % hit limited load-dependence
                                        outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,end).*rate(en_wobuf,:)];
                                    else
                                        outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,min(ni,lldlimit)).*rate(en_wobuf,:)];
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
                                                outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,min(ni,lldlimit)).*rate_k(en_wbuf,:)];
                                            end
                                            outprob = [outprob; ones(size(rate_k(en_wbuf,:),1),1)];
                                        end
                                    end
                                case SchedStrategy.ID_LCFS % move last job in service
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
                                            outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,min(ni,lldlimit)).*rate(en,:)];
                                        end
                                        outprob = [outprob; ones(size(rate(en,:),1),1)];
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
                                            outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,min(ni,lldlimit)).*rate_k(en,:)];
                                        end
                                        outprob = [outprob; ones(size(rate(en,:),1),1)];
                                        space_srv(en_wbuf,Ks(start_svc_class)+kentry) = space_srv(en_wbuf,Ks(start_svc_class)+kentry) - 1;
                                    end
                                case SchedStrategy.ID_LCFSPR % move last job in service
                                    space_srv(en,Ks(class)+k) = space_srv(en,Ks(class)+k) - 1; % record departure
                                    rate(en) = mu{ist}{class}(k)*(phi{ist}{class}(k)).*kir(:,class,k); % assume active
                                    en_wbuf = en & ni>S(ist); %states with jobs in buffer
                                    [~, colfirstnnz] = max( space_buf(en_wbuf,:) ~=0, [], 2 ); % find first nnz column
                                    start_svc_class = space_buf(en_wbuf,colfirstnnz); % job entering service
                                    kentry = space_buf(en_wbuf,colfirstnnz+1); % entry phase of job starting or resuming service
                                    space_buf(en_wbuf,colfirstnnz)=0;% zero popped job
                                    space_buf(en_wbuf,colfirstnnz+1)=0; % zero popped phase
                                    if isempty(start_svc_class)
                                        outspace = [outspace; space_buf(en,:), space_srv(en,:), space_var(en,:)];
                                        if isinf(ni) % hit limited load-dependence
                                            outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,end).*rate(en,:)];
                                        else
                                            outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,min(ni,lldlimit)).*rate(en,:)];
                                        end
                                        outprob = [outprob; ones(size(rate(en,:),1),1)];
                                        return
                                    end
                                    space_srv(en_wbuf,Ks(start_svc_class)+kentry) = space_srv(en_wbuf,Ks(start_svc_class)+kentry) + 1;
                                    % if state is unchanged, still add with rate 0
                                    outspace = [outspace; space_buf(en,:), space_srv(en,:), space_var(en,:)];
                                    if isinf(ni) % hit limited load-dependence
                                        outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,end).*rate(en,:)];
                                    else
                                        outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,min(ni,lldlimit)).*rate(en,:)];
                                    end
                                    outprob = [outprob; ones(size(rate(en,:),1),1)];
                                case SchedStrategy.ID_SIRO
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
                                        outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,min(ni,lldlimit)).*rate(en_wobuf,:)];
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
                                                outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,min(ni,lldlimit)).*rate_k(en_wbuf,:)];
                                            end
                                            outprob = [outprob; ones(size(rate(en_wbuf,:),1),1)];
                                            space_srv_r(en_wbuf,Ks(r)+kentry) = space_srv_r(en_wbuf,Ks(r)+kentry) - 1; % bring job in service
                                        end
                                    end
                                case {SchedStrategy.ID_SEPT,SchedStrategy.ID_LEPT} % move last job in service
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
                                                outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,min(ni,lldlimit)).*rate_k(en,:)];
                                            end
                                            outprob = [outprob; ones(size(rate(en,:),1),1)];
                                        else
                                            outspace = [outspace; space_buf(en,:), space_srv(en,:), space_var(en,:)];
                                            rate_k = rate;
                                            rate_k(en,:) = rate_k(en,:) * pentry(kentry);
                                            if isinf(ni) % hit limited load-dependence
                                                outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,end).*rate_k(en,:)];
                                            else
                                                outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,min(ni,lldlimit)).*rate_k(en,:)];
                                            end
                                            outprob = [outprob; ones(size(rate(en,:),1),1)];
                                        end
                                        space_srv(en_wbuf,Ks(sept_class)+kentry) = space_srv(en_wbuf,Ks(sept_class)+kentry) - 1; % bring job in service
                                    end
                                otherwise
                                    line_error(mfilename,sprintf('Scheduling strategy %s is not supported.', sn.sched{ist}));
                            end
                        end
                    end
                    if isSimulation
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
        case EventType.ID_PHASE
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
                            switch sn.schedid(ist)
                                case SchedStrategy.ID_EXT
                                    rate = proc{ist}{class}{1}(k,kdest); % move next job forward
                                case SchedStrategy.ID_INF
                                    rate = proc{ist}{class}{1}(k,kdest)*kir(:,class,k); % assume active
                                case SchedStrategy.ID_PS
                                    rate = proc{ist}{class}{1}(k,kdest)*kir(:,class,k)./ni(:).*min(ni(:),S(ist)); % assume active
                                case SchedStrategy.ID_DPS
                                    if S(ist) > 1
                                        line_error(mfilename,'Multi-server DPS not supported yet');
                                    end
                                    w_i = sn.schedparam(ist,:);
                                    w_i = w_i / sum(w_i);
                                    rate = proc{ist}{class}{1}(k,kdest)*kir(:,class,k)*w_i(class)./(sum(repmat(w_i,size(nir,1),1)*nir',2)); % assume active
                                case SchedStrategy.ID_GPS
                                    if S(ist) > 1
                                        line_error(mfilename,'Multi-server GPS not supported yet');
                                    end
                                    cir = min(nir,ones(size(nir)));
                                    w_i = sn.schedparam(ist,:); w_i = w_i / sum(w_i);
                                    rate = proc{ist}{class}{1}(k,kdest)*kir(:,class,k)/nir(class)*w_i(class)/(w_i*cir(:)); % assume active

                                case {SchedStrategy.ID_FCFS, SchedStrategy.ID_HOL, SchedStrategy.ID_LCFS, SchedStrategy.ID_LCFSPR, SchedStrategy.ID_SIRO, SchedStrategy.ID_SEPT, SchedStrategy.ID_LEPT}
                                    rate = proc{ist}{class}{1}(k,kdest)*kir(:,class,k); % assume active
                            end
                            % if the class cannot be served locally,
                            % then rate = NaN since mu{i,class}=NaN
                            if isinf(ni) % hit limited load-dependence
                                outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,end).*rate];
                            else
                                outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,min(ni,lldlimit)).*rate];
                            end
                            outspace = [outspace; space_buf_k, space_srv_k, space_var_k];
                            outprob = [outprob; ones(size(rate,1),1)];
                        end
                    end
                end
                if isSimulation
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
elseif sn.isstateful(ind)
    switch sn.nodetype(ind)
        case NodeType.Router
            switch event
                case EventType.ID_ARV
                    space_srv(:,class) = space_srv(:,class) + 1;
                    outspace = [space_srv, space_var]; % buf is empty
                    outrate = -1*ones(size(outspace,1)); % passive action, rate is unspecified
                case EventType.ID_DEP
                    if space_srv(class)>0
                        space_srv(:,class) = space_srv(:,class) - 1;
                        switch sn.routing(ind,class)
                            case RoutingStrategy.ID_RROBIN
                                idx = find(space_var(sum(sn.nvars(ind,1:(R+class)))) == sn.nodeparam{ind}{class}.outlinks);
                                if idx < length(sn.nodeparam{ind}{class}.outlinks)
                                    space_var(sum(sn.nvars(ind,1:(R+class)))) = sn.nodeparam{ind}{class}.outlinks(idx+1);
                                else
                                    space_var(sum(sn.nvars(ind,1:(R+class)))) = sn.nodeparam{ind}{class}.outlinks(1);
                                end
                        end
                        outspace = [space_srv, space_var]; % buf is empty
                        outrate = GlobalConstants.Immediate*ones(size(outspace,1)); % immediate action
                    end
            end
        case NodeType.Cache
            % job arrives in class, then reads and moves into hit or miss
            % class, then departs
            switch event
                case EventType.ID_ARV
                    space_srv(:,class) = space_srv(:,class) + 1;
                    outspace = [space_srv, space_var]; % buf is empty
                    outrate = -1*ones(size(outspace,1)); % passive action, rate is unspecified
                case EventType.ID_DEP
                    if space_srv(class)>0
                        space_srv(:,class) = space_srv(:,class) - 1;
                        switch sn.routing(ind,class)
                            case RoutingStrategy.ID_RROBIN
                                idx = find(space_var(sum(sn.nvars(ind,1:(R+class)))) == sn.nodeparam{ind}{class}.outlinks);
                                if idx < length(sn.nodeparam{ind}{class}.outlinks)
                                    space_var(sum(sn.nvars(ind,1:(R+class)))) = sn.nodeparam{ind}{class}.outlinks(idx+1);
                                else
                                    space_var(sum(sn.nvars(ind,1:(R+class)))) = sn.nodeparam{ind}{class}.outlinks(1);
                                end
                        end
                        outspace = [space_srv, space_var]; % buf is empty
                        outrate = GlobalConstants.Immediate*ones(size(outspace,1)); % immediate action
                    end
                case EventType.ID_READ
                    n = sn.nodeparam{ind}.nitems; % n items
                    m = sn.nodeparam{ind}.itemcap; % capacity
                    ac = sn.nodeparam{ind}.accost; % access cost
                    hitclass = sn.nodeparam{ind}.hitclass;
                    missclass = sn.nodeparam{ind}.missclass;
                    h = length(m);
                    replacement_id = sn.nodeparam{ind}.replacement;
                    if space_srv(class)>0 && sum(space_srv)==1 % is a job of class is in
                        p = sn.nodeparam{ind}.pread{class};
                        en =  space_srv(:,class) > 0;
                        space_srv_k = [];
                        space_var_k = [];
                        outrate = [];
                        if any(en)
                            for e=find(en)'
                                if isSimulation
                                    % pick one item
                                    kset = 1 + max([0,find( rand > cumsum(p) )]);
                                    % pick one entry list
                                    l = 1 + max([0,find( rand > cumsum(ac{class,kset}(1,:)) )]);
                                    outprob = ac{class,kset}(1,l) * p(kset);
                                else
                                    kset = 1:n;
                                end
                                for k=kset % request to item k
                                    space_srv_e = space_srv(e,:);
                                    space_srv_e(class) = space_srv_e(class) - 1;
                                    var = space_var(e,:);
                                    posk = find(k==var);

                                    if isempty(posk) % CACHE MISS, add to list 1
                                        space_srv_e(missclass(class)) = space_srv_e(missclass(class)) + 1;
                                        switch replacement_id
                                            case {ReplacementStrategy.ID_FIFO, ReplacementStrategy.ID_LRU, ReplacementStrategy.ID_SFIFO}
                                                varp = var;
                                                varp(2:m(1)) = var(1:(m(1)-1));
                                                varp(1) = k;
                                                space_srv_k = [space_srv_k; space_srv_e];
                                                space_var_k = [space_var_k; varp];
                                                if isSimulation
                                                    %% no p(k) weighting since that goes in the outprob vec
                                                    outrate(end+1,1) = GlobalConstants.Immediate;
                                                else
                                                    for l=2
                                                        outrate(end+1,1) = ac{class,k}(1,l) * p(k) * GlobalConstants.Immediate;
                                                    end
                                                end
                                            case ReplacementStrategy.ID_RR
                                                if isSimulation
                                                    varp = var;
                                                    r = randi(m(1),1,1);
                                                    varp(r) = k;
                                                    space_srv_k = [space_srv_k; space_srv_e];
                                                    space_var_k = [space_var_k; (varp)];
                                                    outrate(end+1,1) = GlobalConstants.Immediate;
                                                else
                                                    for r=1:m(1) % random position in list 1
                                                        varp = var;
                                                        varp(r) = k;
                                                        space_srv_k = [space_srv_k; space_srv_e];
                                                        space_var_k = [space_var_k; (varp)];
                                                        for l=2
                                                            outrate(end+1,1) = ac{class,k}(1,l) * p(k)/m(1) * GlobalConstants.Immediate;
                                                        end
                                                    end
                                                end
                                        end
                                    elseif posk <= sum(m(1:h-1)) % CACHE HIT in list i < h, move to list i+1
                                        space_srv_e(hitclass(class)) = space_srv_e(hitclass(class)) + 1;
                                        i = min(find(posk <= cumsum(m)));
                                        j = posk - sum(m(1:i-1));

                                        switch replacement_id
                                            case ReplacementStrategy.ID_FIFO
                                                if isSimulation
                                                    varp = var;
                                                    inew = i-1+probchoose(ac{class,k}(i,i:end)/sum(ac{class,k}(i,i:end)))-1; % can choose i
                                                    if inew~=i
                                                        varp(cpos(i,j)) = var(cpos(inew,m(inew)));
                                                        varp(cpos(inew,2):cpos(inew,m(inew))) = var(cpos(inew,1):cpos(inew,m(inew)-1));
                                                        varp(cpos(inew,1)) = k;
                                                    end
                                                    %varp(cpos(i,j)) = var(cpos(i+1,m(i+1)));
                                                    %varp(cpos(i+1,2):cpos(i+1,m(i+1))) = var(cpos(i+1,1):cpos(i+1,m(i+1)-1));
                                                    %varp(cpos(i+1,1)) = k;

                                                    space_srv_k = [space_srv_k; space_srv_e];
                                                    space_var_k = [space_var_k; varp];
                                                    outrate(end+1,1) = GlobalConstants.Immediate;
                                                else
                                                    for inew = i:h
                                                        varp = var;
                                                        varp(cpos(i,j)) = var(cpos(inew,m(inew)));
                                                        varp(cpos(inew,2):cpos(inew,m(inew))) = var(cpos(inew,1):cpos(inew,m(inew)-1));
                                                        varp(cpos(inew,1)) = k;
                                                        space_srv_k = [space_srv_k; space_srv_e];
                                                        space_var_k = [space_var_k; varp];
                                                        outrate(end+1,1) = ac{class,k}(1+i,1+inew) * p(k) * GlobalConstants.Immediate;
                                                    end
                                                end
                                            case ReplacementStrategy.ID_RR
                                                if isSimulation
                                                    inew = i-1+probchoose(ac{class,k}(i,i:end)/sum(ac{class,k}(i,i:end)))-1; % can choose i
                                                    varp = var;
                                                    r = randi(m(inew),1,1);
                                                    varp(cpos(i,j)) = var(cpos(inew,r));
                                                    varp(cpos(inew,r)) = k;
                                                    outprob = outprob * ac{class,k}(1+i,1+inew);
                                                    space_srv_k = [space_srv_k; space_srv_e];
                                                    space_var_k = [space_var_k; varp];
                                                    outrate(end+1,1) = GlobalConstants.Immediate;
                                                else
                                                    for inew = i:h
                                                        for r=1:m(inew) % random position in new list
                                                            varp = var;
                                                            varp(cpos(i,j)) = var(cpos(inew,r));
                                                            varp(cpos(inew,r)) = k;
                                                            space_srv_k = [space_srv_k; space_srv_e];
                                                            space_var_k = [space_var_k; varp];
                                                            outrate(end+1,1) = ac{class,k}(1+i,1+inew) * p(k)/m(inew) * GlobalConstants.Immediate;
                                                        end
                                                    end
                                                end
                                            case {ReplacementStrategy.ID_LRU, ReplacementStrategy.ID_SFIFO}
                                                if isSimulation
                                                    varp = var;
                                                    inew = i-1+ probchoose(ac{class,k}(i,i:end)/sum(ac{class,k}(i,i:end)))-1; % can choose i
                                                    varp(cpos(i,2):cpos(i,j)) = var(cpos(i,1):cpos(i,j-1));
                                                    varp(cpos(i,1)) = var(cpos(inew,m(inew)));
                                                    varp(cpos(inew,2):cpos(inew,m(inew))) = var(cpos(inew,1):cpos(inew,m(inew)-1));
                                                    varp(cpos(inew,1)) = k;
                                                    space_srv_k = [space_srv_k; space_srv_e];
                                                    space_var_k = [space_var_k; varp];
                                                    outrate(end+1,1) = GlobalConstants.Immediate;
                                                else
                                                    for inew = i:h
                                                        varp = var;
                                                        varp(cpos(i,2):cpos(i,j)) = var(cpos(i,1):cpos(i,j-1));
                                                        varp(cpos(i,1)) = var(cpos(inew,m(inew)));
                                                        varp(cpos(inew,2):cpos(inew,m(inew))) = var(cpos(inew,1):cpos(inew,m(inew)-1));
                                                        varp(cpos(inew,1)) = k;
                                                        space_srv_k = [space_srv_k; space_srv_e];
                                                        space_var_k = [space_var_k; varp];
                                                        outrate(end+1,1) = ac{class,k}(1+i,1+inew) * p(k) * GlobalConstants.Immediate;
                                                    end
                                                end
                                        end
                                    else % CACHE HIT in list h
                                        space_srv_e(hitclass(class)) = space_srv_e(hitclass(class)) + 1;
                                        i=h;
                                        j = posk - sum(m(1:i-1));
                                        switch replacement_id
                                            case ReplacementStrategy.ID_RR
                                                space_srv_k = [space_srv_k; space_srv_e];
                                                space_var_k = [space_var_k; (var)];
                                                if isSimulation
                                                    outrate(end+1,1) = GlobalConstants.Immediate;
                                                else
                                                    outrate(end+1,1) = ac{class,k}(1+h,1+h) * p(k) * GlobalConstants.Immediate;
                                                end
                                            case {ReplacementStrategy.ID_FIFO, ReplacementStrategy.ID_SFIFO}
                                                space_srv_k = [space_srv_k; space_srv_e];
                                                space_var_k = [space_var_k; var];
                                                if isSimulation
                                                    outrate(end+1,1) = GlobalConstants.Immediate;
                                                else
                                                    outrate(end+1,1) = ac{class,k}(1+h,1+h) * p(k) * GlobalConstants.Immediate;
                                                end
                                            case ReplacementStrategy.ID_LRU
                                                varp = var;
                                                varp(cpos(h,2):cpos(h,j)) = var(cpos(h,1):cpos(h,j-1));
                                                varp(cpos(h,1)) = var(cpos(h,j));
                                                space_srv_k = [space_srv_k; space_srv_e];
                                                space_var_k = [space_var_k; varp];
                                                if isSimulation
                                                    outrate(end+1,1) = GlobalConstants.Immediate;
                                                else
                                                    outrate(end+1,1) = ac{class,k}(1+h,1+h) * p(k) * GlobalConstants.Immediate;
                                                end
                                        end
                                    end
                                end
                            end
                            % if state is unchanged, still add with rate 0
                            outspace = [space_srv_k, space_var_k];
                        end
                    end
            end
    end % switch nodeType
end
    function pos = cpos(i,j)
        % POS = CPOS(I,J)

        pos = sum(m(1:i-1)) + j;
    end
end
