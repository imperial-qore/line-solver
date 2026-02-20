function [pi,SSq,arvRates,depRates,tranSysState,tranSync,sn]=solver_ssa(sn, init_state, options, eventCache)
% [PI,SSQ,ARVRATES,DEPRATES,TRANSYSSTATE,QN]=SOLVER_SSA(QN,OPTIONS)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

% by default the jobs are all initialized in the first valid state

if ~isfield(options,'seed')
    options.seed = 23000;
end
% Handle parallel computing toolbox gracefully - get worker index
if isMATLABReleaseOlderThan("R2022b")
    % Use labindex for older MATLAB versions
    try
        if ~isempty(getCurrentTask())
            lab_idx = labindex();  %#ok<DLOBIX>
        else
            lab_idx = 1;
        end
    catch
        line_warning(mfilename,'Parallel Computing Toolbox not available or not running in parallel mode. Using labindex = 1.');
        lab_idx = 1;
    end
else
    % Use spmdIndex for R2022b and newer (labindex is deprecated)
    try
        lab_idx = spmdIndex;
        if isempty(lab_idx) || lab_idx == 0
            lab_idx = 1;
        end
    catch
        line_warning(mfilename,'Parallel Computing Toolbox not available or not running in parallel mode. Using labindex = 1.');
        lab_idx = 1;
    end
end
Solver.resetRandomGeneratorSeed(options.seed + lab_idx - 1);

%% generate local state spaces
%nstations = sn.nstations;
nstateful = sn.nstateful;
%init_nserver = sn.nservers; % restore Inf at delay nodes
R = sn.nclasses;
N = sn.njobs';
nnodes = sn.nnodes;
sync = sn.sync;
gsync = sn.gsync;

line_debug('SSA solver starting: nstateful=%d, nclasses=%d, njobs=%s, samples=%d', nstateful, R, mat2str(N), options.samples);
csmask = sn.csmask;

cutoff = options.cutoff;
if isscalar(cutoff)
    cutoff = cutoff * ones(sn.nstations, sn.nclasses);
end

%%
Np = N';
capacityc = zeros(sn.nnodes, sn.nclasses);
original_classcap = sn.classcap; % preserve original classcap for class switching scenarios
for ind=1:sn.nnodes
    if sn.isstation(ind) % place jobs across stations
        ist = sn.nodeToStation(ind);
        %isf = sn.nodeToStateful(ind);
        for r=1:sn.nclasses %cut-off open classes to finite capacity
            c = find(sn.chains(:,r));
            % Check if visits is 0, but also preserve capacity for classes that can
            % receive jobs via class switching (indicated by non-zero original classcap)
            if ~isempty(sn.visits{c}) && sn.visits{c}(ist,r) == 0 && original_classcap(ist,r) == 0
                capacityc(ind,r) = 0;
            elseif ~isempty(sn.proc) && ~isempty(sn.proc{ist}{r}) && any(any(isnan(sn.proc{ist}{r}{1}))) && sn.nodetype(ind) ~= NodeType.Place % disabled (but not Place nodes)
                capacityc(ind,r) = 0;
            else
                if isinf(N(r))
                    capacityc(ind,r) =  min(cutoff(ist,r), sn.classcap(ist,r));
                else
                    capacityc(ind,r) =  sum(sn.njobs(sn.chains(c,:)));
                end
            end
        end
        if isinf(sn.nservers(ist))
            sn.nservers(ist) = sum(capacityc(ind,:));
        end
        sn.cap(ist,:) = sum(capacityc(ind,:));
        sn.classcap(ist,:) = capacityc(ind,:);
    end
end

%%
if any(isinf(Np))
    Np(isinf(Np)) = 0;
end

init_state_hashed = ones(1,nstateful); % pick the first state in init_state{i}

%%
arvRatesSamples = zeros(options.samples,nstateful,R);
depRatesSamples = zeros(options.samples,nstateful,R);
A = length(sync);
G = length(gsync);
samples_collected = 1;
nir = {};
% fill stateCell with initial states
cur_state = cell(nstateful,1); % cell array with current stateful node states
for ind=1:sn.nnodes
    if sn.isstateful(ind)
        isf = sn.nodeToStateful(ind);
        cur_state{isf} = init_state{isf}(init_state_hashed(isf),:);
        if sn.isstation(ind)
            ist = sn.nodeToStation(ind);
            [~,nir{ist}] = State.toMarginal(sn, ind, init_state{isf}(init_state_hashed(isf),:));
            nir{ist} = nir{ist}(:);
        end
    end
end
cur_state_1 = cur_state;
% generate state vector
state = cell2mat(cur_state');
% create function to determine lengths of stateful node states
statelen = cellfun(@length, cur_state);
% data structures to save transient information - pre-allocate for all samples
nSamples = options.samples;
tranSync = zeros(nSamples,1);
tranState = zeros(1+length(state), nSamples);
tranState(1:(1+length(state)),1) = [0, state]';
SSq = zeros(length(cell2mat(nir')), nSamples);
SSq(:,1) = cell2mat(nir');
local = sn.nnodes+1;
last_node_a = 0; % active in the last occurred synchronization
last_node_p = 0; % passive in the last occurred synchronization
for act=1:A
    node_a{act} = sync{act}.active{1}.node;
    node_p{act} = sync{act}.passive{1}.node;
    class_a{act} = sync{act}.active{1}.class;
    class_p{act} = sync{act}.passive{1}.class;
    event_a{act} = sync{act}.active{1}.event;
    event_p{act} = sync{act}.passive{1}.event;
    outprob_a{act} = [];
    outprob_p{act} = [];
end
enabled_next_states = cell(1,A);

%% Start main simulation loop
isSimulation = true; % allow state vector to grow, e.g. for FCFS buffers
samples_collected = 1;
cur_time = 0;
use_inline = true; % true = stable version, false = dev version

try
    while samples_collected < options.samples && cur_time <= options.timespan(2)
        %% This section corresponds to solver_ssa_findenabled in Java
        %% Inlined for performance reasons
        if use_inline
            enabled_sync = []; % row is action label, col1=rate, col2=new state
            enabled_rates = [];
            ctr = 1;
            A = length(sync);
            G = length(gsync);
            for act=1:A
                isf_a = sn.nodeToStateful(node_a{act});

                enabled_next_states{act} = cur_state;
                update_cond_a = true;
                if update_cond_a
                    [enabled_next_states{act}{isf_a}, rate_a{act}, outprob_a{act}, eventCache] =  State.afterEvent(sn, node_a{act}, cur_state{isf_a}, event_a{act}, class_a{act}, isSimulation, eventCache);
                end

                if isempty(enabled_next_states{act}{isf_a}) || isempty(rate_a{act})
                    continue
                end

                for ia=1:size(enabled_next_states{act}{isf_a},1) % for all possible new states, check if they are enabled
                    % if the transition cannot occur
                    if isnan(rate_a{act}(ia)) || rate_a{act}(ia) == 0 % handles degenerate rate values
                        % set the transition with a zero rate so that it is
                        % never selected
                        rate_a{act}(ia) = 1e-38; % ~ zero in 32-bit precision
                    end

                    if enabled_next_states{act}{isf_a}(ia,:) == -1 % hash not found
                        continue
                    end
                    update_cond_p = true; %samples_collected == 1 || ((node_p{act} == last_node_a || node_p{act} == last_node_p)) || isempty(outprob_a{act}) || isempty(outprob_p{act});

                    if rate_a{act}(ia)>0
                        if node_p{act} ~= local
                            if node_p{act} == node_a{act} %self-loop, active and passive are the same
                                isf_p = isf_a;
                                if update_cond_p
                                    [enabled_next_states{act}{isf_p}, ~, outprob_p{act}, eventCache] =  State.afterEvent(sn, node_p{act}, enabled_next_states{act}{isf_p}, event_p{act}, class_p{act}, isSimulation, eventCache);
                                end
                            else % departure
                                isf_p = sn.nodeToStateful(node_p{act});
                                if update_cond_p
                                    [enabled_next_states{act}{isf_p}, ~, outprob_p{act}, eventCache] =  State.afterEvent(sn, node_p{act}, enabled_next_states{act}{isf_p}, event_p{act}, class_p{act}, isSimulation, eventCache);
                                end
                            end
                            if ~isempty(enabled_next_states{act}{isf_p})
                                if sn.isstatedep(node_a{act},3)
                                    prob_sync_p{act} = sync{act}.passive{1}.prob(cur_state, enabled_next_states{act}); %state-dependent
                                else
                                    prob_sync_p{act} = sync{act}.passive{1}.prob;
                                end
                            else
                                prob_sync_p{act} = 0;
                            end
                        end
                        if ~isempty(enabled_next_states{act}{isf_a})
                            if node_p{act} == local
                                prob_sync_p{act} = 1;
                            end
                            if ~isnan(rate_a{act})
                                if all(~cellfun(@isempty,enabled_next_states{act}))
                                    if event_a{act} == EventType.DEP
                                        node_a_sf{act} = isf_a;
                                        node_p_sf{act} = isf_p;
                                        depRatesSamples(samples_collected,node_a_sf{act},class_a{act}) = depRatesSamples(samples_collected,node_a_sf{act},class_a{act}) + outprob_a{act} * outprob_p{act} * rate_a{act}(ia) * prob_sync_p{act};
                                        arvRatesSamples(samples_collected,node_p_sf{act},class_p{act}) = arvRatesSamples(samples_collected,node_p_sf{act},class_p{act}) + outprob_a{act} * outprob_p{act} * rate_a{act}(ia) * prob_sync_p{act};
                                    end
                                    % simulate also self-loops as we need to log them
                                    %if any(~cellfun(@isequal,new_state{act},cur_state))
                                    if node_p{act} < local && ~sn.csmask(class_a{act}, class_p{act}) && sn.nodetype(node_p{act})~=NodeType.Source && (rate_a{act}(ia) * prob_sync_p{act} >0)
                                        line_error(mfilename,sprintf('Error: state-dependent routing at node %d (%s) violates the class switching mask (node %d -> node %d, class %d -> class %d).', node_a{act}, sn.nodenames{node_a{act}}, node_a{act}, node_p{act}, class_a{act}, class_p{act}));
                                    end
                                    enabled_rates(ctr) = rate_a{act}(ia) * prob_sync_p{act};
                                    enabled_sync(ctr) = act;
                                    ctr = ctr + 1;
                                end
                            end
                        end
                    end
                end
            end
            gctr_start = ctr;

            for gact=1:G % event at node ind with global side-effects
                gind = gsync{gact}.active{1}.node; % node index for global event
                [enabled_next_states{A+gact}, outrate, outprob] = State.afterGlobalEvent(sn, gind, cur_state, gsync{gact}, isSimulation);
                for ia=find(outrate .* outprob)
                    enabled_rates(ctr) = outrate(ia) * outprob(ia);
                    enabled_sync(ctr) = A+gact;
                    ctr = ctr + 1;

                    % Record departure/arrival rates for FIRE events at Places
                    if gsync{gact}.active{1}.event == EventType.FIRE
                        mode = gsync{gact}.active{1}.mode;
                        % Get enabling/firing conditions to determine affected classes
                        enabling_m = sn.nodeparam{gind}.enabling{mode};
                        firing_m = sn.nodeparam{gind}.firing{mode};

                        for j=1:length(gsync{gact}.passive)
                            pev = gsync{gact}.passive{j};
                            % Decode linear index to (node, class) - pev.node is a linear index from find() on enabling/firing matrix
                            [pev_node, pev_class] = ind2sub([sn.nnodes, R], pev.node);
                            if pev.event == EventType.PRE
                                % Departure from input Place (consuming tokens)
                                if pev_node <= length(sn.nodeToStateful) && ~isnan(sn.nodeToStateful(pev_node)) && sn.nodeToStateful(pev_node) > 0
                                    ep_isf = sn.nodeToStateful(pev_node);
                                    % Record departures for the specific class from this PRE event
                                    depRatesSamples(samples_collected, ep_isf, pev_class) = ...
                                        depRatesSamples(samples_collected, ep_isf, pev_class) + outrate(ia) * outprob(ia);
                                end
                            elseif pev.event == EventType.POST
                                % Arrival at output Place (producing tokens)
                                if pev_node <= length(sn.nodeToStateful) && ~isnan(sn.nodeToStateful(pev_node)) && sn.nodeToStateful(pev_node) > 0
                                    fp_isf = sn.nodeToStateful(pev_node);
                                    % Record arrivals for the specific class from this POST event
                                    arvRatesSamples(samples_collected, fp_isf, pev_class) = ...
                                        arvRatesSamples(samples_collected, fp_isf, pev_class) + outrate(ia) * outprob(ia);
                                end
                            end
                        end
                    end
                end
            end
        else
            [enabled_next_states,enabled_rates,enabled_sync,gctr_start,depRatesSamples,arvRatesSamples,outprob_a,outprob_p,rate_a, eventCache] = solver_ssa_findenabled(sn,node_a,enabled_next_states,cur_state,outprob_a,event_a,class_a,isSimulation,node_p,local,outprob_p,event_p,class_p,sync,gsync,depRatesSamples,samples_collected,arvRatesSamples,last_node_a,last_node_p,eventCache);
        end
        %% Gillespie direct method
        tot_rate = sum(enabled_rates);
        cum_rate = cumsum(enabled_rates) / tot_rate;
        selected_transition = 1 + max([0,find( rand > cum_rate )]); % select action

        % Update record of last active/passive pair
        if isempty(enabled_sync)
            line_error(mfilename,'SSA simulation entered a deadlock before collecting all samples, no synchronization is enabled.');
        end
        if selected_transition < gctr_start
            % regular event pair
            last_node_a = node_a{enabled_sync(selected_transition)};
            last_node_p = node_p{enabled_sync(selected_transition)};
        else % global event
            last_node_a = NaN;
            last_node_p = NaN;
        end

        %% Update paddings
        % This part is needed to ensure that when the state vector grows the
        % padding of zero is done on the left (e.g., for FCFS buffers)
        for ind=1:sn.nnodes
            if sn.isstation(ind)
                isf = sn.nodeToStateful(ind);
                deltalen = length(cur_state{isf}) - statelen(isf);
                if deltalen>0
                    statelen(isf) = length(cur_state{isf});
                    % here do padding
                    if ind==1
                        shift = 0;
                    else
                        shift = sum(statelen(1:isf-1));
                    end
                    pad = zeros(deltalen, size(tranState,2));
                    tranState = [tranState(1:(shift+1), :); pad ; tranState((shift+1+deltalen):end, :)];
                end
            end
        end

        %% Simulate the time increment
        state = cell2mat(cur_state');
        dt = -(log(rand)/tot_rate);
        cur_time = cur_time + dt;

        %% Save simulation output data
        tranState(1:(1+length(state)),samples_collected) = [dt, state]';
        tranSync(samples_collected,1) = enabled_sync(selected_transition);
        for ind=1:sn.nnodes
            if sn.isstation(ind)
                isf = sn.nodeToStateful(ind);
                ist = sn.nodeToStation(ind);
                [~,nir{ist}] = State.toMarginal(sn, ind, cur_state{isf});
                nir{ist}=nir{ist}(:);
            end
        end
        SSq(:,samples_collected) = cell2mat(nir');

        %% Update current state and sample counter
        cur_state_1 = cur_state;
        cur_state = enabled_next_states{enabled_sync(selected_transition)};

        samples_collected = samples_collected + 1;

        %% Print progress
        print_progress(options,samples_collected);
    end
    % Print newline after progress counter
    if options.verbose
        line_printf('\n');
    end
catch ME
    getReport(ME)
end

% Trim pre-allocated arrays to actual number of samples collected
samples_collected = samples_collected - 1;  % Adjust for the increment at end of loop
tranState = tranState(:, 1:samples_collected);
tranSync = tranSync(1:samples_collected, :);
SSq = SSq(:, 1:samples_collected);

%transient = min([floor(samples_collected/10),1000]); % remove first part of simulation (10% of the samples up to 1000 max)
%transient = 0;
%output = output((transient+1):end,:);
tranState = tranState';


[u,ui,uj] = unique(tranState(:,2:end),'rows');
statesz = cellfun(@length, cur_state_1)';
tranSysState = cell(1,length(cur_state)+1);
tranSysState{1} = cumsum(tranState(:,1));
for j=1:length(statesz)
    tranSysState{1+j} = tranState(:,1+(1+sum(statesz(1:(j-1)))):(1+sum(statesz(1:j))));
end
arvRates = zeros(size(u,1),sn.nstateful,R);
depRates = zeros(size(u,1),sn.nstateful,R);

pi = zeros(1,size(u,1));
for s=1:size(u,1)
    pi(s) = sum(tranState(uj==s,1));
end
SSq = SSq(:,ui)'; % we restrict to unique states in the simulation

for ind=1:sn.nnodes
    if sn.isstateful(ind)
        isf = sn.nodeToStateful(ind);
        if sn.isstation(ind)
            ist = sn.nodeToStation(ind);
            %K = sn.phasessz(ist,:);
            %Ks = sn.phaseshift(ist,:);
        end
        for s=1:size(u,1)
            for r=1:R
                arvRates(s,isf,r) = arvRatesSamples(ui(s),isf,r); % for each unique state, one (any) sample of the rate is enough here
                depRates(s,isf,r) = depRatesSamples(ui(s),isf,r); % for each unique state, one (any) sample of the rate is enough here
            end
        end
    end
end
pi = pi/sum(pi);
%sn.nservers = init_nserver; % restore Inf at delay nodes
end

function print_progress(options,samples_collected)
if options.verbose
    if samples_collected == 1e2
        line_printf(sprintf('\bSSA samples: %6d',samples_collected));
    elseif options.verbose == 2
        if samples_collected == 0
            line_printf(sprintf('\bSSA samples: %6d',samples_collected));
        else
            line_printf(sprintf('\b\b\b\b\b\b%6d',samples_collected));
        end
    elseif mod(samples_collected,1e2)==0 || options.verbose == 2
        line_printf(sprintf('\b\b\b\b\b\b%6d',samples_collected));
    end
end
end