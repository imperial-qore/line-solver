function [pi,SSq,arvRates,depRates,tranSysState,tranSync,sn]=solver_ssa(sn,options)
% [PI,SSQ,ARVRATES,DEPRATES,TRANSYSSTATE,QN]=SOLVER_SSA(QN,OPTIONS)

% Copyright (c) 2012-2024, Imperial College London
% All rights reserved.

% by default the jobs are all initialized in the first valid state

if ~isfield(options,'seed')
    options.seed = 23000;
end
if ~exist('labindex','builtin') || ~exist('labindex','var')
    labindex = 1;
end

Solver.resetRandomGeneratorSeed(options.seed+labindex-1);

%% generate local state spaces
%nstations = sn.nstations;
nstateful = sn.nstateful;
%init_nserver = sn.nservers; % restore Inf at delay nodes
R = sn.nclasses;
N = sn.njobs';
sync = sn.sync;
csmask = sn.csmask;

cutoff = options.cutoff;
if numel(cutoff)==1
    cutoff = cutoff * ones(sn.nstations, sn.nclasses);
end

%%
Np = N';
capacityc = zeros(sn.nnodes, sn.nclasses);
for ind=1:sn.nnodes
    if sn.isstation(ind) % place jobs across stations
        ist = sn.nodeToStation(ind);
        %isf = sn.nodeToStateful(ind);
        for r=1:sn.nclasses %cut-off open classes to finite capacity
            c = find(sn.chains(:,r));
            if ~isempty(sn.visits{c}) && sn.visits{c}(ist,r) == 0
                capacityc(ind,r) = 0;
            elseif ~isempty(sn.proc) && ~isempty(sn.proc{ist}{r}) && any(any(isnan(sn.proc{ist}{r}{1}))) % disabled
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

init_state_hashed = ones(1,nstateful); % pick the first state in space{i}

%%
arvRatesSamples = zeros(options.samples,nstateful,R);
depRatesSamples = zeros(options.samples,nstateful,R);
A = length(sync);
samples_collected = 1;
state = init_state_hashed;
stateCell = cell(nstateful,1);
nir = {};
for ind=1:sn.nnodes
    if sn.isstateful(ind)
        isf = sn.nodeToStateful(ind);
        stateCell{isf} = sn.space{isf}(state(isf),:);
        if sn.isstation(ind)
            ist = sn.nodeToStation(ind);
            [~,nir{ist}] = State.toMarginal(sn, ind, sn.space{isf}(state(isf),:));
            nir{ist} = nir{ist}(:);
        end
    end
end

state = cell2mat(stateCell');
statelen = cellfun(@length, stateCell);
tranSync = zeros(samples_collected,1);
tranState = zeros(1+length(state),samples_collected);
tranState(1:(1+length(state)),1) = [0, state]';
samples_collected = 1;
SSq = cell2mat(nir');
local = sn.nnodes+1;
last_node_a = 0;
last_node_p = 0;
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
newStateCell = cell(1,A);
isSimulation = true; % allow state vector to grow, e.g. for FCFS buffers
cur_time = 0;
while samples_collected < options.samples && cur_time <= options.timespan(2)
    try
    %samples_collected
    ctr = 1;
    enabled_sync = {}; % row is action label, col1=rate, col2=new state
    %enabled_new_state = {};
    enabled_rates = [];
    for act=1:A
        update_cond_a = true; %((node_a{act} == last_node_a || node_a{act} == last_node_p));
        newStateCell{act} = stateCell;
        if update_cond_a || isempty(outprob_a{act})
            isf = sn.nodeToStateful(node_a{act});
            [newStateCell{act}{sn.nodeToStateful(node_a{act})}, rate_a{act}, outprob_a{act}] =  State.afterEvent(sn, node_a{act}, stateCell{isf}, event_a{act}, class_a{act}, isSimulation);
        end
        
        if isempty(newStateCell{act}{sn.nodeToStateful(node_a{act})}) || isempty(rate_a{act})
            continue
        end
        
        for ia=1:size(newStateCell{act}{sn.nodeToStateful(node_a{act})},1) % for all possible new states
            if isnan(rate_a{act}(ia)) || rate_a{act}(ia) == 0 % handles degenerate rate values                                
                if isempty(GlobalConstants.Zero)
                    % the GlobalConstants object does not seem to be accessible in SPMD mode
                    rate_a{act}(ia) = 1e-14; 
                else
                    rate_a{act}(ia) = GlobalConstants.Zero;
                end
            end

            if newStateCell{act}{sn.nodeToStateful(node_a{act})}(ia,:) == -1 % hash not found
                continue
            end
            %update_cond_p = ((node_p{act} == last_node_a || node_p{act} == last_node_p)) || isempty(outprob_p{act});
            update_cond = true; %update_cond_a || update_cond_p;
            
            if rate_a{act}(ia)>0
                if node_p{act} ~= local
                    if node_p{act} == node_a{act} %self-loop
                        if update_cond
                            [newStateCell{act}{sn.nodeToStateful(node_p{act})}, ~, outprob_p{act}] =  State.afterEvent(sn, node_p{act}, newStateCell{act}{sn.nodeToStateful(node_a{act})}, event_p{act}, class_p{act}, isSimulation);
                        end
                    else % departure
                        if update_cond
                            [newStateCell{act}{sn.nodeToStateful(node_p{act})}, ~, outprob_p{act}] =  State.afterEvent(sn, node_p{act}, newStateCell{act}{sn.nodeToStateful(node_p{act})}, event_p{act}, class_p{act}, isSimulation);
                        end
                    end
                    if ~isempty(newStateCell{act}{sn.nodeToStateful(node_p{act})})
                        if sn.isstatedep(node_a{act},3)
                            prob_sync_p{act} = sync{act}.passive{1}.prob(stateCell, newStateCell{act}); %state-dependent
                        else
                            prob_sync_p{act} = sync{act}.passive{1}.prob;
                        end
                    else
                        prob_sync_p{act} = 0;
                    end
                end
                if ~isempty(newStateCell{act}{sn.nodeToStateful(node_a{act})})
                    if node_p{act} == local
                        prob_sync_p{act} = 1; %outprob_a{act}; % was 1
                    end
                    if ~isnan(rate_a{act})
                        if all(~cellfun(@isempty,newStateCell{act}))
                            if event_a{act} == EventType.ID_DEP
                                node_a_sf{act} = sn.nodeToStateful(node_a{act});
                                node_p_sf{act} = sn.nodeToStateful(node_p{act});
                                depRatesSamples(samples_collected,node_a_sf{act},class_a{act}) = depRatesSamples(samples_collected,node_a_sf{act},class_a{act}) + outprob_a{act} * outprob_p{act} * rate_a{act}(ia) * prob_sync_p{act};
                                arvRatesSamples(samples_collected,node_p_sf{act},class_p{act}) = arvRatesSamples(samples_collected,node_p_sf{act},class_p{act}) + outprob_a{act} * outprob_p{act} * rate_a{act}(ia) * prob_sync_p{act};
                            end
                            % simulate also self-loops as we need to log them
                            %if any(~cellfun(@isequal,newStateCell{act},stateCell))
                            if node_p{act} < local && ~csmask(class_a{act}, class_p{act}) && sn.nodetype(node_p{act})~=NodeType.Source && (rate_a{act}(ia) * prob_sync_p{act} >0)
                                line_error(mfilename,sprintf('Error: state-dependent routing at node %d (%s) violates the class switching mask (node %d -> node %d, class %d -> class %d).', node_a{act}, sn.nodenames{node_a{act}}, node_a{act}, node_p{act}, class_a{act}, class_p{act}));
                            end
                            enabled_rates(ctr) = rate_a{act}(ia) * prob_sync_p{act};
                            enabled_sync{ctr} = act;
                            %enabled_new_state{ctr} = new_state;
                            ctr = ctr + 1;
                            %end
                        end
                    end
                end
            end
        end
    end
    
    tot_rate = sum(enabled_rates);
    cum_rate = cumsum(enabled_rates) / tot_rate;
    firing_ctr = 1 + max([0,find( rand > cum_rate )]); % select action
    if isempty(enabled_sync)
        keyboard
        line_error(mfilename,'SSA simulation entered a deadlock before collecting all samples, no synchronization is enabled.');
    end
    last_node_a = node_a{enabled_sync{firing_ctr}};
    last_node_p = node_p{enabled_sync{firing_ctr}};

    % this part is needed to ensure that when the state vector grows the
    % padding of zero is done on the left
    for ind=1:sn.nnodes
        if sn.isstation(ind)
            isf = sn.nodeToStateful(ind);                        
            deltalen = length(stateCell{isf}) > statelen(isf);
            if deltalen>0
                statelen(isf) = length(stateCell{isf});
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
    
    state = cell2mat(stateCell');
    dt = -(log(rand)/tot_rate);
    cur_time = cur_time + dt;

    tranState(1:(1+length(state)),samples_collected) = [dt, state]';
    tranSync(samples_collected,1) = enabled_sync{firing_ctr};
    
    for ind=1:sn.nnodes
        if sn.isstation(ind)
            isf = sn.nodeToStateful(ind);
            ist = sn.nodeToStation(ind);
            [~,nir{ist}] = State.toMarginal(sn, ind, stateCell{isf});
            nir{ist}=nir{ist}(:);
        end
    end
    SSq(:,samples_collected) = cell2mat(nir');
    
    samples_collected = samples_collected + 1;
    stateCell_1 = stateCell;
    stateCell = newStateCell{enabled_sync{firing_ctr}};
    if options.verbose
        if samples_collected == 1e2
            line_printf(sprintf('\b\nSSA samples: %6d\n',samples_collected));
        elseif options.verbose == 2
            if samples_collected == 0
                line_printf(sprintf('\b\nSSA samples: %6d\n',samples_collected));
            else
                line_printf(sprintf('\b\b\b\b\b\b\b%6d\n',samples_collected));
            end
        elseif mod(samples_collected,1e2)==0 || options.verbose == 2
            line_printf(sprintf('\b\b\b\b\b\b\b%6d\n',samples_collected));
        end
    end
    catch ME
        ME.getReport
    end
end

%transient = min([floor(samples_collected/10),1000]); % remove first part of simulation (10% of the samples up to 1000 max)
%transient = 0;
%output = output((transient+1):end,:);
tranState = tranState';


[u,ui,uj] = unique(tranState(:,2:end),'rows');
statesz = cellfun(@length, stateCell_1)';
tranSysState = cell(1,length(stateCell)+1);
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
SSq = SSq(:,ui)';

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
                arvRates(s,isf,r) = arvRatesSamples(ui(s),isf,r); % we just need one sample
                depRates(s,isf,r) = depRatesSamples(ui(s),isf,r); % we just need one sample
            end
        end
    end
end
pi = pi/sum(pi);
%sn.nservers = init_nserver; % restore Inf at delay nodes
end
