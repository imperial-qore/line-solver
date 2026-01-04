function [SSq,SSh,sn] = reachableSpaceGenerator(sn,options)
% [SSQ,SSH,QN] = REACHABLESPACEGENERATOR(QN,OPTIONS)
%
% This differs from spaceGenerator as it is restricted to states reachable
% from the initial state.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

% By default the jobs are all initialized in the first valid state

%nstations = sn.nstations;
nstateful = sn.nstateful;
%init_nserver = sn.nservers; % restore Inf at delay nodes
R = sn.nclasses;
N = sn.njobs';
sync = sn.sync;
csmask = sn.csmask;
stack = {sn.state'};
SSq = [];
A = length(sn.sync);
sync = sn.sync;
isSimulation = false;
local = sn.nnodes+1;

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

space = cell(1,nstateful);
for i=1:nstateful
    space{i} = sn.state{i};
end
SSh = ones(1,nstateful); % initial state
it = 0;
%ih = 1;
stack_index = 1;
maxstatesz = zeros(1,nstateful);
% Q = sparse([]);
% Dfilt = cell(1,A);
% for a=1:A
%     Dfilt{a} = sparse([]);
% end
while 1
    if isempty(stack)
        %Q = ctmc_makeinfgen(Q);
        SSq = []; % initial state
        for i=1:size(SSh,1)
            colctr = 1;
            for j=1:nstateful
                SSq(i,colctr:(colctr+size(space{j},2)-1)) = space{j}(SSh(i,j),:);
                colctr = colctr+size(space{j},2);
            end
        end
        %         for a=1:A
        %             if (size(Dfilt{a},1) < size(SSq,1)) || (size(Dfilt{a},2) < size(SSq,1))
        %                 Dfilt{a}(size(SSq,1),size(SSq,1)) = 0;
        %             end
        %         end
        sn.space = space;
        return
    end
    % pop state
    stateCell = stack{end};
    state = cell2mat(stateCell);
    ih = stack_index(end);
    stack_index(end)=[];
    stack(end)=[];
    statelen = cellfun(@length, stateCell);
    newStateCell = cell(1,A);
    %samples_collected
    ctr = 1;
    enabled_sync = {}; % row is action label, col1=rate, col2=new state
    %enabled_new_state = {};
    enabled_rates = [];
    for act=1:A
        outprob_a{act} = [];
        outprob_p{act} = [];
        update_cond_a = true; %((node_a{act} == last_node_a || node_a{act} == last_node_p));
        newStateCell{act} = stateCell;
        if update_cond_a || isempty(outprob_a{act})
            isf = sn.nodeToStateful(node_a{act});
            [newStateCell{act}{sn.nodeToStateful(node_a{act})}, rate_a{act}, outprob_a{act}] =  State.afterEvent(sn, node_a{act}, stateCell{isf}, event_a{act}, class_a{act}, isSimulation);
            if size(newStateCell{act}{sn.nodeToStateful(node_a{act})},1) == 0
                outprob_a{act} = [];
            end
        end

        if isempty(newStateCell{act}{sn.nodeToStateful(node_a{act})}) || isempty(rate_a{act}) % state not found
            continue
        end

        for ia=1:size(newStateCell{act}{sn.nodeToStateful(node_a{act})},1) % for all possible new states of the active agent
            if isempty(newStateCell{act}{sn.nodeToStateful(node_a{act})})
                continue
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
                            % since self-loop sare multiple instantaneous
                            % state transitions, we save the intermediate
                            % states in the sn.space structure
                            i = node_p{act};
                            nci = newStateCell{act}{sn.nodeToStateful(node_a{act})};
                            for j=1:size(nci,1)
                                if length(nci(j,:)) > size(space{i},2)
                                    space{i} = [zeros(size(space{i},1),length(nci(j,:))-size(space{i},2)),space{i}];
                                elseif length(nci(j,:)) < size(space{i},2)
                                    nci = [zeros(1, size(space{i},2)-length(nci(j,:))),nci(j,:)];
                                end
                                if matchrow(space{i},nci(j,:))==-1
                                    space{i}(end+1,:) = nci(j,:);
                                end
                            end
                            [newStateCell{act}{sn.nodeToStateful(node_p{act})}, ~, outprob_p{act}] =  State.afterEvent(sn, node_p{act}, newStateCell{act}{sn.nodeToStateful(node_a{act})}, event_p{act}, class_p{act}, isSimulation);
                            outprob_a{act} = outprob_p{act};
                        end
                    else % departure from active
                        if update_cond
                            [newStateCell{act}{sn.nodeToStateful(node_p{act})}, ~, outprob_p{act}] =  State.afterEvent(sn, node_p{act}, stateCell{sn.nodeToStateful(node_p{act})}, event_p{act}, class_p{act}, isSimulation);
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
                            if event_a{act} == EventType.DEP
                                node_a_sf{act} = sn.nodeToStateful(node_a{act});
                                node_p_sf{act} = sn.nodeToStateful(node_p{act});
                            end
                            % simulate also self-loops as we need to log them
                            if node_p{act} < local && ~csmask(class_a{act}, class_p{act}) && sn.nodetype(node_p{act})~=NodeType.Source && (rate_a{act}(ia) * prob_sync_p{act} >0)
                                line_error(mfilename,sprintf('Error: state-dependent routing at node %d (%s) violates the class switching mask (node %d -> node %d, class %d -> class %d).', node_a{act}, sn.nodenames{node_a{act}}, node_a{act}, node_p{act}, class_a{act}, class_p{act}));
                            end
                            enabled_rates(ctr) = rate_a{act}(ia) * prob_sync_p{act};
                            enabled_sync{ctr} = act;
                            ctr = ctr + 1; % keep
                        end
                    end
                end
            end
        end
    end

    for firing_ctr=1:length(enabled_rates)
        firing_rate = enabled_rates(firing_ctr);
        last_node_a = node_a{enabled_sync{firing_ctr}};
        last_node_p = node_p{enabled_sync{firing_ctr}};
        act = enabled_sync{firing_ctr};
        netstates = newStateCell{act};
        if enabled_rates(firing_ctr)>0 && ~any(cellfun(@isempty,netstates))
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            nvec = cellfun(@(c) size(c,1)-1, netstates);
            n = pprod(nvec);
            while n>=0
                newstate = [];
                newstatec = {};
                firingprob = 1;
                for i=1:size(netstates,2)
                    if i==last_node_a
                        firingprob = firingprob*outprob_a{enabled_sync{firing_ctr}}(1+n(i));
                    end
                    if i==last_node_p
                        firingprob = firingprob*outprob_p{enabled_sync{firing_ctr}}(1+n(i));
                    end
                    maxstatesz(i) = max(maxstatesz(i),length(netstates{i}(1+n(i),:)));
                    newstate = [newstate,zeros(1,maxstatesz(i)-length(netstates{i}(1+n(i),:))),netstates{i}(1+n(i),:)];
                    newstatec{i} = [zeros(1,maxstatesz(i)-length(netstates{i}(1+n(i),:))),netstates{i}(1+n(i),:)];
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if firingprob>0
                    hashednewstate = zeros(1,nstateful);
                    for i=1:nstateful
                        if length(newstatec{i}) > size(space{i},2)
                            hashednewstate(i) = matchrow([zeros(size(space{i},1),length(newstatec{i})-size(space{i},2)),space{i}], newstatec{i});
                        else
                            hashednewstate(i) = matchrow(space{i}, newstatec{i});
                        end
                    end
                    jh =  matchrow(SSh,hashednewstate);
                    if jh == -1 || (length(newstate)~=length(state) || any(newstate~=state))
                        % store and update jh index
                        if jh == -1
                            for i=1:nstateful
                                if hashednewstate(i) == -1
                                    nci = newstatec{i};
                                    if length(nci) > size(space{i},2)
                                        space{i} = [zeros(size(space{i},1),length(nci)-size(space{i},2)),space{i}];
                                    elseif length(nci) < size(space{i},2)
                                        nci = [zeros(1, size(space{i},2)-length(nci)),nci];
                                    end
                                    space{i}(end+1,:) = nci;
                                    hashednewstate(i) = size(space{i},1);
                                end
                            end
                            SSh = [SSh; hashednewstate];
                            jh = size(SSh,1);
                            stack{end+1} = newstatec;
                            stack_index(end+1) = jh;
                        end
                    end
                    % for some reason the Q generation is buggy so we
                    % disabled it
                    %                     if (size(Q,1)< ih) || (size(Q,2)< jh)
                    %                         Q(ih,jh) = firing_rate;
                    %                         for a=1:A
                    %                             Dfilt{act}(ih,jh) = 0;
                    %                         end
                    %                         Dfilt{act}(ih,jh) = firing_rate;
                    %                     else
                    %                         Q(ih,jh) = Q(ih,jh) + firing_rate;
                    %                         for a=1:A
                    %                             if (size(Dfilt{act},1) < ih) || (size(Dfilt{act},2) < jh)
                    %                                 Dfilt{act}(ih,jh) = 0;
                    %                             end
                    %                         end
                    %                         Dfilt{act}(ih,jh) = Dfilt{act}(ih,jh) + firing_rate;
                    %                     end
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                n = pprod(n,nvec);
            end
        end
    end
    it = it + 1;
end
end
