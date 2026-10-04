function [enabled_next_states,enabled_rates,enabled_sync,gctr_start,depRatesSamples,arvRatesSamples,outprob_a,outprob_p,rate_a,eventCache,startRatesSamples,preemptRatesSamples,enabled_tagS,enabled_tagP] = solver_ssa_findenabled(sn,node_a,enabled_next_states,cur_state,outprob_a,event_a,class_a,isSimulation,node_p,local,outprob_p,event_p,class_p,sync,gsync,depRatesSamples,samples_collected,arvRatesSamples,last_node_a,last_node_p,eventCache,startRatesSamples,preemptRatesSamples)
% The START/PREEMPT accumulators and per-transition tags mirror the inlined
% enumeration in solver_ssa.m: this is the SAME enumeration, kept as a
% function, so tagging only one of the two would make the counters depend on
% the use_inline performance flag.
enabled_sync = []; % row is action label, col1=rate, col2=new state
enabled_rates = [];
enabled_tagS = {};
enabled_tagP = {};
R = sn.nclasses;
ctr = 1;
A = length(sync);
G = length(gsync);
% Global (Whittle) rate scaling: constant within a state, so one evaluation
% serves every transition out of it (see SOLVER_SSA_GDFACTOR).
hasGD = isfield(sn,'gdscaling') && ~isempty(sn.gdscaling);
gdNow = [];
if hasGD
    gdNow = solver_ssa_gdfactor(sn, cur_state);
end
for act=1:A
    isf_a = sn.nodeToStateful(node_a{act});
    % try
    %     isf_p = sn.nodeToStateful(node_p{act});
    %     %update_cond_a = true;
    %     %enabled_next_states{act} = cur_state;        
    %     if isempty(enabled_next_states{act}) || ...
    %             isempty(enabled_next_states{act}{isf_a})|| ...
    %             isempty(enabled_next_states{act}{isf_p}) || ...
    %             length(cur_state) < max(isf_a,isf_p)
    %         % no data or a lot has changed
    %         enabled_next_states{act} = cur_state;
    %         update_cond_a = true;
    %     elseif sn.nodetype(node_p{act}) == NodeType.Cache
    %         enabled_next_states{act} = cur_state;
    %         update_cond_a = true;
    %     elseif (isf_a == isf_p) || (length(cur_state{isf_a}) == length(cur_state_1{isf_a}) && ...
    %             length(cur_state{isf_p}) == length(cur_state_1{isf_p}) && ...
    %             all(cur_state{isf_a}==cur_state_1{isf_a}) && ...
    %             all(cur_state{isf_p}==cur_state_1{isf_p}))
    %         % active is unchanged
    %         enabled_next_states{act}{isf_p} = cur_state{isf_p};
    %         update_cond_a = false;
    %     else
    %         enabled_next_states{act} = cur_state;
    %         update_cond_a = true;
    %     end
    % catch
    %     enabled_next_states{act} = cur_state;
    %     update_cond_a = true;
    % end
        enabled_next_states{act} = cur_state;
        update_cond_a = true;
    % see _kb/06-solver-catalog.md for rationale (SSA immfeed self-loop)
    immfeed_selfloop_act = false;
    if event_a{act}==EventType.DEP && node_p{act}==node_a{act} ...
            && node_a{act}>=1 && node_a{act}<=sn.nnodes && sn.isstation(node_a{act}) ...
            && isfield(sn,'immfeed') && ~isempty(sn.immfeed)
        istA_if = sn.nodeToStation(node_a{act});
        if istA_if>=1 && istA_if<=size(sn.immfeed,1) ...
                && class_p{act}>=1 && class_p{act}<=size(sn.immfeed,2) ...
                && sn.immfeed(istA_if, class_p{act})
            immfeed_selfloop_act = true;
        end
    end
    if update_cond_a
        [enabled_next_states{act}{isf_a}, rate_a{act}, outprob_a{act}, eventCache, start_a, preempt_a] =  State.afterEvent(sn, node_a{act}, cur_state{isf_a}, event_a{act}, class_a{act}, isSimulation, eventCache, [], immfeed_selfloop_act);
    end

    if isempty(enabled_next_states{act}{isf_a}) || isempty(rate_a{act})
        continue
    end

    if hasGD && sn.isstation(node_a{act}) && (event_a{act} == EventType.DEP || event_a{act} == EventType.PHASE)
        rate_a{act} = rate_a{act} * gdNow(sn.nodeToStation(node_a{act}), class_a{act});
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
                        [enabled_next_states{act}{isf_p}, ~, outprob_p{act}, eventCache, start_p, preempt_p] =  State.afterEvent(sn, node_p{act}, enabled_next_states{act}{isf_p}, event_p{act}, class_p{act}, isSimulation, eventCache);
                    end
                else % departure
                    isf_p = sn.nodeToStateful(node_p{act});
                    if update_cond_p
                        [enabled_next_states{act}{isf_p}, ~, outprob_p{act}, eventCache, start_p, preempt_p] =  State.afterEvent(sn, node_p{act}, enabled_next_states{act}{isf_p}, event_p{act}, class_p{act}, isSimulation, eventCache);
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
                        % derived tags of this arc, weighted like enabled_rates
                        % and written for every action, not only departures
                        w_tag = rate_a{act}(ia) * prob_sync_p{act};
                        tagS_here = zeros(0,2);
                        tagP_here = zeros(0,2);
                        if ~isempty(start_a) && ia <= size(start_a,1)
                            startRatesSamples(samples_collected,isf_a,:) = reshape(startRatesSamples(samples_collected,isf_a,:),1,R) + w_tag * start_a(ia,:);
                            preemptRatesSamples(samples_collected,isf_a,:) = reshape(preemptRatesSamples(samples_collected,isf_a,:),1,R) + w_tag * preempt_a(ia,:);
                            tagS_here = [tagS_here; sub_tagrows(isf_a, start_a(ia,:))]; %#ok<AGROW>
                            tagP_here = [tagP_here; sub_tagrows(isf_a, preempt_a(ia,:))]; %#ok<AGROW>
                        end
                        if node_p{act} ~= local && exist('start_p','var') && ~isempty(start_p)
                            startRatesSamples(samples_collected,isf_p,:) = reshape(startRatesSamples(samples_collected,isf_p,:),1,R) + w_tag * start_p(1,:);
                            preemptRatesSamples(samples_collected,isf_p,:) = reshape(preemptRatesSamples(samples_collected,isf_p,:),1,R) + w_tag * preempt_p(1,:);
                            tagS_here = [tagS_here; sub_tagrows(isf_p, start_p(1,:))]; %#ok<AGROW>
                            tagP_here = [tagP_here; sub_tagrows(isf_p, preempt_p(1,:))]; %#ok<AGROW>
                        end
                        enabled_rates(ctr) = w_tag;
                        enabled_sync(ctr) = act;
                        enabled_tagS{ctr} = tagS_here;
                        enabled_tagP{ctr} = tagP_here;
                        ctr = ctr + 1;
                    end
                end
            end
        end
    end
end
gctr_start = ctr;

for gact=1:G % event at node ind with global side-effects
    gind = gsync{gact}.active{1}.node; % get the active node (transition) from the gsync event
    [enabled_next_states{A+gact}, outrate, outprob] = State.afterGlobalEvent(sn, gind, cur_state, gsync{gact}, isSimulation);
    for ia=find(outrate .* outprob)
        enabled_rates(ctr) = outrate(ia) * outprob(ia);
        enabled_sync(ctr) = A+gact;
        % an SPN transition holds no server: no service starts there
        enabled_tagS{ctr} = zeros(0,2);
        enabled_tagP{ctr} = zeros(0,2);
        ctr = ctr + 1;
    end
end
end

function rows = sub_tagrows(isf, tagrow)
% ROWS=SUB_TAGROWS(ISF,TAGROW) one [statefulIndex class] row per tagged job.
rows = zeros(0,2);
for r = find(tagrow(:)' > 0)
    for c = 1:tagrow(r)
        rows(end+1,:) = [isf, r]; %#ok<AGROW>
    end
end
end