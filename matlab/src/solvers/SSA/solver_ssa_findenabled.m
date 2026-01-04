function [enabled_next_states,enabled_rates,enabled_sync,gctr_start,depRatesSamples,arvRatesSamples,outprob_a,outprob_p,rate_a,eventCache] = solver_ssa_findenabled(sn,node_a,enabled_next_states,cur_state,outprob_a,event_a,class_a,isSimulation,node_p,local,outprob_p,event_p,class_p,sync,gsync,depRatesSamples,samples_collected,arvRatesSamples,last_node_a,last_node_p,eventCache)
enabled_sync = []; % row is action label, col1=rate, col2=new state
enabled_rates = [];
ctr = 1;
A = length(sync);
G = length(gsync);
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
    gind = gsync{gact}.active{1}.node; % get the active node (transition) from the gsync event
    [enabled_next_states{A+gact}, outrate, outprob] = State.afterGlobalEvent(sn, gind, cur_state, gsync{gact}, isSimulation);
    for ia=find(outrate .* outprob)
        enabled_rates(ctr) = outrate(ia) * outprob(ia);
        enabled_sync(ctr) = A+gact;
        ctr = ctr + 1;
    end
end
end