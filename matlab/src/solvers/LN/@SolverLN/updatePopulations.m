function updatePopulations(self, it)
% Update the populations in each layer to address interlocking
lqn = self.lqn;
ilscaling = ones(lqn.nhosts+ lqn.ntasks, lqn.nhosts+ lqn.ntasks); % interlock scaling factors

% we now define from where we count the multiplicity of a caller tasks
% self.njobsorig(caller_tidx,idx) is the multiplicity in the upper layer
% of a caller to the server in model idx

call_mult_count = self.njobsorig;

for h=1:lqn.nhosts
    minremote = Inf;
    for hops = 1:self.nlayers
        hidx = h;
        ilscaling(hidx) = 1.0;
        % the following are remote (indirect) callers that are certain to be
        % callers of task t, hence if they have multiplicity m then task t
        % cannot have as a matter of fact multiplicity more than m
        callers = lqn.tasksof{hidx};
        caller_conn_components = lqn.conntasks(callers-lqn.tshift);
        multcallers = sum(call_mult_count(callers,hidx));
        indirect_callers = find(self.ptaskcallers_step{hops}(hidx,:));
        multremote = 0;
        for remidx=indirect_callers(:)'
            % first we consider the update where the remote caller is an infinite server
            % but since the ref task has finite multiplicity it is treated
            % as normal
            if lqn.sched(remidx) == SchedStrategy.INF && ~lqn.isref(remidx) 
                % first we consider the update where the remote caller is an infinite server but we do not scale populations if a ref task is a direct client
                multremote = Inf; % do not apply interlock correction
            else
                % now we multiply the probability that a request to hidx
                % originates from remidx
                multremote = multremote + self.ptaskcallers_step{hops}(hidx,remidx)*call_mult_count(remidx,hidx);
            end
        end
        if multcallers > multremote && multremote > GlobalConstants.CoarseTol && ~isinf(multremote) && multremote < minremote
            minremote = multremote;
            % [multcallers, multremote]
            % we spread the scaling proportionally to the direct
            % caller probabilities
            caller_spreading_ratio = self.ptaskcallers(hidx,callers); % this a probability vector so no further renormalization is needed
            for u=unique(caller_conn_components)
                caller_spreading_ratio(caller_conn_components==u) = caller_spreading_ratio(caller_conn_components==u)/sum(caller_spreading_ratio(caller_conn_components==u));
            end
            for c=callers
                ilscaling(c,hidx) = min(1, multremote / multcallers .* caller_spreading_ratio(c==callers));
            end
        end
    end
end

maxhops = self.nlayers; % max backward hops in interlocking correction
%maxhops = 1; % max backward hops in interlocking correction
for t=1:lqn.ntasks
    minremote = Inf;
    for hops = 1:maxhops
        tidx = lqn.tshift + t;
        if ~lqn.isref(tidx)
            % the following are remote (indirect) callers that are certain to be
            % callers of task t, hence if they have multiplicity m then task t
            % cannot have as a matter of fact multiplicity more than m
            [calling_idx, called_entries] = find(lqn.iscaller(:, lqn.entriesof{tidx})); %#ok<ASGLU>
            callers = intersect(lqn.tshift+(1:lqn.ntasks), unique(calling_idx)');
            caller_conn_components = lqn.conntasks(callers-lqn.tshift);
            multcallers = sum(call_mult_count(callers,tidx));
            indirect_callers = find(self.ptaskcallers_step{hops}(tidx,:)); % caller at step hops from the node
            multremote = 0;
            for remidx=indirect_callers(:)'
                if lqn.sched(remidx) == SchedStrategy.INF % first we consider the update where t is an infinite server
                    multremote = Inf; % do not apply interlock correction
                else
                    multremote = multremote + self.ptaskcallers_step{hops}(tidx,remidx)*call_mult_count(remidx,tidx);
                end
            end
            if multcallers > multremote && multremote > GlobalConstants.CoarseTol && ~isinf(multremote)  && multremote < minremote
                minremote = multremote;
                % [multcallers, multremote]
                % we spread the scaling proportionally to the direct
                % caller probabilities
                caller_spreading_ratio = self.ptaskcallers(tidx,callers); % this a probability vector so no further renormalization is needed
                for u=unique(caller_conn_components)
                    caller_spreading_ratio(caller_conn_components==u) = caller_spreading_ratio(caller_conn_components==u)/sum(caller_spreading_ratio(caller_conn_components==u));
                end
                for c=callers
                    ilscaling(c,tidx) = min(1, multremote / multcallers .* caller_spreading_ratio(c==callers));
                end
            end
        end
    end
end
self.ilscaling = ilscaling;
self.njobs = self.njobsorig .* ilscaling;
end
