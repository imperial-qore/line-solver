function updateThinkTimes(self, it)
% Update the think times of all callers at iteration it. The method handles
% differently the case where a caller is a ref task than the case where the
% caller is a queueing station. A coarse heuristic is used when one or more
% callers are themselves infinite servers.

% create local variable due to MATLAB's slow access to self properties
lqn = self.lqn;
idxhash = self.idxhash;
results = self.results;

% main code starts here
if size(lqn.iscaller,2) > 0 % ignore models without callers
    torder = 1:(lqn.ntasks); % set sequential order to update the tasks
    % solve all task models
    for t = torder
        tidx = lqn.tshift + t;
        tidx_thinktime = lqn.think{tidx}.getMean; % user specified think time
        %if ~lqn.isref(tidx) && ~isnan(idxhash(tidx)) % update tasks ignore ref tasks and empty tasks
        if ~isnan(self.idxhash(tidx)) % this skips all REF tasks
            % obtain total self.tput of task t
            % mean throughput of task t in the model where it is a server, summed across replicas
            njobs = max(self.njobs(tidx,:)); % we use njobs to adapt to interlocking corrections
            self.tput(tidx) = lqn.repl(tidx)*sum(results{end,idxhash(tidx)}.TN(self.ensemble{idxhash(tidx)}.attribute.serverIdx,:),2);
            if lqn.sched(tidx) == SchedStrategy.INF % first we consider the update where t is an infinite server
                % obtain total self.utilization of task t
                self.util(tidx) = sum(results{end,idxhash(tidx)}.UN(self.ensemble{idxhash(tidx)}.attribute.serverIdx,:),2);
                % key think time update formula for LQNs, this accounts for the fact that in LINE infinite server self.utilization is dimensionally a mean number of jobs
                self.thinkt(tidx) = max(GlobalConstants.Zero, (njobs-self.util(tidx)) / self.tput(tidx) - tidx_thinktime);
            else % otherwise we consider the case where t is a regular queueing station (other than an infinite server)
                self.util(tidx) = sum(results{end,idxhash(tidx)}.UN(self.ensemble{idxhash(tidx)}.attribute.serverIdx,:),2); % self.utilization of t as a server
                % key think time update formula for LQNs, this accounts that in LINE self.utilization is scaled in [0,1] for all queueing stations irrespectively of the number of servers
                self.thinkt(tidx) = max(GlobalConstants.Zero, njobs*abs(1-self.util(tidx)) / self.tput(tidx) - tidx_thinktime);
            end
            % Recover from Inf/NaN: snap back to previous iteration's value
            if it > 1 && ~isnan(self.thinkt_prev(tidx))
                if isinf(self.thinkt(tidx)) || isnan(self.thinkt(tidx))
                    self.thinkt(tidx) = self.thinkt_prev(tidx);
                end
            end
            % Apply under-relaxation to think time if enabled
            omega = self.relax_omega;
            if omega < 1.0 && it > 1 && ~isnan(self.thinkt_prev(tidx))
                rawT = self.thinkt(tidx);
                prevT = self.thinkt_prev(tidx);
                % If recovering from crash (prev much larger than raw), snap to raw
                if prevT > 10 * rawT && rawT > GlobalConstants.FineTol
                    self.thinkt_prev(tidx) = rawT; % reset prev to allow recovery
                end
                self.thinkt(tidx) = omega * self.thinkt(tidx) + (1 - omega) * self.thinkt_prev(tidx);
            end
            self.thinkt_prev(tidx) = self.thinkt(tidx);
            self.thinktproc{tidx} = Exp.fitMean(self.thinkt(tidx) + tidx_thinktime);
        else % ref task or forwarding target (no task layer)
            % Check if this is a forwarding target task
            isFwdTarget = false;
            if ~lqn.isref(tidx)
                for eidx_fwd = lqn.entriesof{tidx}
                    for cidx_fwd = 1:lqn.ncalls
                        if lqn.calltype(cidx_fwd) == CallType.FWD && lqn.callpair(cidx_fwd, 2) == eidx_fwd
                            isFwdTarget = true;
                            source_eidx = lqn.callpair(cidx_fwd, 1);
                            source_tidx = lqn.parent(source_eidx);
                            fwd_prob = lqn.callproc{cidx_fwd}.getMean();
                            break;
                        end
                    end
                    if isFwdTarget; break; end
                end
            end
            if isFwdTarget
                % see _kb/06-solver-catalog.md (LN section) for rationale
                njobs = max(self.njobs(tidx,:));
                tidx_thinktime = lqn.think{tidx}.getMean;
                arrival_rate = self.tput(source_tidx) * fwd_prob;
                if arrival_rate > GlobalConstants.FineTol && njobs > 0
                    self.tput(tidx) = arrival_rate;
                    % Subtract the processor response time for the target's
                    % activities (already computed by updateMetricsDefault)
                    target_eidx = lqn.callpair(cidx_fwd, 2);
                    host_residt = 0;
                    for aidx_fwd = lqn.actsof{target_eidx}
                        host_residt = host_residt + self.residt(aidx_fwd);
                    end
                    self.thinkt(tidx) = max(GlobalConstants.Zero, njobs / arrival_rate - host_residt - tidx_thinktime);
                else
                    % Source throughput not yet available; use large think time
                    self.thinkt(tidx) = 1000;
                end
                % Apply under-relaxation
                omega = self.relax_omega;
                if omega < 1.0 && it > 1 && ~isnan(self.thinkt_prev(tidx))
                    self.thinkt(tidx) = omega * self.thinkt(tidx) + (1 - omega) * self.thinkt_prev(tidx);
                end
                self.thinkt_prev(tidx) = self.thinkt(tidx);
                self.thinktproc{tidx} = Exp.fitMean(self.thinkt(tidx) + tidx_thinktime);
            else
                self.thinkt(tidx) = GlobalConstants.FineTol;
                self.thinktproc{tidx} = Immediate();
            end
        end
    end
end
end
