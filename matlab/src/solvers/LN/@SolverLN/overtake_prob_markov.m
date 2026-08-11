function prOt = overtake_prob_markov(self, eidx)
% OVERTAKE_PROB_MARKOV Overtaking probability via the LQNS phased-server Markov chain.
%
% PROT = OVERTAKE_PROB_MARKOV(SELF, EIDX) computes the probability that a new
% arrival to server entry EIDX finds the server busy in phase 2 (post-reply
% processing), using the LQNS V6 slice/overtake Markov chain (Layer 1) rather
% than the reduced 3-state CTMC in OVERTAKE_PROB.
%
% This is the input-mapping layer (Layer 1c): it derives the per-client-phase
% slice parameters (nSlices, host residence, calls to the server task, calls to
% other tasks, delay at other tasks) from the LayeredNetworkStruct and the
% current per-layer MVA residence times, then delegates the chain solution to
% LQN_OVERTAKE_MARKOV. Contributions of all synchronous caller entries are
% summed and truncated to 1 (cf. LQNS Markov_Phased_Server::PrOT_e).
%
% Reference: Franks & Woodside, "Effectiveness of early replies in client-server
% systems", Perf. Eval. 36 (1999). See [[overtaking-markov-port]].

    lqn = self.lqn;
    prOt = 0.0;

    % Server phase-2 residence (x_j for the tested server phase j=2).
    xj = self.servt_ph2(eidx);
    if ~(xj > GlobalConstants.FineTol)
        return;
    end
    server_tidx = lqn.parent(eidx);

    % Synchronous caller activities into this server entry.
    caller_acts = full(lqn.callpair(lqn.callpair(:,2) == eidx, 1));
    if isempty(caller_acts)
        return;
    end

    % Group caller activities by their owning client entry.
    caller_entries = [];
    for ci = 1:numel(caller_acts)
        aidx = caller_acts(ci);
        ceidx = local_entry_of_activity(lqn, aidx);
        if ceidx > 0 && ~any(caller_entries == ceidx)
            caller_entries(end+1) = ceidx; %#ok<AGROW>
        end
    end

    for ceidx = caller_entries
        ctidx = lqn.parent(ceidx);
        acts = lqn.actsof{ceidx};
        if isempty(acts)
            continue;
        end

        % Maximum client phase (LINE activities carry phase 1 or 2).
        maxPhaseA = 1;
        for aidx = acts
            a = aidx - lqn.ashift;
            if a >= 1 && a <= lqn.nacts
                maxPhaseA = max(maxPhaseA, lqn.actphase(a));
            end
        end

        % clientPhases rows p=0..maxPhaseA: [nSlices service y_ij y_ik t_k].
        nStates = maxPhaseA + 1;
        clientPhases = zeros(nStates, 5);
        clientPhases(1,1) = 1.0;                          % think slice: nSlices=1
        clientPhases(1,2) = local_think_time(lqn, ctidx); % service = client think time

        y_aj = zeros(1, maxPhaseA + 1);
        for p = 1:maxPhaseA
            nSlices = 1.0; service = 0.0;
            y_ij = 0.0; y_ik = 0.0; tk_num = 0.0;
            for aidx = acts
                a = aidx - lqn.ashift;
                if a < 1 || a > lqn.nacts || lqn.actphase(a) ~= p
                    continue;
                end
                service = service + self.servt(aidx);     % host residence of the phase
                for c = lqn.callsof{aidx}
                    if lqn.calltype(c) ~= CallType.SYNC
                        continue;
                    end
                    y = lqn.callproc_mean(c);
                    if y == 0, continue; end
                    nSlices = nSlices + y;
                    dst_tidx = lqn.parent(lqn.callpair(c,2));
                    if dst_tidx == server_tidx
                        y_ij = y_ij + y;
                    else
                        y_ik = y_ik + y;
                        tk_num = tk_num + y * self.callresidt(c);  % rendezvous delay
                    end
                end
            end
            if y_ik > 0.0, t_k = tk_num / y_ik; else, t_k = 0.0; end
            clientPhases(p+1,:) = [nSlices, service, y_ij, y_ik, t_k];
            y_aj(p+1) = y_ij;
            y_aj(1)   = y_aj(1) + y_ij;
        end

        if y_aj(1) == 0.0
            continue;   % this client does not call the server task
        end

        % Client entry visit probability (1 for reference/sole entry).
        prVisit = 1.0;
        if self.tput(ctidx) > GlobalConstants.FineTol && self.tput(ceidx) > GlobalConstants.FineTol
            prVisit = self.tput(ceidx) / self.tput(ctidx);
        end

        prOt = prOt + lqn_overtake_markov(clientPhases, prVisit, xj, y_aj);
    end

    prOt = max(0.0, min(1.0, prOt));
end

% ---------------------------------------------------------------------------
function ceidx = local_entry_of_activity(lqn, aidx)
% Absolute entry index owning activity aidx (via actsof), else -1.
    ceidx = -1;
    for e = 1:lqn.nentries
        eabs = lqn.eshift + e;
        if any(lqn.actsof{eabs} == aidx)
            ceidx = eabs;
            return;
        end
    end
end

function z = local_think_time(lqn, tidx)
% Client task think time (0 when unspecified).
    z = 0.0;
    if isfield(lqn, 'think_mean') && numel(lqn.think_mean) >= tidx
        v = lqn.think_mean(tidx);
        if isfinite(v) && v > 0, z = v; end
    end
end
