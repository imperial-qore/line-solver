function prot = overtake_prob_markov(self, eidx)
% OVERTAKE_PROB_MARKOV Overtaking probability from the phased-server jump chain.
%
% PROT = OVERTAKE_PROB_MARKOV(SELF, EIDX) computes the probability that a new
% arrival to server entry EIDX finds the server busy in phase 2 (post-reply
% processing), from the jump chain of Franks (1999), Sec. 5.4, rather than
% from the reduced 3-state CTMC in OVERTAKE_PROB.
%
% This routine maps the LayeredNetworkStruct and the current per-layer MVA
% residence times onto the chain parameters: the slice count of Eq. (3.1),
% the host residence of each client phase, the calls to the server task and
% to every other task, and the delay incurred at those other tasks. The
% chain itself is solved by LQN_OVERTAKE_MARKOV. Every synchronous caller
% entry contributes a term of Eq. (5.6) and the total is truncated to 1.
%
% Reference: G. Franks, "Performance Analysis of Distributed Server
% Systems", PhD thesis, Carleton University, 1999, Sec. 5.4; published as
% G. Franks and M. Woodside, "Effectiveness of early replies in
% client-server systems", Perform. Eval. 36 (1999) 165-183.

    lqn = self.lqn;
    prot = 0.0;

    % Residence s_jx of the server phase under test (x = 2).
    srvresid = self.servt_ph2(eidx);
    if ~(srvresid > GlobalConstants.FineTol)
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
        maxphase = 1;
        for aidx = acts
            a = aidx - lqn.ashift;
            if a >= 1 && a <= lqn.nacts
                maxphase = max(maxphase, lqn.actphase(a));
            end
        end

        % phasetab rows p=0..maxphase: [nslices service y_ij y_ik t_k].
        nphases = maxphase + 1;
        phasetab = zeros(nphases, 5);
        phasetab(1,1) = 1.0;                          % think slice: nslices=1
        phasetab(1,2) = local_think_time(lqn, ctidx); % service = client think time

        ycalls = zeros(1, maxphase + 1);
        for p = 1:maxphase
            nslices = 1.0; service = 0.0;
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
                    nslices = nslices + y;
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
            phasetab(p+1,:) = [nslices, service, y_ij, y_ik, t_k];
            ycalls(p+1) = y_ij;
            ycalls(1)   = ycalls(1) + y_ij;
        end

        if ycalls(1) == 0.0
            continue;   % this client does not call the server task
        end

        % Client entry visit probability (1 for reference/sole entry).
        prvisit = 1.0;
        if self.tput(ctidx) > GlobalConstants.FineTol && self.tput(ceidx) > GlobalConstants.FineTol
            prvisit = self.tput(ceidx) / self.tput(ctidx);
        end

        prot = prot + lqn_overtake_markov(phasetab, prvisit, srvresid, ycalls);
    end

    prot = max(0.0, min(1.0, prot));
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
